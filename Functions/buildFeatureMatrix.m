function [X, feature_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, waveforms, acgs, sr_wf, acg_bin_size, acg_lag)
% BUILDFEATUREMATRIX  Build a classifier-compatible feature matrix from raw waveforms and ACGs.
%
% Converts raw waveforms and ACGs (from DeePhys units or any external source)
% into the feature format used by the CellTypeClassifier. All data is harmonized
% to ctc.Parameters.Harmonization regardless of its origin.
%
% Waveform pipeline:
%   1. Resample to reference rate (if sr_wf differs from DeePhys recording rate)
%   2. Trim / zero-pad to reference sample count
%   3. Interpolate to WaveformTargetSamplingRate with makima spline
%   4. Normalise each waveform by max(abs)
%   5. Circular-shift to align troughs to population mean position
%   6. Trim 5×interp_factor samples from each edge
%
% ACG pipeline:
%   1. Validate bin count against acg_bin_size, acg_lag, and Harmonization params
%   2. Normalise each column by its maximum
%
% Feature order: FullACG1…N, Waveform1…M  (matches DeePhys convention)
%
% INPUTS:
%   ctc          - CellTypeClassifier (provides Harmonization parameters)
%   waveforms    - (N_samples × N_units) raw, unnormalised waveforms
%                  Waveforms at a different rate are resampled automatically.
%   acgs         - (N_bins × N_units) autocorrelograms (central peak NOT removed)
%                  Must satisfy: N_bins == round(2*acg_lag/acg_bin_size) + 1
%   sr_wf        - waveform sampling rate in Hz
%   acg_bin_size - ACG bin size in seconds   (default: ctc.Parameters.Harmonization.ACGBinSize)
%   acg_lag      - one-sided ACG lag in seconds (default: ctc.Parameters.Harmonization.ACGLag)
%
% OUTPUTS:
%   X             - (N_units × N_features) raw (un-z-scored) feature matrix
%   feature_names - (1 × N_features) string column names
%   aligned_wf    - (N_wf_samples × N_units) processed waveforms (upsampled, aligned, trimmed)
%   norm_acgs     - (N_bins × N_units) normalised ACGs

arguments
    ctc          CellTypeClassifier
    waveforms    (:,:) double
    acgs         (:,:) double
    sr_wf        (1,1) double
    acg_bin_size (1,1) double = ctc.Parameters.Harmonization.ACGBinSize
    acg_lag      (1,1) double = ctc.Parameters.Harmonization.ACGLag
end

ph = ctc.Parameters.Harmonization;

% ── Reference waveform parameters ───────────────────────────────────────────
if ~isempty(ctc.UnitList)
    ref_unit = ctc.UnitList(1);
    sr_ref   = ref_unit.MEArecording.RecordingInfo.SamplingRate;
    n_ref    = length(ref_unit.ReferenceWaveform);
    fprintf('Waveform reference from ctc.UnitList: sr=%g Hz, %i samples\n', sr_ref, n_ref);
else
    sr_ref = ctc.Parameters.External.WaveformSamplingRate;
    n_ref  = ctc.Parameters.External.WaveformNSamples;
    fprintf('Waveform reference from ctc.Parameters.External: sr=%g Hz, %i samples\n', sr_ref, n_ref);
end

% ── Input orientation: enforce (N_samples × N_units) and (N_bins × N_units) ─
% The expected waveform sample count is known from the reference sampling rate
% and the standard Kilosort template width (82 samples at 30 kHz). Using this
% avoids the ambiguous heuristic when N_units > N_samples.
[r_wf, c_wf]   = size(waveforms);
[r_acg, c_acg] = size(acgs);

expected_wf_samples = n_ref;  % from ctc.Parameters or the DeePhys default
if r_wf ~= expected_wf_samples && c_wf == expected_wf_samples
    waveforms = waveforms'; [~, c_wf] = size(waveforms);
elseif r_wf ~= expected_wf_samples && c_wf ~= expected_wf_samples
    % Fallback: use the heuristic but warn
    warning('buildFeatureMatrix:orient', ...
        'Neither waveform dimension matches expected %i samples — using size heuristic.', expected_wf_samples);
    if c_wf > r_wf, waveforms = waveforms'; [~, c_wf] = size(waveforms); end
end
% ACG orientation: use the expected bin count
n_bins_expected_orient = round(2 * acg_lag / acg_bin_size) + 1;
if r_acg ~= n_bins_expected_orient && c_acg == n_bins_expected_orient
    acgs = acgs'; [~, c_acg] = size(acgs);
elseif r_acg ~= n_bins_expected_orient && c_acg ~= n_bins_expected_orient
    warning('buildFeatureMatrix:orient', ...
        'Neither ACG dimension matches expected %i bins — using size heuristic.', n_bins_expected_orient);
    if c_acg > r_acg, acgs = acgs'; [~, c_acg] = size(acgs); end
end

assert(c_wf == c_acg, ...
    'Inconsistent number of units: waveforms has %i, acgs has %i (after auto-orient).', ...
    c_wf, c_acg);
N_units = c_wf;

% ── Validate ACG parameters ──────────────────────────────────────────────────
n_bins_in       = size(acgs, 1);
n_bins_expected = round(2 * acg_lag / acg_bin_size) + 1;

assert(n_bins_in == n_bins_expected, ...
    ['ACG bin count mismatch: got %i bins but acg_bin_size=%.4f s, acg_lag=%.4f s implies %i bins.\n' ...
     'Check acg_bin_size and acg_lag arguments.'], ...
    n_bins_in, acg_bin_size, acg_lag, n_bins_expected);

n_bins_harm = round(2 * ph.ACGLag / ph.ACGBinSize) + 1;
assert(n_bins_in == n_bins_harm, ...
    ['ACG bin count (%i) does not match Harmonization params (%i bins).\n' ...
     'Harmonization: ACGBinSize=%.4f s, ACGLag=%.4f s → %i bins.\n' ...
     'Recompute your ACGs with the matching parameters or update ctc.Parameters.Harmonization.'], ...
    n_bins_in, n_bins_harm, ph.ACGBinSize, ph.ACGLag, n_bins_harm);

% ── Resample waveforms to reference rate ─────────────────────────────────────
if sr_wf ~= sr_ref
    [p, q]    = rat(sr_ref / sr_wf, 1e-6);
    waveforms = resample(double(waveforms), p, q);
    fprintf('Resampled waveforms %g Hz → %g Hz (factor %i/%i)\n', sr_wf, sr_ref, p, q);
end

% Trim / zero-pad to exactly n_ref samples
n_got = size(waveforms, 1);
if n_got > n_ref
    waveforms = waveforms(1:n_ref, :);
elseif n_got < n_ref
    waveforms = [waveforms; zeros(n_ref - n_got, N_units)];
end

% ── Waveform alignment ───────────────────────────────────────────────────────
%  1. Interpolate to WaveformTargetSamplingRate with makima spline
%  2. Normalise each waveform by max(abs)
%  3. Circular-shift to align troughs to population mean position
%  4. Trim 5×interp_factor samples from each edge
target_sr         = ph.WaveformTargetSamplingRate;
interp_factor_raw = target_sr / sr_ref;
interp_factor     = round(interp_factor_raw);
if abs(interp_factor_raw - interp_factor) > 1e-9
    warning('buildFeatureMatrix:nonIntegerFactor', ...
        ['WaveformTargetSamplingRate (%g Hz) is not an integer multiple of the reference rate (%g Hz). ' ...
        'Rounding factor from %.6f to %i (effective: %g Hz).'], ...
        target_sr, sr_ref, interp_factor_raw, interp_factor, sr_ref * interp_factor);
end

ref_wf = double(waveforms);
ref_wf(sum(ref_wf, 2) == 0, :) = [];    % drop all-zero rows

[~, trough_idx] = min(ref_wf, [], 1);
peak_idx        = mean(trough_idx);
max_offset      = round(peak_idx / 2);
n_samp          = size(ref_wf, 1);
buffer_n        = n_samp + 2 * max_offset;

x  = max_offset : (n_samp + max_offset - 1);
xq = linspace(1, buffer_n, buffer_n * interp_factor);

interp_wf = interp1(x, ref_wf, xq, 'makima');
interp_wf = interp_wf( (max_offset*interp_factor + 1) : ((buffer_n-max_offset-1)*interp_factor), :);
interp_wf = interp_wf ./ max(abs(interp_wf), [], 1);

[~, min_idx] = min(interp_wf, [], 1);
shifts = round(peak_idx * interp_factor) - min_idx;
for u = 1:N_units
    interp_wf(:, u) = circshift(interp_wf(:, u), shifts(u));
end

aligned_wf = interp_wf((5*interp_factor) : (end - 5*interp_factor), :);

% ── ACG normalisation ────────────────────────────────────────────────────────
norm_acgs = acgs ./ max(acgs, [], 1);
norm_acgs(~isfinite(norm_acgs)) = 0;

% ── Assemble feature matrix (FullACG first, Waveform second) ─────────────────
acg_names    = "FullACG"  + (1:size(norm_acgs, 1));
wf_names     = "Waveform" + (1:size(aligned_wf, 1));
feature_names = [acg_names, wf_names];

X = [norm_acgs', aligned_wf'];

fprintf('Feature matrix: %i units × %i features (%i ACG + %i waveform)\n', ...
    N_units, size(X, 2), size(norm_acgs, 1), size(aligned_wf, 1));
end
