function [X, feature_names, aligned_wf, norm_acgs, wf_window] = buildFeatureMatrix(ctc, waveforms, acgs, sr_wf, acg_bin_size, acg_lag)
% BUILDFEATUREMATRIX  Build a classifier-compatible feature matrix from raw waveforms and ACGs.
%
% Converts raw waveforms and ACGs (from DeePhys units or any external source)
% into the feature format used by the CellTypeClassifier. All data is harmonized
% to ctc.Parameters.Harmonization regardless of its origin.
%
% Waveform pipeline (per unit):
%   1. Find trough (minimum) in raw waveform
%   2. Build time axis relative to trough (trough = 0 ms)
%   3. Interpolate onto fixed output grid [−PreTrough, +PostTrough] at TargetSR
%   4. Normalise by max(abs)
%
% Output waveform dimensions are determined entirely by Harmonization params
% (WaveformPreTrough, WaveformPostTrough, WaveformTargetSamplingRate) and
% are independent of input sampling rate or sample count.
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

% ── Harmonized output parameters (independent of input source) ─────────────
% All waveforms are resampled to target_sr and trimmed to a fixed time
% window around the trough. This guarantees identical output dimensions
% regardless of input sampling rate or sample count.
target_sr     = ph.WaveformTargetSamplingRate;
pre_trough_ms  = ph.WaveformPreTrough;   % ms before trough (positive number, e.g. 1.0)
post_trough_ms = ph.WaveformPostTrough;  % ms after trough  (positive number, e.g. 0.5)
n_out_pre      = round(pre_trough_ms  * target_sr / 1000);  % samples before trough at target rate
n_out_post     = round(post_trough_ms * target_sr / 1000);  % samples after  trough at target rate
n_out          = n_out_pre + n_out_post + 1;                 % total output samples
fprintf('Waveform harmonization: [-%g, +%g] ms @ %g Hz → %d samples\n', ...
    pre_trough_ms, post_trough_ms, target_sr, n_out);

% ── Input orientation: enforce (N_samples × N_units) and (N_bins × N_units) ─
% Use expected ACG bin count as the reliable anchor; waveform orientation
% follows from matching unit count with ACGs.
[r_wf, c_wf]   = size(waveforms);
[r_acg, c_acg] = size(acgs);

n_bins_expected = round(2 * acg_lag / acg_bin_size) + 1;
if r_acg ~= n_bins_expected && c_acg == n_bins_expected
    acgs = acgs'; [r_acg, c_acg] = size(acgs);
elseif r_acg ~= n_bins_expected && c_acg ~= n_bins_expected
    warning('buildFeatureMatrix:orient', ...
        'Neither ACG dimension matches expected %i bins — using size heuristic.', n_bins_expected);
    if c_acg > r_acg, acgs = acgs'; [r_acg, c_acg] = size(acgs); end
end

% Orient waveforms to match ACG unit count
if c_wf == c_acg
    % Already correct orientation
elseif r_wf == c_acg
    waveforms = waveforms'; [r_wf, c_wf] = size(waveforms);
else
    % Neither dimension matches — try smaller-dim-is-samples heuristic
    if c_wf > r_wf, waveforms = waveforms'; [r_wf, c_wf] = size(waveforms); end
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

% ── Waveform harmonization ────────────────────────────────────────────────────
% Goal: produce aligned waveforms at target_sr covering a fixed time window
%       around the trough. Output dimensions depend on edge_mode:
%   "zero" / "edge" → fixed (n_out × N_units) using full requested window
%   "trim"          → (n_trim × N_units) shrunk to the shortest available span
%
% Strategy per unit:
%   1. Find trough (minimum) sample index
%   2. Build time axis relative to trough (trough = 0 ms)
%   3. Interpolate onto the output time grid at target_sr
%   4. Handle out-of-range samples according to edge_mode
%   5. Normalise by max(abs)

wf_raw    = double(waveforms);
dt_in     = 1000 / sr_wf;                                  % input sample spacing in ms
dt_out    = 1000 / target_sr;                               % output sample spacing in ms
edge_mode = ph.WaveformEdgeMode;                            % "zero" | "edge" | "trim"

% ── Pass 1 (trim mode only): find tightest available window across all units
if edge_mode == "trim"
    global_pre_ms  = pre_trough_ms;
    global_post_ms = post_trough_ms;
    for u = 1:N_units
        wf_u = wf_raw(:, u);
        if all(wf_u == 0); continue; end
        [~, trough_samp] = min(wf_u);
        avail_pre_ms  = (trough_samp - 1) * dt_in;         % ms before trough
        avail_post_ms = (size(wf_u,1) - trough_samp) * dt_in; % ms after trough
        global_pre_ms  = min(global_pre_ms,  avail_pre_ms);
        global_post_ms = min(global_post_ms, avail_post_ms);
    end
    % Recompute output grid with trimmed window
    n_out_pre  = round(global_pre_ms  * target_sr / 1000);
    n_out_post = round(global_post_ms * target_sr / 1000);
    n_out      = n_out_pre + n_out_post + 1;
    fprintf('Trim mode: window narrowed to [-%g, +%g] ms → %d samples\n', ...
        global_pre_ms, global_post_ms, n_out);
end

t_out = (-n_out_pre : n_out_post)' * dt_out;               % (n_out × 1) output time grid in ms
aligned_wf = zeros(n_out, N_units);
n_truncated = 0;

for u = 1:N_units
    wf_u = wf_raw(:, u);

    % Skip all-zero waveforms
    if all(wf_u == 0)
        continue
    end

    % Find trough (minimum) — assumes negative-polarity waveforms (Maxwell HDMEA).
    % Other acquisition systems may require polarity inversion before this step.
    [~, trough_samp] = min(wf_u);

    % Build input time axis relative to trough (ms)
    n_in   = size(wf_u, 1);
    t_in   = ((1:n_in)' - trough_samp) * dt_in;            % trough at t=0

    % Determine which output samples fall within the input range
    t_min_avail = t_in(1);
    t_max_avail = t_in(end);
    in_range = (t_out >= t_min_avail) & (t_out <= t_max_avail);

    if sum(in_range) < 3
        warning('buildFeatureMatrix:shortWaveform', ...
            'Unit %d: input waveform too short (%.2f ms available, [%.1f, %.1f] ms requested) — zeroed.', ...
            u, t_max_avail - t_min_avail, -pre_trough_ms, post_trough_ms);
        continue
    end

    % Interpolate within available range
    wf_interp = zeros(n_out, 1);
    wf_interp(in_range) = interp1(t_in, wf_u, t_out(in_range), 'makima');

    % Fill out-of-range samples according to edge_mode
    if ~all(in_range) && edge_mode == "edge"
        % Hold boundary values (nearest-neighbor extrapolation)
        first_valid = find(in_range, 1, 'first');
        last_valid  = find(in_range, 1, 'last');
        wf_interp(1:first_valid-1)   = wf_interp(first_valid);
        wf_interp(last_valid+1:end)  = wf_interp(last_valid);
    end
    % "zero" mode: wf_interp was initialised to zeros, out-of-range stays 0
    % "trim" mode: in_range should be all true after pass 1

    if ~all(in_range)
        n_truncated = n_truncated + 1;
    end

    % Normalise by max(abs)
    amp = max(abs(wf_interp));
    if amp > 0
        wf_interp = wf_interp / amp;
    end

    aligned_wf(:, u) = wf_interp;
end

if n_truncated > 0
    fprintf('Waveform edge handling (%s): %d/%d units had incomplete coverage\n', ...
        edge_mode, n_truncated, N_units);
end
fprintf('Waveform alignment: %d units, %d input samples @ %g Hz → %d output samples @ %g Hz\n', ...
    N_units, size(wf_raw, 1), sr_wf, n_out, target_sr);

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

% ── Return actual waveform window (for cross-dataset harmonization) ─────────
wf_window.pre_trough_ms  = n_out_pre  * dt_out;
wf_window.post_trough_ms = n_out_post * dt_out;
wf_window.n_samples      = n_out;
wf_window.target_sr      = target_sr;
end
