function [cls, info] = acgSubtypeLegacy(acg, bin_size, opts)
% ACGSUBTYPELEGACY  Rule-based ACG subtype of the AdvSci 2026 paper (Fig. 6), kept unchanged
% for reproducibility.
%
%   [cls, info] = acgSubtypeLegacy(acg, bin_size)
%
% Reproduces the classification loop of EI_Paper_Figures/scripts/acg_osc_analysis_main.m. The
% exponential envelope fit of the original is omitted because it does not enter the class.
%
% INPUTS:
%   acg       - (N_units x N_bins) two-sided ACGs including the zero-lag bin (odd N_bins);
%               the paper used 2 ms bins over +-2 s, normalised to the maximum
%   bin_size  - bin width in seconds (paper: 0.002)
%   opts      - name-value overrides of the published thresholds (see arguments block)
%
% OUTPUTS:
%   cls  - (N_units x 1) string: "noisy" | "noisy_weak" | "oscillatory" | "bursty"
%          ("harmonic" units are merged into "oscillatory", as in the paper)
%   info - (N_units x 1) struct array: RI, PeakFreqFFT, PeakFreqACG, CV, IsHarmonic, IsBursty
%
% Known limitations (see spikeRhythmicity for the statistical alternative):
%   - RI is the peak of the squared ACG spectrum over its mean in 2-60 Hz, without an
%     aperiodic (1/f) model; slowly decaying ACGs peak at the 2 Hz band edge.
%   - No null model: thresholds are not tied to a false-positive rate.
%   - ACGs reflect network rhythms as well as the unit's own firing.
arguments
    acg      double
    bin_size (1,1) double = 0.002
    opts.MinFreq            (1,1) double = 2
    opts.MaxFreq            (1,1) double = 60
    opts.MinPeaks           (1,1) double = 3
    opts.MaxCV              (1,1) double = 0.4
    opts.RIOscThresh        (1,1) double = 50
    opts.RIWeakThresh       (1,1) double = 20
    opts.BurstLagThresh     (1,1) double = 100     % ms
    opts.HarmonicTol        (1,1) double = 0.15
    opts.HarmonicRatioThresh(1,1) double = 0.25
    opts.BurstProminenceMult(1,1) double = 3
    opts.BurstTailPct       (1,1) double = 0.20
end
if isvector(acg), acg = acg(:)'; end
n_units = size(acg, 1);
n_bins  = size(acg, 2);
lags     = ((1:n_bins) - ceil(n_bins / 2)) * bin_size;
zi       = lags == 0;
lags_nz  = lags(~zi);
pos      = lags_nz > 0;
lags_pos = lags_nz(pos) * 1000;                                  % ms

cls  = strings(n_units, 1);
info = repmat(struct('RI', NaN, 'PeakFreqFFT', NaN, 'PeakFreqACG', NaN, 'CV', NaN, ...
    'IsHarmonic', false, 'IsBursty', false), n_units, 1);
for i = 1:n_units
    acg_full = acg(i, ~zi);
    acg_pos  = acg_full(pos);
    acg_smooth = smoothdata(acg_pos - mean(acg_pos), 'gaussian', 5);
    [~, locs] = findpeaks(acg_smooth, 'MinPeakProminence', std(acg_smooth), ...
        'MinPeakDistance', round(0.01 / bin_size));

    % peak spacing of the ACG
    cv = NaN; f_acg = NaN;
    if numel(locs) >= opts.MinPeaks
        intervals = diff(lags_pos(locs)) / 1000;
        f_acg = 1 / median(intervals);
        cv    = std(intervals) / mean(intervals);
    end
    valid_acg = ~isnan(f_acg) && f_acg >= opts.MinFreq && f_acg <= opts.MaxFreq && cv <= opts.MaxCV;

    % harmonics among all pairwise peak intervals
    is_harmonic = false;
    if numel(locs) >= 3 && ~isnan(f_acg)
        pt = lags_pos(locs);
        [a, b] = find(triu(true(numel(pt)), 1));
        all_f = 1 ./ ((pt(b) - pt(a)) / 1000);
        all_f = all_f(all_f >= opts.MinFreq & all_f <= opts.MaxFreq);
        if numel(all_f) >= 4
            for ratio = [2, 3]
                f_h = f_acg * ratio;
                if sum(abs(all_f - f_h) <= opts.HarmonicTol * f_h) / numel(all_f) >= opts.HarmonicRatioThresh ...
                        && f_h <= opts.MaxFreq
                    is_harmonic = true;
                    break
                end
            end
        end
    end

    % rhythmicity index: peak of the squared ACG spectrum over the mean of the other bins
    x  = acg_full - mean(acg_full);
    n  = numel(x);
    pw = abs(fft(x)).^2;  pw = pw(1:floor(n / 2) + 1);
    fq = (0:floor(n / 2)) / (n * bin_size);
    band = fq >= opts.MinFreq & fq <= opts.MaxFreq;
    pb = pw(band); fb = fq(band);
    [pk, ip] = max(pb);
    other = pb; other(ip) = [];
    ri = pk / mean(other);

    % burst shoulder: first ACG peak before BurstLagThresh, prominent against the ACG tail
    is_bursty = false;
    if ~isempty(locs)
        first_lag = lags_pos(locs(1));
        n_tail = max(1, round(opts.BurstTailPct * numel(acg_pos)));
        tail   = acg_pos(end - n_tail + 1:end);
        if first_lag < opts.BurstLagThresh
            tail_mean  = mean(tail);
            prominence = max(acg_pos(lags_pos <= opts.BurstLagThresh)) - tail_mean;
            is_bursty  = (tail_mean + prominence) >= tail_mean + opts.BurstProminenceMult * std(tail);   % original form
        end
    end

    if     is_bursty && (isnan(f_acg) || ri < opts.RIWeakThresh), c = "bursty";
    elseif ri < opts.RIWeakThresh,                                c = "noisy";
    elseif ri < opts.RIOscThresh && ~valid_acg,                   c = "noisy_weak";
    elseif valid_acg,                                             c = "oscillatory";   % incl. harmonic
    elseif ri >= opts.RIWeakThresh,                               c = "noisy_weak";
    else,                                                         c = "noisy";         % RI undefined (empty ACG)
    end
    cls(i) = c;
    info(i) = struct('RI', ri, 'PeakFreqFFT', fb(ip), 'PeakFreqACG', f_acg, 'CV', cv, ...
        'IsHarmonic', is_harmonic, 'IsBursty', is_bursty);
end
end
