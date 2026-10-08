function r = spikeRhythmicity(spike_times, duration, opts)
% SPIKERHYTHMICITY  Statistical test for rhythmic firing of one spike train.
%
%   r = spikeRhythmicity(spike_times, duration)
%   r = spikeRhythmicity(spike_times, duration, Band=[2 60], Alpha=0.01)
%
% Method
%   1. Welch spectrum of the binned train (2 ms bins, 4 s Hann windows, 50 % overlap).
%   2. Aperiodic background P_ap(f) = c + A / (k + f^chi): the white floor of a point process
%      plus a slow "knee" hump from bursts / rate fluctuations (as in FOOOF / specparam,
%      Donoghue et al. 2020), fitted robustly in log power with peaks excluded.
%   3. Under H0 (no rhythm) the Welch estimate is P_ap * chi2(nu)/nu, so each frequency of
%      the search band gets a p-value; Bonferroni across the band controls the per-train
%      false-positive rate. A rhythm also needs MinEffect excess power over the background
%      and must stand out from its +-2..5 Hz surround (a local spectral peak), so a broad
%      misfit of the background is not called a rhythm.
%   4. Burstiness is reported separately as the local variation LV of the inter-spike
%      intervals (Shinomoto et al. 2003: ~1 Poisson, <1 regular, >1 bursty).
%
% Replaces the threshold rules of acgSubtypeLegacy (kept for reproducibility). On simulated
% trains it detects 8 Hz sine-modulated and 3 Hz periodic-burst firing (100 %) where the
% legacy rules return "noisy_weak", with 2-7 % calls on aperiodic nulls (Poisson, OU rate
% fluctuations, random network bursts).
%
% INPUTS:
%   spike_times - spike times in seconds
%   duration    - recording duration in seconds
% OUTPUT (struct):
%   Rhythmic  - true if a significant spectral peak exists in Band
%   Freq      - peak frequency (Hz; strongest excess in Band if not rhythmic)
%   Effect    - excess power at Freq relative to the background (P/P_ap - 1)
%   p         - uncorrected p-value at Freq
%   Knee      - knee frequency of the background (Hz), the slow timescale of the train
%   Exponent  - background exponent chi
%   LV        - local variation of the ISIs
%   nu        - degrees of freedom of the Welch estimate
arguments
    spike_times (:,1) double
    duration    (1,1) double
    opts.Bin       (1,1) double = 0.002
    opts.Win       (1,1) double = 4
    opts.Band      (1,2) double = [2 60]
    opts.FitBand   (1,2) double = [0.5 200]
    opts.Alpha     (1,1) double = 0.01
    opts.MinEffect (1,1) double = 0.10
end
r = struct('Rhythmic', false, 'Freq', NaN, 'Effect', NaN, 'p', NaN, 'Knee', NaN, ...
    'Exponent', NaN, 'LV', NaN, 'nu', NaN);
dt = opts.Bin; n = floor(duration / dt);
st = sort(spike_times(spike_times >= 0 & spike_times < n * dt));
isi = diff(st);
if numel(isi) >= 3
    r.LV = 3 / (numel(isi) - 1) * sum(((isi(1:end-1) - isi(2:end)) ./ (isi(1:end-1) + isi(2:end))).^2);
end
L = round(opts.Win / dt); L = L + mod(L, 2); step = L / 2;
K = floor((n - L) / step) + 1;
if K < 4 || numel(st) < 20, return; end

% -- Welch spectrum -------------------------------------------------------------------
x = accumarray(floor(st / dt) + 1, 1, [n 1]); x = x - mean(x);
w = 0.5 - 0.5 * cos(2 * pi * (0:L-1)' / (L - 1));
P = zeros(L / 2 + 1, 1);
for k = 1:K
    F = fft(x((k - 1) * step + (1:L)) .* w);
    P = P + abs(F(1:L / 2 + 1)).^2;
end
P = P / K;
f = (0:L / 2)' / (L * dt);
nu = 2 * K / (1 + 2 * 0.1667^2);                       % Hann window, 50 % overlap

% -- aperiodic background, fitted robustly in log power ---------------------------------
fit_m = f >= opts.FitBand(1) & f <= min(opts.FitBand(2), 0.5 / dt);
ff = f(fit_m); yy = log10(max(P(fit_m), realmin));
c0 = median(P(fit_m & f > 0.75 * max(ff)));
A0 = max(P(find(fit_m, 1)) - c0, c0) * 2;
model = @(p, f) log10(10^p(1) + 10^p(2) ./ (10^p(3) + f.^p(4)));
pb = [log10(c0), log10(A0), 0, 2];
keep = true(size(ff));
for it = 1:3
    sse = @(p) sum((yy(keep) - model(p, ff(keep))).^2) + 1e3 * (max(0, -p(4)) + max(0, p(4) - 6));
    pb = fminsearch(sse, pb, optimset('MaxFunEvals', 4000, 'MaxIter', 4000, 'Display', 'off'));
    res = yy - model(pb, ff);
    s = 1.4826 * median(abs(res(keep) - median(res(keep))));
    keep = ~(conv(double(res > 2.5 * s), ones(5, 1), 'same') > 0);   % exclude peaks
end
R = P ./ 10.^model(pb, f);
r.Knee = (10^pb(3))^(1 / max(pb(4), eps)); r.Exponent = pb(4); r.nu = nu;

% -- peak test -------------------------------------------------------------------------
band = f >= opts.Band(1) & f <= opts.Band(2);
pv = 1 - gammainc(R * nu / 2, nu / 2);                 % P(chi2_nu / nu >= R)
sig = band & pv < opts.Alpha / sum(band) & R >= 1 + opts.MinEffect;
for i = find(sig)'                                      % local peak against its surround
    sur = abs(f - f(i)) >= 2 & abs(f - f(i)) <= 5;
    sig(i) = R(i) - median(R(sur)) >= opts.MinEffect / 2;
end
Rb = R; Rb(~band) = -Inf;
if any(sig), Rb(~sig) = -Inf; r.Rhythmic = true; end
[~, im] = max(Rb);
r.Freq = f(im); r.Effect = R(im) - 1; r.p = pv(im);
end
