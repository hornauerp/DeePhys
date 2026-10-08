function R = burstRhythm(spike_times, spike_units, groups, bursts, opts)
% BURSTRHYTHM  Fast structure inside network bursts: burst-locked sub-peaks and shared
% fluctuations that are not locked to the burst, each against its own null.
%
%   R = burstRhythm(spike_times, spike_units, groups, bursts)
%   R = burstRhythm(..., Pairs={"E", "I"}, NSplits=10)
%
% 1. Locked (what a burst-averaged trace shows, e.g. the E peaks of a "triphasic" burst). Per
%    group, the population rate (BinSize, Hz per unit) is averaged over bursts aligned at T_peak;
%    its fast part is the average minus a Gaussian smooth (SmoothSigma). Split-half
%    reliability: correlation of the fast parts of the odd- and the even-numbered bursts inside
%    LockedWindow (Spearman-Brown corrected in Reliability); P = share of circular shifts of one
%    half (20 ms .. half the window) with at least that correlation. Freq = 1 / lag of the
%    highest autocorrelation peak of the fast part within FreqLags. A window shorter than
%    2 x FreqLags(1) cannot hold two peaks: a core-only window returns the burst width as the
%    frequency.
% 2. Shared, not locked. For every burst lasting >= SharedWindow, the counts of two disjoint unit
%    sets in [Onset, Onset + SharedWindow] (BinSize), minus their own Gaussian smooth, Hann
%    taper: cross-spectrum of the real pairing minus the mean over NShift pairings of set A of
%    burst k with set B of burst k + j (shift predictor: keeps whatever is locked to burst
%    onset, removes co-fluctuations that are not), normalised to coherency. Single-unit
%    structure (refractoriness, intrinsic bursts) does not enter. Pairs: two random halves of a
%    group (NSplits draws, averaged) or two groups (e.g. E and I; the phase of the coherency at
%    the peak gives the lag of the second relative to the first, Lag > 0 = second lags). Peak =
%    largest local maximum in Band, as height above its 3-8 Hz surround; P = share of shifted
%    pairings, each tested against the others, with at least that height.
%
% INPUTS:
%   spike_times, spike_units - (S x 1)
%   groups - struct of (1 x U) logical unit masks, e.g. struct("E", lab == 1, "I", lab == 2)
%   bursts - struct with Onset, Offset, T_peak (detectNetworkBursts)
% OPTIONS:
%   BinSize (0.002 s), SmoothSigma (0.025 s), LockedWindow ([-0.1 0.3] s re T_peak),
%   FreqLags ([0.025 0.2] s), MinBursts (10), MinUnits (3), SharedWindow (0.25 s),
%   Pairs ({} = random halves of every group; or a cell of {groupA, groupB} / {group} rows),
%   NSplits (5), NShift (20), Band ([5 60] Hz), NFFT (1000 points), Seed (1)
% OUTPUT (struct):
%   Locked - table per group: Group, NUnits, NBursts, SplitR, Reliability, P, Freq
%   Shared - table per pair: Pair, NBursts, PeakFreq, PeakExcess, P, Lag, Coh5_15, Coh15_30,
%            Coh30_60 (mean real coherency in the band, shift-corrected)
%   Traces - struct per group: Time, Mean, Fast (burst-averaged rate and its fast part)
%   Spectra - struct: Freq (Hz) and per pair the shift-corrected coherency
%
% See also burstEIMetrics, detectNetworkBursts, spikeRhythmicity.

arguments
    spike_times (:,1) double
    spike_units (:,1) double
    groups      struct
    bursts      struct
    opts.BinSize      (1,1) double = 0.002
    opts.SmoothSigma  (1,1) double = 0.025
    opts.LockedWindow (1,2) double = [-0.1 0.3]
    opts.FreqLags     (1,2) double = [0.025 0.2]
    opts.MinBursts    (1,1) double = 10
    opts.MinUnits     (1,1) double = 3
    opts.SharedWindow (1,1) double = 0.25
    opts.Pairs        cell = {}
    opts.NSplits      (1,1) double = 5
    opts.NShift       (1,1) double = 20
    opts.Band         (1,2) double = [5 60]
    opts.NFFT         (1,1) double = 1000
    opts.Seed         (1,1) double = 1
end
dt = opts.BinSize;
sg = opts.SmoothSigma / dt; g = exp(-0.5 * ((-ceil(4 * sg):ceil(4 * sg)) / sg).^2); g = g(:) / sum(g);
gnames = string(fieldnames(groups))';
[ss, o] = sort(spike_times); uu = spike_units(o);
rs = RandStream('twister', 'Seed', opts.Seed);

% -- 1. locked -----------------------------------------------------------------------------
pk = bursts.T_peak(:)'; pk = pk(isfinite(pk));
W = round(max(abs(opts.LockedWindow)) / dt) + round(4 * sg);
tw = (-W:W) * dt; win = tw >= opts.LockedWindow(1) & tw <= opts.LockedWindow(2);
rows = cell(0, 7); R.Traces = struct();
for gname = gnames
    sel = groups.(gname)(:)';
    if sum(sel) < opts.MinUnits || numel(pk) < opts.MinBursts
        rows(end+1, :) = {gname, sum(sel), numel(pk), NaN, NaN, NaN, NaN}; %#ok<AGROW>
        continue
    end
    m = sel(uu); C = countsAround(ss(m), pk, -W * dt, (2 * W + 1), dt) / sum(sel) / dt;   % bursts x bins
    fast = @(x) x(:) - conv(x(:), g, 'same') ./ conv(ones(numel(x), 1), g, 'same');
    A1 = fast(mean(C(1:2:end, :), 1)); A2 = fast(mean(C(2:2:end, :), 1)); Af = fast(mean(C, 1));
    a1 = A1(win); a2 = A2(win); r = corr(a1, a2);
    nw = numel(a1); sh = unique(round(linspace(round(0.02 / dt), max(round(nw / 2), round(0.02 / dt)), 100)));
    rn = arrayfun(@(k) corr(a1, circshift(a2, k)), sh);
    x = Af(win) - mean(Af(win)); ac = xcorr(x, 'coeff'); ac = ac(numel(x) + 1:end);
    lag = (1:numel(ac))' * dt; okl = lag >= opts.FreqLags(1) & lag <= opts.FreqLags(2);
    fq = NaN;
    if sum(okl) >= 3
        [pv, pl] = findpeaks(ac(okl));
        if ~isempty(pl), [~, j] = max(pv); lg = lag(okl); fq = 1 / lg(pl(j)); end
    end
    rows(end+1, :) = {gname, sum(sel), numel(pk), r, 2 * r / (1 + r), mean(rn >= r), fq}; %#ok<AGROW>
    R.Traces.(gname) = struct('Time', tw, 'Mean', mean(C, 1), 'Fast', Af');
end
R.Locked = cell2table(rows, 'VariableNames', {'Group', 'NUnits', 'NBursts', 'SplitR', 'Reliability', 'P', 'Freq'});

% -- 2. shared, not locked --------------------------------------------------------------------
n = round(opts.SharedWindow / dt); nf = opts.NFFT;
f = (0:floor(nf / 2))' / (nf * dt); band = f >= opts.Band(1) & f <= opts.Band(2);
on = bursts.Onset(:)'; off = bursts.Offset(:)';
on = on(isfinite(on) & isfinite(off) & off - on >= opts.SharedWindow);
nb = numel(on); NS = min(opts.NShift, nb - 1);
ec = conv(ones(n, 1), g, 'same'); hw = 0.5 - 0.5 * cos(2 * pi * (0:n - 1)' / (n - 1));
spec = @(X) fft((X - conv2(X, g, 'same') ./ ec) .* hw, nf);
pairs = opts.Pairs;
if isempty(pairs), pairs = num2cell(gnames(:)); end
rows = cell(0, 9); R.Spectra = struct('Freq', f);
for p = 1:numel(pairs)
    pr = string(pairs{p}); name = strjoin(pr, "-"); if isscalar(pr), name = pr + " halves"; end
    if nb < opts.MinBursts || NS < 2 || any(arrayfun(@(q) sum(groups.(q)) < opts.MinUnits * (1 + isscalar(pr)), pr))
        rows(end+1, :) = {name, nb, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        continue
    end
    if isscalar(pr)
        idx = find(groups.(pr)); splits = cell(opts.NSplits, 2);
        for s = 1:opts.NSplits
            h = false(size(groups.(pr))); h(idx(randperm(rs, numel(idx), floor(numel(idx) / 2)))) = true;
            splits(s, :) = {h, groups.(pr) & ~h};
        end
    else
        splits = {groups.(pr(1)), groups.(pr(2))};
    end
    Creal = zeros(numel(f), 1); Csh = zeros(numel(f), NS);
    for s = 1:size(splits, 1)
        FA = spec(countsAround(ss(splits{s, 1}(uu)), on, 0, n, dt)'); FB = spec(countsAround(ss(splits{s, 2}(uu)), on, 0, n, dt)');
        FA = FA(1:numel(f), :); FB = FB(1:numel(f), :);
        nrm = sqrt(mean(abs(FA).^2, 2) .* mean(abs(FB).^2, 2)); nrm(nrm == 0) = Inf;
        Creal = Creal + mean(FA .* conj(FB), 2) ./ nrm / size(splits, 1);
        for j = 1:NS
            Csh(:, j) = Csh(:, j) + mean(FA .* conj(FB(:, circshift(1:nb, j))), 2) ./ nrm / size(splits, 1);
        end
    end
    D = Creal - mean(Csh, 2);
    [pf, pe] = peakStat(real(D), f, band);
    pn = arrayfun(@(j) nthOut(@() peakStat(real(Csh(:, j) - mean(Csh(:, setdiff(1:NS, j)), 2)), f, band), 2), 1:NS);
    lagv = NaN; if isfinite(pf), lagv = angle(D(f == pf)) / (2 * pi * pf); end
    bandMean = @(a, b) mean(real(D(f >= a & f <= b)));
    rows(end+1, :) = {name, nb, pf, pe, mean(pn >= pe), lagv, bandMean(5, 15), bandMean(15, 30), bandMean(30, 60)}; %#ok<AGROW>
    R.Spectra.(matlab.lang.makeValidName(name)) = D;
end
R.Shared = cell2table(rows, 'VariableNames', {'Pair', 'NBursts', 'PeakFreq', 'PeakExcess', 'P', 'Lag', 'Coh5_15', 'Coh15_30', 'Coh30_60'});
end

function C = countsAround(t, t0, start, nbins, dt)
% (numel(t0) x nbins) spike counts of the sorted times t in nbins bins of dt from t0 + start
C = zeros(numel(t0), nbins);
for k = 1:numel(t0)
    a = t0(k) + start; i0 = find(t >= a, 1); if isempty(i0), continue; end
    i1 = find(t < a + nbins * dt, 1, 'last'); if isempty(i1) || i1 < i0, continue; end
    C(k, :) = accumarray(min(floor((t(i0:i1) - a) / dt) + 1, nbins), 1, [nbins 1])';
end
end

function [pf, pe] = peakStat(r, f, band)
% strongest local maximum of r (3-point smoothed) in band, as height above its 3-8 Hz surround
r = movmean(r, 3); pf = NaN; best = -inf;
for j = find(band)'
    if j == 1 || j == numel(r) || r(j) < r(j - 1) || r(j) < r(j + 1), continue; end
    sur = abs(f - f(j)) >= 3 & abs(f - f(j)) <= 8;
    e = r(j) - mean(r(sur));
    if e > best, best = e; pf = f(j); end
end
pe = max(best, 0);
end

function v = nthOut(fn, k)
% k-th output of fn()
out = cell(1, k); [out{:}] = fn(); v = out{k};
end
