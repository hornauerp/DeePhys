function [bursts, info] = detectNetworkBursts(spike_times, spike_units, duration, opts)
% DETECTNETWORKBURSTS  Network bursts by ISI_N, with N and the threshold chosen jointly against
% surrogate data.
%
%   [bursts, info] = detectNetworkBursts(spike_times, spike_units, duration)
%   [...] = detectNetworkBursts(..., TonicMethod="none", NSurrogates=50)
%   [...] = detectNetworkBursts(..., N=20, Threshold=0.1)        % plain ISI_N, fixed parameters
%
% 1. Uncoupled units are removed (TonicMethod): units whose firing does not follow the rest of
%    the network add spikes between bursts and blur the ISI_N distribution.
%      "coupling" - population coupling (Okun et al. 2015): correlation of the unit's binned
%                   rate with the summed rate of all other units, z-scored against circular
%                   shifts of the unit's own train; removed if z < MinZ
%      "gini"     - Gini coefficient of the unit's binned counts against Poisson firing at the
%                   unit's rate (z vs Poisson trains with the same count); bins hold
%                   GiniSpikesPerBin spikes on average (GiniBin "adaptive") or are fixed;
%                   removed if z < MinZ (not more clustered than Poisson)
%      "regular"  - same Gini, removed only if z < -MinZ (significantly more regular than
%                   Poisson, i.e. tonic); Poisson-like units are kept
%      "coupling+regular" - remove if uncoupled or more regular than Poisson
%      "either" / "both" - remove if either / both of "coupling" and "gini" flag the unit
%      "legacy"   - removeTonic (high-rate outliers with low CV(ISI); needs FiringRates, CVISI)
%      "none"
%    Both metrics are always computed and returned in info.Units.
% 2. For every N (log grid) the log10 ISI_N of the merged train of the kept units is split at
%    every threshold t. The split quality G(N, t) is the per-window log-likelihood gain of two
%    Gaussians, one per side, over one Gaussian (minimum-error thresholding, Kittler &
%    Illingworth 1986): unlike Otsu's criterion it is not dominated by the background when burst
%    windows are rare. Surrogates in which every kept unit's train is circularly shifted by its
%    own random offset (rate, ISI distribution and single-unit bursts kept, co-activation
%    destroyed) give the significance z(N, t) = (G - mean over surrogates of max_t G) / SD.
% 3. Levels. Distinct maxima of z (prominence - height above the highest saddle to a higher
%    maximum - >= max(MinZ, LevelProminence x height)) are levels of grouping, e.g. network
%    bursts and the slow activity epochs or superbursts that group them. No maximum is "the"
%    burst scale by itself: G grows towards coarse groupings, z is highest at the fine end of a
%    burst ridge, and slow epochs are often the highest maximum. Along the ridge of one level the
%    same events appear at many (N, t): small N / t include the smaller events too, large N / t
%    only the large ones. Candidates: the best t at every N within each level, if significant
%    (z >= MinZ and G - surrogate >= MinEffect). Each gives bursts: runs of consecutive spikes
%    linked by a window of N spikes spanning <= t (ISI_N rule of Bakkum et al. 2014) that
%    recruit >= MinBurstUnits distinct units and significantly more than chance: MinZ SD above
%    the units expected to fire at least once anyway within the burst's core (middle 80 % of its
%    spikes) at their rates outside the bursts (Poisson: in a long window of dense activity many
%    units fire once). A candidate needs >= 5 such bursts.
% 4. Burst level and operating point. Network bursts are continuous above-background activity
%    separated by silence of at least MinQuietGap (0.2 s). Contiguity of a candidate = (share of
%    bursts without such a silence inside) x (share of gaps with one) x (share of the strong
%    sustained active episodes it covers, time-weighted). Chance share of a candidate: bursts
%    the background alone gives at the same (N, t) - the spikes outside its bursts, on a time
%    axis with the bursts cut out, each unit circularly shifted - over its burst count
%    (coincidences of a few units at small N / t).
%    Each level is represented by the candidate of smallest chance share among those within 0.05
%    of its highest contiguity (the same bursts at different (N, t)); a level more than 25 %
%    chance is dropped. The burst level is the finest level whose representative has no silence
%    inside most bursts, silence in most gaps and covers most strong activity (else the most
%    contiguous level); finer structure with active gaps is sub-events, coarser structure groups
%    bursts. Within the burst level the operating point is the finest candidate that finds
%    >= 90 % of the representative's bursts, splits <= 10 % of them and has a chance share
%    <= MaxChance (or, in weakly bursting cultures where no scale gets there, <= the
%    representative's): the same bursts plus all smaller events with >= MinBurstUnits units - one
%    rule in every culture. Also reported in
%    info.Scales: sub-events (the most significant finer candidate that splits most bursts) and
%    other levels that group the bursts, each only if the gaps it bridges between the finer
%    events are the short ones of a bimodal gap distribution (minimum-error gain of that split vs
%    exponential gaps of the same mean, GapZ >= MinZ). MainScale "score": the burst level is the
%    one of the highest score.
%    Superbursts: bursts inside a superburst are separated by silence relative to their own
%    peaks, not to the background of the whole recording, so ISI_N merges them. A main burst
%    whose population rate (kept units, 5 ms bins, 25 ms smoothing) stays below baseline +
%    SplitLevel x (its peak - baseline) for >= MinQuietGap inside its span is cut in the middle
%    of each such stretch; pieces that recruit as above are bursts, the rest (e.g. the
%    afterdischarge of a few units) is dropped. Bursts cut from one event share its index in
%    bursts.Superburst; info.Superbursts holds the superburst spans.
%    Gate: no valid level = no network bursts (info.Gate = false, no bursts returned).
%    Reported spans (T_start / T_end, info.ScaleBursts): first and last spike of each ISI_N run (of
%    each piece for bursts cut from a superburst) -
%    which spikes belong to a burst. They are not timing measures: the first spike can be a
%    background spike up to t before the burst, and t differs between cultures. Timing and shape
%    (Onset, Offset, Duration, Rise, Fall, FWHM, T_peak) come from burstShape at ShapeLevel
%    (default 0.25 of each burst's peak above baseline, population rate of all given spikes), so
%    every burst is measured at the same point of its own profile. For the scale selection
%    only, episodes of sustained above-background activity (kept units, 5 ms bins, 25 ms
%    smoothing, background + 2 SD) are used.
%
% INPUTS:
%   spike_times - (S x 1) spike times (s); spike_units - (S x 1) unit of each spike
%   duration    - recording duration (s)
% OPTIONS:
%   TonicMethod ("gini") with TonicMinRate ("median": units below the recording's median rate are
%   never removed - low-rate units carry burst spikes and look Poisson-like whatever their
%   timing). Without the floor, "gini" can return no bursts in weakly bursting recordings; with
%   it the detector agrees with the legacy detector. Known limit: with quasi-periodic bursting,
%   units firing ~1 spike per burst can look regular at 3 spikes per bin and are removed.
%   TonicMinRate also takes a rate in Hz or 0 (no floor).
%   FiringRates / CVISI (per unit 1..max(spike_units), "legacy" only), NSurrogates (20), NValues
%   ([] = 30 log-spaced values from MinBurstUnits to 20 x kept units), LogBin (0.02 decades),
%   MinZ (3), MinEffect (0.05 nats: absolute floor on G - surrogate; sharp bursts that are rare
%   in a dense background stay below it - 0 makes them candidates but can add small events of a
%   few units), MinClassFraction (0.001 and >= 50 windows on each side of a split: rare bursts
%   in dense background), CouplingBin (0.05 s), GiniBin
%   ("adaptive": GiniSpikesPerBin (3) spikes per bin on average, per unit; or seconds),
%   MinUnitSpikes (20: units with fewer spikes are kept and not tested), MinBurstUnits ([] =
%   max(5, 10 % of the active units, >= MinUnitSpikes spikes)), MainScale ("contiguity" |
%   "score"), MinQuietGap (0.2 s), LevelProminence (0.25), MaxChance (0.02), SplitLevel (0.5:
%   dips below half the peak for >= MinQuietGap separate the bursts of a superburst), ShapeLevel
%   (0.25: onset / offset at this fraction of each burst's peak above baseline; see burstShape),
%   N / Threshold (fixed: plain ISI_N, no search), Seed (1)
% OUTPUTS:
%   bursts - struct as computeBursts: T_start, T_end (1 x B, s: first / last spike), best_N,
%            best_ISI_N, Method; timing and shape from burstShape: T_peak, Onset, Offset,
%            Duration (Offset - Onset), Rise, Fall, FWHM, PeakRate (spikes/s), ShapeLevel; and the
%            size of each burst between Onset and Offset: NUnits (distinct kept units), NSpikes
%            (spikes of kept units) - e.g. to keep only large bursts; Superburst (index into
%            info.Superbursts, 0 = not part of a superburst)
%   info   - struct: Gate, Scales (table, main bursts first, then sub-events and grouping levels:
%            N, Threshold, Gain, Score, Z, NBursts, MedianUnits, RecruitedShare, Contiguity, Level,
%            Inside, Contains, GapZ), ScaleBursts (cell of [T_start; T_end]), Levels (table: Level,
%            PeakN, PeakThreshold, PeakZ, Prominence, RefN, RefThreshold, RefNBursts, ChanceShare, Contiguity,
%            NoQuietInside, QuietGaps, Coverage, Main), RejectedScales (with Level and Reason
%            "units" / "too few bursts" / "chance" / "no grouping"), MinBurstUnits, TonicMinRate,
%            Superbursts (2 x K: first / last spike), Units (tonic metrics, Removed),
%            Landscape (N, LogT, Gain, GainSurrogateMean, NullMaxMean, NullMaxSD)
arguments
    spike_times (:,1) double
    spike_units (:,1) double
    duration (1,1) double {mustBePositive}
    opts.TonicMethod (1,1) string {mustBeMember(opts.TonicMethod, ...
        ["coupling", "gini", "regular", "coupling+regular", "either", "both", "legacy", "none"])} = "gini"
    opts.FiringRates double = []
    opts.CVISI double = []
    opts.NSurrogates (1,1) double = 20
    opts.NValues double = []
    opts.LogBin (1,1) double = 0.02
    opts.MinZ (1,1) double = 3
    opts.MinEffect (1,1) double = 0.05
    opts.MinClassFraction (1,1) double = 0.001
    opts.CouplingBin (1,1) double = 0.05
    opts.GiniBin = "adaptive"
    opts.GiniSpikesPerBin (1,1) double = 3
    opts.MinUnitSpikes (1,1) double = 20
    opts.TonicMinRate = "median"
    opts.MinBurstUnits double = []
    opts.MainScale (1,1) string {mustBeMember(opts.MainScale, ["contiguity", "score"])} = "contiguity"
    opts.MinQuietGap (1,1) double = 0.2
    opts.LevelProminence (1,1) double = 0.25
    opts.MaxChance (1,1) double = 0.02
    opts.ShapeLevel (1,1) double {mustBeInRange(opts.ShapeLevel, 0, 1, "exclusive")} = 0.25
    opts.SplitLevel (1,1) double {mustBeInRange(opts.SplitLevel, 0, 1, "exclusive")} = 0.5
    opts.N double = []
    opts.Threshold double = []
    opts.Seed (1,1) double = 1
end
rs = RandStream('twister', 'Seed', opts.Seed);
[spike_times, o] = sort(spike_times); spike_units = spike_units(o);
bursts = struct('T_start', zeros(1, 0), 'T_end', zeros(1, 0), 'T_peak', zeros(1, 0), 'Onset', zeros(1, 0), ...
    'Offset', zeros(1, 0), 'Duration', zeros(1, 0), 'Rise', zeros(1, 0), 'Fall', zeros(1, 0), 'FWHM', zeros(1, 0), ...
    'PeakRate', zeros(1, 0), 'ShapeLevel', opts.ShapeLevel, 'NUnits', zeros(1, 0), 'NSpikes', zeros(1, 0), 'Superburst', zeros(1, 0), 'best_N', NaN, 'best_ISI_N', NaN, ...
    'Method', "ISINJoint");
info = struct('Gate', false, 'Scales', table(), 'ScaleBursts', {{}}, 'ScalePeaks', {{}}, 'RejectedScales', table(), ...
    'MinBurstUnits', NaN, 'Units', table(), 'Landscape', struct(), 'Levels', table(), 'Superbursts', zeros(2, 0));

% -- fixed parameters: plain ISI_N on all spikes ------------------------------------------------
if ~isempty(opts.N) || ~isempty(opts.Threshold)
    assert(~isempty(opts.N) && ~isempty(opts.Threshold), 'detectNetworkBursts:fixed', ...
        'Give both N and Threshold for fixed-parameter ISI_N.');
    [bursts.T_start, bursts.T_end] = isinBursts(spike_times, opts.N, opts.Threshold, NaN);   % first / last spike
    bursts.best_N = opts.N; bursts.best_ISI_N = opts.Threshold; bursts.Method = "ISIN";
    info.Gate = NaN;
    return
end

% -- 1. uncoupled / tonic units -------------------------------------------------------------------
U = unitMetrics(spike_times, spike_units, duration, opts, rs);
switch opts.TonicMethod
    case "coupling", rm = U.CouplingZ < opts.MinZ;
    case "gini",     rm = U.GiniZ < opts.MinZ;
    case "regular",  rm = U.GiniZ < -opts.MinZ;
    case "coupling+regular", rm = U.CouplingZ < opts.MinZ | U.GiniZ < -opts.MinZ;
    case "either",   rm = U.CouplingZ < opts.MinZ | U.GiniZ < opts.MinZ;
    case "both",     rm = U.CouplingZ < opts.MinZ & U.GiniZ < opts.MinZ;
    case "legacy"
        assert(~isempty(opts.FiringRates) && ~isempty(opts.CVISI), 'detectNetworkBursts:legacy', ...
            'TonicMethod "legacy" needs FiringRates and CVISI.');
        fr = opts.FiringRates(:); cv = opts.CVISI(:);
        tonic_ids = find((isoutlier(fr, "ThresholdFactor", 1) & fr > median(fr)) & cv < mean(cv));
        rm = ismember(U.Unit, tonic_ids);
    otherwise,       rm = false(height(U), 1);
end
% units below the rate floor are never removed: low-rate units are rarely tonic, and removing
% them loses more burst spikes than background spikes
floor_hz = opts.TonicMinRate;
if (isstring(floor_hz) || ischar(floor_hz)) && string(floor_hz) == "median"
    floor_hz = median(U.Rate(U.NSpikes >= opts.MinUnitSpikes));
end
info.TonicMinRate = floor_hz;
U.Removed = rm & U.NSpikes >= opts.MinUnitSpikes & U.Rate >= floor_hz;
info.Units = U;
keep_units = U.Unit(~U.Removed);
k = ismember(spike_units, keep_units);
t = spike_times(k); su = spike_units(k);
n_kept = numel(unique(su));
min_units = opts.MinBurstUnits;
% a network burst recruits >= 10 % of the active units of the recording (not only of the kept
% ones: when coupling removal leaves few units, a handful of co-firing units is no network burst)
if isempty(min_units), min_units = max(5, ceil(0.1 * sum(U.NSpikes >= opts.MinUnitSpikes))); end
info.MinBurstUnits = min_units;
if n_kept < min_units || numel(t) < 100
    return
end

% -- 2. G(N, t) of the data and of the surrogates -------------------------------------------------
Nv = opts.NValues;
if isempty(Nv)
    % from the minimum recruitment up to 20 spikes per kept unit (spans bursts of sub-events)
    n_lo = max(3, min_units);
    Nv = unique(round(logspace(log10(n_lo), log10(max(2 * n_lo, min(20 * n_kept, numel(t) / 20))), 30)));
end
edges = (floor(log10(1e-4) / opts.LogBin):ceil(log10(duration) / opts.LogBin)) * opts.LogBin;
G = gainMap(t, Nv, edges, opts.MinClassFraction);
S = opts.NSurrogates;
null_max = nan(S, numel(Nv)); g_sur = zeros(size(G));
[~, ~, iu] = unique(su);
for s = 1:S
    shift = duration * rand(rs, max(iu), 1);
    ts = sort(mod(t + shift(iu), duration));
    Gs = gainMap(ts, Nv, edges, opts.MinClassFraction);
    null_max(s, :) = max(Gs, [], 2, 'omitnan')';
    g_sur = g_sur + fillmissing(Gs, 'constant', 0) / S;
end
mu0 = mean(null_max, 1, 'omitnan'); sd0 = max(std(null_max, 0, 1, 'omitnan'), 1e-6);
centers = edges(1:end-1) + opts.LogBin / 2;
info.Landscape = struct('N', Nv, 'LogT', centers, 'Gain', G, 'GainSurrogateMean', g_sur, ...
    'NullMaxMean', mu0, 'NullMaxSD', sd0);

% -- 3. levels: distinct maxima of the significance map --------------------------------------------
% z(N, t) = (G - mean of the surrogates' best G at N) / SD. Its distinct maxima (prominence, i.e.
% height above the highest saddle to a higher maximum, >= max(MinZ, LevelProminence x height))
% are levels of grouping, e.g. bursts and the slow epochs or superbursts that group them. Along
% the ridge of one level the same events appear at many (N, t): small N / t include the smaller
% events too, large N / t only the large ones. Candidates: the best t at every N within each
% level (significant: z >= MinZ and G - surrogate >= MinEffect).
Z = (G - mu0') ./ sd0';
Zs = movmean(Z, 3, 2, 'omitnan'); Zs(isnan(G)) = NaN;
sig = Z >= opts.MinZ & G - mu0' >= opts.MinEffect;
[lab, lvl_peak, lvl_prom] = levelBasins(Zs, opts.MinZ, opts.LevelProminence);
nl = numel(lvl_peak);
cand = zeros(0, 3);                                         % [N index, t index, level]
for L = 1:nl
    for n = 1:numel(Nv)
        cols = find(lab(n, :) == L & sig(n, :));
        if isempty(cols), continue; end
        [~, m] = max(Zs(n, cols)); cand(end+1, :) = [n, cols(m), L]; %#ok<AGROW>
    end
end
info.Levels = table();
if isempty(cand)
    return
end
% bursts per candidate: recruitment, >= 5 bursts, and the contiguity with the population rate
[rr, cc, lvl] = popRate(t, duration, 0.005, 5);
minq = max(1, round(opts.MinQuietGap / 0.005));
nc = size(cand, 1);
thr = 10 .^ edges(cand(:, 2) + 1)';                         % upper edge of the split bin
cz = Z(sub2ind(size(Z), cand(:, 1), cand(:, 2))); cg = G(sub2ind(size(G), cand(:, 1), cand(:, 2)));
cs = cg - mu0(cand(:, 1))';
valid = false(nc, 1); Braw = cell(nc, 1); Bref = Braw; Pk = Braw; nbur = zeros(nc, 1); mu = nan(nc, 1);
mx = mu; q4 = nan(nc, 4); rej = {};
prev = prevSame(iu);                                        % previous spike of the same unit
for c = 1:nc
    N = Nv(cand(c, 1)); x = log10(max(t(N:end) - t(1:end-N+1), 1e-6)); lt = log10(thr(c));
    edge_isi = [];
    if any(x <= lt) && any(x > lt), edge_isi = 10 ^ ((mean(x(x <= lt)) + mean(x(x > lt))) / 2) / (N - 1); end
    [b0, b1, i0, i1] = isinBursts(t, N, thr(c), edge_isi);
    % recruitment: >= MinBurstUnits distinct units, and significantly more (MinZ SD) than expected
    % to fire anyway within the burst's core (middle 80 % of its spikes) at the units' rates
    % outside the bursts - in a long window of dense activity many units fire once by chance
    nu_b = runUnits(prev, i0, i1);
    ok_b = recruited(nu_b, unitRates(t, iu, i0, i1, duration), coreDuration(t, i0, i1), min_units, opts.MinZ);
    mu(c) = median(nu_b); mx(c) = mean(ok_b);
    if ~any(ok_b)
        rej(end+1, :) = {N, thr(c), cs(c), cz(c), numel(b0), mu(c), cand(c, 3), "units"}; %#ok<AGROW>
        continue
    end
    B = [b0(ok_b); b1(ok_b)];
    if size(B, 2) < 5                       % a burst scale needs >= 5 bursts (slow epochs: 1-4)
        rej(end+1, :) = {N, thr(c), cs(c), cz(c), size(B, 2), mu(c), cand(c, 3), "too few bursts"}; %#ok<AGROW>
        continue
    end
    valid(c) = true; Braw{c} = B; nbur(c) = size(B, 2);
    [Bref{c}, Pk{c}] = refineEdges(t, B, duration, 0.005, 5);
    [q4(c, 1), q4(c, 2), q4(c, 3), q4(c, 4)] = contiguity(rr, cc, lvl, Bref{c}, minq, 2 * 5);
end
if ~isempty(rej)
    info.RejectedScales = cell2table(rej, 'VariableNames', {'N', 'Threshold', 'Score', 'Z', 'NBursts', 'MedianUnits', 'Level', 'Reason'});
end
if ~any(valid)
    return
end

% -- 4. the burst level, its operating point, and the other levels --------------------------------
% Every level is represented by its candidate of highest contiguity (bursts = continuous
% above-background activity separated by >= MinQuietGap of silence; see contiguity). The burst
% level (MainScale "contiguity") is the finest level whose representative has no silence inside
% most events, silence in most gaps and covers most sustained activity: finer levels with active
% gaps are sub-events, coarser ones group bursts (superbursts, activity epochs). If none
% qualifies, the level of highest contiguity. MainScale "score": the level of the highest score.
% Chance: the background alone - the spikes outside a candidate's bursts, on a time axis with the
% bursts cut out, each unit circularly shifted - gives the bursts expected by chance at the same
% (N, t) (coincidences of a few units at small N / t; cutting out spurious events does not hide
% them, the shifted background forms them again at the same rate). Each level is represented by
% the candidate with the smallest chance share among those within 0.05 of its highest
% contiguity (the same bursts at different (N, t)); a level more than 25 % chance is dropped
% (synthetic weak bursts: 6-12 %).
% Operating point: the finest candidate of the burst level that finds >= 90 % of the
% representative's bursts, splits <= 10 % of them (two of its bursts in one) and whose chance
% share is <= MaxChance or, in weakly bursting cultures where no scale gets there, <= the
% representative's: the same bursts plus the smaller events that recruit >= MinBurstUnits units -
% one rule in all cultures.
chance = nan(nc, 1);
getChance = @(c) chanceBursts(backgroundTrains(t, iu, Bref{c}, 0.05, duration, 5, rs), Nv(cand(c, 1)), thr(c), min_units, opts.MinZ);
rep = zeros(nl, 1); tol = nan(nl, 1);
for L = 1:nl
    k = find(valid & cand(:, 3) == L);
    if isempty(k), continue; end
    k = k(q4(k, 1) >= max(q4(k, 1)) - 0.05);
    for c = k', chance(c) = getChance(c); end
    [~, o] = sortrows([chance(k) ./ nbur(k), -q4(k, 1)]); c = k(o(1));
    if chance(c) <= 0.25 * nbur(c), rep(L) = c; tol(L) = max(opts.MaxChance, chance(c) / nbur(c)); end
end
lv_ok = find(rep > 0);
if ~isempty(lv_ok)
    if opts.MainScale == "score"
        k = lv_ok; [~, m] = max(cs(rep(k))); main_L = k(m);
    else
        qual = lv_ok(q4(rep(lv_ok), 2) >= 0.5 & q4(rep(lv_ok), 3) >= 0.5 & q4(rep(lv_ok), 4) >= 0.5);
        if ~isempty(qual)
            [~, m] = min(thr(rep(qual))); main_L = qual(m);
        else
            [~, m] = max(q4(rep(lv_ok), 1)); main_L = lv_ok(m);
        end
    end
    R = Braw{rep(main_L)};
    hits = @(A, Bx) sum(A(1, :)' <= Bx(2, :) & A(2, :)' >= Bx(1, :), 2);   % per A burst: # Bx overlapping
    k = find(valid & cand(:, 3) == main_L);
    [~, o] = sortrows([thr(k), Nv(cand(k, 1))']); k = k(o);
    op = rep(main_L);
    for c = k'
        if thr(c) >= thr(op), break; end
        h = hits(R, Braw{c});
        if mean(h >= 1) < 0.9 || mean(h >= 2) > 0.1, continue; end
        if isnan(chance(c)), chance(c) = getChance(c); end
        if chance(c) <= tol(main_L) * nbur(c), op = c; break; end
    end
end
tested = find(~isnan(chance) & ~(chance <= tol(cand(:, 3)) .* nbur));
if ~isempty(tested)
    extra = table(Nv(cand(tested, 1))', thr(tested), cs(tested), cz(tested), nbur(tested), mu(tested), cand(tested, 3), ...
        repmat("chance", numel(tested), 1), 'VariableNames', {'N', 'Threshold', 'Score', 'Z', 'NBursts', 'MedianUnits', 'Level', 'Reason'});
    if isempty(info.RejectedScales), info.RejectedScales = extra; else, info.RejectedScales = [info.RejectedScales; extra]; end
end
if isempty(lv_ok)
    return
end
% sub-events: the most significant finer candidate (any level) that splits most main bursts
% (kept below only if the gaps it leaves inside the bursts are the short ones of a bimodal
% gap distribution)
sub = [];
kf = find(valid & thr < thr(op));
if ~isempty(kf)
    kf = kf(arrayfun(@(c) mean(hits(Braw{op}, Braw{c}) >= 2), kf) >= 0.5);
    if ~isempty(kf), [~, m] = max(cz(kf)); sub = kf(m); end
end
% the other levels relative to the main bursts (Inside: sub-events; Contains: groups of bursts);
% a level whose representative finds the main bursts or the sub-events again is not reported
% separately
others = rep(setdiff(lv_ok, main_L));
others = others(~arrayfun(@(c) sameEvents(Braw{c}, Braw{op}), others));
if ~isempty(sub), others = others(others ~= sub & ~arrayfun(@(c) sameEvents(Braw{c}, Braw{sub}), others)); end
sc = [op; sub; others];
nsc = numel(sc);
sb = Bref(sc); pks = Pk(sc); nb = nbur(sc);
inside_main = nan(nsc, 1); contains_main = nan(nsc, 1); gap_z = nan(nsc, 1); keep_sc = true(nsc, 1);
within = @(a, b) arrayfun(@(i) any(a(1, i) >= b(1, :) - 1e-9 & a(2, i) <= b(2, :) + 1e-9), 1:size(a, 2));
for c = 2:nsc
    inside_main(c) = mean(within(sb{c}, sb{1}));
    contains_main(c) = mean(within(sb{1}, sb{c}));
    % a level that groups (or splits) the main bursts is real only if the gaps it bridges between
    % the finer bursts are the short ones of a bimodal gap distribution (exponential gaps of the
    % same mean = no grouping); large windows otherwise just merge neighbours
    if max(inside_main(c), contains_main(c)) >= 0.5
        if contains_main(c) >= inside_main(c), fine = sb{1}; coarse = sb{c}; else, fine = sb{c}; coarse = sb{1}; end
        [gz, gg] = splitGapTest(fine, coarse, rs);
        gap_z(c) = gz;
        keep_sc(c) = gz >= opts.MinZ && gg >= opts.MinEffect;
    elseif max(mean(any(sb{c}(1, :)' <= sb{1}(2, :) & sb{c}(2, :)' >= sb{1}(1, :), 2)), ...
               mean(any(sb{1}(1, :)' <= sb{c}(2, :) & sb{1}(2, :)' >= sb{c}(1, :), 2))) >= 0.2
        keep_sc(c) = false;            % partial overlap: fragments of the main events, not a burst type
    end
end
info.Gate = true;
pn = mod(lvl_peak - 1, size(Zs, 1)) + 1; pt = floor((lvl_peak - 1) / size(Zs, 1)) + 1;
li = (1:nl)'; rp = rep; rp(rp == 0) = NaN; okr = ~isnan(rp);
refN = nan(nl, 1); refT = refN; refQ = nan(nl, 4); refNB = refN;
refN(okr) = Nv(cand(rp(okr), 1)); refT(okr) = thr(rp(okr)); refQ(okr, :) = q4(rp(okr), :); refNB(okr) = nbur(rp(okr)); refF = nan(nl, 1); refF(okr) = chance(rp(okr)) ./ nbur(rp(okr));
info.Levels = table(li, Nv(pn)', 10 .^ (centers(pt)' + opts.LogBin / 2), Zs(lvl_peak(:)), lvl_prom(:), refN, refT, refNB, refF, ...
    refQ(:, 1), refQ(:, 2), refQ(:, 3), refQ(:, 4), li == main_L, ...
    'VariableNames', {'Level', 'PeakN', 'PeakThreshold', 'PeakZ', 'Prominence', 'RefN', 'RefThreshold', 'RefNBursts', 'ChanceShare', ...
    'Contiguity', 'NoQuietInside', 'QuietGaps', 'Coverage', 'Main'});
info.Levels(isnan(info.Levels.RefN), :) = [];
if any(~keep_sc)
    r = find(~keep_sc);
    extra = table(Nv(cand(sc(r), 1))', thr(sc(r)), cs(sc(r)), cz(sc(r)), nb(r), mu(sc(r)), cand(sc(r), 3), ...
        repmat("no grouping", numel(r), 1), 'VariableNames', {'N', 'Threshold', 'Score', 'Z', 'NBursts', 'MedianUnits', 'Level', 'Reason'});
    if isempty(info.RejectedScales), info.RejectedScales = extra; else, info.RejectedScales = [info.RejectedScales; extra]; end
end
sc = sc(keep_sc); sb = sb(keep_sc); pks = pks(keep_sc); nb = nb(keep_sc);
inside_main = inside_main(keep_sc); contains_main = contains_main(keep_sc); gap_z = gap_z(keep_sc);
info.Scales = table(Nv(cand(sc, 1))', thr(sc), cg(sc), cs(sc), cz(sc), nb, mu(sc), mx(sc), q4(sc, 1), cand(sc, 3), ...
    inside_main, contains_main, gap_z, 'VariableNames', ...
    {'N', 'Threshold', 'Gain', 'Score', 'Z', 'NBursts', 'MedianUnits', 'RecruitedShare', 'Contiguity', 'Level', ...
    'Inside', 'Contains', 'GapZ'});
% Reported spans: first and last spike of each ISI_N run (no trimming; the rate-based edges
% above only serve the scale selection) that contains an accepted burst. Timing and shape:
% burstShape at ShapeLevel, the same point of every burst's own rate profile.
raw = cell(numel(sc), 1);
for c = 1:numel(sc)
    [r0, r1] = isinBursts(t, Nv(cand(sc(c), 1)), thr(sc(c)), NaN);
    A = Braw{sc(c)}; ok_r = any(r0(:) <= A(2, :) & r1(:) >= A(1, :), 2)';
    raw{c} = [r0(ok_r); r1(ok_r)];
end
% Superbursts: bursts inside a superburst are separated by silence relative to their own peaks,
% not to the background of the whole recording. A main burst whose population rate (kept units)
% stays below baseline + ShapeLevel x (its peak - baseline) for >= MinQuietGap is split there;
% pieces that recruit >= MinBurstUnits units beyond chance are bursts, the rest (e.g. the
% afterdischarge of a few units) is dropped. Split bursts keep the index of their superburst.
[raw{1}, sb_of, info.Superbursts] = splitBursts(raw{1}, t, iu, prev, rr, cc, duration, minq, opts.SplitLevel, min_units, opts.MinZ);
info.Scales.NBursts(1) = size(raw{1}, 2);
info.ScaleBursts = raw;
info.ScalePeaks = pks(:);                       % rate maxima of the selection episodes
bursts.T_start = raw{1}(1, :); bursts.T_end = raw{1}(2, :); bursts.Superburst = sb_of;
S = burstShape(spike_times, bursts.T_start, bursts.T_end, duration, Level=opts.ShapeLevel);
bursts.T_peak = S.Peak'; bursts.Onset = S.Onset'; bursts.Offset = S.Offset'; bursts.Duration = S.Duration';
bursts.Rise = S.Rise'; bursts.Fall = S.Fall'; bursts.FWHM = S.FWHM'; bursts.PeakRate = S.PeakRate';
bursts.ShapeLevel = opts.ShapeLevel;
% size of each main burst: kept units and their spikes between Onset and Offset (the first-to-
% last-spike span can include background spikes up to t on either side)
[bursts.NUnits, bursts.NSpikes] = deal(nan(1, numel(bursts.T_start)));
for kb = find(isfinite(bursts.Onset) & isfinite(bursts.Offset))
    j0 = searchSorted(t, bursts.Onset(kb), false); j1 = searchSorted(t, bursts.Offset(kb), true) - 1;
    bursts.NUnits(kb) = numel(unique(su(j0:j1))); bursts.NSpikes(kb) = j1 - j0 + 1;
end
bursts.best_N = Nv(cand(sc(1), 1)); bursts.best_ISI_N = thr(sc(1));
end

% =================================================================================================
function [lab, peak, prom] = levelBasins(A, minz, frac)
% Levels of the map A (NaN = outside): maxima whose prominence - height above the highest saddle
% to a higher maximum - is >= max(minz, frac x height), found by flooding from the top
% (8-neighbourhood, union-find). Lower bumps merge into the higher component at their saddle;
% every cell belongs to the level whose component it joined. lab: level index per cell (0 =
% none); peak: linear index of each level's maximum; prom: its prominence.
[r, c] = size(A); v = A(:); ok = isfinite(v);
[~, order] = sort(v, 'descend'); order = order(ok(order));
parent = zeros(numel(v), 1); pk = zeros(numel(v), 1); closed = false(numel(v), 1); pr = nan(numel(v), 1);
for idx = order'
    [i, j] = ind2sub([r c], idx); roots = [];
    for di = -1:1
        for dj = -1:1
            ii = i + di; jj = j + dj;
            if (di == 0 && dj == 0) || ii < 1 || ii > r || jj < 1 || jj > c, continue; end
            nbr = sub2ind([r c], ii, jj);
            if parent(nbr) == 0, continue; end
            while parent(nbr) ~= nbr, nbr = parent(nbr); end
            roots(end+1) = nbr; %#ok<AGROW>
        end
    end
    roots = unique(roots);
    if isempty(roots), parent(idx) = idx; pk(idx) = idx; continue; end
    [~, k] = max(v(pk(roots))); top = roots(k);
    for q = roots(:)'
        if q == top || closed(q), continue; end
        h = v(pk(q));
        if h - v(idx) >= max(minz, frac * h)
            closed(q) = true; pr(q) = h - v(idx);       % a level of its own
        else
            parent(q) = top;                            % a bump on the flank of a higher maximum
        end
    end
    parent(idx) = top;
end
root = zeros(numel(v), 1);
for idx = find(parent > 0)'
    q = idx; while parent(q) ~= q, q = parent(q); end
    root(idx) = q;
end
rts = unique(root(root > 0));
rts = rts(v(pk(rts)) >= minz);
for q = rts(~closed(rts))'                              % the global maximum, isolated islands
    pr(q) = v(pk(q)) - min(v(root == q));
end
rts = rts(pr(rts) >= max(minz, frac * v(pk(rts))));
[~, o] = sort(v(pk(rts)), 'descend'); rts = rts(o);
lab = zeros(r, c);
for L = 1:numel(rts), lab(root == rts(L)) = L; end
peak = pk(rts); prom = pr(rts);
end

function prev = prevSame(iu)
% index of the previous spike of the same unit (0 if none); iu = unit index per spike, time order
prev = zeros(numel(iu), 1); last = zeros(max(iu), 1);
for i = 1:numel(iu)
    prev(i) = last(iu(i)); last(iu(i)) = i;
end
end

function n = runUnits(prev, i0, i1)
% distinct units in each run of spikes i0(k):i1(k): spikes whose unit has not fired earlier in it
n = arrayfun(@(k) sum(prev(i0(k):i1(k)) < i0(k)), 1:numel(i0));
end

function sur = backgroundTrains(t, iu, B, pad, duration, K, rs)
% K surrogate trains of the background: the spikes outside the bursts B (+-pad), on a time axis
% with the bursts cut out, each unit circularly shifted by its own random offset
sur = {};
lo = max(B(1, :) - pad, 0); hi = min(B(2, :) + pad, duration);
[lo, o] = sort(lo); hi = hi(o);
L = lo(1); H = hi(1);
for i = 2:numel(lo)                                      % merge overlapping intervals
    if lo(i) <= H(end), H(end) = max(H(end), hi(i)); else, L(end+1) = lo(i); H(end+1) = hi(i); end %#ok<AGROW>
end
cut = [0, cumsum(H - L)];                                % time cut out before each interval
T = duration - cut(end);
if T <= 0, return; end
ns = discretize(t, [-inf, L, inf]) - 1;                  % intervals starting at or before each spike
inside = ns >= 1;
inside(inside) = t(inside) <= H(ns(inside))';
tc = t(~inside) - cut(ns(~inside) + 1)'; uk = iu(~inside);
sur = cell(K, 1);
for s = 1:K
    shift = T * rand(rs, max(iu), 1);
    [ts, o] = sort(mod(tc + shift(uk), T));
    sur{s} = struct('t', ts, 'prev', prevSame(uk(o)), 'u', uk(o), 'T', T);
end
end

function rate = unitRates(t, u, i0, i1, duration)
% rate of every unit outside the runs of spikes i0(k):i1(k) (its spikes there / time there)
M = numel(t); d = zeros(M + 1, 1);
d(i0) = d(i0) + 1; d(i1 + 1) = d(i1 + 1) - 1;
in = cumsum(d(1:M)) > 0;
Tout = max(duration - sum(t(i1) - t(i0)), 1);
rate = accumarray(u(~in), 1, [max(u) 1]) / Tout;
end

function ok = recruited(nu, rate, T, min_units, minz)
% bursts recruiting >= min_units distinct units, and minz SD more than the units expected to fire
% at least once by chance within their core duration T (Poisson at the given rates: each unit
% fires with p = 1 - exp(-rate x T); mean sum p, variance sum p (1 - p))
p = 1 - exp(-rate(:) * reshape(T, 1, []));
E = sum(p, 1); V = sum(p .* (1 - p), 1);
ok = reshape(nu, 1, []) >= min_units & reshape(nu, 1, []) - E >= minz * sqrt(V);
end

function T = coreDuration(t, i0, i1)
% time spanned by the middle 80 % of the spikes of each run i0(k):i1(k)
n = i1 - i0;
T = reshape(t(i0 + round(0.9 * n)) - t(i0 + round(0.1 * n)), 1, []);
end

function [B2, sb_of, SB] = splitBursts(B, t, iu, prev, r, c, duration, minq, level, min_units, minz)
% Split bursts at silences relative to their own peak: runs of >= minq bins of the rate r below
% base + level x (peak - base) inside a burst (not touching its span edges), base = median rate
% outside all bursts (+-0.1 s); the cut is at the middle of the run and the span edges stay.
% Pieces that recruit >= min_units units beyond chance (unit rates outside the bursts) are
% bursts; the rest is dropped. A burst split into >= 2 is a superburst (SB: 2 x K spans; sb_of: superburst of
% each burst, 0 = none).
nb = size(B, 2); B2 = zeros(2, 0); sb_of = zeros(1, 0); SB = zeros(2, 0);
if nb == 0, return; end
i0 = arrayfun(@(k) searchSorted(t, B(1, k), false), 1:nb);
i1 = arrayfun(@(k) searchSorted(t, B(2, k), true) - 1, 1:nb);
rate = unitRates(t, iu, i0, i1, duration);
inb = false(size(r));
for k = 1:nb, inb(c >= B(1, k) - 0.1 & c <= B(2, k) + 0.1) = true; end
base = median(r(~inb)); if isempty(base) || isnan(base), base = 0; end
for k = 1:nb
    w = find(c >= B(1, k) & c <= B(2, k));
    P = B(:, k);
    if numel(w) > minq
        th = base + level * (max(r(w)) - base);
        e = diff([0, r(w) < th, 0]); s = find(e == 1); f = find(e == -1) - 1;
        long = (f - s + 1) >= minq & s > 1 & f < numel(w);
        if any(long)
            cut = (c(w(s(long))) + c(w(f(long)))) / 2;
            lo = [B(1, k), cut]; hi = [cut, B(2, k)];
            P = zeros(2, 0);
            for q = find(hi > lo)
                j0 = searchSorted(t, lo(q), false); j1 = searchSorted(t, hi(q), true) - 1;
                if j1 <= j0, continue; end
                if recruited(runUnits(prev, j0, j1), rate, coreDuration(t, j0, j1), min_units, minz)
                    P(:, end+1) = [t(j0); t(j1)]; %#ok<AGROW>
                end
            end
            if isempty(P), P = B(:, k); end
        end
    end
    if size(P, 2) >= 2
        SB(:, end+1) = [P(1, 1); P(2, end)]; sb_of = [sb_of, repmat(size(SB, 2), 1, size(P, 2))]; %#ok<AGROW>
    else
        sb_of(end+1) = 0; %#ok<AGROW>
    end
    B2 = [B2, P]; %#ok<AGROW>
end
end

function n = chanceBursts(sur, N, thr, min_units, minz)
% mean number of ISI_N bursts (N, thr; >= min_units units beyond chance) in the surrogate trains
% (0 if none)
n = 0;
for s = 1:numel(sur)
    [~, ~, i0, i1] = isinBursts(sur{s}.t, N, thr, NaN);
    n = n + sum(recruited(runUnits(sur{s}.prev, i0, i1), unitRates(sur{s}.t, sur{s}.u, i0, i1, sur{s}.T), ...
        coreDuration(sur{s}.t, i0, i1), min_units, minz));
end
n = n / max(numel(sur), 1);
end

function [r, c, level] = popRate(t, duration, bin, sm)
% smoothed population rate and the background level (mean + 3 SD of the bins below the 80th
% percentile, i.e. outside bursts as long as they fill < 20 % of the time)
eb = 0:bin:duration + bin;
r = movmean(histcounts(t, eb), sm) / bin;
c = eb(1:end-1) + bin / 2;
bg = r(r <= prctile(r, 80));
level = mean(bg) + 3 * std(bg);
end

function [q, qin, qgap, cov] = contiguity(r, c, level, B, minquiet, minactive)
% How well bursts B match continuous above-background activity separated by silence: q = qin x
% qgap x cov, with qin the share of bursts without a quiet stretch >= minquiet bins inside, qgap
% the share of gaps containing one, cov the share of sustained active episodes (>= minactive
% bins above level) overlapping a burst (time-weighted)
nb = size(B, 2); q = 0; qin = 0; qgap = 0; cov = 0;
if nb == 0, return; end
bin = c(2) - c(1);
e = diff([0, r <= level, 0]); s = find(e == 1); f = find(e == -1) - 1;
long = (f - s + 1) >= minquiet; qs = c(s(long)); qf = c(f(long));
in_q = arrayfun(@(i) any(qs > B(1, i) & qf < B(2, i)), 1:nb);
qin = mean(~in_q); qgap = 1; cov = 1;
if nb > 1
    gap_q = arrayfun(@(j) any(min(qf, B(1, j+1)) - max(qs, B(2, j)) >= (minquiet - 1) * bin), 1:nb-1);
    qgap = mean(gap_q);
end
e = diff([0, r > level, 0]); s = find(e == 1); f = find(e == -1) - 1;
act = (f - s + 1) >= minactive;
if any(act)
    % strong episodes only (peak >= 25 % of the 90th percentile of episode peaks: in a nearly
    % silent background a few co-firing units exceed the level too), time-weighted
    pk = arrayfun(@(k) max(r(s(k):f(k))), 1:numel(s));
    act = act & pk >= 0.25 * prctile(pk(act), 90);
    as = c(s(act)); af = c(f(act)); w = f(act) - s(act) + 1;
    covered = arrayfun(@(k) any(as(k) <= B(2, :) & af(k) >= B(1, :)), 1:numel(as));
    cov = sum(w(covered)) / sum(w);
end
q = qin * qgap * cov;
end

function [z, gain] = splitGapTest(fine, coarse, rs)
% The coarser scale bridges some gaps between consecutive finer bursts (both inside one coarse
% burst) and not others. It reflects a real grouping only if that split separates short from
% long gaps: minimum-error gain of this split of the log gaps, z-scored against the best split
% of exponential gaps with the same mean (no grouping)
z = NaN; gain = NaN;
gaps = fine(1, 2:end) - fine(2, 1:end-1);
idx = nan(1, size(fine, 2));
for i = 1:size(fine, 2)
    h = find(fine(1, i) >= coarse(1, :) - 1e-9 & fine(2, i) <= coarse(2, :) + 1e-9, 1);
    if ~isempty(h), idx(i) = h; end
end
bridged = idx(1:end-1) == idx(2:end) & ~isnan(idx(1:end-1));
ok = gaps > 0; gaps = log10(gaps(ok)); bridged = bridged(ok);
if sum(bridged) < 3 || sum(~bridged) < 3 || mean(gaps(bridged)) >= mean(gaps(~bridged)), return; end
gain = splitGain(gaps, bridged);
K = 100; g0 = zeros(K, 1); m = mean(10 .^ gaps);
for k = 1:K
    g0(k) = kiGain(log10(-m * log(rand(rs, numel(gaps), 1))));
end
z = (gain - mean(g0)) / max(std(g0), 1e-6);
gain = gain - mean(g0);
end

function g = splitGain(x, cls)
% log-likelihood gain per sample of two Gaussians (classes cls / ~cls) over one
x = x(:); cls = logical(cls(:)); n = numel(x);
P1 = mean(cls); P2 = 1 - P1;
v1 = max(var(x(cls), 1), 1e-6); v2 = max(var(x(~cls), 1), 1e-6); vT = max(var(x, 1), 1e-6);
g = 0.5 * log(vT) - (P1 * 0.5 * log(v1) + P2 * 0.5 * log(v2)) + P1 * log(P1) + P2 * log(P2);
if n < 6, g = NaN; end
end

function g = kiGain(x)
% best two-Gaussian minimum-error split of the samples x: log-likelihood gain per sample over one
% Gaussian (each side keeps >= max(3, 5 %) of the samples)
x = sort(x(:)); n = numel(x); m = max(3, ceil(0.05 * n));
s1 = cumsum(x); q1 = cumsum(x.^2); k = (m:n-m)';
if isempty(k), g = NaN; return; end
P1 = k / n; P2 = 1 - P1;
v1 = max(q1(k) ./ k - (s1(k) ./ k).^2, 1e-6);
v2 = max((q1(end) - q1(k)) ./ (n - k) - ((s1(end) - s1(k)) ./ (n - k)).^2, 1e-6);
vT = max(q1(end) / n - (s1(end) / n)^2, 1e-6);
g = max(0.5 * log(vT) - (P1 .* 0.5 .* log(v1) + P2 .* 0.5 .* log(v2)) + P1 .* log(P1) + P2 .* log(P2));
end

function [B, pk] = refineEdges(t, B, duration, bin, sm)
% Burst edges from the population rate of the kept units (bin s, sm-bin moving mean): each burst
% becomes the sustained episode above background + 2 SD (runs >= 2 x sm bins) that overlaps its
% ISI_N core (+-50 ms), extended outwards to the ends of those runs (e.g. afterdischarges) but
% never into the neighbouring bursts. Only the edges move: dips inside a burst (sub-events) do
% not split it. Background = the rate outside all bursts (padded by 0.1 s). Bursts without such
% a run keep their ISI_N edges. pk: time of the rate maximum of each burst.
nbr = size(B, 2); pk = nan(1, nbr);
if nbr == 0, return; end
eb = 0:bin:duration + bin;
r = movmean(histcounts(t, eb), sm) / bin;
c = eb(1:end-1) + bin / 2;
inb = false(size(r));
for i = 1:nbr, inb(c >= B(1, i) - 0.1 & c <= B(2, i) + 0.1) = true; end
level = Inf;
if sum(~inb) >= 10, level = mean(r(~inb)) + 2 * std(r(~inb)); end
e = diff([0, r > level, 0]); rs0 = find(e == 1); rs1 = find(e == -1) - 1;
long = (rs1 - rs0 + 1) >= 2 * sm; rs0 = rs0(long); rs1 = rs1(long);
pad = 0.05; B0 = B;
for i = 1:nbr
    lo = -Inf; if i > 1, lo = B0(2, i-1); end
    hi = Inf; if i < nbr, hi = B0(1, i+1); end
    hit = find(c(rs1) >= B0(1, i) - pad & c(rs0) <= B0(2, i) + pad);
    if ~isempty(hit)
        B(1, i) = max(eb(rs0(hit(1))), lo); B(2, i) = min(eb(rs1(hit(end)) + 1), hi);
    end
    kk = find(c >= B(1, i) & c <= B(2, i));
    if isempty(kk), [~, kk] = min(abs(c - mean(B(:, i)))); end
    [~, m] = max(r(kk)); pk(i) = c(kk(m));
end
end

function tf = sameEvents(A, B)
% the same bursts: similar counts and >= 80 % of each set overlapping the other
na = size(A, 2); nb = size(B, 2);
tf = false;
if min(na, nb) / max(na, nb) < 0.8, return; end
ov = @(X, Y) mean(any(X(1, :)' <= Y(2, :) & X(2, :)' >= Y(1, :), 2));
tf = ov(A, B) >= 0.8 && ov(B, A) >= 0.8;
end

function G = gainMap(t, Nv, edges, min_frac)
% Minimum-error thresholding of log10 ISI_N at each bin edge: per-window log-likelihood gain of
% two Gaussians (one per side of the split) over one; NaN where a side holds < min_frac
nb = numel(edges) - 1; c = (edges(1:end-1) + edges(2:end)) / 2; w = edges(2) - edges(1);
vmin = w^2 / 12;
G = nan(numel(Nv), nb);
for i = 1:numel(Nv)
    N = Nv(i);
    if numel(t) <= N, continue; end
    x = log10(max(t(N:end) - t(1:end-N+1), 10^edges(1)));
    h = histcounts(x, edges); H = sum(h);
    if H == 0, continue; end
    p = h / H;
    P1 = cumsum(p); m1 = cumsum(p .* c); q1 = cumsum(p .* c.^2);
    P2 = 1 - P1; m2 = m1(end) - m1; q2 = q1(end) - q1;
    v1 = max(q1 ./ max(P1, eps) - (m1 ./ max(P1, eps)).^2, vmin);
    v2 = max(q2 ./ max(P2, eps) - (m2 ./ max(P2, eps)).^2, vmin);
    vT = max(q1(end) - m1(end)^2, vmin);
    g = 0.5 * log(vT) - (P1 .* 0.5 .* log(v1) + P2 .* 0.5 .* log(v2)) ...
        + P1 .* log(max(P1, eps)) + P2 .* log(max(P2, eps));
    m = max(min_frac, 50 / H);               % each side >= min_frac of the windows and >= 50 windows
    g(P1 < m | P2 < m) = NaN;
    G(i, :) = g;
end
end

function [T_start, T_end, i0, i1] = isinBursts(t, N, thr, edge_isi)
% Runs of consecutive spikes linked by at least one window of N spikes spanning <= thr: link i
% (spike i to i+1) is in a burst if a qualifying window j <= i < j+N-1 contains both spikes.
% i0 / i1: index of the first / last spike of each burst.
% Edges: a window that ends inside a burst can start with background spikes up to thr earlier.
% Edge spikes are dropped while the kw = min(N, 10) spikes from the edge inwards span more than
% (kw - 1) x edge_isi (default thr / (N - 1); [] = default, NaN = no trimming), so bursts start
% and end where the local rate crosses 1 / edge_isi.
if nargin < 4 || isempty(edge_isi), edge_isi = thr / (N - 1); end
t = t(:); M = numel(t);
T_start = zeros(1, 0); T_end = T_start; i0 = T_start; i1 = T_start;
if M < max(N, 2), return; end
js = find(t(N:end) - t(1:end-N+1) <= thr);
if isempty(js), return; end
d = zeros(M, 1);
d(js) = d(js) + 1; d(js + N - 1) = d(js + N - 1) - 1;
link = cumsum(d(1:M-1)) > 0;
e = diff([0; link; 0]);
i0 = find(e == 1)'; i1 = (find(e == -1) - 1)' + 1;
if isfinite(edge_isi) && N > 1
    kw = min(N, 10); span = (kw - 1) * edge_isi;
    for k = 1:numel(i0)
        while i1(k) - i0(k) >= kw && t(i0(k) + kw - 1) - t(i0(k)) > span, i0(k) = i0(k) + 1; end
        while i1(k) - i0(k) >= kw && t(i1(k)) - t(i1(k) - kw + 1) > span, i1(k) = i1(k) - 1; end
    end
end
T_start = t(i0)'; T_end = t(i1)';
end

function i = searchSorted(x, v, after)
% first index of sorted x with x >= v (after = false) or x > v (after = true); numel(x) + 1 if none
lo = 1; hi = numel(x) + 1;
while lo < hi
    m = floor((lo + hi) / 2);
    if x(m) < v || (after && x(m) == v), lo = m + 1; else, hi = m; end
end
i = lo;
end

function U = unitMetrics(t, su, duration, opts, rs)
% Population coupling and Gini excess per unit, each z-scored against its own surrogates
units = unique(su); nu = numel(units); K = opts.NSurrogates;
[~, iu] = ismember(su, units);
n_spk = accumarray(iu, 1, [nu 1]);
% population coupling
eb = 0:opts.CouplingBin:duration + opts.CouplingBin;
X = zeros(nu, numel(eb) - 1);
b = discretize(t, eb);
ok = ~isnan(b);
X(:) = accumarray(sub2ind(size(X), iu(ok), b(ok)), 1, [numel(X) 1]);
pop = sum(X, 1);
cz = nan(nu, 1); cc = cz;
for u = 1:nu
    if n_spk(u) < opts.MinUnitSpikes, continue; end
    x = X(u, :); other = pop - x;
    cc(u) = corrFast(x, other);
    cn = zeros(K, 1);
    for kk = 1:K
        cn(kk) = corrFast(circshift(x, randi(rs, numel(x))), other);
    end
    cz(u) = (cc(u) - mean(cn)) / max(std(cn), 1e-6);
end
% Gini of binned counts against Poisson firing at the unit's rate: GiniRatio = Gini / expected
% Gini of Poisson counts with the same mean (exp(-2 mu) (I0(2 mu) + I1(2 mu))); GiniZ against
% homogeneous Poisson trains with the same spike count. Regular (tonic) units fall below 1,
% bursty units above; sparse units (< ~1 spike per bin) stay near 1 whatever their timing.
gbin = opts.GiniBin;
g = nan(nu, 1); gz = g; gx = g;
for u = 1:nu
    if n_spk(u) < opts.MinUnitSpikes, continue; end
    if (isstring(gbin) || ischar(gbin)) && string(gbin) == "adaptive", bu = opts.GiniSpikesPerBin / (n_spk(u) / duration); else, bu = gbin; end
    gb = 0:bu:duration + bu;
    x = histcounts(t(iu == u), gb);
    g(u) = giniCounts(x);
    mu = mean(x);
    gx(u) = g(u) / (besseli(0, 2 * mu, 1) + besseli(1, 2 * mu, 1));
    gn = zeros(K, 1);
    for kk = 1:K
        gn(kk) = giniCounts(histcounts(duration * rand(rs, n_spk(u), 1), gb));
    end
    gz(u) = (g(u) - mean(gn)) / max(std(gn), 1e-6);
end
U = table(units(:), n_spk, n_spk / duration, cc, cz, g, gx, gz, 'VariableNames', ...
    {'Unit', 'NSpikes', 'Rate', 'Coupling', 'CouplingZ', 'Gini', 'GiniRatio', 'GiniZ'});
end

function r = corrFast(a, b)
a = a - mean(a); b = b - mean(b);
r = sum(a .* b) / max(sqrt(sum(a.^2) * sum(b.^2)), eps);
end

function G = giniCounts(x)
x = sort(x(:)); n = numel(x); s = sum(x);
if s == 0, G = NaN; return; end
G = 2 * sum((1:n)' .* x) / (n * s) - (n + 1) / n;
end
