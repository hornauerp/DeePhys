function R = burstEIMetrics(spike_times, spike_units, labels, bursts, duration, opts)
% BURSTEIMETRICS  E / I population activity around network bursts, against a rate-matched null.
%
%   R = burstEIMetrics(spike_times, spike_units, labels, bursts, duration)
%   R = burstEIMetrics(..., NNull=500, PreWindow=0.3, Seed=2)
%
% Every unit's firing is summarised once, so any population (and every null draw) is an
% average over its units:
%   epoch rates (Hz) from each burst's own onset and offset -
%     In      [Onset, Offset] of every burst
%     Pre     PreWindow before each onset (outside every burst)
%     Post    PostWindow after each offset (outside every burst and every Pre window)
%     Between outside [Onset - Margin, Offset + Margin] of every burst (the interburst reference)
%   burst-aligned mean rates (Hz, BinSize bins over Window) aligned at Onset, T_peak and Offset.
% Population rates are means over units, I / E ratios are taken after averaging (bins without E
% spikes are not dropped), and folds are relative to the population's own Between rate.
%
% Metrics (R.Metrics, one row each):
%   IFoldIn, EFoldIn, IEContrastIn     In / Between of I and of E, and their ratio (I/E lower
%                                      during bursts: IEContrastIn < 1)
%   IFoldPre, EFoldPre, IEContrastPre  the same for the window before onset
%   IFoldPost, EFoldPost, IEContrastPost  the same after offset (terminal I activity)
%   ILeadTime   time of the I maximum of the onset-aligned trace in LeadWindow (s, re onset)
%   IEMaxTime   time of the I/E maximum of the onset-aligned traces in LeadWindow (s, re onset)
%   IEMinTime   time of the I/E minimum of the peak-aligned traces in MinWindow (s, re peak)
%   ITerminal   maximum of the offset-aligned I trace in TerminalWindow / I Between rate
%   ITerminalVsIn  the same / I In rate (a rebound after the in-burst silence)
%   IAfterPeak  after-peak of I: height of the largest maximum of the peak-aligned I trace in
%               AfterWindow above the lowest point between the burst peak and it, / I Between
%               rate (0 if I only decays after the peak); IAfterPeakTime its time (s, re peak)
%   EAfterPeak  the same for E (an after-peak in both = network reverberation, not I alone)
%   IAfterRel, EAfterRel  after-peak height / the range (max - min) of the population's peak-
%               aligned trace, 0..1, comparable between I and E
%   AntiPhase   bursts inside a superburst after its first (bursts.Superburst):
%               (I in the gaps between them / I inside them) / (same for E); NaN without
%               >= MinInnerBursts such bursts
% Traces are smoothed over SmoothBins bins before times and maxima are read. Folds are NaN when a
% population never fires between bursts or the Between epoch is shorter than MinBetween (frequent
% bursts leave no interburst reference).
%
% Null: NNull random "I" sets with the same number of units in each firing-rate stratum
% (NullStrata quantiles of the labelled units' rates), drawn from all labelled units; the other
% labelled units are E. A label that only reflects firing rate gives values inside the null.
%
% INPUTS:
%   spike_times, spike_units - (S x 1) spike times (s) and unit indices 1..U
%   labels   - (1 x U) 1 = excitatory, 2 = inhibitory, NaN = left out of both populations
%   bursts   - struct with Onset, Offset, T_peak (1 x B, s); optional Superburst (1 x B, 0 = none)
%              as returned by detectNetworkBursts
%   duration - recording duration (s)
% OPTIONS:
%   PreWindow (0.3), PostWindow (0.3), Margin (0.5) s; BinSize (0.01) s; Window ([-1.5 1.5]) s;
%   SmoothBins (3); LeadWindow ([-0.5 0.05]); MinWindow ([-0.2 0.3]); TerminalWindow
%   ([-0.05 0.2]); AfterWindow ([0.02 0.3]); MinBetween (10 s); NNull (200), NullStrata (5), Seed (1); MinUnits (3) per population;
%   MinBursts (5); MinInnerBursts (5)
% OUTPUT (struct):
%   Valid    - false if too few bursts or units (metrics NaN)
%   Metrics  - table: Metric, Value, NullMedian, NullLow, NullHigh (2.5 / 97.5 %), PLower
%              (share of null draws <= Value), PUpper (share >= Value)
%   Rates    - table per population (I, E): In, Pre, Post, Between (Hz per unit)
%   Traces   - struct per alignment (Onset, Peak, Offset): Time, I, E, IE, IENullLow, IENullHigh
%   NBursts, NI, NE, NInnerBursts, Units (table: Unit, Label, Rate, In, Pre, Post, Between)
%
% See also detectNetworkBursts, burstUnitTiming, burstRhythm, EIAnalyzer.

arguments
    spike_times (:,1) double
    spike_units (:,1) double
    labels      (1,:) double
    bursts      struct
    duration    (1,1) double
    opts.PreWindow      (1,1) double = 0.3
    opts.PostWindow     (1,1) double = 0.3
    opts.Margin         (1,1) double = 0.5
    opts.BinSize        (1,1) double = 0.01
    opts.Window         (1,2) double = [-1.5 1.5]
    opts.SmoothBins     (1,1) double = 3
    opts.LeadWindow     (1,2) double = [-0.5 0.05]
    opts.MinWindow      (1,2) double = [-0.2 0.3]
    opts.TerminalWindow (1,2) double = [-0.05 0.2]
    opts.AfterWindow    (1,2) double = [0.02 0.3]
    opts.MinBetween     (1,1) double = 10
    opts.NNull         (1,1) double = 200
    opts.NullStrata     (1,1) double = 5
    opts.Seed           (1,1) double = 1
    opts.MinUnits       (1,1) double = 3
    opts.MinBursts      (1,1) double = 5
    opts.MinInnerBursts (1,1) double = 5
end

names = ["IFoldIn", "EFoldIn", "IEContrastIn", "IFoldPre", "EFoldPre", "IEContrastPre", ...
    "IFoldPost", "EFoldPost", "IEContrastPost", "ILeadTime", "IEMaxTime", "IEMinTime", ...
    "ITerminal", "ITerminalVsIn", "IAfterPeak", "IAfterPeakTime", "EAfterPeak", "IAfterRel", "EAfterRel", "AntiPhase"];
U = numel(labels);
isI = labels == 2; isE = labels == 1; lab = isI | isE;
on = bursts.Onset(:)'; off = bursts.Offset(:)'; pk = bursts.T_peak(:)';
ok = isfinite(on) & isfinite(off) & isfinite(pk) & off >= on;
sb = zeros(size(on)); if isfield(bursts, 'Superburst') && ~isempty(bursts.Superburst), sb = bursts.Superburst(:)'; end
on = on(ok); off = off(ok); pk = pk(ok); sb = sb(ok);
R = emptyResult(names, sum(ok), sum(isI), sum(isE));
if numel(on) < opts.MinBursts || sum(isI) < opts.MinUnits || sum(isE) < opts.MinUnits, return; end
keep = spike_units >= 1 & spike_units <= U & spike_times >= 0 & spike_times < duration;
st = spike_times(keep); su = spike_units(keep);

% -- per-unit epoch rates on a 1 ms grid --------------------------------------------------
dt = 1e-3; ng = ceil(duration / dt);
inM = gridMask([on; off], dt, ng);
preM = gridMask([on - opts.PreWindow; on], dt, ng) & ~inM;
postM = gridMask([off; off + opts.PostWindow], dt, ng) & ~inM & ~preM;
betM = ~gridMask([on - opts.Margin; off + opts.Margin], dt, ng);
gi = min(floor(st / dt) + 1, ng);
epochs = {inM, preM, postM, betM};
rates = zeros(U, 4);
for e = 1:4
    m = epochs{e};
    rates(:, e) = accumarray(su(m(gi)), 1, [U 1]) / max(sum(m) * dt, eps);
end
if sum(betM) * dt < opts.MinBetween, rates(:, 4) = NaN; end     % too little interburst time
rate_all = accumarray(su, 1, [U 1]) / duration;

% -- per-unit aligned mean traces ---------------------------------------------------------
edges = opts.Window(1):opts.BinSize:opts.Window(2);
time = edges(1:end-1) + opts.BinSize / 2;
tr.Onset = alignedRates(st, su, on, edges, U);
tr.Peak = alignedRates(st, su, pk, edges, U);
tr.Offset = alignedRates(st, su, off, edges, U);

% -- superbursts: gaps between consecutive bursts of one superburst ------------------------
inner = [false, sb(2:end) > 0 & sb(2:end) == sb(1:end-1)];
gapC = []; insC = [];
if sum(inner) >= opts.MinInnerBursts
    k = find(inner); gaps = [off(k - 1); on(k)]; gaps = gaps(:, gaps(2, :) > gaps(1, :));
    gM = gridMask(gaps, dt, ng); iM = gridMask([on(k); off(k)], dt, ng);
    gapC = accumarray(su(gM(gi)), 1, [U 1]) / max(sum(gM) * dt, eps);
    insC = accumarray(su(iM(gi)), 1, [U 1]) / max(sum(iM) * dt, eps);
end

P = struct('time', time, 'lead', time >= opts.LeadWindow(1) & time <= opts.LeadWindow(2), ...
    'minw', time >= opts.MinWindow(1) & time <= opts.MinWindow(2), ...
    'term', time >= opts.TerminalWindow(1) & time <= opts.TerminalWindow(2), 'sm', opts.SmoothBins, ...
    'after', time >= opts.AfterWindow(1) & time <= opts.AfterWindow(2), 'zero', find(time >= 0, 1));
[val, trv] = metricsOf(isI, isE, rates, tr, gapC, insC, P);

% -- rate-matched null ----------------------------------------------------------------------
rs = RandStream('twister', 'Seed', opts.Seed);
li = find(lab);
q = prctile(rate_all(li), linspace(0, 100, opts.NullStrata + 1)); q = q(:)';
stratum = discretize(rate_all(li), unique([-inf, q(2:end-1), inf]));
nullv = nan(opts.NNull, numel(names)); nullie = struct('Onset', [], 'Peak', [], 'Offset', []);
for f = ["Onset", "Peak", "Offset"], nullie.(f) = nan(opts.NNull, numel(time)); end
for r = 1:opts.NNull
    sI = false(1, U);
    for s = unique(stratum)'
        pool = li(stratum == s); n_s = sum(isI(pool));
        if n_s > 0, sI(pool(randperm(rs, numel(pool), min(n_s, numel(pool))))) = true; end
    end
    [nullv(r, :), tn] = metricsOf(sI, lab & ~sI, rates, tr, gapC, insC, P);
    for f = ["Onset", "Peak", "Offset"], nullie.(f)(r, :) = tn.(f); end
end

% -- output -----------------------------------------------------------------------------------
R.Valid = true;
pl = mean(nullv <= val, 1); pu = mean(nullv >= val, 1);
pl(isnan(val)) = NaN; pu(isnan(val)) = NaN;
R.Metrics = table(names(:), val(:), median(nullv, 1, 'omitnan')', prctile(nullv, 2.5, 1)', ...
    prctile(nullv, 97.5, 1)', pl(:), pu(:), 'VariableNames', ...
    {'Metric', 'Value', 'NullMedian', 'NullLow', 'NullHigh', 'PLower', 'PUpper'});
pop = @(sel) mean(rates(sel, :), 1);
R.Rates = array2table([pop(isI); pop(isE)], 'VariableNames', {'In', 'Pre', 'Post', 'Between'}, ...
    'RowNames', {'I', 'E'});
for f = ["Onset", "Peak", "Offset"]
    R.Traces.(f) = struct('Time', time, 'I', mean(tr.(f)(isI, :), 1), 'E', mean(tr.(f)(isE, :), 1), ...
        'IE', trv.(f), 'IENullLow', prctile(nullie.(f), 2.5, 1), 'IENullHigh', prctile(nullie.(f), 97.5, 1));
end
R.NInnerBursts = sum(inner);
lb = labels(:); lb(~lab) = NaN;
R.Units = table((1:U)', lb, rate_all, rates(:, 1), rates(:, 2), rates(:, 3), rates(:, 4), ...
    'VariableNames', {'Unit', 'Label', 'Rate', 'In', 'Pre', 'Post', 'Between'});
end

% =============================================================================================
function [v, ie] = metricsOf(isI, isE, rates, tr, gapC, insC, P)
% all metrics of one I / E split; ie = smoothed I/E traces per alignment
rI = mean(rates(isI, :), 1); rE = mean(rates(isE, :), 1);         % In, Pre, Post, Between
fI = dv(rI(1:3), rI(4)); fE = dv(rE(1:3), rE(4));            % NaN without a Between rate
sm = @(x) movmean(x, P.sm);
for f = ["Onset", "Peak", "Offset"]
    I.(f) = sm(mean(tr.(f)(isI, :), 1)); E.(f) = sm(mean(tr.(f)(isE, :), 1));
    ie.(f) = I.(f) ./ max(E.(f), eps);
end
t = P.time;
[~, k] = max(I.Onset .* P.lead - ~P.lead); lead = t(k);
[~, k] = max(ie.Onset .* P.lead - ~P.lead); iemax = t(k);
x = ie.Peak; x(~P.minw) = inf; [~, k] = min(x); iemin = t(k);
term = max(I.Offset(P.term));
[ia, ta] = afterPeak(I.Peak, P); ea = afterPeak(E.Peak, P);
anti = NaN;
if ~isempty(gapC)
    anti = dv(dv(mean(gapC(isI)), mean(insC(isI))), dv(mean(gapC(isE)), mean(insC(isE))));
end
v = [fI(1), fE(1), dv(fI(1), fE(1)), fI(2), fE(2), dv(fI(2), fE(2)), ...
    fI(3), fE(3), dv(fI(3), fE(3)), lead, iemax, iemin, ...
    dv(term, rI(4)), dv(term, rI(1)), dv(ia, rI(4)), ta, dv(ea, rE(4)), ...
    dv(ia, max(I.Peak) - min(I.Peak)), dv(ea, max(E.Peak) - min(E.Peak)), anti];
end

function r = dv(a, b)
% a ./ b, NaN where b is not positive or not finite (undefined fold)
r = a ./ b; bad = ~(b > 0) | ~isfinite(b);
if isscalar(b), if bad, r(:) = NaN; end, else, r(bad) = NaN; end
end

function [h, t] = afterPeak(x, P)
% largest rise of x within P.after above its running minimum since the burst peak (time 0)
k = P.zero:find(P.after, 1, 'last');
rise = x(k) - cummin(x(k)); rise(~P.after(k)) = -inf;
[h, j] = max(rise); t = P.time(k(j));
if h <= 0, h = 0; t = NaN; end
end

function M = alignedRates(st, su, t0, edges, U)
% (U x T) mean rate (Hz) of each unit around the times t0, bins given by edges (relative)
[ss, o] = sort(st); uu = su(o);
nb = numel(edges) - 1; bin = edges(2) - edges(1);
C = zeros(U, nb);
for k = 1:numel(t0)
    i0 = searchFirst(ss, t0(k) + edges(1)); i1 = searchFirst(ss, t0(k) + edges(end)) - 1;
    if i1 < i0, continue; end
    b = min(floor((ss(i0:i1) - t0(k) - edges(1)) / bin) + 1, nb);
    C = C + accumarray([uu(i0:i1), b], 1, [U nb]);
end
M = C / (numel(t0) * bin);
end

function i = searchFirst(x, v)
% first index with x(i) >= v (numel(x) + 1 if none); x sorted ascending
lo = 1; hi = numel(x) + 1;
while lo < hi
    mid = floor((lo + hi) / 2);
    if x(mid) < v, lo = mid + 1; else, hi = mid; end
end
i = lo;
end

function m = gridMask(iv, dt, ng)
% logical grid (ng bins of dt) covering the union of the intervals iv (2 x n)
d = zeros(ng + 1, 1);
a = max(floor(iv(1, :) / dt) + 1, 1); b = min(floor(iv(2, :) / dt) + 1, ng + 1);
v = b > a;
d = d + accumarray(a(v)', 1, [ng + 1 1]) - accumarray(b(v)', 1, [ng + 1 1]);
m = cumsum(d(1:ng)) > 0;
end

function R = emptyResult(names, nb, ni, ne)
n = numel(names);
R = struct('Valid', false, 'Metrics', table(names(:), nan(n, 1), nan(n, 1), nan(n, 1), nan(n, 1), ...
    nan(n, 1), nan(n, 1), 'VariableNames', {'Metric', 'Value', 'NullMedian', 'NullLow', 'NullHigh', 'PLower', 'PUpper'}), ...
    'Rates', table(), 'Traces', struct(), 'NBursts', nb, 'NI', ni, 'NE', ne, 'NInnerBursts', 0, 'Units', table());
end
