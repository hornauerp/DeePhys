function [network_table, unit_table] = computeCellTypeBurstEI(spike_times, spike_units, labels, ...
    bursts, duration, opts)
% COMPUTECELLTYPEBURSTEI  E / I burst features of one recording, against a rate-matched null.
%
%   [network_table, unit_table] = computeCellTypeBurstEI(spike_times, spike_units, labels, ...
%       bursts, duration)
%   [...] = computeCellTypeBurstEI(..., NNull=200, NullStrata=5, Seed=1, MinUnits=3, MinBursts=5)
%
% Bursts are timed by their onset and offset (detectNetworkBursts), not by the first and last
% spike of a run: a run's first spike can be a background spike before the burst, which would
% give the population of that spike a lead time of exactly 0 by construction.
%
% Network features (one row):
%   BurstEI_<m>, BurstEI_<m>_Excess  every burstEIMetrics metric m (epoch folds, I/E contrasts,
%                                    lead / terminal / after-peak measures) and its value minus
%                                    the median of random "I" sets matched in firing rate
%   BurstUnit<x>_E, BurstUnit<x>_I   mean over E / I units of burstUnitTiming's per-unit x
%                                    (Latency, Participation, PreRatio, InRatio, PostRatio)
%   BurstUnit<x>_IminusE, BurstUnit<x>_IminusE_Excess  I mean minus E mean, and that difference
%                                    minus its median over the same rate-matched random I sets
% Per-unit measures (participation, latency) depend on firing rate: faster units take part in
% more bursts and fire earlier by chance. The _Excess columns are what the labels add beyond
% firing rate; a label that only reflects rate gives values near 0.
%
% Unit features (one row per unit 1..U, label-free): BurstUnitLatency, BurstUnitParticipation,
%   BurstUnitPreRatio, BurstUnitInRatio, BurstUnitPostRatio (see burstUnitTiming).
%
% Too few bursts (MinBursts) or units of either type (MinUnits) give NaN network features.
%
% INPUTS:
%   spike_times, spike_units - (S x 1) spike times (s) and unit indices 1..U
%   labels   - (1 x U) 1 = excitatory, 2 = inhibitory, NaN = left out of both populations
%   bursts   - struct with Onset, Offset, T_peak (detectNetworkBursts)
%   duration - recording duration (s)
%
% See also burstEIMetrics, burstUnitTiming, cellTypeBurstFeatureNames.

arguments
    spike_times (:,1) double
    spike_units (:,1) double
    labels      (1,:) double
    bursts      struct
    duration    (1,1) double
    opts.NNull      (1,1) double = 200
    opts.NullStrata (1,1) double = 5
    opts.Seed       (1,1) double = 1
    opts.MinUnits   (1,1) double = 3
    opts.MinBursts  (1,1) double = 5
end

U = numel(labels);
[net_names, unit_names, metric_names] = cellTypeBurstFeatureNames();
timing_vars = erase(unit_names, "BurstUnit");
net = nan(1, numel(net_names));

T = burstUnitTiming(spike_times, spike_units, bursts, duration, NUnits=U);
unit_table = array2table(T{:, cellstr(timing_vars)}, 'VariableNames', cellstr(unit_names));

R = burstEIMetrics(spike_times, spike_units, labels, bursts, duration, NNull=opts.NNull, ...
    NullStrata=opts.NullStrata, Seed=opts.Seed, MinUnits=opts.MinUnits, MinBursts=opts.MinBursts);
if R.Valid
    [~, loc] = ismember(metric_names, string(R.Metrics.Metric));
    val = R.Metrics.Value(loc)'; nul = R.Metrics.NullMedian(loc)';
    net(ismember(net_names, "BurstEI_" + metric_names)) = val;
    net(ismember(net_names, "BurstEI_" + metric_names + "_Excess")) = val - nul;

    % Per-unit timing by type, against the same kind of null as burstEIMetrics: random I sets
    % with as many units per firing-rate stratum as the I population has
    isI = labels(:) == 2; isE = labels(:) == 1; lab = isI | isE;
    li = find(lab);
    rate = T.Rate;
    q = prctile(rate(li), linspace(0, 100, opts.NullStrata + 1));
    stratum = discretize(rate(li), unique([-inf, q(2:end-1), inf]));
    rs = RandStream('twister', 'Seed', opts.Seed);
    X = T{:, cellstr(timing_vars)};
    null_diff = nan(opts.NNull, numel(timing_vars));
    for r = 1:opts.NNull
        sI = false(U, 1);
        for s = unique(stratum)'
            pool = li(stratum == s); n_s = sum(isI(pool));
            if n_s > 0, sI(pool(randperm(rs, numel(pool), min(n_s, numel(pool))))) = true; end
        end
        null_diff(r, :) = mean(X(sI, :), 1, 'omitnan') - mean(X(lab & ~sI, :), 1, 'omitnan');
    end
    mE = mean(X(isE, :), 1, 'omitnan'); mI = mean(X(isI, :), 1, 'omitnan');
    d = mI - mE;
    for k = 1:numel(timing_vars)
        base = "BurstUnit" + timing_vars(k);
        net(net_names == base + "_E") = mE(k);
        net(net_names == base + "_I") = mI(k);
        net(net_names == base + "_IminusE") = d(k);
        net(net_names == base + "_IminusE_Excess") = d(k) - median(null_diff(:, k), 'omitnan');
    end
end
network_table = array2table(net, 'VariableNames', cellstr(net_names));
end
