function T = burstUnitTiming(spike_times, spike_units, bursts, duration, opts)
% BURSTUNITTIMING  Each unit's firing relative to network burst onsets.
%
%   T = burstUnitTiming(spike_times, spike_units, bursts, duration)
%   T = burstUnitTiming(..., PreWindow=0.5, NUnits=120)
%
% Per unit, over the bursts with a finite onset and offset:
%   Latency       median over bursts of the first spike in [Onset - PreWindow, Offset],
%                 relative to Onset (s; negative = fires before the burst starts); bursts
%                 without such a spike do not count
%   Participation share of bursts with at least one spike in [Onset, Offset]
%   PreRatio      rate in [Onset - PreWindow, Onset) (outside every burst) / the unit's mean rate
%   InRatio       rate in [Onset, Offset] / mean rate
%   PostRatio     rate in (Offset, Offset + PostWindow] (outside bursts and Pre windows) / mean rate
% These are defined per unit, so labels can be applied afterwards, soft or partial, and E and I
% compared within firing-rate strata.
%
% INPUTS:
%   spike_times, spike_units - (S x 1); bursts - struct with Onset, Offset (detectNetworkBursts)
%   duration - recording duration (s)
% OPTIONS:
%   PreWindow (0.3 s), PostWindow (0.3 s), NUnits ([] = max(spike_units))
% OUTPUT:
%   T - table, one row per unit 1..NUnits: Unit, NSpikes, Rate, Latency, Participation,
%       PreRatio, InRatio, PostRatio
%
% See also burstEIMetrics, detectNetworkBursts, EIAnalyzer.computeUnitTiming.

arguments
    spike_times (:,1) double
    spike_units (:,1) double
    bursts      struct
    duration    (1,1) double
    opts.PreWindow  (1,1) double = 0.3
    opts.PostWindow (1,1) double = 0.3
    opts.NUnits     double = []
end
U = opts.NUnits; if isempty(U), U = max([spike_units; 0]); end
on = bursts.Onset(:)'; off = bursts.Offset(:)';
ok = isfinite(on) & isfinite(off) & off >= on; on = on(ok); off = off(ok); nb = numel(on);
keep = spike_units >= 1 & spike_units <= U; st = spike_times(keep); su = spike_units(keep);
nsp = accumarray(su, 1, [U 1]); rate = nsp / duration;
lat = nan(U, 1); part = nan(U, 1); pre = nan(U, 1); inr = nan(U, 1); post = nan(U, 1);
if nb > 0
    % first spike of each unit in each burst window
    [ss, o] = sort(st); uu = su(o);
    first = nan(U, nb); hit = false(U, nb);
    for k = 1:nb
        i0 = find(ss >= on(k) - opts.PreWindow, 1); if isempty(i0), continue; end
        i1 = find(ss <= off(k), 1, 'last'); if isempty(i1) || i1 < i0, continue; end
        t = ss(i0:i1); u = uu(i0:i1);
        [uq, ia] = unique(u, 'first');                % spikes are time-sorted: first spike per unit
        first(uq, k) = t(ia) - on(k);
        hit(unique(u(t >= on(k))), k) = true;
    end
    lat = median(first, 2, 'omitnan');
    part = mean(hit, 2);
    dt = 1e-3; ng = ceil(duration / dt); gi = min(floor(st / dt) + 1, ng);
    inM = gridMask([on; off], dt, ng);
    preM = gridMask([on - opts.PreWindow; on], dt, ng) & ~inM;
    postM = gridMask([off; off + opts.PostWindow], dt, ng) & ~inM & ~preM;
    r = @(m) accumarray(su(m(gi)), 1, [U 1]) / max(sum(m) * dt, eps) ./ max(rate, eps);
    pre = r(preM); inr = r(inM); post = r(postM);
    none = nsp == 0; pre(none) = NaN; inr(none) = NaN; post(none) = NaN;
end
T = table((1:U)', nsp, rate, lat, part, pre, inr, post, 'VariableNames', ...
    {'Unit', 'NSpikes', 'Rate', 'Latency', 'Participation', 'PreRatio', 'InRatio', 'PostRatio'});
end

function m = gridMask(iv, dt, ng)
% logical grid (ng bins of dt) covering the union of the intervals iv (2 x n)
d = zeros(ng + 1, 1);
a = max(floor(iv(1, :) / dt) + 1, 1); b = min(floor(iv(2, :) / dt) + 1, ng + 1);
v = b > a;
d = d + accumarray(a(v)', 1, [ng + 1 1]) - accumarray(b(v)', 1, [ng + 1 1]);
m = cumsum(d(1:ng)) > 0;
end
