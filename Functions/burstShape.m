function S = burstShape(spike_times, T_start, T_end, duration, opts)
% BURSTSHAPE  Burst timing and shape from the population rate, the same way for every burst.
%
%   S = burstShape(spike_times, T_start, T_end, duration)
%   S = burstShape(..., Level=0.2)
%
% A detector says which spikes form a burst (T_start / T_end, e.g. first and last spike of an
% ISI_N run). Where a burst starts and ends in time, and its shape, are measured here on the
% smoothed population rate, at the same point of every burst's own profile:
%   rate      all given spikes, Bin s bins, Gaussian smoothing (SD Sigma s)
%   baseline  median rate outside all bursts (+-0.1 s)
%   peak      rate maximum within [T_start - Pad, T_end + Pad]
%   onset / offset  crossings of baseline + Level x (peak - baseline), walking out from the
%             first / last sub-event (linear interpolation): sub-events are the stretches above
%             that level that contain the peak or overlap [T_start, T_end] and reach at least
%             twice the level, so dips inside a burst do not cut it while noise flicker around
%             the level on a ramp is ignored; a span that misses part of the ramp still gets the
%             whole ramp. Searched only up to the rate minimum in the gap between neighbouring
%             spans (between the peaks if spans overlap), so bursts never overlap and the silent
%             dips between the sub-bursts of a nested burst do not cut it.
%   FWHM      width of the peak at baseline + 0.5 x (peak - baseline) (local to the peak)
% Because the levels are relative to each burst's own amplitude, onset and offset sit at the
% same point of the ramp for small and large bursts, and hardly depend on where the detector put
% T_start / T_end. Edges based on spikes (first / last spike) or absolute rate thresholds mix
% bursts measured from the baseline with bursts measured from part-way up the ramp.
%
% INPUTS:
%   spike_times - (S x 1) spike times (s) of the units to use (e.g. all QC units)
%   T_start, T_end - (1 x B) burst spans from the detector
%   duration    - recording duration (s)
% OPTIONS:
%   Level (0.25: 20-30 % is robust to background noise yet keeps most of the ramp),
%   Bin (0.002 s), Sigma (0.005 s), Pad (0.05 s)
% OUTPUT: table with Peak (time, s), PeakRate (spikes/s), Baseline, Onset, Offset,
%   Duration (Offset - Onset), Rise (Peak - Onset), Fall (Offset - Peak), FWHM
arguments
    spike_times (:,1) double
    T_start (1,:) double
    T_end (1,:) double
    duration (1,1) double {mustBePositive}
    opts.Level (1,1) double {mustBeInRange(opts.Level, 0, 1, 'exclusive')} = 0.25
    opts.Bin (1,1) double = 0.002
    opts.Sigma (1,1) double = 0.005
    opts.Pad (1,1) double = 0.05
end
nb = numel(T_start);
names = {'Peak', 'PeakRate', 'Baseline', 'Onset', 'Offset', 'Duration', 'Rise', 'Fall', 'FWHM'};
if nb == 0
    S = array2table(zeros(0, numel(names)), 'VariableNames', names);
    return
end
[T_start, o] = sort(T_start); T_end = T_end(o);
eb = 0:opts.Bin:duration + opts.Bin;
r = histcounts(spike_times, eb) / opts.Bin;
k = -ceil(4 * opts.Sigma / opts.Bin):ceil(4 * opts.Sigma / opts.Bin);
g = exp(-0.5 * (k * opts.Bin / opts.Sigma).^2);
r = conv(r, g / sum(g), 'same');
c = eb(1:end-1) + opts.Bin / 2;
inb = false(size(r));
for i = 1:nb, inb(c >= T_start(i) - 0.1 & c <= T_end(i) + 0.1) = true; end
base = median(r(~inb));
if isempty(base) || isnan(base), base = prctile(r, 20); end
pk = zeros(1, nb);
for i = 1:nb
    w = find(c >= T_start(i) - opts.Pad & c <= T_end(i) + opts.Pad);
    if isempty(w), [~, w] = min(abs(c - (T_start(i) + T_end(i)) / 2)); end
    [~, m] = max(r(w)); pk(i) = w(m);
end
% boundary between neighbouring bursts: the rate minimum in the gap between their spans (the
% first minimum after a peak can be a dip between sub-bursts of a nested burst); if the spans
% overlap, the minimum between the two peaks
bnd = zeros(1, max(nb - 1, 0));
for i = 1:nb - 1
    w = find(c >= T_end(i) & c <= T_start(i+1));
    if isempty(w), w = min(pk(i), pk(i+1)):max(pk(i), pk(i+1)); end
    [~, m] = min(r(w)); bnd(i) = w(m);
end
[on, off, fl, fr] = deal(nan(1, nb));
for i = 1:nb
    lo = 1;  if i > 1,  lo = min(bnd(i-1), pk(i)); end
    hi = numel(r); if i < nb, hi = max(bnd(i), pk(i)); end
    amp = r(pk(i)) - base;
    if amp <= 0, continue; end
    % sub-events: stretches above the level within [lo, hi] that contain the peak or overlap the
    % span and rise to at least twice the level (noise flicker around the level on a ramp does
    % not); onset / offset by walking out from the maxima of the first / last sub-event
    lev = base + opts.Level * amp; lev2 = base + min(1, 2 * opts.Level) * amp;
    seg = lo:hi; e = diff([0, r(seg) > lev, 0]);
    ra = seg(find(e == 1)); rb = seg(find(e == -1) - 1);
    mx = zeros(size(ra)); mi = zeros(size(ra));
    for q = 1:numel(ra), [mx(q), m] = max(r(ra(q):rb(q))); mi(q) = ra(q) + m - 1; end
    keep = ((ra <= pk(i) & rb >= pk(i)) | (c(rb) >= T_start(i) & c(ra) <= T_end(i))) & mx >= lev2;
    first = min([mi(keep), pk(i)]); last = max([mi(keep), pk(i)]);
    on(i)  = crossing(r, c, first, lo, lev, -1);
    off(i) = crossing(r, c, last, hi, lev, +1);
    fl(i)  = crossing(r, c, pk(i), lo, base + 0.5 * amp, -1);
    fr(i)  = crossing(r, c, pk(i), hi, base + 0.5 * amp, +1);
end
tp = c(pk);
S = table(tp(:), r(pk)', repmat(base, nb, 1), on(:), off(:), off(:) - on(:), tp(:) - on(:), ...
    off(:) - tp(:), fr(:) - fl(:), 'VariableNames', names);
S(o, :) = S;                                    % back to the input order
end

function t = crossing(r, c, from, limit, level, dir)
% time where r falls to <= level walking from 'from' towards 'limit' (dir -1 left, +1 right),
% linearly interpolated between bins; the limit itself if it never does
j = from;
while j ~= limit && r(j) > level, j = j + dir; end
if r(j) > level, t = c(j); return; end
j2 = j - dir;                                   % last bin above the level
f = (r(j2) - level) / max(r(j2) - r(j), eps);
t = c(j2) + dir * f * (c(2) - c(1));
end
