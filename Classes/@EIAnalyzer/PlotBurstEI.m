function fig = PlotBurstEI(eia, r, opts)
% PLOTBURSTEI  E and I rates around the network bursts of recording r (index into BurstEI),
% aligned at onset, peak and offset, and I/E with the 95 % band of the rate-matched null.
%
%   fig = eia.PlotBurstEI(r)
%   fig = eia.PlotBurstEI(r, XLim=[-0.5 0.5])

arguments
    eia EIAnalyzer
    r   (1,1) double
    opts.XLim (1,2) double = [-0.6 0.9]
end
assert(~isempty(eia.BurstEI) && r <= numel(eia.BurstEI), 'EIAnalyzer:noBurstEI', 'Run computeBurstEI() first.');
R = eia.BurstEI(r).Result;
assert(isfield(R, 'Valid') && R.Valid, 'EIAnalyzer:invalid', 'Recording %d has too few bursts or units.', r);
fs = eia.Parameters.Plotting.FontSize;
fig = figure('Color', 'w', 'Position', [100 100 1100 520]);
try, theme(fig, 'light'); catch, end   % R2025a+ default theme can be dark
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact');
title(tl, sprintf('%s: %d bursts, %d I / %d E units', eia.BurstEI(r).Culture, R.NBursts, R.NI, R.NE), ...
    'FontSize', fs + 1, 'Interpreter', 'none');
al = ["Onset", "Peak", "Offset"];
for k = 1:3
    t = R.Traces.(al(k));
    ax = nexttile(tl, k); hold(ax, 'on');
    plot(ax, t.Time, t.E / max(t.E), 'Color', [0.2 0.4 0.9], 'LineWidth', 1.2);
    plot(ax, t.Time, t.I / max(t.I), 'Color', [0.9 0.2 0.2], 'LineWidth', 1.2);
    xline(ax, 0, 'k:'); xlim(ax, opts.XLim); ylabel(ax, 'rate / max'); title(ax, "aligned at " + lower(al(k)));
    if k == 1, legend(ax, ["E", "I"], 'Location', 'northwest', 'Box', 'off'); end
    set(ax, 'FontSize', fs);
    ax = nexttile(tl, 3 + k); hold(ax, 'on');
    fill(ax, [t.Time, fliplr(t.Time)], [t.IENullLow, fliplr(t.IENullHigh)], [0.85 0.85 0.85], 'EdgeColor', 'none');
    plot(ax, t.Time, t.IE, 'k', 'LineWidth', 1.2);
    xline(ax, 0, 'k:'); xlim(ax, opts.XLim); ylabel(ax, 'I / E'); xlabel(ax, "s from " + lower(al(k)));
    set(ax, 'FontSize', fs);
end
end
