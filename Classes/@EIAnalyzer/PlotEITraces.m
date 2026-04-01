function ax = PlotEITraces(eia, culture_idx, opts)
% PLOTEITRACES  Stacked plot of E/I firing rates and total network activity.
%
% Top panel:    normalised inhibitory and excitatory firing rate traces.
% Bottom panel: total network firing rate in kHz.
%
% INPUTS:
%   eia          - EIAnalyzer object
%   culture_idx  - Index of the culture to plot (default 1)
%
% OPTIONS (name-value):
%   StackCutout  - bin index range to plot; default = all bins
%   Parent       - tiledlayout or figure to plot into; creates a new figure if not provided
%
% OUTPUTS:
%   ax           - (1 x 2) axes handle array [top, bottom]

arguments
    eia          EIAnalyzer
    culture_idx  (1,1) double = 1
    opts.StackCutout (1,:) double = []
    opts.Parent  = []
end

assert(~isempty(eia.Activity), 'Run computeActivity() before plotting');

act = eia.Activity(culture_idx);
p   = eia.Parameters.Activity;

stack_cutout = opts.StackCutout;
if isempty(stack_cutout)
    stack_cutout = 1:numel(act.total);
end

max_I = max(act.I);  if max_I == 0; max_I = 1; end
max_E = max(act.E);  if max_E == 0; max_E = 1; end

if isempty(opts.Parent)
    f  = figure('Color', 'w');
    tl = tiledlayout(f, 2, 1, 'Padding', 'tight', 'TileSpacing', 'tight');
else
    tl = opts.Parent;
    f  = tl.Parent;
end

ax(1) = nexttile(tl);
hold(ax(1), 'on');
plot(ax(1), act.I(stack_cutout) / max_I, 'r', 'LineWidth', 0.8);
plot(ax(1), act.E(stack_cutout) / max_E, 'b', 'LineWidth', 0.8);
ax(1).XAxis.Visible = 'off';
box(ax(1), 'off'); axis(ax(1), 'tight');
leg = legend(ax(1), ["Inhibitory", "Excitatory"], 'Box', 'off', 'Location', 'northwest');
leg.Title.String  = 'Putative neuron type';
leg.ItemTokenSize = [5, 5];
ylabel(ax(1), 'Normalised firing rate');

ax(2) = nexttile(tl);
plot(ax(2), (act.total(stack_cutout) / p.BinSize) / 1000, 'k', 'LineWidth', 0.8);
box(ax(2), 'off'); axis(ax(2), 'tight');
ylabel(ax(2), 'Network firing rate (kHz)');
xlabel(ax(2), 'Time (s)');
t_bins = 0 : 100 : numel(stack_cutout);
xticks(ax(2), t_bins);
xticklabels(ax(2), t_bins * p.BinSize);

setAllFontSizes(f, eia.Parameters.Plotting.FontSize);
end
