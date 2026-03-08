function PlotEITraces(eia, culture_idx, stack_cutout)
% PLOTEITRACES  Stacked plot of E/I firing rates and total network activity.
%
% Top panel: normalised inhibitory and excitatory firing rate traces.
% Bottom panel: total network firing rate in kHz.
%
% INPUTS:
%   eia          - EIAnalyzer object
%   culture_idx  - Index of the culture to plot (default 1)
%   stack_cutout - (optional) bin index range to plot; default = all bins

arguments
    eia          EIAnalyzer
    culture_idx  (1,1) double = 1
    stack_cutout (1,:) double = []
end

assert(~isempty(eia.Activity), 'Run computeActivity() before plotting');

act = eia.Activity(culture_idx);
p   = eia.Parameters.Activity;

if isempty(stack_cutout)
    stack_cutout = 1:length(act.total);
end

max_I = max(act.I);  if max_I == 0; max_I = 1; end
max_E = max(act.E);  if max_E == 0; max_E = 1; end

f = figure('Color', 'w');
tiledlayout(2, 1, 'Padding', 'tight', 'TileSpacing', 'tight');

ax(1) = nexttile;
hold on;
plot(act.I(stack_cutout) / max_I, 'r', 'LineWidth', 0.8);
plot(act.E(stack_cutout) / max_E, 'b', 'LineWidth', 0.8);
ax(1).XAxis.Visible = 'off';
box off; axis tight;
leg = legend(["Inhibitory", "Excitatory"], ...
    'Box', 'off', 'Location', 'northwest');
leg.Title.String    = 'Putative neuron type';
leg.ItemTokenSize   = [5, 5];
ylabel('Normalised firing rate');

ax(2) = nexttile;
plot((act.total(stack_cutout) / p.BinSize) / 1000, 'k', 'LineWidth', 0.8);
box off; axis tight;
ylabel('Network firing rate (kHz)');
xlabel('Time (s)');
t_bins = 0 : 100 : length(stack_cutout);
xticks(t_bins);
xticklabels(t_bins * p.BinSize);

setAllFontSizes(f, 7);
end
