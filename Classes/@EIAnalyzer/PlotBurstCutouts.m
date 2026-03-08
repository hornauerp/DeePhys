function PlotBurstCutouts(eia, culture_idx)
% PLOTBURSTCUTOUTS  Plot mean E and I burst profiles and their time derivatives.
%
% Top panel: normalised mean excitatory and inhibitory firing rates during bursts.
% Bottom panel: time derivatives (gradient) of both traces.
%
% INPUTS:
%   eia          - EIAnalyzer object
%   culture_idx  - Index of the culture to plot (default 1)

arguments
    eia         EIAnalyzer
    culture_idx (1,1) double = 1
end

assert(~isempty(eia.BurstCutouts), 'Run extractBurstCutouts() before plotting');

bc = eia.BurstCutouts(culture_idx);
if isempty(bc.total)
    warning('No burst cutouts found for culture %i — skipping plot', culture_idx);
    return
end

mean_i = mean(bc.i_mat, 1, 'omitnan');
mean_e = mean(bc.e_mat, 1, 'omitnan');
max_i  = max(mean_i);  if max_i == 0; max_i = 1; end
max_e  = max(mean_e);  if max_e == 0; max_e = 1; end

figure('Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
plot(mean_e / max_e, 'b', 'LineWidth', 1);
plot(mean_i / max_i, 'r', 'LineWidth', 1);
axis tight; box off;
legend(["Excitatory", "Inhibitory"], 'Box', 'off', 'Location', 'northwest');
ylabel('Normalised firing rate');
title(sprintf('Culture %i — mean burst cutout', culture_idx), 'FontWeight', 'normal');

nexttile;
h = 2;
hold on;
plot(gradient(mean_e / max_e, h), 'b', 'LineWidth', 1);
plot(gradient(mean_i / max_i, h), 'r', 'LineWidth', 1);
axis tight; box off;
legend(["Exc \partialFR/\partialt", "Inh \partialFR/\partialt"], 'Box', 'off', 'Location', 'northwest');
ylabel('\partialFR/\partialt');
xlabel('Time (bins)');
end
