function ax = PlotBurstCutouts(eia, culture_idx, opts)
% PLOTBURSTCUTOUTS  Plot mean E and I burst profiles and their time derivatives.
%
% Plots the accepted burst cutouts (after outlier rejection). If the data has
% been clustered into multiple groups (k>1), each cluster is plotted as a
% separate panel row.
%
% Top row:    normalised mean excitatory and inhibitory firing rates during bursts.
% Bottom row: time derivatives (gradient) of both traces.
%
% INPUTS:
%   eia          - EIAnalyzer object
%   culture_idx  - Index of the culture to plot (default 1)
%
% OPTIONS (name-value):
%   Parent       - tiledlayout or figure to plot into; creates a new figure if not provided
%
% OUTPUTS:
%   ax           - axes handle array (2 x K, where K is the number of clusters)

arguments
    eia         EIAnalyzer
    culture_idx (1,1) double = 1
    opts.Parent = []
end

assert(~isempty(eia.BurstCutouts), 'Run extractBurstCutouts() before plotting');
assert(culture_idx >= 1 && culture_idx <= numel(eia.BurstCutouts), ...
    'culture_idx %d out of range [1, %d]', culture_idx, numel(eia.BurstCutouts));

bc = eia.BurstCutouts(culture_idx);
if isempty(bc.accepted.total)
    warning('EIAnalyzer:PlotBurstCutouts', ...
        'No accepted burst cutouts found for culture %i — skipping plot', culture_idx);
    ax = gobjects(0);
    return
end

clusters  = bc.clusters;
n_clusters = numel(clusters);

if isempty(opts.Parent)
    f = figure('Color', 'w');
    tl = tiledlayout(f, 2, n_clusters, 'TileSpacing', 'compact', 'Padding', 'compact');
else
    tl = opts.Parent;
    f  = tl.Parent;
end

ax = gobjects(2, n_clusters);

for k = 1:n_clusters
    cl = clusters(k);
    if isempty(cl.total); continue; end

    mean_i = mean(cl.i_mat, 1, 'omitnan');
    mean_e = mean(cl.e_mat, 1, 'omitnan');
    max_i  = max(mean_i);  if max_i == 0; max_i = 1; end
    max_e  = max(mean_e);  if max_e == 0; max_e = 1; end

    cluster_label = '';
    if n_clusters > 1
        cluster_label = sprintf(' — cluster %d/%d (n=%d)', k, n_clusters, size(cl.total, 1));
    end

    % Top: normalised firing rates
    ax(1,k) = nexttile(tl, k);
    hold(ax(1,k), 'on');
    plot(ax(1,k), mean_e / max_e, 'b', 'LineWidth', 1);
    plot(ax(1,k), mean_i / max_i, 'r', 'LineWidth', 1);
    axis(ax(1,k), 'tight'); box(ax(1,k), 'off');
    legend(ax(1,k), ["Excitatory", "Inhibitory"], 'Box', 'off', 'Location', 'northwest');
    ylabel(ax(1,k), 'Normalised firing rate');
    title(ax(1,k), sprintf('Culture %i%s', culture_idx, cluster_label), 'FontWeight', 'normal');

    % Bottom: time derivatives
    ax(2,k) = nexttile(tl, k + n_clusters);
    hold(ax(2,k), 'on');
    plot(ax(2,k), gradient(mean_e / max_e, 2), 'b', 'LineWidth', 1);
    plot(ax(2,k), gradient(mean_i / max_i, 2), 'r', 'LineWidth', 1);
    axis(ax(2,k), 'tight'); box(ax(2,k), 'off');
    legend(ax(2,k), ["Exc \partial FR / \partial t", "Inh \partial FR / \partial t"], ...
        'Box', 'off', 'Location', 'northwest');
    ylabel(ax(2,k), '\partial FR / \partial t');
    xlabel(ax(2,k), 'Time (bins)');
end

setAllFontSizes(f, eia.Parameters.Plotting.FontSize);
end
