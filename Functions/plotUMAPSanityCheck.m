function fig = plotUMAPSanityCheck(ctc)
% PLOTUMAPSANITYCHECK  Unsupervised embedding coloured by training seeds and by final labels.
%
%   fig = plotUMAPSanityCheck(ctc)
%
% Tile 1: drug-responsive candidates (orange), inhibitory seeds (red) and excitatory
%         seeds (blue) on the unsupervised embedding of all unique units.
% Tile 2: the same embedding coloured by the predicted labels (NaN = grey).
%
% The embedding is for inspection only; labels come from the unit graph. Call after
% ctc.classifyUnits() with Parameters.UMAP.ComputeEmbedding = true (the default).

assert(~isempty(ctc.UnitLabels), 'Run classifyUnits() before plotUMAPSanityCheck().');
assert(isstruct(ctc.Reduction) && isfield(ctc.Reduction, 'Unsupervised') && ...
    ~isempty(ctc.Reduction.Unsupervised), ...
    'No embedding: set Parameters.UMAP.ComputeEmbedding = true and rerun generateTrainLabels().');

r   = ctc.Reduction.Unsupervised;           % (n_unique x NDims)
rep = ctc.NormalizedFeatures.unique_to_rep; % representative row of each unique unit
labels = ctc.UnitLabels(rep);
resp   = ctc.ResponsiveUnitIdx(rep);
tl_ids = ctc.TrainLabels.sorted_train_ids;
tl_y   = ctc.TrainLabels.sorted_y_train;

exc_color = [0.22, 0.49, 0.87];
inh_color = [0.87, 0.28, 0.20];
gray      = [0.70, 0.70, 0.70];

fig = figure;   % follows the MATLAB theme (light or dark)
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Embedding sanity check', 'FontWeight', 'normal');

nexttile(tl); hold on
scatter(r(:, 1), r(:, 2), 5, gray, 'filled', 'MarkerFaceAlpha', 0.25, ...
    'DisplayName', sprintf('All units (n=%d)', size(r, 1)));
scatter(r(resp, 1), r(resp, 2), 10, [0.95, 0.62, 0.10], 'filled', 'MarkerFaceAlpha', 0.7, ...
    'DisplayName', sprintf('Responders (n=%d)', sum(resp)));
seed_exc = tl_ids(tl_y == 1); seed_inh = tl_ids(tl_y == 2);
scatter(r(seed_exc, 1), r(seed_exc, 2), 14, exc_color, 'filled', ...
    'DisplayName', sprintf('Excitatory seeds (n=%d)', numel(seed_exc)));
scatter(r(seed_inh, 1), r(seed_inh, 2), 14, inh_color, 'filled', ...
    'DisplayName', sprintf('Inhibitory seeds (n=%d)', numel(seed_inh)));
xlabel('UMAP 1'); ylabel('UMAP 2'); title('Training seeds', 'FontWeight', 'normal');
legend('Location', 'best', 'FontSize', 7, 'Box', 'off'); box off; hold off

nexttile(tl); hold on
m_unc = isnan(labels);
if any(m_unc)
    scatter(r(m_unc, 1), r(m_unc, 2), 5, gray, 'filled', 'MarkerFaceAlpha', 0.3, ...
        'DisplayName', sprintf('Unlabelled (n=%d)', sum(m_unc)));
end
scatter(r(labels == 1, 1), r(labels == 1, 2), 6, exc_color, 'filled', 'MarkerFaceAlpha', 0.5, ...
    'DisplayName', sprintf('Excitatory (n=%d)', sum(labels == 1)));
scatter(r(labels == 2, 1), r(labels == 2, 2), 6, inh_color, 'filled', 'MarkerFaceAlpha', 0.6, ...
    'DisplayName', sprintf('Inhibitory (n=%d)', sum(labels == 2)));
xlabel('UMAP 1'); ylabel('UMAP 2'); title('Predicted labels', 'FontWeight', 'normal');
legend('Location', 'best', 'FontSize', 7, 'Box', 'off'); box off; hold off
end
