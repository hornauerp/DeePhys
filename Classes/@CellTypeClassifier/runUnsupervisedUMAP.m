function [reduction, n_neighbors] = runUnsupervisedUMAP(ctc, X_all, p_umap, reuse_existing)
% RUNUNSUPERVISEDMAP  Compute or reuse the unsupervised UMAP embedding for all unique units.
%
% Sets ctc.UMAP and ctc.Reduction.Unsupervised when a new embedding is computed.
% Called by generateTrainLabels; also callable directly after changing Parameters.
%
% INPUTS:
%   ctc            - CellTypeClassifier handle
%   X_all          - (N x F) normalized feature matrix for all unique units
%   p_umap         - ctc.Parameters.UMAP struct
%   reuse_existing - (logical) if true, skip recompute when a valid embedding exists

N_all = size(X_all, 1);
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_all)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_all);
else
    n_neighbors = p_umap.NNeighbors;
end

umap_ready = ~isempty(ctc.UMAP) && ~isempty(ctc.Reduction) ...
             && isfield(ctc.Reduction, 'Unsupervised') ...
             && ~isempty(ctc.Reduction.Unsupervised) ...
             && size(ctc.Reduction.Unsupervised, 1) == N_all;

if reuse_existing && umap_ready
    reduction = ctc.Reduction.Unsupervised;
    fprintf('runUnsupervisedUMAP: reusing existing embedding (N=%d). Pass reuse_existing=false to force rerun.\n', N_all);
    return
end

[reduction, umap_model, ~, ~] = run_umap(X_all, ...
    'n_components',  p_umap.NDims, ...
    'n_neighbors',   n_neighbors, ...
    'min_dist',      p_umap.MinDist, ...
    'spread',        p_umap.Spread, ...
    'method', 'Java', 'verbose', 'none');
ctc.UMAP = umap_model;
if isempty(ctc.Reduction); ctc.Reduction = struct(); end
ctc.Reduction.Unsupervised = reduction;
fprintf('Unsupervised UMAP: %d units, %d features, %d dimensions\n', ...
    N_all, size(X_all, 2), p_umap.NDims);
end
