function [reduction, n_neighbors] = runUnsupervisedUMAP(ctc, X_all, p_umap, reuse_existing)
% RUNUNSUPERVISEDMAP  Build the unit graph and the unsupervised UMAP embedding for all unique units.
%
% Sets ctc.Graph (Parameters.Graph.Method = "exact": deterministic kNN graph, see
% buildKNNGraph) and, when an embedding is computed, ctc.UMAP and ctc.Reduction.Unsupervised.
% Called by generateTrainLabels; also callable directly after changing Parameters.
%
% INPUTS:
%   ctc            - CellTypeClassifier handle
%   X_all          - (N x F) normalized feature matrix for all unique units
%   p_umap         - ctc.Parameters.UMAP struct
%   reuse_existing - (logical) if true, skip recompute when a valid graph/embedding exists
%
% OUTPUT:
%   reduction      - (N x NDims) embedding, or [] when UMAP.ComputeEmbedding is false

N_all = size(X_all, 1);
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_all)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_all);
else
    n_neighbors = p_umap.NNeighbors;
end
graph_method = "umap";
if isfield(ctc.Parameters, 'Graph') && isfield(ctc.Parameters.Graph, 'Method')
    graph_method = lower(string(ctc.Parameters.Graph.Method));
end
want_embedding = ~isfield(p_umap, 'ComputeEmbedding') || p_umap.ComputeEmbedding || graph_method == "umap";

umap_ready = ~isempty(ctc.UMAP) && ~isempty(ctc.Reduction) ...
             && isfield(ctc.Reduction, 'Unsupervised') ...
             && ~isempty(ctc.Reduction.Unsupervised) ...
             && size(ctc.Reduction.Unsupervised, 1) == N_all;

% Same group weighting as optimizeUnsupervisedUMAP, so a recomputed embedding (e.g. one per
% ensemble seed) matches the optimised one and the ACG bins do not dominate distances.
if ~isempty(ctc.NormalizationParams) && isfield(ctc.NormalizationParams, 'feat_names_trimmed')
    X_all = applyFeatureGroupWeights(X_all, ctc.NormalizationParams.feat_names_trimmed, p_umap);
end

% The exact graph is rebuilt every time (seconds; never stale). "Reuse" applies to the
% expensive UMAP embedding only.
if graph_method == "exact"
    t0 = tic;
    ctc.Graph = struct('W', buildKNNGraph(X_all, n_neighbors), 'k', n_neighbors, 'Method', "exact");
    fprintf('Unit graph: exact %d-NN fuzzy graph on %d units (%.1f s)\n', n_neighbors, N_all, toc(t0));
elseif graph_method == "wnn"
    t0 = tic;
    is_acg = acgColumns(ctc.NormalizationParams.feat_names_trimmed);
    [W, winfo] = buildWNNGraph(X_all, is_acg, n_neighbors);
    ctc.Graph = struct('W', W, 'k', n_neighbors, 'Method', "wnn", 'ModalityWeights', winfo.weights);
    fprintf(['Unit graph: weighted-nearest-neighbour %d-NN graph on %d units (%.1f s); ' ...
             'median modality weight ACG %.2f, waveform %.2f\n'], n_neighbors, N_all, toc(t0), ...
        median(winfo.weights(:, 1)), median(winfo.weights(:, 2)));
else
    ctc.Graph = struct();
end

if reuse_existing && umap_ready
    reduction = ctc.Reduction.Unsupervised;
    fprintf('runUnsupervisedUMAP: reusing existing embedding (N=%d). Pass reuse_existing=false to force rerun.\n', N_all);
    return
end

reduction = [];
if want_embedding
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
end
