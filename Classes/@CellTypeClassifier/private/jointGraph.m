function [G, reduction] = jointGraph(ctc, X, n_neighbors, feat_names)
% JOINTGRAPH  Unit graph (and optional embedding) for a combined DeePhys + external feature set.
%
% Uses the same construction as the main pipeline: the exact kNN fuzzy graph when
% Parameters.Graph.Method is "exact" (deterministic), otherwise the UMAP model's graph.
% The embedding is computed when UMAP.ComputeEmbedding is true (or required by "umap").
    p_umap = ctc.Parameters.UMAP;
    method = "umap";
    if isfield(ctc.Parameters, 'Graph') && isfield(ctc.Parameters.Graph, 'Method')
        method = lower(string(ctc.Parameters.Graph.Method));
    end
    want_embedding = ~isfield(p_umap, 'ComputeEmbedding') || p_umap.ComputeEmbedding || method == "umap";
    N = size(X, 1);
    reduction = [];
    if want_embedding
        rng(ctc.Parameters.RNGSeed, 'twister');
        [reduction, umap_model, ~, ~] = run_umap(X, ...
            'n_components', p_umap.NDims, 'n_neighbors', n_neighbors, ...
            'min_dist', p_umap.MinDist, 'spread', p_umap.Spread, ...
            'method', 'Java', 'verbose', 'none');
    end
    if method == "exact"
        G = buildKNNGraph(X, n_neighbors);
    elseif method == "wnn"
        G = buildWNNGraph(X, acgColumns(feat_names), n_neighbors);   % feat_names: columns of X
    else
        G = sparse(umap_model.head, umap_model.tail, umap_model.weights, N, N);
        G = (G + G') / 2;
    end
end
