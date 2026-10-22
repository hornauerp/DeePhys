function key = embeddingKey(X, n_neighbors, min_dist, spread, n_dims)
% EMBEDDINGKEY  Fingerprint of the input of an unsupervised UMAP embedding: size and row /
%   column sums of the (weighted) feature matrix plus the UMAP settings. An embedding is
%   reused only when the key matches, so changed features or settings recompute it.
    key = paramHash(struct('Size', size(X), 'RowSum', sum(X, 2)', 'ColSum', sum(X, 1), ...
        'NNeighbors', n_neighbors, 'MinDist', min_dist, 'Spread', spread, 'NDims', n_dims));
end
