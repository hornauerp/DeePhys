function W = buildKNNGraph(X, k)
% BUILDKNNGRAPH  Deterministic UMAP-style fuzzy kNN graph with exact nearest neighbours.
%
%   W = buildKNNGraph(X, k)
%
% Same construction as UMAP's fuzzy simplicial set (McInnes et al. 2018), but with an
% exact neighbour search instead of the randomised approximate one, so the graph (and
% everything computed on it: communities, label propagation) is identical across runs:
%   rho_i   = distance to the nearest neighbour at a non-zero distance
%   sigma_i : sum_j exp(-max(0, d_ij - rho_i) / sigma_i) = log2(k)   (binary search)
%   a_ij    = exp(-max(0, d_ij - rho_i) / sigma_i)
%   W       = A + A' - A .* A'                                         (fuzzy union)
%
% INPUTS:
%   X - (N x F) feature matrix (rows = units)
%   k - number of neighbours (excluding the unit itself)
% OUTPUT:
%   W - (N x N) sparse symmetric weight matrix with zero diagonal
    N = size(X, 1);
    k = min(k, N - 1);
    if k < 1
        W = sparse(N, N);
        return
    end
    [idx, d] = knnsearch(X, X, 'K', k + 1);
    % Drop each unit's own entry (normally the first; with duplicate rows it may not be)
    self = idx == (1:N)';
    no_self = ~any(self, 2);
    self(no_self, end) = true;            % no self match: drop the farthest neighbour instead
    idx = reshape(idx', [], 1); d = reshape(d', [], 1); keep = ~reshape(self', [], 1);
    idx = reshape(idx(keep), k, N)';
    d   = reshape(d(keep), k, N)';

    % rho: distance to the nearest neighbour at a non-zero distance
    dpos = d; dpos(dpos <= 0) = Inf;
    rho = min(dpos, [], 2); rho(isinf(rho)) = 0;

    % sigma: binary search per unit (vectorised), as UMAP's smooth_knn_dist
    target = log2(k);
    lo = zeros(N, 1); hi = inf(N, 1); sigma = ones(N, 1);
    dr = max(d - rho, 0);
    for it = 1:64
        psum = sum(exp(-dr ./ sigma), 2);
        too_big = psum > target;
        hi(too_big) = sigma(too_big);
        lo(~too_big) = sigma(~too_big);
        fin = isfinite(hi);
        sigma(fin) = (lo(fin) + hi(fin)) / 2;
        sigma(~fin) = sigma(~fin) * 2;
    end
    mean_d = mean(d, 2);
    sigma = max(sigma, 1e-3 * mean_d);   % UMAP's MIN_K_DIST_SCALE floor
    sigma(sigma <= 0) = eps;

    a = exp(-dr ./ sigma);
    A = sparse(repmat((1:N)', 1, k), idx, a, N, N);
    W = A + A' - A .* A';
    W = W - spdiags(diag(W), 0, N, N);
end
