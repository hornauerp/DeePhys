function [M, Q, info] = louvainBest(G, resolution, restarts, rng_seed, algorithm)
% LOUVAINBEST  Best-modularity partition over several restarts (Louvain or Leiden).
%
%   [M, Q, info] = louvainBest(G, resolution, restarts, rng_seed)
%   [M, Q, info] = louvainBest(G, resolution, restarts, rng_seed, "leiden")
%
% The graph is fixed, so the restarts only differ in the algorithm's random node order;
% with enough of them the best-Q partition no longer depends on rng_seed.
% info.n_best: number of restarts that reached the best partition (a stability indicator).
%
% INPUTS:
%   G          - (N x N) symmetric weight matrix (sparse or full)
%   resolution - modularity resolution gamma
%   restarts   - number of runs
%   rng_seed   - base seed (run r uses rng_seed + r)
%   algorithm  - "louvain" (BCT community_louvain, full matrix; default) or "leiden"
%                (leidenCommunities: connected communities, sparse graph)
    if nargin < 5 || isempty(algorithm), algorithm = "louvain"; end
    algorithm = lower(string(algorithm));
    assert(ismember(algorithm, ["louvain", "leiden"]), 'CellTypeClassifier:communityAlgorithm', ...
        'Community.Algorithm must be "louvain" or "leiden", not "%s".', algorithm);
    N = size(G, 1);
    if algorithm == "louvain"
        W = full(G);   % BCT community_louvain requires a full matrix
    else
        W = sparse(G);
    end
    Q = -Inf; M = ones(N, 1);
    parts = cell(1, restarts); qs = -inf(1, restarts);
    for r = 1:restarts
        rng(rng_seed + r, 'twister');
        try
            if algorithm == "louvain"
                [M_r, Q_r] = community_louvain(W, resolution);
            else
                [M_r, Q_r] = leidenCommunities(W, resolution);
            end
            parts{r} = M_r(:); qs(r) = Q_r;
            if Q_r > Q, Q = Q_r; M = M_r(:); end
        catch ME
            warning('CellTypeClassifier:louvainFailed', '%s restart %d/%d failed: %s', ...
                algorithm, r, restarts, ME.message);
        end
    end
    % Agreement of every restart with the best partition (adjusted Rand index). Exact
    % identity is too strict for thousands of nodes (a few boundary units always move).
    ari = cellfun(@(p) adjRand(p, M), parts(~cellfun(@isempty, parts)));
    info = struct('n_comm', max(M), 'restarts', restarts, 'Q_all', qs, ...
        'ari_median', median(ari), 'n_best', sum(ari > 0.999), 'algorithm', algorithm);
    info.partitions = parts(~cellfun(@isempty, parts));   % all restarts (consensus use)
end

function a = adjRand(p, q)
% Adjusted Rand index between two partitions (label vectors)
    [~, ~, p] = unique(p(:)); [~, ~, q] = unique(q(:));
    n = numel(p);
    N = accumarray([p q], 1);
    c2 = @(x) x .* (x - 1) / 2;
    s_ij = sum(c2(N(:))); s_a = sum(c2(sum(N, 2))); s_b = sum(c2(sum(N, 1)));
    e = s_a * s_b / c2(n);
    a = (s_ij - e) / ((s_a + s_b) / 2 - e);
end
