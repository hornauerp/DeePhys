function [M, Q] = leidenCommunities(A, gamma, opts)
% LEIDENCOMMUNITIES  Leiden community detection for modularity (Traag, Waltman & van Eck 2019).
%
%   [M, Q] = leidenCommunities(A)
%   [M, Q] = leidenCommunities(A, gamma, Randomness=0.01, MaxPasses=10, Initial=M0)
%
% Optimises the same quality as BCT community_louvain(A, gamma):
%   Q = 1/(2m) * sum_ij (A_ij - gamma k_i k_j / (2m)) delta(c_i, c_j).
% Unlike Louvain, Leiden refines each community into well-connected subcommunities
% before aggregating, so every community is guaranteed to be connected (Louvain can
% return arbitrarily badly connected communities), and it repeats passes from the
% previous partition until the quality no longer improves. It works on a sparse graph
% (Louvain in BCT needs a full N x N matrix).
%
% INPUTS:
%   A     - (n x n) symmetric non-negative weight matrix (sparse or full); the diagonal
%           is ignored
%   gamma - resolution (default 1)
%   Randomness - theta of the refinement step (default 0.01, as in leidenalg)
%   MaxPasses  - maximum number of Leiden passes (default 10)
%   Initial    - (n x 1) starting partition (default: singletons)
%
% OUTPUTS:
%   M - (n x 1) community labels 1..K (largest community first)
%   Q - modularity of M
%
% Uses the current random number generator (set rng for reproducible results).

    arguments
        A
        gamma (1,1) double = 1
        opts.Randomness (1,1) double {mustBePositive} = 0.01
        opts.MaxPasses (1,1) double {mustBePositive} = 10
        opts.Initial = []
    end
    A = sparse(double(A));
    n = size(A, 1);
    A(1:n+1:end) = 0;
    if ~isequal(A, A')
        warning('leidenCommunities:asymmetric', 'A is not symmetric; using (A + A'') / 2.');
        A = (A + A') / 2;
    end
    assert(all(nonzeros(A) >= 0), 'leidenCommunities:negative', 'Weights must be non-negative.');
    k  = full(sum(A, 2));
    m2 = sum(k);
    if m2 == 0
        M = (1:n)'; Q = 0;
        return
    end

    if isempty(opts.Initial)
        M = (1:n)';
    else
        [~, ~, M] = unique(opts.Initial(:));
    end
    Q = modularityOf(A, k, m2, gamma, M);
    for pass = 1:opts.MaxPasses
        M_new = leidenPass(A, k, m2, gamma, M, opts.Randomness);
        Q_new = modularityOf(A, k, m2, gamma, M_new);
        if Q_new <= Q + 1e-12 && pass > 1
            break
        end
        if Q_new >= Q
            M = M_new; Q = Q_new;
        end
    end
    % labels 1..K, largest community first
    [~, ~, M] = unique(M);
    [~, order] = sort(accumarray(M, 1), 'descend');
    new_label = zeros(numel(order), 1);
    new_label(order) = 1:numel(order);
    M = new_label(M);
end


function M = leidenPass(A, k, m2, gamma, M0, theta)
% One Leiden pass: local moving, refinement, aggregation, until no aggregation is left.
    G = A; kk = k;
    node_map = (1:size(A, 1))';      % original node -> node of the current graph
    [~, ~, P] = unique(M0);
    while true
        P = moveNodesFast(G, kk, m2, gamma, P);
        n_nodes = size(G, 1);
        if max(P) == n_nodes          % every node its own community
            break
        end
        R = refinePartition(G, kk, m2, gamma, P, theta);
        n_ref = max(R);
        if n_ref == n_nodes           % nothing merged: aggregation would not shrink the graph
            break
        end
        S  = sparse(1:n_nodes, R, 1, n_nodes, n_ref);     % node -> refined community
        G  = S' * G * S;                                   % aggregate (internal weight on diagonal)
        kk = S' * kk;
        P_agg = zeros(n_ref, 1);
        P_agg(R) = P;                 % each refined community lies inside one community of P
        [~, ~, P] = unique(P_agg);
        node_map = R(node_map);
    end
    M = P(node_map);
end


function P = moveNodesFast(G, kk, m2, gamma, P)
% Queue-based local moving: a node is revisited only when a neighbour changed community.
    n = size(G, 1);
    [nbr, wts] = adjacency(G);
    K  = accumarray(P, kk, [n 1]);    % community totals (labels <= n)
    sz = accumarray(P, 1, [n 1]);
    free = find(sz == 0);             % stack of empty labels
    queue = randperm(n)'; head = 1; tail = n; in_queue = true(n, 1);
    count = n;
    while count > 0
        v = queue(head); head = mod(head, n) + 1; count = count - 1; in_queue(v) = false;
        nb = nbr{v};
        if isempty(nb), continue; end
        w = wts{v}; kv = kk(v); c_old = P(v);
        [uc, ~, ic] = unique(P(nb));
        wc = accumarray(ic, w);
        K(c_old) = K(c_old) - kv;
        gains = wc - gamma * kv * K(uc) / m2;
        own = find(uc == c_old, 1);
        if isempty(own)
            g_own = -gamma * kv * K(c_old) / m2;
        else
            g_own = gains(own);
        end
        [g_best, ib] = max(gains);
        c_new = c_old;
        if g_best > max(g_own, 0) + 1e-12
            c_new = uc(ib);
        elseif g_own < -1e-12 && sz(c_old) > 1 && ~isempty(free)   % alone is better
            c_new = free(end); free(end) = [];
        end
        K(c_new) = K(c_new) + kv;
        if c_new ~= c_old
            P(v) = c_new;
            sz(c_old) = sz(c_old) - 1; sz(c_new) = sz(c_new) + 1;
            if sz(c_old) == 0, free(end + 1) = c_old; end %#ok<AGROW>
            add = nb(P(nb) ~= c_new & ~in_queue(nb));
            for a = add(:)'
                tail = mod(tail, n) + 1; queue(tail) = a; count = count + 1;
            end
            in_queue(add) = true;
        end
    end
    [~, ~, P] = unique(P);
end


function R = refinePartition(G, kk, m2, gamma, P, theta)
% Split each community of P into well-connected subcommunities: start from singletons
% and merge well-connected singletons into well-connected subcommunities of the same
% community, choosing randomly with probability ~ exp(gain / theta) among those that do
% not decrease the quality.
    n = size(G, 1);
    [nbr, wts] = adjacency(G);
    [i, j, w] = find(G);
    off = i ~= j;
    same = off & P(i) == P(j);
    wS = accumarray(i(same), w(same), [n 1]);          % weight from v to the rest of its community
    KS = accumarray(P, kk);                               % community totals
    R  = (1:n)';
    KR = kk;                                              % refined community totals
    szR = ones(n, 1);
    ext = wS;                                             % E(C, S - C) of refined communities
    well = wS >= gamma * kk .* (KS(P) - kk) / m2;
    for v = randperm(n)
        if ~well(v) || szR(R(v)) ~= 1, continue; end
        s = P(v);
        nb = nbr{v}; wv = wts{v};
        in_s = P(nb) == s;
        if ~any(in_s), continue; end
        [uc, ~, ic] = unique(R(nb(in_s)));
        wc = accumarray(ic, wv(in_s));
        gains = wc - gamma * kk(v) * KR(uc) / m2;
        ok = gains >= 0 & ext(uc) >= gamma * KR(uc) .* (KS(s) - KR(uc)) / m2;
        if ~any(ok), continue; end
        g = gains(ok); cand = uc(ok); wcand = wc(ok);
        p = exp((g - max(g)) / theta);
        pick = find(rand * sum(p) <= cumsum(p), 1);
        c = cand(pick);
        old = R(v);
        R(v) = c;
        KR(c) = KR(c) + kk(v); KR(old) = 0;
        szR(c) = szR(c) + 1;   szR(old) = 0;
        ext(c) = ext(c) + wS(v) - 2 * wcand(pick);
    end
    [~, ~, R] = unique(R);
end


function [nbr, wts] = adjacency(G)
% Neighbour lists (self-loops excluded) of a symmetric sparse matrix.
    n = size(G, 1);
    [i, j, w] = find(G);
    keep = i ~= j;
    i = i(keep); j = j(keep); w = w(keep);
    counts = accumarray(j, 1, [n 1]);      % find returns entries sorted by column
    nbr = mat2cell(i, counts, 1);
    wts = mat2cell(w, counts, 1);
end


function Q = modularityOf(A, k, m2, gamma, M)
    [~, ~, M] = unique(M);
    [i, j, w] = find(A);
    same = M(i) == M(j);
    in_c = accumarray(M(i(same)), w(same), [max(M) 1]);
    K_c  = accumarray(M, k, [max(M) 1]);
    Q = sum(in_c / m2 - gamma * (K_c / m2) .^ 2);
end
