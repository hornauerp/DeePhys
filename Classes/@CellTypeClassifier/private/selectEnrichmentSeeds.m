function det = selectEnrichmentSeeds(G, responsive_local, subset_global_idx, resp_unique, ...
        eligible_global, p_comm, rng_seed, strata_global)
% SELECTENRICHMENTSEEDS  Training seeds from communities statistically enriched / depleted in responders.
%
% Replaces the community heuristics (EnrichmentFactor, InhibitoryCommunityRelThresh,
% PuritySigmaThreshold, CommunityFallbackThreshold) and the iforest fallback with one test:
%   1. Louvain on the fixed graph, LouvainRestarts times.
%   2. Per community of each restart, a hypergeometric test of its responder count among
%      eligible units (units in cultures where responders can occur), in both directions,
%      stratified by culture (Community.StratifyEnrichment) so that a culture's responder
%      rate does not make its units' communities look enriched; Benjamini-Hochberg FDR at
%      Community.SeedFDR across communities.
%   3. Consensus over restarts (the partitions can be degenerate, i.e. restarts need not
%      agree): inhibitory seeds are responders in an enriched community in >= half of the
%      restarts; excitatory seeds come from non-responders in a depleted community in >= half,
%      as many as inhibitory seeds, allocated by community size (best partition) and taken
%      where the graph neighbourhood holds the fewest responders (deterministic; no embedding).
% Errors when no responder qualifies: the responder calls then carry no structure in the
% feature graph and any training set would be arbitrary.
%
% INPUTS:
%   G                - (N_all x N_all) symmetric graph (CellTypeClassifier.graphMatrix)
%   responsive_local - (1 x R) local (subset) indices of responsive units
%   subset_global_idx- (1 x N_subset) global unique-unit index per subset row
%   resp_unique      - (1 x N_all) logical, responsive per unique unit
%   eligible_global  - (1 x N_all) logical, unit belongs to a culture where responders occur
%   p_comm           - Parameters.Community (LouvainResolution, LouvainRestarts, SeedFDR)
%   rng_seed         - Parameters.RNGSeed
%   strata_global    - (1 x N_all) culture per unique unit: the strata of the enrichment test
%                      when p_comm.StratifyEnrichment is true ([] or omitted = one stratum)
    N_all = size(G, 1);
    fdr = 0.05;
    if isfield(p_comm, 'SeedFDR'), fdr = p_comm.SeedFDR; end
    if nargin < 8 || ~isfield(p_comm, 'StratifyEnrichment') || ~p_comm.StratifyEnrichment
        strata_global = [];
    end

    algorithm = "louvain";
    if isfield(p_comm, 'Algorithm'), algorithm = p_comm.Algorithm; end
    [M, Q, linfo] = louvainBest(G, p_comm.LouvainResolution, p_comm.LouvainRestarts, rng_seed, algorithm);
    n_comm = max(M);
    fprintf('%s: %d communities, Q=%.3f (resolution %.2f; restarts agree with the best partition: median ARI %.3f over %d)\n', ...
        linfo.algorithm, n_comm, Q, p_comm.LouvainResolution, linfo.ari_median, linfo.restarts);

    in_subset = false(1, N_all); in_subset(subset_global_idx) = true;
    elig = eligible_global(:)' & in_subset;
    [T, inh_comm, exc_comm] = communityEnrichmentTest(M, elig, resp_unique(:)', fdr, strata_global);
    tested = T.NEligible' > 0;
    fprintf(['Responder enrichment (%d responders among %d eligible units, FDR %.2f): ' ...
             '%d enriched, %d depleted, %d neutral communities\n'], sum(T.NResponsive), sum(T.NEligible), fdr, ...
        numel(inh_comm), numel(exc_comm), sum(tested) - numel(inh_comm) - numel(exc_comm));
    for c = inh_comm
        fprintf('   enriched community %d: %d units, %d responders (%.1f expected), q=%.2g\n', ...
            c, T.NEligible(c), T.NResponsive(c), T.Expected(c), T.qEnriched(c));
    end
    % -- consensus over all Louvain restarts --------------------------------------------
    % Partitions of these graphs can be degenerate (restarts need not agree), so seeds
    % taken from one partition change with small input changes. A unit is an inhibitory seed
    % when it is a responder in an enriched community in at least half of the restarts, and
    % joins the excitatory pool when it sits in a depleted community in at least half.
    parts = linfo.partitions; R = numel(parts);
    inh_vote = zeros(1, N_all); exc_vote = zeros(1, N_all);
    for r = 1:R
        Mr = parts{r};
        [~, ic, ec] = communityEnrichmentTest(Mr, elig, resp_unique(:)', fdr, strata_global);
        inh_vote = inh_vote + ismember(Mr', ic);
        exc_vote = exc_vote + ismember(Mr', ec);
    end
    inh_freq = inh_vote / R; exc_freq = exc_vote / R;

    % -- inhibitory seeds -----------------------------------------------------------------
    resp_global = subset_global_idx(responsive_local);
    in_inh = inh_freq(resp_global) >= 0.5;
    clean_global = resp_global(in_inh);
    if isempty(clean_global)
        error('CellTypeClassifier:noEnrichedCommunity', ...
            ['No responder lies in a significantly enriched community (FDR %.2f) in at least half of ' ...
             'the %d Louvain restarts. The responder calls are not reflected in the feature graph, so ' ...
             'no reliable training set exists. Check the responder calls, the features, or NNeighbors.'], fdr, R);
    end

    % -- excitatory seeds -----------------------------------------------------------------
    % Responder density of each unit's graph neighbourhood (weighted). Only responders of the
    % training subset count: in leave-one-culture-out evaluation the held-out culture's
    % responders must not decide which seeds are kept.
    deg = full(sum(G, 2))'; deg(deg == 0) = 1;
    dens = full(double(resp_unique(:)' & in_subset) * G) ./ deg;
    pool_mask = elig & ~resp_unique(:)' & exc_freq >= 0.5;
    if ~any(pool_mask)
        warning('CellTypeClassifier:noDepletedCommunity', ...
            ['No unit lies in a significantly depleted community in at least half of the restarts; ' ...
             'excitatory seeds are taken from non-responders outside enriched communities.']);
        pool_mask = elig & ~resp_unique(:)' & inh_freq < 0.5;
    end
    exc_comm = unique(M(pool_mask))';           % stratify the pool by the best partition
    n_target = numel(clean_global);
    pool_c = accumarray(M(pool_mask), 1, [n_comm 1])';
    alloc = zeros(1, n_comm);
    alloc(exc_comm) = pool_c(exc_comm) / max(sum(pool_c(exc_comm)), 1) * n_target;
    alloc = floor(alloc);
    short = n_target - sum(alloc);                 % distribute the rounding remainder
    [~, order] = sort(pool_c .* ismember(1:n_comm, exc_comm) - alloc, 'descend');
    for c = order(1:min(short, numel(order)))
        alloc(c) = alloc(c) + 1;
    end
    ce_global = zeros(1, 0);
    for c = exc_comm
        cand = find(pool_mask & M' == c);
        if isempty(cand) || alloc(c) == 0, continue; end
        [~, o] = sortrows([dens(cand)', cand']);   % lowest density first, ties by index
        ce_global = [ce_global, cand(o(1:min(alloc(c), numel(cand))))]; %#ok<AGROW>
    end

    % -- balance: the propagation compares per-seed absorption, i.e. assumes equal class
    % priors at the seeds. With a small excitatory pool (e.g. pure-culture labels after a rate
    % filter) the surplus inhibitory seeds would make each excitatory seed weigh more and bias
    % the labels towards excitatory. Keep as many inhibitory seeds as excitatory ones,
    % allocated by community and taken where the neighbourhood holds the most responders.
    n_inh_all = numel(clean_global);
    if ~isempty(ce_global) && numel(ce_global) < n_inh_all
        n_keep = numel(ce_global);
        Mi = M(clean_global)'; comms = unique(Mi);
        cnt = accumarray(Mi(:), 1, [n_comm 1])';
        alloc_i = floor(cnt / n_inh_all * n_keep);
        [~, order] = sort(cnt / n_inh_all * n_keep - alloc_i, 'descend');
        for c = order(1:(n_keep - sum(alloc_i))), alloc_i(c) = alloc_i(c) + 1; end
        keep = zeros(1, 0);
        for c = comms
            cand = clean_global(Mi == c);
            [~, o] = sortrows([-dens(cand)', cand']);   % highest responder density first
            keep = [keep, cand(o(1:min(alloc_i(c), numel(cand))))]; %#ok<AGROW>
        end
        clean_global = sort(keep);
    end
    fprintf(['Training seeds (consensus of %d restarts): %d inhibitory (of %d in enriched communities, ' ...
             '%d responders), %d excitatory from a pool of %d\n'], R, numel(clean_global), n_inh_all, ...
        numel(resp_global), numel(ce_global), sum(pool_mask));

    [~, clean_local] = ismember(clean_global, subset_global_idx);
    [~, ce_local]    = ismember(ce_global, subset_global_idx);
    det = struct( ...
        'responsive_local',             responsive_local, ...
        'clean_resp_local',             clean_local(clean_local > 0), ...
        'counterexample_idx_local',     ce_local(ce_local > 0), ...
        'tf_resp_outlier',              ~in_inh(:)', ...
        'inh_centroid',                 [], ...
        'resp_geom_mask',               logical([]), ...
        'resp_community_outlier_mask',  ~in_inh(:)', ...
        'resp_purity_outlier_mask',     false(1, sum(in_inh)), ...
        'ce_outlier_mask',              logical([]), ...
        'ce_distances_from_inh',        [], ...
        'community_ids',                M, ...
        'inh_comm_ids',                 inh_comm(:)', ...
        'Q_modularity',                 Q, ...
        'use_community_path',           true, ...
        'n_ce_pool',                    sum(pool_mask), ...
        'success',                      true, ...
        'enrichment',                   T, ...
        'louvain',                      linfo);
end
