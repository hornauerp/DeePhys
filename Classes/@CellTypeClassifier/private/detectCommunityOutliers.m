function det = detectCommunityOutliers(X_feat, responsive_local, subset_responsive, ...
        subset_global_idx, reduction, N_all, n_neighbors, resp_unique, ...
        p_comm, umap_model, rng_seed)
% DETECTCOMMUNITYOUTLIERS  Louvain community detection + graph purity filter + k-medoids CE.
%
% Returns det.success=false when community structure does not separate responsive
% units (triggers fallback to iforest in the orchestrator).
%
% INPUTS:
%   X_feat            - (N_subset x F) normalized feature matrix for training subset
%   responsive_local  - (1 x R) local indices of responsive units in X_feat
%   subset_responsive - (N_subset x 1) logical
%   subset_global_idx - (1 x N_subset) global unique-unit index per subset row
%   reduction         - (N_all x D) full UMAP embedding
%   N_all             - total number of unique units
%   n_neighbors       - resolved NNeighbors (used for graph k and purity filter)
%   resp_unique       - (1 x N_all) logical, true = responsive in full population
%   p_comm            - ctc.Parameters.Community
%   umap_model        - ctc.UMAP (must have .head/.tail/.weights fields)
%   rng_seed          - ctc.Parameters.RNGSeed (for Louvain restarts and k-medoids)

assert(~isempty(umap_model) && ~isempty(umap_model.head), ...
    ['ctc.UMAP graph data is empty. ' ...
     'Run optimizeUnsupervisedUMAP() or generateTrainLabels() first, ' ...
     'or set Parameters.OutlierDetection.Method = "iforest" to skip community detection.']);

G_umap  = sparse(umap_model.head, umap_model.tail, umap_model.weights, N_all, N_all);
G_umap  = (G_umap + G_umap') / 2;
k_graph = min(n_neighbors, N_all - 1);

% -- Louvain community detection (best of LouvainRestarts runs) ---------------
W_lou  = full(G_umap);
best_Q = -Inf;
best_M = ones(N_all, 1);
for r_lou = 1:p_comm.LouvainRestarts
    rng(rng_seed + r_lou, 'twister');
    try
        [M_r, Q_r] = community_louvain(W_lou, p_comm.LouvainResolution);
        if Q_r > best_Q; best_Q = Q_r; best_M = M_r(:); end
    catch ME
        warning('CellTypeClassifier:louvainFailed', ...
            'Louvain restart %d/%d failed: %s', r_lou, p_comm.LouvainRestarts, ME.message);
    end
end
community_ids = best_M;
Q_modularity  = best_Q;
n_comm        = max(best_M);
fprintf('Louvain: %d communities, Q=%.3f (resolution=%.2f, %d restarts)\n', ...
    n_comm, best_Q, p_comm.LouvainResolution, p_comm.LouvainRestarts);

% -- Identify inhibitory communities (enriched in responsive units) -----------
resp_global_idx_all = subset_global_idx(responsive_local);
resp_global_all     = false(1, N_all);
resp_global_all(resp_global_idx_all) = true;

comm_resp_frac = arrayfun(@(c) ...
    sum(resp_global_all(best_M == c)) / max(sum(best_M == c), 1), 1:n_comm);
max_frac   = max(comm_resp_frac);
p_resp_gen = sum(resp_unique) / N_all;
inh_by_rel = comm_resp_frac >= max_frac * p_comm.InhibitoryCommunityRelThresh;
inh_by_abs = comm_resp_frac >= p_comm.EnrichmentFactor * p_resp_gen;
inh_comm_ids = find(inh_by_rel & inh_by_abs);
rel_thr = max_frac * p_comm.InhibitoryCommunityRelThresh;
abs_thr = p_comm.EnrichmentFactor * p_resp_gen;
fprintf('Inhibitory communities: %d/%d (rel_thr=%.3f [%.2f×max=%.3f], abs_thr=%.3f [%.1f×p_resp=%.3f])\n', ...
    numel(inh_comm_ids), n_comm, rel_thr, p_comm.InhibitoryCommunityRelThresh, max_frac, ...
    abs_thr, p_comm.EnrichmentFactor, p_resp_gen);

in_inh_comm_global = ismember(best_M', inh_comm_ids);
frac_resp_in_inh   = sum(resp_global_all & in_inh_comm_global) / max(sum(resp_global_all), 1);
fprintf('Responsive units in inhibitory communities: %.1f%%\n', 100*frac_resp_in_inh);

% -- Fallback check -----------------------------------------------------------
if frac_resp_in_inh < p_comm.CommunityFallbackThreshold
    warning('CellTypeClassifier:communityFallback', ...
        ['Only %.0f%% of responsive units in inhibitory communities ' ...
         '(threshold %.0f%%).\nFalling back to iforest path.\n' ...
         'Consider re-running optimizeUnsupervisedUMAP(), adjusting ' ...
         'Parameters.Community settings, or setting ' ...
         'Parameters.OutlierDetection.Method = "iforest".'], ...
        100*frac_resp_in_inh, 100*p_comm.CommunityFallbackThreshold);
    det = emptyDet(responsive_local);
    det.community_ids = community_ids;
    det.inh_comm_ids  = inh_comm_ids(:)';
    det.Q_modularity  = Q_modularity;
    det.success       = false;
    return
end

% -- Community outlier removal ------------------------------------------------
in_inh_resp = ismember(best_M(resp_global_idx_all)', inh_comm_ids);
resp_community_outlier_mask = ~in_inh_resp;
resp_in_inh_global = resp_global_idx_all(in_inh_resp);
fprintf('Community outlier removal: %d/%d responsive units excluded\n', ...
    sum(resp_community_outlier_mask), numel(resp_global_idx_all));

[resp_clean_global, resp_purity_outlier_mask] = graphPurityOutlier( ...
    reduction, resp_in_inh_global, resp_global_all, k_graph, ...
    p_comm.PuritySigmaThreshold);
fprintf('Graph purity outlier removal: %d/%d remaining responsive excluded\n', ...
    sum(resp_purity_outlier_mask), numel(resp_in_inh_global));

if isempty(resp_clean_global)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'All responsive candidates removed by outlier detection — using community filter only.');
    resp_clean_global        = resp_in_inh_global;
    resp_purity_outlier_mask = false(1, numel(resp_in_inh_global));
end

[~, loc] = ismember(resp_clean_global, subset_global_idx);
clean_resp_local = loc(loc > 0);

tf_resp_outlier = resp_community_outlier_mask(:)';
pass_comm_idx   = find(~resp_community_outlier_mask);
tf_resp_outlier(pass_comm_idx(resp_purity_outlier_mask)) = true;

inh_centroid = mean(X_feat(clean_resp_local, :), 1);

% -- k-medoids CE selection from non-inhibitory communities ------------------
non_resp_global  = subset_global_idx(~subset_responsive);
in_ce_comm       = ~ismember(best_M(non_resp_global)', inh_comm_ids);
ce_pool_global   = non_resp_global(in_ce_comm);
ce_pool_comm_ids = best_M(ce_pool_global)';
ce_communities_u = unique(ce_pool_comm_ids);
n_ce_pool        = numel(ce_pool_global);
fprintf('CE pool: %d non-responsive units in %d communities\n', ...
    n_ce_pool, numel(ce_communities_u));

n_clean_resp = numel(clean_resp_local);
k_alloc      = zeros(1, numel(ce_communities_u));
for ci = 1:numel(ce_communities_u)
    k_alloc(ci) = max(p_comm.MinCEPerCommunity, ...
        round(n_clean_resp * sum(ce_pool_comm_ids == ce_communities_u(ci)) / n_ce_pool));
end
k_diff = sum(k_alloc) - n_clean_resp;
if k_diff ~= 0
    [~, largest] = max(k_alloc);
    k_alloc(largest) = max(p_comm.MinCEPerCommunity, k_alloc(largest) - k_diff);
end

counterexample_global = zeros(1, 0);
for ci = 1:numel(ce_communities_u)
    c        = ce_communities_u(ci);
    c_mask   = ce_pool_comm_ids == c;
    c_global = ce_pool_global(c_mask);
    k_c      = min(k_alloc(ci), sum(c_mask));
    if k_c < 1; continue; end
    umap_c = reduction(c_global, :);
    rng(rng_seed, 'twister');
    [~, ~, ~, ~, midx] = kmedoids(umap_c, k_c, 'Distance', 'sqeuclidean', 'Replicates', 3);
    counterexample_global = [counterexample_global, c_global(midx(:))]; %#ok<AGROW>
end

[~, ce_loc] = ismember(counterexample_global, subset_global_idx);
counterexample_idx_local = ce_loc(ce_loc > 0);
fprintf('CE selection: %d medoids (target: %d = n_clean_resp)\n', ...
    numel(counterexample_idx_local), n_clean_resp);

det = struct( ...
    'responsive_local',             responsive_local, ...
    'clean_resp_local',             clean_resp_local, ...
    'counterexample_idx_local',     counterexample_idx_local, ...
    'tf_resp_outlier',              tf_resp_outlier, ...
    'inh_centroid',                 inh_centroid, ...
    'resp_geom_mask',               logical([]), ...
    'resp_community_outlier_mask',  resp_community_outlier_mask(:)', ...
    'resp_purity_outlier_mask',     resp_purity_outlier_mask(:)', ...
    'ce_outlier_mask',              logical([]), ...
    'ce_distances_from_inh',        [], ...
    'community_ids',                community_ids, ...
    'inh_comm_ids',                 inh_comm_ids(:)', ...
    'Q_modularity',                 Q_modularity, ...
    'use_community_path',           true, ...
    'n_ce_pool',                    n_ce_pool, ...
    'success',                      true);
end

function det = emptyDet(responsive_local)
    det = struct( ...
        'responsive_local',             responsive_local, ...
        'clean_resp_local',             responsive_local, ...
        'counterexample_idx_local',     zeros(1, 0), ...
        'tf_resp_outlier',              false(1, numel(responsive_local)), ...
        'inh_centroid',                 [], ...
        'resp_geom_mask',               logical([]), ...
        'resp_community_outlier_mask',  logical([]), ...
        'resp_purity_outlier_mask',     logical([]), ...
        'ce_outlier_mask',              logical([]), ...
        'ce_distances_from_inh',        [], ...
        'community_ids',                [], ...
        'inh_comm_ids',                 [], ...
        'Q_modularity',                 NaN, ...
        'use_community_path',           false, ...
        'n_ce_pool',                    0, ...
        'success',                      false);
end
