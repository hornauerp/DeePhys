function det = selectIforestCE(X_feat, responsive_local, subset_responsive, ...
        subset_global_idx, reduction, p_outlr)
% SELECTIFORESTCE  Isolation forest outlier detection + farthest-point CE selection.
%
% INPUTS:
%   X_feat            - (N_subset x F) normalized feature matrix for training subset
%   responsive_local  - (1 x R) local indices of responsive units in X_feat
%   subset_responsive - (N_subset x 1) logical
%   subset_global_idx - (1 x N_subset) global unique-unit index per subset row
%   reduction         - (N_all x D) full UMAP embedding
%   p_outlr           - ctc.Parameters.OutlierDetection

outlier_domain = lower(string(p_outlr.Domain));
contamination  = p_outlr.ContaminationFraction;
gmm_sep        = p_outlr.AutoGMMSeparation;

% -- Clean inhibitory candidates with iforest ---------------------------------
inh_input = prepareIforestInput(X_feat, responsive_local, ...
    reduction, subset_global_idx, outlier_domain, p_outlr);
[tf_resp_outlier, resp_info] = runOutlierDetection(inh_input, contamination, gmm_sep);
fprintf('Iforest (inh, %s): %d/%d flagged%s\n', outlier_domain, ...
    sum(tf_resp_outlier), numel(tf_resp_outlier), resp_info.detail);

clean_resp_local = responsive_local(~tf_resp_outlier(:)');
if isempty(clean_resp_local)
    clean_resp_local = responsive_local;
end

% -- Geometric consistency filter ---------------------------------------------
inh_centroid = mean(X_feat(clean_resp_local, :), 1);
exc_centroid = mean(X_feat(~subset_responsive, :), 1);
d_to_inh     = pdist2(inh_centroid, X_feat(clean_resp_local, :), 'correlation');
d_to_exc     = pdist2(exc_centroid, X_feat(clean_resp_local, :), 'correlation');
geom_ok      = d_to_inh(:) < d_to_exc(:);
resp_geom_mask   = ~geom_ok(:)';
clean_resp_local = clean_resp_local(geom_ok);
if isempty(clean_resp_local)
    clean_resp_local = responsive_local;
    resp_geom_mask   = false(1, numel(clean_resp_local));
    inh_centroid     = mean(X_feat(clean_resp_local, :), 1);
end

% -- Farthest-from-inhibitory CE selection ------------------------------------
non_responsive_local = find(~subset_responsive);
ce_ratio    = p_outlr.CounterexampleRatio;
n_pick      = min(ce_ratio * numel(clean_resp_local), numel(non_responsive_local));
pool_dists  = pdist2(inh_centroid, X_feat(non_responsive_local, :), 'correlation');
dist_threshold = prctile(pdist2(inh_centroid, X_feat(clean_resp_local, :), 'correlation'), ...
    p_outlr.CounterexampleDistancePercentile);
far_pool = non_responsive_local(pool_dists(:) >= dist_threshold);
if numel(far_pool) < n_pick; far_pool = non_responsive_local; end
n_ce_pool     = numel(far_pool);
perm_order    = randperm(numel(far_pool));
sampled_local = far_pool(perm_order(1:min(n_pick, numel(far_pool))));

% -- Clean CE candidates with iforest -----------------------------------------
ce_input = prepareIforestInput(X_feat, sampled_local, ...
    reduction, subset_global_idx, outlier_domain, p_outlr);
[tf_ce_outlier, ~] = runOutlierDetection(ce_input, contamination, gmm_sep);
counterexample_idx_local = sampled_local(~tf_ce_outlier(:)');
if isempty(counterexample_idx_local)
    counterexample_idx_local = sampled_local;
end

det = struct( ...
    'responsive_local',             responsive_local, ...
    'clean_resp_local',             clean_resp_local, ...
    'counterexample_idx_local',     counterexample_idx_local, ...
    'tf_resp_outlier',              tf_resp_outlier(:)', ...
    'inh_centroid',                 inh_centroid, ...
    'resp_geom_mask',               resp_geom_mask, ...
    'resp_community_outlier_mask',  logical([]), ...
    'resp_purity_outlier_mask',     logical([]), ...
    'ce_outlier_mask',              tf_ce_outlier(:)', ...
    'ce_distances_from_inh',        pool_dists, ...
    'community_ids',                [], ...
    'inh_comm_ids',                 [], ...
    'Q_modularity',                 NaN, ...
    'use_community_path',           false, ...
    'n_ce_pool',                    n_ce_pool, ...
    'success',                      true);
end
