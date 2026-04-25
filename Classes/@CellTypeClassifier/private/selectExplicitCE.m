function det = selectExplicitCE(X_feat, responsive_local, subset_responsive, ...
        subset_global_idx, ce_local_all, reduction, p_outlr)
% SELECTEXPLICITCE  Explicit counterexample path (metadata labels for both classes).
%
% Cleans inhibitory candidates (responsive units) with iforest + geometric filter.
% Cleans excitatory CE units with iforest.
%
% INPUTS:
%   X_feat            - (N_subset x F) normalized feature matrix for training subset
%   responsive_local  - (1 x R) local indices of responsive units in X_feat
%   subset_responsive - (N_subset x 1) logical, true = responsive
%   subset_global_idx - (1 x N_subset) global unique-unit index per subset row
%   ce_local_all      - (1 x C) local indices of explicit CE units in X_feat
%   reduction         - (N_all x D) full UMAP embedding
%   p_outlr           - ctc.Parameters.OutlierDetection

outlier_domain = lower(string(p_outlr.Domain));
contamination  = p_outlr.ContaminationFraction;
gmm_sep        = p_outlr.AutoGMMSeparation;

% -- Clean inhibitory candidates --
inh_input = prepareIforestInput(X_feat, responsive_local, ...
    reduction, subset_global_idx, outlier_domain, p_outlr);
[tf_resp_outlier, resp_info] = runOutlierDetection(inh_input, contamination, gmm_sep);
fprintf('Isolation forest (inh, %s, %s): %d/%d flagged (%.0f%%)%s\n', ...
    outlier_domain, resp_info.mode, sum(tf_resp_outlier), ...
    numel(tf_resp_outlier), sum(tf_resp_outlier)/numel(tf_resp_outlier)*100, ...
    resp_info.detail);

clean_resp_local = responsive_local(~tf_resp_outlier(:)');
if isempty(clean_resp_local)
    clean_resp_local = responsive_local;
    tf_resp_outlier  = false(size(tf_resp_outlier));
end

% -- Geometric consistency filter --
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

% -- Clean explicit CE units --
ce_input = prepareIforestInput(X_feat, ce_local_all, ...
    reduction, subset_global_idx, outlier_domain, p_outlr);
[tf_ce_outlier, ~] = runOutlierDetection(ce_input, contamination, gmm_sep);
counterexample_idx_local = ce_local_all(~tf_ce_outlier(:)');
if isempty(counterexample_idx_local)
    counterexample_idx_local = ce_local_all;
end

ce_distances_from_inh = [];
if ~isempty(inh_centroid) && ~isempty(ce_local_all)
    ce_distances_from_inh = pdist2(inh_centroid, X_feat(ce_local_all, :), 'correlation');
end
fprintf('Metadata counterexamples: %d labeled, %d after outlier cleaning\n', ...
    numel(ce_local_all), numel(counterexample_idx_local));

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
    'ce_distances_from_inh',        ce_distances_from_inh, ...
    'community_ids',                [], ...
    'inh_comm_ids',                 [], ...
    'Q_modularity',                 NaN, ...
    'use_community_path',           false, ...
    'n_ce_pool',                    numel(ce_local_all), ...
    'success',                      true);
end
