function labels = assembleTrainLabels(det, subset_global_idx, nf, p_train, ...
        unit_ids_table, responsive_strength)
% ASSEMBLETRAINLABELS  Map local detection results to global unique-unit indices.
%
% Builds the labels struct consumed by classifyUnits.
%
% INPUTS:
%   det               - detection result struct from selectExplicit/Iforest/Community/NoCE
%   subset_global_idx - (1 x N_subset) global unique-unit index per subset row
%   nf                - NormalizedFeatures struct (n_unique, all_to_unique, n_units_full, unique_ud)
%   p_train           - ctc.Parameters.TrainLabels
%   unit_ids_table    - (N_full x 1) string UnitID from FeatureStore.UnitTable
%   responsive_strength - (1 x N_full) continuous response strength, or []

resp_label    = p_train.ResponsiveClassLabel;
counter_label = 3 - resp_label;

in_train_id = subset_global_idx(det.clean_resp_local);
ex_train_id = subset_global_idx(det.counterexample_idx_local);

train_ids = [ex_train_id, in_train_id];
y_train   = [counter_label * ones(1, numel(ex_train_id)), ...
             resp_label    * ones(1, numel(in_train_id))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids = sorted_train_ids;
labels.sorted_y_train   = y_train(sort_idx);
labels.umap_train_idx   = false(1, nf.n_unique);
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx    = ~labels.umap_train_idx;

outlier_local = det.responsive_local(det.tf_resp_outlier(:)');
labels.outlier_global_idx              = subset_global_idx(outlier_local);
labels.excitatory_candidate_global_idx = subset_global_idx(det.counterexample_idx_local);

labels.n_unique      = nf.n_unique;
labels.all_to_unique = nf.all_to_unique;
labels.n_units_full  = nf.n_units_full;
labels.has_explicit_ce = det.has_explicit_ce;

% Diagnostic fields
labels.resp_outlier_mask           = det.tf_resp_outlier(:)';
labels.resp_geom_mask              = det.resp_geom_mask;
labels.resp_community_outlier_mask = det.resp_community_outlier_mask;
labels.resp_purity_outlier_mask    = det.resp_purity_outlier_mask;
labels.community_ids               = det.community_ids;
labels.inh_comm_ids                = det.inh_comm_ids;
labels.Q_modularity                = det.Q_modularity;
labels.use_community_path          = det.use_community_path;
labels.ce_outlier_mask             = det.ce_outlier_mask;
labels.ce_distances_from_inh_centroid = det.ce_distances_from_inh;
labels.inh_centroid                = det.inh_centroid;

labels.filter_counts = struct( ...
    'n_responsive_in',    numel(det.responsive_local), ...
    'n_after_outlier',    numel(det.responsive_local) - sum(det.tf_resp_outlier), ...
    'n_after_geometric',  numel(det.clean_resp_local), ...
    'n_ce_in',            det.n_ce_pool, ...
    'n_ce_after_outlier', numel(det.counterexample_idx_local));

% Inhibitory strength for training units
if ~isempty(responsive_strength)
    in_train_ud_ids = string({nf.unique_ud(in_train_id).UnitID});
    [~, fs_locs]    = ismember(in_train_ud_ids, unit_ids_table);
    labels.inhibitory_strength = responsive_strength(fs_locs);
else
    labels.inhibitory_strength = ones(1, numel(in_train_id));
end
end
