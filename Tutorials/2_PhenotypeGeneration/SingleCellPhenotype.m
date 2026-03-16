% SingleCellPhenotype.m
%
% Use extracted features at the single-cell (unit) level to find clusters
% corresponding to different cell types, then classify and compare them.
%
% This extends the WaveMAP approach by incorporating activity, regularity,
% and Catch22 features alongside waveform shape.
%
% WHAT CHANGED (2026 refactor):
%   - New features: CV2, LV, LvR, FanoFactor, HalfWidth enrich the
%     single-cell feature space.
%   - UMAP defaults: min_dist=0.1, spread=1.0, n_neighbors=15 (was 1/5/100).
%   - Louvain clustering uses the BCT community_louvain implementation.
%   - Random forest: proper feature subsampling + MinLeafSize=5 regularisation.
%   - cv_grouping defaults to "recording" to prevent cross-recording leakage.
%   - For dedicated cell-type classification using pharmacological data,
%     see the CellTypeClassification tutorial (uses @CellTypeClassifier).

%% 0 - Setup
addpath(genpath("/path/to/DeePhys")) % CHANGE

%% 1 - Load data
root_path = "/path/to/your/data/"; % CHANGE
path_logic = {'2406*', 'T002*', 'Network', 'w*', 'sorter_output'};
path_list = generate_sorting_path_list(root_path, path_logic);

rm_con = true;
mearec_array = recording_array_from_single_files(path_list, rm_con);
mearec_array = remove_low_unit_recordings(mearec_array, 1);

%% 2 - Create RecordingGroup
rg_params.Selection.Inclusion = {{'DIV', 14}, {'RecordingDate', "240610"}};
rg_params.Selection.Exclusion = {};
pattern_rg = RecordingGroup(mearec_array, rg_params);

%% 3 - Unit-level dimensionality reduction
dr_level = "Unit";
dr_method = "UMAP";
n_dims = 2;

unit_features = ["all"];
feature_names = [];
network_features = [];          % not relevant for unit-level
useClustered = false;
normalization_var = ["PlatingDate"];
grouping_var = [];
grouping_values = [];
normalization = [];
tolerance = 1;

pattern_rg.reduceDimensionality(dr_level, dr_method, n_dims, ...
    unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% 4 - Visualise by metadata
grouping_var_plot = ["Mutation"];
figure("Color", "w");
[cluster_idx, group_labels_comb] = pattern_rg.plot_true_clusters(dr_level, dr_method, grouping_var_plot);

% Split across subplots (cleaner for many units)
pattern_rg.plot_dimensionality_reduction( ...
    pattern_rg.DimensionalityReduction.Unit.UMAP.Reduction, ...
    cluster_idx, group_labels_comb);

%% 5 - Generate clusters
% Louvain typically works best for single-cell clustering. Resolution
% parameter controls granularity (smaller = fewer clusters).
clust_method = "louvain";
clust_param = 1;                % resolution parameter
pattern_rg.clusterByFeatures(dr_method, dr_level, clust_method, clust_param);

n_clust = max(pattern_rg.Clustering.Unit.(clust_method).Index);
fprintf('Found %i clusters\n', n_clust)

%% 6 - Visualise clusters
reduction = pattern_rg.DimensionalityReduction.Unit.(dr_method).Reduction;
cluster_idx = pattern_rg.Clustering.Unit.(clust_method).Index;

plot_centroid = false;
nodeSz = 5;
mapSz = 500;
sigma = 10;
cmap = othercolor('Set16', n_clust);

figure("Color", "w")
ax = gca;
[reduction, cmap] = RecordingGroup.plot_cluster_outlines(reduction, ...
    cluster_idx, ax, plot_centroid, nodeSz, mapSz, sigma, cmap);

%% 7 - Cluster densities
pooling_vals = {};
n_bins = 50;
smoothing_factor = 0.5;

[value_array, group_labels_comb] = pattern_rg.plot_cluster_densities( ...
    dr_level, dr_method, cluster_idx, pooling_vals, n_bins, smoothing_factor);

%% 8 - Cluster distribution changes across conditions
metadata_var = "Mutation";
[value_array, group_labels_comb] = pattern_rg.plot_cluster_densities( ...
    dr_level, dr_method, metadata_var, pooling_vals, n_bins, smoothing_factor);

figure("Color", "w")
pattern_rg.plot_cluster_shifts(value_array, group_labels_comb)

%% 9 - Feature heatmap by cluster
feature_names_heatmap = [];     % empty = all features
pattern_rg.plot_unit_cluster_heatmap(cluster_idx, feature_names_heatmap);

% Box plots for selected features
feature_names_box = ["FiringRate", "CVInterSpikeInterval", "CV2InterSpikeInterval", "HalfWidth"];
pattern_rg.plot_unit_cluster_features(cluster_idx, feature_names_box, cmap);

%% 10 - Cluster proportions across conditions
metadata_var = "Mutation";
colors = othercolor('Spectral9', n_clust);
figure('Color', 'w');
pattern_ratios = pattern_rg.plot_cluster_proportions(metadata_var, {}, clust_method, colors);

%% 11 - WaveMAP (waveform-only)
% To run the original WaveMAP approach using only waveform shape:
wavemap_rg = RecordingGroup(mearec_array, rg_params);

wavemap_rg.reduceDimensionality("Unit", "UMAP", 2, ...
    "ReferenceWaveform", [], [], false, ["PlatingDate"], [], [], [], 0);
wavemap_rg.clusterByFeatures("UMAP", "Unit", "louvain", 1);
wavemap_rg.plot_cluster_waveforms("louvain");

%% 12 - Supervised classification of clusters
cluster_idx = pattern_rg.Clustering.Unit.louvain.Index;
feature_groups = "all";
N_hyper = 0;
K_fold = 4;                     % dont use LOO for 1000s of units

result = pattern_rg.classifyClusteredUnits(cluster_idx, feature_groups, ...
    [], [], [], [], [], N_hyper, K_fold, 0);

[metrics_special, train_acc, avg_pred_imp, conf_mat] = ...
    pattern_rg.assessClassifier(result);

figure('Color', 'w');
confusionchart(conf_mat)

%% 13 - Feature importance
[sorted_pI, sort_idx] = sort(avg_pred_imp, 'descend');
figure('Color', 'w');
bar(sorted_pI);
xticks(1:length(avg_pred_imp));
xticklabels(result(1).Features(sort_idx));
box off; ylabel('Feature importance');
setAllFontSizes(gcf, 7)

%% 14 - Apply classifier to unseen data
% Train on one condition, predict another.
pooling_vals = {};
metadata_var = "Mutation";
[group_idx, group_labels] = pattern_rg.combineMetadataIndices("Unit", metadata_var, pooling_vals);

test_idx = group_idx == 3;      % pick one condition as test
train_idx = ~test_idx;
train_cluster_idx = cluster_idx(train_idx);
clf = [];                       % retrain from scratch (or pass a trained clf)

result = pattern_rg.applyClusteredUnitClassifier(train_cluster_idx, ...
    train_idx, test_idx, feature_groups, clf, [], [], [], [], [], N_hyper, 0);

figure("Color", "w")
confusionchart(cluster_idx(test_idx), result.Y_pred);

%% 15 - Removing a cluster (USE WITH CARE)
% THIS REMOVES UNITS FROM MEMORY. Reload data to get them back.
%
% calc_feat = true;  % recalculate clustered features
% concatenated = false;
% wavemap_rg.assignUnitClusterIdx(clust_method, calc_feat, concatenated);
% wavemap_rg.removeUnitsByCluster(clust_method, 2, concatenated);
%
% % Update spike times and features after removal
% for r = 1:length(wavemap_rg.Recordings)
%     wavemap_rg.Recordings(r).updateSpikeTimes();
%     wavemap_rg.Recordings(r).aggregateSingleCellFeatures();
% end
