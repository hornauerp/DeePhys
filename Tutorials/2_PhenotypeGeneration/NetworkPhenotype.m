% NetworkPhenotype.m
%
% Evaluate extracted features at the recording/network level to generate
% phenotypes. Includes unsupervised dimensionality reduction, supervised
% classification, feature importance, and feature trajectory plots.
%
% WHAT CHANGED (2026 refactor):
%   - Cross-validation now defaults to recording-level grouping (GroupedCV)
%     to prevent data leakage from the same recording.
%   - Random forest uses proper feature subsampling (NumVariablesToSample
%     = 'auto' instead of 'all') and regularisation (MinLeafSize=5).
%   - NaN imputation uses per-column median instead of zero-fill.
%   - UMAP defaults changed: min_dist=0.1, spread=1.0, n_neighbors=15.
%     These are now configurable via umap_min_dist/spread/n_neighbors args.
%   - StatisticalTest class auto-selects rank-sum vs t-test, uses
%     Lilliefors (not KS) for normality, and Nakagawa-Schielzeth R2 for
%     mixed-effects models.
%   - NormalizationPipeline prevents data leakage (fit on train, transform
%     test separately).
%
% OLD APPROACH (pre-refactor):
%   - cv_grouping defaulted to "none" (random unit-level splits).
%   - RF was actually bagging (all features per tree).
%   - Missing features were zero-imputed.
%   - UMAP used hardcoded min_dist=1, spread=5, n_neighbors=100.

%% 0 - Setup
addpath(genpath("/path/to/DeePhys")) % CHANGE

%% 1 - Load data
root_path = "/path/to/your/data/"; % CHANGE
path_logic = {'2406*', 'T002*', 'Network', 'w*', 'sorter_output'};
path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n", length(path_list))

% Load in parallel (removes CCG data by default to save memory)
rm_con = true;
mearec_array = recording_array_from_single_files(path_list, rm_con);

% Optional: remove recordings with too few units
unit_threshold = 1;
mearec_array = remove_low_unit_recordings(mearec_array, unit_threshold);

%% 2 - Create RecordingGroup
% Filter to a subset of recordings using metadata-based inclusion/exclusion.
rg_params.Selection.Inclusion = {{'DIV', 14}, {'RecordingDate', "240610"}};
rg_params.Selection.Exclusion = {};
pattern_rg = RecordingGroup(mearec_array, rg_params);

% RecordingGroup auto-detects cultures (grouped by ChipID + PlatingDate).
% Recordings within a culture are sorted by DIV if RecordingDate is set.

%% 3 - Unsupervised dimensionality reduction
dr_level = "Recording";
dr_method = "UMAP"; % "UMAP", "tSNE", "PCA"
n_dims = 2;

unit_features = ["all"];
network_features = ["Catch22", "Regularity", "GraphFeatures"];
feature_names = [];             % empty = use all from selected groups
useClustered = false;
normalization_var = ["PlatingDate"]; % z-score within each batch
grouping_var = [];              % set to "DIV" for longitudinal data
grouping_values = nan;          % nan = use all available
normalization = [];             % "baseline" or "scaled" for longitudinal
tolerance = 0;

pattern_rg.reduceDimensionality(dr_level, dr_method, n_dims, ...
    unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% 4 - Visualise reduction
cluster_var = ["Mutation"]; % colour by this metadata variable

nodeSz = 20;
mapSz = 500;
sigma = 10;
cmap = othercolor('Set16', 6); % adjust for number of conditions

figure("Color", "w");
[cluster_idx, group_labels_comb] = pattern_rg.plot_true_clusters( ...
    dr_level, dr_method, cluster_var, cmap=cmap, nodeSz=nodeSz, ...
    mapSz=mapSz, sigma=sigma);

%% 5 - Optional: pooling conditions
% Combine multiple conditions into groups without a dedicated metadata field.
% Three levels of nesting: {per_metadata_var{per_group{conditions_in_group}}}
pooling_vals = {{{"C1", "C4"}, {"C2", "C3", "C5", "C6"}}};

figure("Color", "w");
[cluster_idx, group_labels_comb] = pattern_rg.plot_true_clusters( ...
    dr_level, dr_method, cluster_var, cmap=cmap, nodeSz=nodeSz, ...
    mapSz=mapSz, sigma=sigma, pooling_vals=pooling_vals);

%% 6 - Cluster purity
clust_method = "spectral"; % "kmeans", "hierarchical", "spectral", "gmm", "louvain"
cluster_purity = pattern_rg.calculateClusterPurity(dr_method, clust_method, cluster_var);
fprintf("Clustering by '%s': purity = %.3f\n", cluster_var, cluster_purity)

%% 7 - Supervised classification
% Build a separate RecordingGroup for classification (different filters).
rg_params.Selection.Inclusion = {};
rg_params.Selection.Exclusion = {{'Mutation', "WT"}};
clf_rg = RecordingGroup(mearec_array, rg_params);

clf_level = "Recording";
alg = 'rf';                     % 'svm', 'knn', 'cnb', 'rf'
clf_var = "Mutation";
pooling_vals = {};

network_features_clf = ["Regularity", "GraphFeatures"];
unit_features_clf = ["WaveformFeatures", "ActivityFeatures", "RegularityFeatures"];
feature_names = [];
useClustered = false;
grouping_var = [];
grouping_values = nan;
normalization = [];
normalization_var = [];
N_hyper = 0;                    % hyperparameter search iterations (0 = skip)
K_fold = -1;                    % -1 = leave-one-out CV
tolerance = 0;

% NOTE: cv_grouping now defaults to "recording" (was "none" previously).
% This prevents units from the same recording leaking between train/test.
result = clf_rg.classifyByFeatureGroups(clf_level, alg, clf_var, ...
    pooling_vals, network_features_clf, unit_features_clf, feature_names, ...
    useClustered, grouping_var, grouping_values, normalization, ...
    normalization_var, N_hyper, K_fold, tolerance);

%% 8 - Evaluate classifier
[metrics_special, train_acc, avg_pred_imp, conf_mat] = clf_rg.assessClassifier(result);
fprintf("F1 score: %.3f\n", metrics_special.F1Score)

figure('Color', 'w');
confusionchart(conf_mat)

%% 9 - Feature importance
[sorted_pI, sort_idx] = sort(avg_pred_imp, 'descend');

figure('Color', 'w');
bar(sorted_pI);
xticks(1:length(avg_pred_imp));
xticklabels(result(1).Features(sort_idx));
box off; ylabel('Feature importance');
setAllFontSizes(gcf, 7)

%% 10 - Feature trajectories
% Useful for longitudinal data. For single-timepoint data, it plots
% one point per condition with error bars.
plt_feature_names = result(1).Features(sort_idx(1:5));

grouping_var = "DIV";
comp_var = "Mutation";
[feature_table, mean_mat, sd_mat, group_labels_comb] = ...
    clf_rg.plot_feature_trajectories(clf_level, grouping_var, ...
    grouping_values, network_features_clf, unit_features_clf, ...
    normalization, plt_feature_names, comp_var, pooling_vals);

%% 11 - Feature heatmap
% Compare two conditions across features. Best with longitudinal data.
pooling_vals = {{{"C1", "C4"}, {"C2", "C3", "C5", "C6"}}};
color_lim = 3;

[color_mat, features] = clf_rg.plot_feature_heatmap(clf_level, ...
    grouping_var, grouping_values, network_features_clf, ...
    unit_features_clf, feature_names, normalization, comp_var, ...
    pooling_vals, useClustered, tolerance, color_lim);
