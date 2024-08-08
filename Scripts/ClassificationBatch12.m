addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'DeePhysS*','19*','*','w*','sorted'};
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data to save RAM
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Source','FCDI'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',12}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
batch12group = RecordingGroup(rec_array, rg_params);

%% Find most discriminative features
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_test_val = []; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [7:7:35]; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = "PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = -1; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;

save_path = "/home/phornauer/Git/DeePhys/Data/Figure3";

%% Single cell features
clf_network_features = [];
clf_unit_features = "all";
feature_names = batch12group.returnFeatureNames("UnitFeatures");
sc_train_acc = nan(size(feature_names{1}));
sc_test_acc = sc_train_acc;
sc_score = sc_train_acc;
sc_pred_imp = nan(length(feature_names{1}),length(clf_grouping_values));

%% Main loop for single cell features
for fn = 1:length(feature_names{1})
    clf_feature_names = feature_names{1}(fn);
    clf_result = batch12group.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_feature_names, clf_useClustered,...
        clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
    
    [sc_train_acc(fn), sc_test_acc(fn), sc_score(fn), sc_pred_imp(fn,:)] = batch12group.assessClassifier(clf_result);
%     save(fullfile(save_path,'single_cell_features.mat'),'sc_train_acc', 'sc_test_acc', 'sc_score', 'sc_pred_imp','feature_names')
    fprintf('Finished feature %i/%i\n', fn, length(feature_names{1}))
end

%% Run MLM for single cell feature
[sc_mlm_result_array, sc_p_G, sc_p_GxT, sc_features] = runMLM(batch12group, clf_network_features, clf_unit_features, feature_names{1}, clf_useClustered, clf_grouping_var, clf_classification_var);
save(fullfile(save_path, 'single_cell_p_vals.mat'), 'sc_mlm_result_array','sc_p_G','sc_p_GxT','sc_features')

%% Timelines single cell
pooling_vals = {};
[sc_feature_table, sc_mean_mat, sc_sd_mat, sc_group_labels_comb] = batch12group.plot_feature_trajectories(clf_level, clf_grouping_var, clf_grouping_values, clf_network_features, ...
    clf_unit_features, clf_normalization, feature_names{1}, clf_classification_var, pooling_vals, clf_useClustered, clf_tolerance);
save(fullfile(save_path, 'single_cell_values.mat'),'sc_feature_table', 'sc_mean_mat', 'sc_sd_mat', 'sc_group_labels_comb')

%% Heatmap single cell
color_lim = 3;
[sc_color_mat, sc_features] = batch12group.plot_feature_heatmap(clf_level, clf_grouping_var, clf_grouping_values, clf_network_features, clf_unit_features, clf_normalization, feature_names{1}, clf_classification_var, clf_classification_test_val,...
    clf_useClustered, clf_tolerance, color_lim);
save(fullfile(save_path, 'single_cell_heatmap.mat'),'sc_color_mat', 'sc_features')

%% Network features
clf_network_features = "all";
clf_unit_features = [];

feature_names = batch12group.returnFeatureNames("NetworkFeatures");
rm_id = contains(feature_names{1},["Motif","Rich"]);
feature_names{1}(rm_id) = [];
nw_train_acc = nan(size(feature_names{1}));
nw_test_acc = nw_train_acc;
nw_score = nw_train_acc;
nw_pred_imp = cell(size(feature_names{1}));

%% Main loop for network features
for fn = 1:length(feature_names{1})
    clf_feature_names = feature_names{1}(fn);
    clf_result = batch12group.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_feature_names, clf_useClustered,...
        clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
    
    [nw_train_acc(fn), nw_test_acc(fn), nw_score(fn), nw_pred_imp(fn,:)] = batch12group.assessClassifier(clf_result);
%     save(fullfile(save_path,'network_features.mat'),'nw_train_acc', 'nw_test_acc', 'nw_score', 'nw_pred_imp','feature_names')
    fprintf('Finished feature %i/%i\n', fn, length(feature_names{1}))
end

%% Run MLM for single cell feature
[nw_mlm_result_array, nw_p_G, nw_p_GxT, nw_features] = runMLM(batch12group, clf_network_features, clf_unit_features, feature_names{1}, clf_useClustered, clf_grouping_var, clf_classification_var);
save(fullfile(save_path, 'network_p_vals.mat'), 'nw_mlm_result_array','nw_p_G','nw_p_GxT','nw_features')

%% Timelines network
pooling_vals = {};
[nw_feature_table, nw_mean_mat, nw_sd_mat, nw_group_labels_comb] = batch12group.plot_feature_trajectories(clf_level, clf_grouping_var, clf_grouping_values, clf_network_features, ...
    clf_unit_features, clf_normalization, feature_names{1}, clf_classification_var, pooling_vals, clf_useClustered, clf_tolerance);
save(fullfile(save_path, 'network_values.mat'),'nw_feature_table', 'nw_mean_mat', 'nw_sd_mat', 'nw_group_labels_comb')

%% Heatmap network
color_lim = 3;
[nw_color_mat, nw_features] = batch12group.plot_feature_heatmap(clf_level, clf_grouping_var, clf_grouping_values, clf_network_features, clf_unit_features, feature_names{1}, clf_normalization, clf_classification_var, clf_classification_test_val,...
    clf_useClustered, clf_tolerance, color_lim);
save(fullfile(save_path, 'network_heatmap.mat'),'nw_color_mat', 'nw_features')

%% Check recordings unsupervised
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "RecordingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [14:7:35]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

plot_var = ["Mutation"];
clust_method = "kmeans";

%% Single cell features
feature_names = sc_features_sorted;
sc_reduction = batch12group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

figure("Color","w");
[sc_true_idx, sc_group_labels_comb] = batch12group.plot_true_clusters(dr_level, dr_method, plot_var);
sc_assigned_idx = batch12group.clusterByFeatures(dr_method,dr_level,clust_method);
sc_clust_coef = sum(sc_true_idx == sc_assigned_idx) / length(sc_true_idx);
sc_clust_coef = max([sc_clust_coef, 1 - sc_clust_coef]);
save(fullfile(save_path, 'single_cell_clustering.mat'),'sc_reduction', 'sc_true_idx', 'sc_clust_coef', 'sc_group_labels_comb')

%% Network features
feature_names = nw_features_sorted;
nw_reduction = batch12group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

figure("Color","w");
[nw_true_idx, nw_group_labels_comb] = batch12group.plot_true_clusters(dr_level, dr_method, plot_var);
nw_assigned_idx = batch12group.clusterByFeatures(dr_method,dr_level,clust_method);
nw_clust_coef = sum(nw_true_idx == nw_assigned_idx) / length(nw_true_idx);
nw_clust_coef = max([nw_clust_coef, 1 - nw_clust_coef]);
save(fullfile(save_path, 'network_clustering.mat'),'nw_reduction', 'nw_true_idx', 'nw_clust_coef', 'nw_group_labels_comb')

%% Combined features
feature_names = [sc_features_sorted, nw_features_sorted];
comb_reduction = batch12group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

figure("Color","w");
[comb_true_idx, comb_group_labels_comb] = batch12group.plot_true_clusters(dr_level, dr_method, plot_var);
comb_assigned_idx = batch12group.clusterByFeatures(dr_method,dr_level,clust_method);
comb_clust_coef = sum(comb_true_idx == comb_assigned_idx) / length(comb_true_idx);
comb_clust_coef = max([comb_clust_coef, 1 - comb_clust_coef]);
save(fullfile(save_path, 'comb_clustering.mat'),'comb_reduction', 'comb_true_idx', 'comb_clust_coef', 'comb_group_labels_comb')

%% Feature groups x grouping values
pooling_vals = {};
grouping_value_list = [6:7:34,nan];
[train_acc, test_acc, avg_score] = batch12group.classifyByFeatureGroupsAndGroupingVar(clf_level, clf_alg, clf_classification_var, pooling_vals, clf_useClustered,...
                clf_grouping_var, grouping_value_list, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
save(fullfile(save_path, 'feature_group_heatmap.mat'), 'train_acc', 'test_acc', 'avg_score', 'feature_groups')
            
%% Age regression
level = "Recording"; %Unit or Recording
alg = "rf"; %rf, svm, cnb, knn
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset 
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on untreated, test on treated)
pooling_vals = {};
network_features = ["all"];
unit_features = ["all"];
useClustered = false;
normalization_var = "PlatingDate";
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = -1; % number of K-fold CV

%% Age regression WT
rg_params.Selection.Inclusion = {{'Mutation','WT'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',12},{'PlatingDate',200121}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_rg = RecordingGroup(rec_array, rg_params);

wt_age_result = wt_rg.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);

%% Age regression A53T
rg_params.Selection.Inclusion = {{'Mutation','A53T'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',12},{'PlatingDate',200121}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
a53t_rg = RecordingGroup(rec_array, rg_params);

a53t_age_result = a53t_rg.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);

%% Prepare for saving
wt_age = rmfield(wt_age_result, 'Mdl');
wt_age = rmfield(wt_age, 'objects');

a53t_age = rmfield(a53t_age_result, 'Mdl');
a53t_age = rmfield(a53t_age, 'objects');

save(fullfile(save_path,'base_age_prediction.mat'),'wt_age','a53t_age','age_feature_names')

%% Check for mutation clusters
plot_var = ["Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = batch12group.plot_true_clusters(dr_level, dr_method, plot_var);

%% Check units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform"];["ReferenceWaveform","ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "PlatingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = 28; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
batch12group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Clusters in one plot
grouping_var = ["DIV","Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = batch12group.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Clusters split across subplots
batch12group.plot_dimensionality_reduction(batch12group.DimensionalityReduction.Unit.UMAP.Reduction, cluster_idx,group_labels_comb);

%% Plot cluster densities
grouping_var = ["DIV","Mutation"];
n_bins = 100;
value_array = batch12group.plot_cluster_densities(dr_level,dr_method,grouping_var,n_bins);

%% Plot cluster shifts
figure('Color','w');
batch12group.plot_cluster_shifts(value_array, group_labels_comb);

%% Generate clusters
clust_method = "louvain";
batch12group.clusterByFeatures(dr_method,dr_level,clust_method,0.2);

%% Plot waveforms
% figure("Color","w");
batch12group.plot_cluster_waveforms(clust_method);

%% Classification
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_test_val = []; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_network_features = ["Regularity","Burst","Catch22"];%["Regularity","Burst","Catch22","GraphFeatures"]
clf_unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_feature_names = "MeanInterSpikeInterval";
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:35]; %Selected values corresponding to grouping_var
clf_normalization = "baseline"; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = [];%"PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = -1; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;
clf_result = batch12group.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_feature_names, clf_useClustered,... 
                                    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
                                
[train_acc, test_acc, avg_score] = batch12group.assessClassifier(clf_result);

%% Age regression
level = "Recording"; %Unit or Recording
alg = "rf"; %rf, svm, cnb, knn
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset 
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on untreated, test on treated)
pooling_vals = [];
network_features = ["Regularity","Burst","Catch22"];
unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures","Catch22"];
useClustered = false;
normalization_var = "PlatingDate";
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = 5; % number of K-fold CV

age_result = batch12group.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);

%%
y_pred_wt = vertcat(age_result.Y_pred);
y_test_wt = vertcat(age_result.Y_test);

%%
y_pred_a53t = vertcat(age_result.Y_pred);
y_test_a53t = vertcat(age_result.Y_test);
y_pred = [y_pred_wt; y_pred_a53t];
y_test = [y_test_wt; y_test_a53t];
color_idx = [zeros(size(y_pred_wt)); ones(size(y_pred_a53t))];

min_diff = min(diff(unique(y_test)));
x_vals = y_test/min_diff;

figure('Color','w');
b = boxchart(x_vals, y_pred,'GroupByColor',color_idx,'MarkerSize',1);
hold on
p = plot([0 max(x_vals)*1.2],[0 max(y_test)*1.2],'k--'); p.Color(4) = 0.3;

%% Age regression plot
regression_var = "DIV";
color_var = "Mutation";
color_order = [];
batch12group.plot_regression_results(regression_var, color_var, color_order)

%% Feature trajectories
level = "Recording";
grouping_var = "DIV";
grouping_values = [7:7:35];
network_features = ["Graph"];
unit_features = [];["ActivityFeatures"];
normalization = [];
feature_names = [];
comp_var = "Mutation"; %Compare groups by this metadata field
pooling_vals = {};
useClustered = false;
tolerance = 1;

batch12group.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance)

%% Feature heatmap
level = "Recording";
grouping_var = "DIV";
grouping_values = [7:7:35];
network_features = ["all"];
unit_features = ["all"];
normalization = [];
comp_var = "Mutation"; %Compare groups by this metadata field
comp_val = "WT";
useClustered = false;
tolerance = 1;
color_lim = 3; %set maximum fold change to be visualized (clim capped at this value)

batch12group.plot_feature_heatmap(level, grouping_var, grouping_values, network_features, unit_features, normalization, comp_var, comp_val, useClustered, tolerance, color_lim)
