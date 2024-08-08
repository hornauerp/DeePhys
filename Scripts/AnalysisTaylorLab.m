addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
% path_logic = {'SCR*','*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
% sorting_path_list = generate_sorting_path_list(root_path, path_logic);
path_logic = {'230123','How_*','M*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "day_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for Taylor cells
rg_params.Selection.Inclusion = {{'Source','Taylor'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',9,12},{'Treatment',"ASO"},{'Mutation',"GBA"},{'Concentration',16,66}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
taylor = RecordingGroup(rec_array, rg_params);

%% Plot representative scatter plots
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"WT"},{'Treatment',"Untreated"}};
wt_taylor = RecordingGroup(rec_array, rg_params);
wt_units = arrayfun(@(x) length(x.Units),wt_taylor.Recordings);

rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"SNCA"},{'Treatment',"Untreated"}};
snca_taylor = RecordingGroup(rec_array, rg_params);
snca_units = arrayfun(@(x) length(x.Units),snca_taylor.Recordings);

%% Scatter plot
time_cutout = [0 120];
figure('Color','w');
subplot(2,1,1)
wt_taylor.Recordings(end-1).PlotNetworkScatter(time_cutout,blue)
subplot(2,1,2)
snca_taylor.Recordings(end).PlotNetworkScatter(time_cutout,red)

%% Check recordings unsupervised
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "RecordingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [21:7:35]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
feature_names = [];
sc_reduction = taylor.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Plot UMAP
plot_var = ["Mutation","Concentration"];
plot_centroid = false;
nodeSz = 20;
mapSz = 300;
sigma = mapSz/60;
cmap = othercolor('RdBu9',8);
cmap = cmap([1:3,8:-1:6],:);
pooling_vals = {{"WT","SNCA"},{0, 33, 100}};
figure("Color","w");
[sc_true_idx, sc_group_labels_comb] = taylor.plot_true_clusters(dr_level, dr_method, plot_var, pooling_vals, plot_centroid, nodeSz, mapSz, sigma, cmap);

%% Classification between mutations
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_pooling_val = {}; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_network_features = ["all"];["Regularity","Burst","Catch22"];%,"GraphFeatures"];%["Regularity","Burst","Catch22","GraphFeatures"]
clf_unit_features = ["all"];["ActivityFeatures","RegularityFeatures","WaveformFeatures","Catch22"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_feature_names = [];
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [21:7:35]; %Selected values corresponding to grouping_var
clf_normalization = "RecordingDate"; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = []; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = -1; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;
clf_result = taylor.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_pooling_val, clf_network_features, clf_unit_features, clf_feature_names, clf_useClustered,... 
                                    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
                                
[train_acc, test_acc, avg_score] = taylor.assessClassifier(clf_result);

%% Classification between mutations
clf_classification_var = "Treatment";
clf_result = taylor.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_pooling_val, clf_network_features, clf_unit_features, clf_feature_names, clf_useClustered,... 
                                    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
                                
[train_acc, test_acc, avg_score] = taylor.assessClassifier(clf_result);

%% Concentration regression (does not work well)
level = "Recording"; %Unit or Recording
regression_var = "Concentration";
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset 
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on untreated, test on treated)
pooling_vals = {};
network_features = ["all"];
unit_features = ["all"];
useClustered = false;
normalization_var = "PlatingDate";
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = -1; % number of K-fold CV

taylor.regressionByFeatureGroups(level, regression_var, stratification_var, stratification_values, pooling_vals, grouping_var, grouping_values, ...
    network_features, unit_features, useClustered, normalization_var, normalization, tolerance, N_hyper, K_fold);
taylor.plot_regression_results(regression_var);

%% Feature heatmap
level = "Recording";
grouping_var = "DIV";
grouping_values = [14:7:35];
network_features = ["Regularity","Burst","Catch22"];%["all"];
unit_features = [];["all"];
feature_names = [];
normalization = [];
comp_var = "Mutation"; %Compare groups by this metadata field
pooling_vals = {{"WT","SNCA"}};
useClustered = false;
tolerance = 0;
color_lim = 3; %set maximum fold change to be visualized (clim capped at this value)

taylor.plot_feature_heatmap(level, grouping_var, grouping_values, network_features, unit_features, feature_names, normalization, comp_var, pooling_vals, useClustered, tolerance, color_lim);

%% Feature trajectories
level = "Recording";
grouping_var = "DIV";
grouping_values = [14:7:35];
network_features = ["Regularity","Burst","Catch22"];
unit_features = [];["ActivityFeatures"];
normalization = [];
feature_names = [];
comp_var = "Mutation"; %Compare groups by this metadata field
pooling_vals = {};
useClustered = false;
tolerance = 0;

taylor.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance);

%% Check units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform","ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "PlatingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = []; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
taylor.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Clusters in one plot
grouping_var = ["DIV","Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = taylor.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Clusters split across subplots
taylor.plot_dimensionality_reduction(taylor.DimensionalityReduction.Unit.UMAP.Reduction, cluster_idx,group_labels_comb);

%% Plot cluster densities
grouping_var = ["DIV","Mutation"];
n_bins = 100;
value_array = taylor.plot_cluster_densities(dr_level,dr_method,grouping_var,n_bins);

%% Plot cluster shifts
figure('Color','w');
taylor.plot_cluster_shifts(value_array, group_labels_comb);

%% Generate clusters
clust_method = "louvain";
taylor.clusterByFeatures(dr_method,dr_level,clust_method,0.2);

%% Plot waveforms
% figure("Color","w");
colors = taylor.plot_cluster_waveforms(clust_method);
                
%% Cluster proportions
taylor.assignUnitClusterIdx("louvain",true,false);

separation_var = "ShortName";
pooling_vals = {};%{{"WT",["A53T","SNCA"]}};
% colors = flipud(colors);
taylor.plot_cluster_proportions(separation_var, pooling_vals, clust_method, colors);
setAllFontSizes(gcf,6.5)

%% Trajectories across unit clusters
level = "Recording";
grouping_var = "DIV";
grouping_values  = [14:7:35];
network_features = [];
unit_features = "RegularityFeatures";
normalization = []; %"baseline" or "scaled"
feature_names = [];
comp_var = "Mutation";
pooling_vals = {};
useClustered = true;
tolerance = 0;
taylor.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization,...
                feature_names, comp_var, pooling_vals, useClustered, tolerance, colors);

%% Pooled FCDI + Taylor
rg_params.Selection.Inclusion = {{'PlatingDate', 190308, 190903, 200121, 221209},{'DIV',20,21,27,28,34,35}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Treatment',"ASO","LNA"},{'Mutation',"GBA"},{'Concentration',16,66}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
pooled = RecordingGroup(rec_array, rg_params);

%% Feature groups x grouping values 
save_path = "/home/phornauer/Git/DeePhys/Data/Figure5";
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_pooling_val = {{"WT",["A53T","SNCA"]}}; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
grouping_value_list = [nan]; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = "RecordingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = -1; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;
[train_acc, test_acc, avg_score] = pooled.classifyByFeatureGroupsAndGroupingVar(clf_level, clf_alg, clf_classification_var, clf_pooling_val, clf_useClustered,...
                clf_grouping_var, grouping_value_list, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
save(fullfile(save_path, 'feature_group_heatmap.mat'), 'train_acc', 'test_acc', 'avg_score', 'feature_groups')

%% Classification between mutations
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_pooling_val = {{"WT",["A53T","SNCA"]}}; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_network_features = ["all"];%,"GraphFeatures"];%["Regularity","Burst","Catch22","GraphFeatures"]
clf_unit_features = ["all"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_feature_names = [];
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [21:7:35]; %Selected values corresponding to grouping_var
clf_normalization = "RecordingDate"; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = []; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = -1; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;
clf_result = pooled.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_pooling_val, clf_network_features, clf_unit_features, clf_feature_names, clf_useClustered,...
    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);

%% Pooled FCDI + Taylor
rg_params.Selection.Inclusion = {{'PlatingDate', 190308, 190903, 200121, 221209},{'DIV',20, 21}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',9, 12},{'Treatment',"ASO","LNA"},{'Mutation',"GBA"},{'Concentration',16,66}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
pooled = RecordingGroup(rec_array, rg_params);

%% Check pooled units 
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform"];["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "PlatingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = []; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
pooled.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Clusters in one plot
grouping_var = ["ShortName"];
figure("Color","w");
[cluster_idx, group_labels_comb] = pooled.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Clusters split across subplots
pooled.plot_dimensionality_reduction(pooled.DimensionalityReduction.Unit.UMAP.Reduction, cluster_idx,group_labels_comb);

%% Plot cluster densities
grouping_var = ["DIV","Mutation"];
n_bins = 100;
value_array = pooled.plot_cluster_densities(dr_level,dr_method,grouping_var,n_bins);

%% Plot cluster shifts
figure('Color','w');
pooled.plot_cluster_shifts(value_array, group_labels_comb);

%% Generate clusters
clust_method = "spectral";
pooled.clusterByFeatures(dr_method,dr_level,clust_method,4);

%% Plot waveforms
% figure("Color","w");
colors = pooled.plot_cluster_waveforms(clust_method);

%%
pooled.assignUnitClusterIdx("louvain",true,false);

separation_var = "ShortName";
pooling_vals = {};%{{"WT",["A53T","SNCA"]}};
colors = flipud(colors);
pooled.plot_cluster_proportions(separation_var, pooling_vals, clust_method, colors);
setAllFontSizes(gcf,6.5)

