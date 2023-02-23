addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'DeePhysS*','*','*','w*'};
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Source','FCDI'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',12},{'Treatment',"ASO"},{'Mutation',"LRRK2"},{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
batch12group = RecordingGroup(rec_array, rg_params);

%% Check unsupervised recordings
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = [];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "RecordingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [14:7:28]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = "PlatingDate"; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
batch12group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for mutation clusters
grouping_var = ["Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = batch12group.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Check units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = []; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
batch12group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
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
clust_method = "spectral";
batch12group.clusterByFeatures(dr_method,dr_level,clust_method,4);

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
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:28]; %Selected values corresponding to grouping_var
clf_normalization = "baseline"; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = [];%"PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = -1; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;
clf_result = batch12group.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_useClustered,... 
                                    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
                                
[train_acc, test_acc, avg_score] = batch12group.assessClassifier(clf_result);

%% Age regression
level = "Recording"; %Unit or Recording
alg = "rf"; %rf, svm, cnb, knn
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset (e.g. train on wildtype, test on mutation)
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data
network_features = "all";
unit_features = "all";
useClustered = false;
normalization_var = [];
N_hyper (1,1) = 0; %If >0 do hyperparameter optimization
K_fold (1,1) = -1; % number of K-fold CV

age_result = batch12group.predictAge(level, alg, stratification_var, stratification_values, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);

y_pred = vertcat(age_result.Y_pred);
y_test = vertcat(age_result.Y_test);
age_objects = horzcat(age_result.objects);
[age_group_idx, age_group] = batch12group.combineMetadataIndices(age_objects,"Mutation");
age_pred{1} = y_pred(age_group_idx == 2);
age_pred{2} = y_pred(age_group_idx == 1);
age{1} = y_test(age_group_idx == 2);
age{2} = y_test(age_group_idx == 1);
mean(abs(y_pred - y_test))

%%
level = "Recording";
grouping_var = "DIV";
grouping_values = [7:7:28];
network_features = ["Burst","Regularity"];
unit_features = ["ActivityFeatures"];
normalization = [];
feature_names = "NetworkRegularityFrequency";
comp_var = "Mutation"; %Compare groups by this metadata field
useClustered = false;
tolerance = 1;

batch12group.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, useClustered, tolerance)

%%
level = "Recording";
grouping_var = "DIV";
grouping_values = [7:7:35];
network_features = ["GraphFeatures"];%["Burst","Regularity","Catch22"];
unit_features = [];%["WaveformFeatures"];
normalization = [];
comp_var = "Mutation"; %Compare groups by this metadata field
comp_val = "WT";
useClustered = false;
tolerance = 1;
color_lim = 3; %set maximum fold change to be visualized (clim capped at this value)

batch12group.plot_feature_heatmap(level, grouping_var, grouping_values, network_features, unit_features, normalization, comp_var, comp_val, useClustered, tolerance, color_lim)
