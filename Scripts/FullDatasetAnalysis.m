addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'SCR*','*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
sorting_path_list_new = generate_sorting_path_list(root_path, path_logic);
path_logic = {'DeePhysS*','*','*','w*','sorted'};
sorting_path_list_old = generate_sorting_path_list(root_path, path_logic);
sorting_path_list = [sorting_path_list_new sorting_path_list_old];
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings  
full_rec_array = recording_array_from_single_files(sorting_path_list);

%% Remove MEArecordings without units 
min_N_units = 20;
N_units = arrayfun(@(x) length(full_rec_array(x).Units),1:length(full_rec_array));
rec_array_filtered = full_rec_array(N_units >= min_N_units);
fprintf("Kept %i out of %i units\n", length(rec_array_filtered),length(full_rec_array))

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Source','Taylor'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
%,{'Mutation','WT'} 
rg_params.Selection.Exclusion = {{'DIV',12},{'Treatment',"ASO"},{'Mutation',"LRRK2","GBA"}}; %Cell array of cell arrays with fieldname + value
full_rg = RecordingGroup(full_rec_array, rg_params);

%% Check unsupervised 
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform","ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["Regularity","Burst","Catch22","GraphFeatures"];%"all"; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = [];%"PlatingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = [];%"DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [28:7:35]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
full_rg.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for LNA clusters
grouping_var = ["ShortName","Treatment"];
figure("Color","w");
[cluster_idx, group_labels_comb] = full_rg.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Check units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = []; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
full_rg.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Clusters in one plot
grouping_var = ["Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = full_rg.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Clusters split across subplots
full_rg.plot_dimensionality_reduction(full_rg.DimensionalityReduction.Unit.UMAP.Reduction, cluster_idx,group_labels_comb);

%% Plot cluster densities
grouping_var = ["Mutation"];
n_bins = 100;
value_array = full_rg.plot_cluster_densities(dr_level,dr_method,grouping_var,n_bins);

%% Plot cluster shifts
figure('Color','w');
full_rg.plot_cluster_shifts(value_array, group_labels_comb);

%% Generate clusters
clust_method = "spectral";
full_rg.clusterByFeatures(dr_method,dr_level,clust_method);

%% Plot waveforms
% figure("Color","w");
full_rg.plot_cluster_waveforms(clust_method);

%% Classification
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_test_val = "WT"; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_network_features = ["Regularity","Burst","Catch22"];%,"GraphFeatures"];%["Regularity","Burst","Catch22","GraphFeatures"]
clf_unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures","Catch22"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:35]; %Selected values corresponding to grouping_var
clf_normalization = "RecordingDate"; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = []; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = -1; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;
clf_result = full_rg.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_useClustered,... 
                                    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
                                
[train_acc, test_acc, avg_score] = full_rg.assessClassifier(clf_result);

%% Train on first two batches and test the third batch, assess treatment effect
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_val = []; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_test_var = "PlatingDate";
clf_test_val = [200121, 221207, 221208];
clf_network_features = ["Regularity","Burst", "Catch22"];%,"GraphFeatures"]
clf_unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:28]; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = "PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_tolerance = 1;
applied_result = full_rg.applyClassifier(clf_level, clf_alg, clf_classification_var, clf_classification_val, clf_test_var, clf_test_val, clf_network_features,...
    clf_unit_features, clf_useClustered, clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_tolerance);

assessment_var = "Treatment";
assessment_val = "LNA";
accuracy_mat = full_rg.assessAppliedClassifier(applied_result, assessment_var, assessment_val);

%% Check LNA concentration
for i = 1:length(full_rg.Recordings)
    if full_rg.Recordings(i).Metadata.Treatment == "ntLNA"
        full_rg.Recordings(i).Metadata.Concentration = 0;
    end
end

level = "Recording"; %Unit or Recording
regression_var = "Concentration";
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset (e.g. train on wildtype, test on mutation)
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data
grouping_var = "DIV";
grouping_values = [21:7:35];
network_features = ["Regularity","Burst", "Catch22"];
unit_features = [];%["ActivityFeatures","RegularityFeatures","WaveformFeatures"];
useClustered = false;
normalization_var = [];
normalization = [];
tolerance = 1;
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = 5; % number of K-fold CV
lna_result = full_rg.regressionByFeatures(level, regression_var, stratification_var, stratification_values, grouping_var, grouping_values, ...
    network_features, unit_features, useClustered, normalization_var, normalization, tolerance, N_hyper, K_fold);
y_pred = vertcat(lna_result.Y_pred);
y_test = vertcat(lna_result.Y_test);
reg_objects = [lna_result.objects];
[mut_idx, mut_labels] = full_rg.returnMetadataArray(reg_objects,"Mutation");

full_rg.plot_regression_results(regression_var);

%%
