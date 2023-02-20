addpath(genpath("/home/phornauer/Git/DeePhys"))
%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'SCR*','*','*','w*'}; %Each cell corresponds to one subdirectory
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings  
rec_array = recording_array_from_single_files(sorting_path_list);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Source','FCDI'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
%,{'Mutation','WT'}
rg_params.Selection.Exclusion = {{'DIV',7,12,14},{'Treatment',"ASO"},{'Mutation',"LRRK2"}}; %Cell array of cell arrays with fieldname + value
lna_group = RecordingGroup(rec_array_filtered, rg_params);

%% Check unsupervised 
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["WaveformFeatures","ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%"all"; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
lna_group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for LNA clusters
grouping_var = ["Treatment"];
figure("Color","w");
[cluster_idx, group_labels_comb] = lna_group.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Classification
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %rf, svm,
clf_classification_var = "Treatment"; %Metadata field that determines y_test/y_train
clf_classification_test_val = "LNA"; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_network_features = "all";
clf_unit_features = "all";
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = nan; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = [];%"PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = 25; % number of K-fold CV
clf_tolerance = 1;
result = lna_group.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_useClustered,... 
                                    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
                                
%% Evaluation
y_pred = vertcat(result.Y_pred);
y_test = vertcat(result.Y_test);

accuracy = sum(y_pred == y_test)/length(y_pred)

%% Score
test_object_array = [result.objects];
if iscell(test_object_array)
    test_object_array = cellfun(@(x) x(1), test_object_array);
end

scores = vertcat(result.scores);
metadata = [test_object_array.Metadata];
treatment = [metadata.Treatment];
lna_concentration = [metadata.Concentration];
lna_concentration(treatment == "ntLNA") = 0;

mutation = [metadata.Mutation];
[r,p] = corrcoef(lna_concentration,scores(:,1))