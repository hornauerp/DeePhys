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
%,{'Mutation','WT'}
rg_params.Selection.Exclusion = {{'DIV',12},{'Treatment',"ASO"},{'Mutation',"LRRK2"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
batch123group = RecordingGroup(rec_array, rg_params);

%% Check unsupervised 
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ActivityFeatures","RegularityFeatures"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["Regularity","Burst"];%["Regularity","Burst","Catch22"];%"all"; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = [];%["PlatingDate"]; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [14:7:28]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
batch123group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for mutation clusters
grouping_var = ["Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = batch123group.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Train on first two batches and test the third batch, assess treatment effect
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_val = []; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_test_var = "PlatingDate";
clf_test_val = 200121;
clf_network_features = ["Regularity","Burst"];%,"GraphFeatures"]
clf_unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:28]; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = "PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_tolerance = 1;
result = batch123group.applyClassifier(clf_level, clf_alg, clf_classification_var, clf_classification_val, clf_test_var, clf_test_val, clf_network_features,...
    clf_unit_features, clf_useClustered, clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_tolerance);

assessment_var = "Treatment";
assessment_val = "LNA";
accuracy_mat = batch123group.assessAppliedClassifier(result, assessment_var, assessment_val);