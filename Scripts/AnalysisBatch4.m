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
rg_params.Selection.Exclusion = {{'DIV',12},{'Treatment',"ASO"},{'Mutation',"LRRK2"}}; %Cell array of cell arrays with fieldname + value
batch4_group = RecordingGroup(rec_array, rg_params);

%% Check unsupervised 
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["WaveformFeatures","ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["Regularity","Burst","Catch22"];%"all"; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [14:7:28]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
batch4_group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for LNA clusters
grouping_var = ["Mutation","Treatment"];
figure("Color","w");
[cluster_idx, group_labels_comb] = batch4_group.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Classification
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_test_val = "WT"; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_network_features = ["Regularity","Burst","Catch22"];%["Regularity","Burst","Catch22","GraphFeatures"]
clf_unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:28]; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = []; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_K_fold = 5; % number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;
clf_result = batch4_group.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_useClustered,... 
                                    clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
                                
[train_acc, test_acc, avg_score] = batch4_group.assessClassifier(clf_result);

%% Check if LRRK2 is assessed as WT or A53T
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_val = []; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_test_var = "Mutation";
clf_test_val = "LRRK2";
clf_network_features = ["Regularity","Burst","Catch22"];%,"GraphFeatures"]
clf_unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:28]; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = "PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_tolerance = 1;
result = batch4_group.applyClassifier(clf_level, clf_alg, clf_classification_var, clf_classification_val, clf_test_var, clf_test_val, clf_network_features,...
    clf_unit_features, clf_useClustered, clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_tolerance);

%% Check LNA concentration
for i = 1:length(batch4_group.Recordings)
    if batch4_group.Recordings(i).Metadata.Treatment == "ntLNA"
        batch4_group.Recordings(i).Metadata.Concentration = 0;
    end
end

level = "Recording"; %Unit or Recording
regression_var = "Concentration";
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset (e.g. train on wildtype, test on mutation)
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data
grouping_var = "DIV";
grouping_values = [14:7:28];
network_features = [];"all";
unit_features = "all";
useClustered = false;
normalization_var = [];
normalization = [];
tolerance = 1;
N_hyper (1,1) = 0; %If >0 do hyperparameter optimization
K_fold (1,1) = 5; % number of K-fold CV
lna_result = batch4_group.regressionByFeatures(level, regression_var, stratification_var, stratification_values, grouping_var, grouping_values, ...
    network_features, unit_features, useClustered, normalization_var, normalization, tolerance, N_hyper, K_fold);
y_pred = vertcat(lna_result.Y_pred);
y_test = vertcat(lna_result.Y_test);
reg_objects = [lna_result.objects];
[mut_idx, mut_labels] = batch4_group.returnMetadataArray(reg_objects,"Mutation");

figure; boxchart(y_test/10,y_pred,'BoxWidth',1)
%%
figure;
% plot(y_test, y_pred,'k.')
plot(y_test(mut_idx==1),y_pred(mut_idx==1),'b.')
hold on
plot(y_test(mut_idx==2),y_pred(mut_idx==2),'r.')
%%
