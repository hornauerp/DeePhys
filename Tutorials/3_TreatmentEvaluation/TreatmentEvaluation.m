%% Treatment Evaluation
% Here, we demonstrate how you can evaluate treatments, either within one
% experiment or by including previous batches of the same cell lines.

%% Add the toolbox to the path
addpath(genpath("/home/phornauer/Git/DeePhys")) %Change

%% Data Loading
% First, we need to load in the objects that were returned from the feature extraction script (1_FeatureExtraction). 
% If you have a larger number of recordings/sortings it is sensible to store the objects (MEArec) individually. 
% In that case, we use the loading approach described here to speed things up. Otherwise, you can also just load them manually/individually or - if possible - the array as a whole.

%% Similar to what we did in the feature extraction, we try to automatically generate a list of paths containing the objects so that we can then parallelize the loading process. 
% Here, you need to specify your own 'root_path'  (common path trunk) and 'path_logic' (starting from the 'root_path').
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'DeePhysS*','19*','*','w*','sorted'}; %Change
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Now, we can go ahead and load the data (this might take a while depending on the number of sortings)
rm_con = true; %Remove connectivity data to save memory
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);


%% Select only the drug experiment for starters
% Here, we will look only at the experiment on the chronic LNA application 

rg_params.Selection.Inclusion = {{'PlatingDate',200121}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
lna = RecordingGroup(rec_array, rg_params);

%% Dimensionality Reduction for high-level info
% First, we check for any obvious treatment effects by looking at an
% aggregated plot, i.e. a dimensionality reduction
% This code only looks at single-cell features, but can be adjusted to your
% liking
% We exclude the first week, as the treatment effect appears to start after
% two weeks

dr_level = "Recording";                 %"Unit" or "Recording"
dr_method = "UMAP";                     %"UMAP", "tSNE", "PCA"
n_dims = 2;                             %Number of output dimensions
unit_features = "all";                  %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22" or simply "all"
network_features = [];                  %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"] or simply "all"
feature_names = [];                     %Here you can specify feature names to include in the analysis, e.g. features that proved to be predictive in prior experiments
useClustered = false;                   %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = [];                 %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV";                   %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [13:7:27];            %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = [];                     %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0;                          %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

sc_reduction = lna.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Plot Dimensionality Reduction
plot_vars = ["Mutation","Treatment"]; %Color coding variables, (almost) any metadata field name can be used

figure('Color','w');
[sc_true_idx, sc_group_labels_comb] = lna.plot_true_clusters(dr_level, dr_method, plot_vars);

%% Train on first two batches and test the third batch
clf_level = "Recording"; %Unit or Recording
clf_alg = "rf"; %"rf", "svm", "cnb"
clf_classification_var = "Mutation"; %Metadata field that determines y_test/y_train
clf_classification_val = []; %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_test_var = "PlatingDate"; %Metadata field that determines train/test set
clf_test_val = 200121; %Values corresponding to clf_test_var that 
clf_network_features = ["Regularity","Burst"];%,"GraphFeatures"]
clf_unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
clf_useClustered = false;
clf_grouping_var = "DIV"; %Metadata field that groups recordings to cultures
clf_grouping_values = [14:7:28]; %Selected values corresponding to grouping_var
clf_normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = "PlatingDate"; % normalization by each value of normalization_var
clf_N_hyper = 0; %If >0 do hyperparameter optimization
clf_tolerance = 1;
result = lna.applyClassifier(clf_level, clf_alg, clf_classification_var, clf_classification_val, clf_test_var, clf_test_val, clf_network_features,...
    clf_unit_features, clf_useClustered, clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_tolerance);

%% Assess performance on unknown batch
% Now we can check if 1. the performance of the clasifier on the control cultures of the
% unknown batch is good and only then can we 2. assess the effect of the
% treatment. 
assessment_var = "Treatment"; %Metadata field that contains the different conditions
pooling_vals = {{"LNA"},{"ntLNA","Untreated"}}; %Conditions that should be compared; if empty, then all individual conditions will be compared
accuracy_mat = lna.assessAppliedClassifier(result, assessment_var, pooling_vals);
% This returns a N_celllines x N_treatments matrix (length of pooling vals) containing the accuracy
% values of the corresponding cellline x treatment combination

%% Finally, we can check if the age prediction is affected
% This might give insight about developmental alterations due to the
% treatment
% Here, we only compare WT cultures
rg_params.Selection.Inclusion = {{'Mutation',"WT"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV', 34}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt = RecordingGroup(rec_array, rg_params);

%% Age prediction WT 
level = "Recording"; %Unit or Recording
alg = "rf"; %rf, svm, cnb, knn
stratification_var = "Treatment"; %Specify the variable by which to split training and test dataset 
stratification_values = "Untreated"; %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on untreated, test on treated)
pooling_vals = [];
network_features = ["Regularity","Burst","Catch22"];
unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures","Catch22"];
useClustered = false;
normalization_var = "PlatingDate";
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = 5; % number of K-fold CV

wt_age_result = wt.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);

%% Age regression plot
regression_var = "DIV"; %Metadata field containing the regressed variable
color_var = "Mutation";
color_order = [];
wt.plot_regression_results(regression_var, color_var, color_order)