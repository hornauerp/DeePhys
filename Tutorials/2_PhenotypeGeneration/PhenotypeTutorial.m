%% Phenotype Generation
% Here, we demonstrate how you can systematically compare two or more cell lines and find metrics that are most reliable in discriminating between them. 
% This example can be easily adapted by changing any of the parameters.

%% Add the toolbox to the path
addpath(genpath("/home/phornauer/Git/DeePhys")) %Change

%% Data Loading
% First, we need to load in the objects that were returned from the feature extraction script (1_FeatureExtraction). 
% If you have a larger number of recordings/sortings it is sensible to store the objects (MEArec) individually. 
% In that case, we use the loading approach described here to speed things up. Otherwise, you can also just load them manually/individually or - if possible - the whole array as a whole.

%% Similar to what we did in the feature extraction, we try to automatically generate a list of paths containing the objects so that we can then parallelize the loading process. 
% Here, you need to specify your own 'root_path'  (common path trunk) and 'path_logic' (starting from the 'root_path').
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'DeePhysS*','19*','*','w*','sorted'}; %Change
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Now, we can go ahead and load the data (this might take a while depending on the number of sortings)
rm_con = true; %Remove connectivity data to save RAM
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Finally, you might want to look only at a specific subset of the data, or exclude individual recordings for various reasons. 
% Hereby you can use any metadata information that you provided during the feature extraction. You can also check again which information is available to you:
rec_array(1).Metadata

%% Any of these fields can now be used to include/exclude recordings:
rg_params.Selection.Inclusion = {{'Source','FCDI'}};    %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',12}};           %Cell array of cell arrays with fieldname + value
rg = RecordingGroup(rec_array, rg_params);

% Now we will continue using just the RecordingGroup object rg. We will start by building classifiers that will distinguish between the cell lines that we are interested in, here we use WT and A53T. 
% A lot of different customizations/normalization techniques are available and can be adjusted to the specific use case:

%% Set classification parameters
clf_level = "Recording";                        %Level of classification "Unit" or "Recording"
clf_alg = "rf";                                 %Classifier algorithm "rf", "svm", "cnb"
clf_classification_var = "Mutation";            %Metadata field that determines y_test/y_train // used for stratificatoin
clf_classification_test_val = [];               %Values corresponding to classification_var that is used as the test data (e.g. LNA)
clf_useClustered = false;                       %Use features that were calculated based on single-cell clusters
clf_grouping_var = "DIV";                       %Metadata field that groups recordings to cultures
clf_grouping_values = 7:7:35;                   %Selected values corresponding to grouping_var
clf_normalization = [];                         %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
clf_normalization_var = "PlatingDate";          %Normalization by each value of normalization_var (e.g. PlatingDate normalizes each batch and reduces batch-to-batch variability)
clf_N_hyper = 0;                                %If >0 do hyperparameter optimization
clf_K_fold = -1;                                %Number of K-fold CV // -1 corresponds to leave-one-out-CV
clf_tolerance = 1;                              %Tolerance of clf_grouping_values i.e. grouping_value +- tolerance
% (e.g. grouping_value = 6 and tolerance = 1 will consider values 5-7 to be
% identical; usefule for recordings that were performed on two subsequent
% days)

%% Preallocation of some variables
% For this tutorial, we will only look at the Unit/Single-cell features, but this can be easily extended:

save_path = "/home/phornauer/Git/DeePhys/Data/Figure3";
clf_network_features = [];
clf_unit_features = "all";                                                  %Unit features to use
feature_names = batch12group.returnFeatureNames("UnitFeatures");            %This returns a list of features that were extracted
sc_train_acc = nan(size(feature_names{1}));                                 %Training accuracy (use to optimize classifier if necessary)
sc_test_acc = sc_train_acc;                                                 %Test accuracy (what we are mostly interested in)
sc_score = sc_train_acc;                                                    %Classifer score ("certainty" of the classification)
sc_pred_imp = nan(length(feature_names{1}),length(clf_grouping_values));    %Predictor importance (importance of individual time points)

%% Here we loop over all features and build classifiers for each one individually and store accuracy and predictor importance values:
% Main loop for single cell features

for fn = 1:length(feature_names{1})
    clf_feature_names = feature_names{1}(fn);
    clf_result = rg.classifyByFeatureGroups(clf_level, clf_alg, clf_classification_var, clf_classification_test_val, clf_network_features, clf_unit_features, clf_feature_names, clf_useClustered,...
        clf_grouping_var, clf_grouping_values, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
    
    [sc_train_acc(fn), sc_test_acc(fn), sc_score(fn), sc_pred_imp(fn,:)] = rg.assessClassifier(clf_result);
    
    save(fullfile(save_path,'single_cell_features.mat'),'sc_train_acc', 'sc_test_acc', 'sc_score', 'sc_pred_imp','feature_names')
    fprintf('Finished feature %i/%i\n', fn, length(feature_names{1}))
end

%% And compare our results to a mixed linear model, assessing statistical differences between the cell lines and cell line x development interactions:
% sc_p_G are the p-values for a MLM comparing between groups
% sc_p_G is the corresponding p-value for a group-by-time interaction

[sc_mlm_result_array, sc_p_G, sc_p_GxT, sc_features] = runMLM(rg, clf_network_features, clf_unit_features, feature_names{1}, clf_useClustered, clf_grouping_var, clf_classification_var);
save(fullfile(save_path, 'single_cell_p_vals.mat'), 'sc_mlm_result_array','sc_p_G','sc_p_GxT','sc_features')

%% Let's get some visual representation of the data:
pooling_vals = {}; % Not necessary here, but can be used to pool conditions (e.g. different cell lines, treatments)

%% Feature trajectories over development
% Change variables to only plot subset of the data

[sc_feature_table, sc_mean_mat, sc_sd_mat, sc_group_labels_comb] = rg.plot_feature_trajectories(clf_level, clf_grouping_var, clf_grouping_values, clf_network_features, ...
    clf_unit_features, clf_normalization, feature_names{1}, clf_classification_var, pooling_vals, clf_useClustered, clf_tolerance);
save(fullfile(save_path, 'single_cell_values.mat'),'sc_feature_table', 'sc_mean_mat', 'sc_sd_mat', 'sc_group_labels_comb')

%% Feature heatmaps
color_lim = 3; %Indicating color limits for the heatmaps (roughly equivalent to matlab 'clim') 
[sc_color_mat, sc_features] = rg.plot_feature_heatmap(clf_level, clf_grouping_var, clf_grouping_values, clf_network_features, clf_unit_features, feature_names{1}, clf_normalization, clf_classification_var, clf_classification_test_val,...
    clf_useClustered, clf_tolerance, color_lim);
save(fullfile(save_path, 'single_cell_heatmap.mat'),'sc_color_mat', 'sc_features')

%% Dimensionality reduction
dr_level = "Recording";                                                                     %"Unit" or "Recording"
dr_method = "UMAP";                                                                         %"UMAP", "tSNE", "PCA"
n_dims = 2;                                                                                 %Number of output dimensions
unit_features = ["ReferenceWaveform","ActivityFeatures","RegularityFeatures","Catch22"];    %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["Regularity","Burst","Catch22","GraphFeatures"];                        %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false;                                                                       %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "PlatingDate";                                                          %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV";                                                                       %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = 14:7:35;                                                                  %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = [];                                                                         %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1;                                                                              %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

rg.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

                
%% Now we can plot it
plot_var = "Mutation"; %Color indicator, can be (almost) any metadata field)
figure("Color","w");
batch12group.plot_true_clusters(dr_level, dr_method, plot_var);

%% Classify feature groups x grouping values
% We might also be interested in having a higher level overview across the
% development. Here, we do not classify by individual features, but by
% feature groups and developmental time points:
pooling_vals = {};
grouping_value_list = [6:7:34,nan]; %We add 'nan', which includes all possible grouping values
[train_acc, test_acc, avg_score] = rg.classifyByFeatureGroupsAndGroupingVar(clf_level, clf_alg, clf_classification_var, pooling_vals, clf_useClustered,...
                clf_grouping_var, grouping_value_list, clf_normalization, clf_normalization_var, clf_N_hyper, clf_K_fold, clf_tolerance);
save(fullfile(save_path, 'feature_group_heatmap.mat'), 'train_acc', 'test_acc', 'avg_score', 'feature_groups')

%% Age regression parameters 
level = "Recording"; %Unit or Recording
alg = "rf"; %rf, svm, cnb, knn
stratification_var = "Mutation"; %Specify the variable by which to split training and test dataset 
stratification_values = []; %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on untreated, test on treated)
pooling_vals = {};
network_features = "all";
unit_features = "all";
useClustered = false;
normalization_var = "PlatingDate";
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = -1; % number of K-fold CV

%% Age regression for WT cultures
rg_params.Selection.Inclusion = {{'Mutation','WT'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',12},{'PlatingDate',200121}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_rg = RecordingGroup(rec_array, rg_params);

wt_age_result = wt_rg.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);