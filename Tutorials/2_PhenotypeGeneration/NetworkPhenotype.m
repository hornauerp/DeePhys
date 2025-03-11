% This tutorial shows how to evaluate the extracted features based on some
% metadata condition (e.g. cell lines, mutations) to generate phenotyes.
% Specify the DeePhys root path
addpath(genpath("/home/phornauer/Git/DeePhys"))

% Or if your working directory is 2_FeatureExtraction we can also use a
% relative path
addpath(genpath("../../Functions"))
addpath(genpath("../../Classes"))
addpath(genpath("../../Toolboxes"))

%% Find relevant MEArecording objects
% We use a similar approach as for the feature extraction. If this is not
% applicable to your own data, you can specify the path_list any other way
% as well. 

root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/phornauer/iNeurons_dataset/network/"; %Root path
path_logic = {'2406*','T002*','Network','w*','sorter_output'}; %Variable parts
path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load data
% We do this in parallel to speed things up. By default we also dont load
% the CCG matrix between units, as this uses a lot of memory. 

rm_con = true; %If true, removes CCG data to save memory
mearec_array = recording_array_from_single_files(path_list, rm_con);

% The number of objects retrieved may be lower than the length of the path
% list, as sortings might have failed or provided too few units to do the
% feature extraction.

%% Optional: Remove networks that are too small
% As MEA data can be very heterogeneous, you will often time find that
% networks with too few units behave drastically different. You can remove
% them easily here. However, use with caution as this might bias your
% results!
unit_threshold = 1; % Minimum number of units
mearec_array = remove_low_unit_recordings(mearec_array, unit_threshold);

%% Filter for relevant dataset
% If we want to restrict the analysis to a certain subset, you can specify
% inclusion/exclusion criteria here. The filter only works with metadata
% names, e.g. mearec_array(1).Metadata
rg_params.Selection.Inclusion = {{'DIV',14},{'RecordingDate',"240610"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
pattern_rg = RecordingGroup(mearec_array, rg_params);

% This is our primary analysis object. Running an analysis several time
% will overwrite the previous results. So if you want to compare different
% parameters make sure to create several of these objects.

% If we look at the pattern_rg object now, it only contains the parameters
% and the recordings/units that we used to instantiate it. Additionally, it
% automatically detects recordings that were plated on the same day (same
% PlatingDate in the metadata) and have the same chipID. If the
% RecordingDate is provided, this will sort them by their DIV.

%% Check unsupervised network phenotypes
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions

% Here, we specify the input features. This can be a list of the feature
% names or simply "all". Specifying "ReferenceWaveform" basically runs the
% WaveMAP algorithm.
unit_features = ["all"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22" or "all"

% Network features are only relevant for dr_level = Recording
network_features =["Catch22","Regularity","GraphFeatures"]; % ["Regularity","Burst","Catch22","GraphFeatures"] or "all"

% If you already know which features you want to look at you can specify
% them here directly. You find the names by looking within the associated
% feature group, e.g. mearec_array(1).UnitFeatures.Catch22

% If used, it's best to set both unit_features and network_features to
% "all"
feature_names = [];

% Advanced feature, if you have previously run a unit clustering. This will
% use unit features for each unit cluster instead of one representative
% value for the entire network
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand

% Very useful when you have several experiments/batches that typically
% behave slightly different. Depending on the Metadata variable specified
% here, this will z-score feature for each group before further processing
normalization_var = ["PlatingDate"]; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)

% This variable will group recordings into cultures based on the metadata
% variable provided (e.g. "Timepoint" or "DIV"). This is important if you
% analyze a drug response or the developmental trajectory. In the second
% case, you would e.g. specify "DIV" to track the network across the
% development.
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked 

% This is only relevant if grouping_var is not empty.
% If you only want to group only certain DIVs (e.g., you want to exclude very
% early time points), then you can specify those here. Otherwise, leave it
% at nan.
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)

% This is only relevant if grouping_var is not empty.
% Here, you can specify an additional normalization, which might be
% relevant for drug experimens where you want to compare to the baseline.
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]

% Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
% THIS MIGHT OCCASIONALLY RESULT IN ERRORS!
tolerance = 0; 

% Now, we run the actual algorithm
pattern_rg.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Evaluating clustering results
% The results are stored in the pattern_rg object under
% DimensionalityReduction. Results are then sorted by dr_level and
% algorithm. Here e.g. pattern_rg.DimensionalityReduction.Recording.UMAP

% Of course we would like to see what we just did. We specify the variable
% by which we want to color the results and plot the low-dimensional
% representation.

cluster_var = ["Patterning"]; % Can also be a list of several metadata vars

nodeSz = 20; % Scatter point size
mapSz = 500; % Size of the plot
sigma = 10; % Spread of the background coloring
cmap = othercolor('Set16',6); % Adjust for the number of conditions you want to plot simultaneously

figure("Color","w");
[cluster_idx, group_labels_comb] = pattern_rg.plot_true_clusters(dr_level, dr_method, cluster_var, ...
    cmap=cmap,nodeSz=nodeSz, mapSz=mapSz, sigma=sigma);
% Legend labels are of course only relevant for this specific dataset
legend(["CHIR","RA","BMP4","CHIR+BMP4","RA+BMP4","CTRL"])

%% Optional: The pooling_vars variable
% Sometimes we might want to groups as a conditions, without a specific
% metadata variable being available to do so. E.g., maybe we want to
% compare conditions with CHIR to those without CHIR. This is possible
% using the pooling_vals variable. This is a cell array of cell arrays,
% for each metadata variable that each represent a grouping of conditions:
pooling_vals = {{{"C1","C4"},{"C2","C3","C5","C6"}}};
% It looks weird, but we really need 3 nested cell arrays :)
% number of metadata variables * number of groups * conditions in group

figure("Color","w");
[cluster_idx, group_labels_comb] = pattern_rg.plot_true_clusters(dr_level, dr_method, cluster_var, ...
    cmap=cmap,nodeSz=nodeSz, mapSz=mapSz, sigma=sigma,pooling_vals=pooling_vals);
legend(["ChIR","Non-CHIR"])

% The pooling_vars variable is available for a number of functions that
% need a grouping/clustering of samples by some metadata variable. It works
% exactly the same as here, so we will not mention it from here (as it is
% not really relevant for this example anyways).

%% 
% We would also like to put a number on how well our conditions separate.
% To this end we can calculate the cluster purity (ranges from -1 (worst) to 1 (best))

clust_method = "spectral"; % "kmeans", "hierarchical", "spectral", "gmm", "louvain"
cluster_purity = pattern_rg.calculateClusterPurity(dr_method, clust_method, cluster_var);
fprintf("Clustering by '%s' results in a purity of %.3f\n",cluster_var, cluster_purity)

% You can find the clustering results under pattern_rg.Clustering. If you
% want to make your own comparisons, you can get the true labels either
% from the cluster_idx output of the plotting function, or by using:
[c_idx, c_labels] = pattern_rg.combineMetadataIndices(pattern_rg.Recordings,cluster_var);
% Again, cluster_var can be a list, so you can get as granular as you want.

%% Supervised analysis
% If the unsupervised methods are promising, we might be interested in
% understanding which features drive those differences or if the aren't, we
% might want to try if supervised methods work better. To this end, we will
% generate a classifier. This function takes many arguments, however, most
% are identical to the unsupervised approach.

rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings

% We exclude the control group here, as we are interested in the
% differences between the patterned conditions
rg_params.Selection.Exclusion = {{'Patterning',"C6"}};% , %Cell array of cell arrays with fieldname + value
clf_rg = RecordingGroup(mearec_array, rg_params);

clf_level = "Recording";

% Classification algorithm to use. By default we use rf (random forest), as
% it typically performs best and allows for a straightworward inference of
% feature importances.
alg = 'rf'; %'svm', 'knn', 'cnb', 'rf'
clf_var = "Patterning";
pooling_vals = {};

% Since we're interested in the more intuitive features, we dont include
% the catch22 features here. 
network_features_clf = ["Regularity","GraphFeatures"];
unit_features_clf = ["WaveformFeatures","ActivityFeatures","RegularityFeatures"];

% If you have an a priori hypothesis about driving factors or you did some
% sort of feature selection procedure, you can specify those here.
feature_names = [];

% See above.
useClustered = false;
grouping_var = [];
grouping_values = nan;
normalization = [];
normalization_var = [];

% Indicates the number of iterations for the hyperparameter optimization, 0
% skips it. This increases the runtime A LOT, so I would only run this if
% you are not satisfied with the default hyperparameter performance.
N_hyper = 0;

% Indicates the k-fold split for cross-validation. -1 performs
% leave-one-out cross-validation. Since we typically work with small
% datasets for a machine learning setting, -1 is the default. This value
% directly correlates with the runtime, so increasing it will drastically
% increase performance. 
K_fold = -1;
tolerance = 0;

% Runs the actual algorithm
result = clf_rg.classifyByFeatureGroups(clf_level, alg, clf_var, pooling_vals, network_features_clf, ...
    unit_features_clf, feature_names, useClustered, grouping_var, grouping_values, normalization, ...
    normalization_var, N_hyper, K_fold, tolerance);

% You can also find the results in the clf_rg obejct under Classification. 

%% Evaluate classification results
% For a basic assessment we can run 
[metrics_special, train_acc, avg_pred_imp, conf_mat] = clf_rg.assessClassifier(result);
% All metrics_special are informative, but the F1score typically gives a
% good balanced view (0 to 1).

figure('Color','w');
confusionchart(conf_mat, ["CHIR","RA","BMP4","CHIR+BMP4","RA+BMP4"])

%% Find important features
% If the performance of the classifier is good enough, we can further look
% for the important features that underlie it. 

% We sort them for visual clarity
[sorted_pI, sort_idx] = sort(avg_pred_imp,'descend');

figure('Color','w');
bar(sorted_pI); xticks(1:length(avg_pred_imp)); xticklabels(result(1).Features(sort_idx)); box off
ylabel('Feature importance'); setAllFontSizes(gcf,7)

% A next step might now be to select some of the top performing features
% (ideally from different feature groups), and narrow down the number of
% features needed to perform a good classification. This makes the
% classifier more robust to outliers and should make it more generalizable.
% This optimal classifier can then be used to e.g. evaluate treatment approaches
% for a disease line. While we do not cover this in the tutorial here, as
% it doesnt apply to our dataset, you can look at the 

%% Look at those features
plt_feature_names = [result(1).Features(sort_idx(1:5))];

% Normally, this would actually show a line plot with errorbars, so e need 
% to specify some metadata variable here that connects the data points, 
% although it is not really needed here (suboptimally coded).
grouping_var = "DIV"; 
comp_var = "Patterning";
[feature_table, mean_mat, sd_mat, group_labels_comb] = clf_rg.plot_feature_trajectories(clf_level, grouping_var, ...
    grouping_values, network_features_clf, unit_features_clf, normalization, plt_feature_names, comp_var, pooling_vals);

%% Comparing two conditions in a heatmap
% To get a more aggregate view, one can also use a heatmap to show the
% development in feature differences over time (see original DeePhys
% paper). While this is not really applicable here, we can make up our own
% use case here, by comparing again CHIR vs non-CHIR recordings
pooling_vals = {{{"C1","C4"},{"C2","C3","C5","C6"}}};

% The color limit shows the how extreme the value differences can be before
% they are cutoff. This prevents one feature that is very different from
% skewing the entire heatmap.
color_lim = 3;

[color_mat, features] = clf_rg.plot_feature_heatmap(clf_level, grouping_var, grouping_values, network_features_clf, ...
    unit_features_clf, feature_names, normalization, comp_var, pooling_vals, useClustered, tolerance, color_lim);

% As you see in the x-axis, this is not an ideal visualization as we only
% have one DIV available. This works analogously if you have several
% recordings across the development.


% That is it for the core analysis pipeline for network analyses, check also the X
% tutorial for more utility functions and more advanced options :)