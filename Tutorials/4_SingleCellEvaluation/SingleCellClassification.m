%% Single Cell Evaluation
% The other tutorials always reduced the complexity of different cell
% populations down to one representative value for the whole culture. Here,
% we will show how we can resolve differences on the single-cell level.

%% Add the toolbox to the path
addpath(genpath("/home/phornauer/Git/DeePhys")) %Change

%% Data Loading
% First, we need to load in the objects that were returned from the feature extraction script (1_FeatureExtraction). 
% If you have a larger number of recordings/sortings it is sensible to store the objects (MEArec) individually. 
% In that case, we use the loading approach described here to speed things up. Otherwise, you can also just load them manually/individually or - if possible - the array as a whole.

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
% path_logic = {'SCR*','*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
path_logic = {'230123','How_*','M*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "day_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for training data
% These recordings will be used to detect single-cell clusters and train
% the classifiers. We would therefore recommend to use control cultures
% that you assume to have as many of your cell types as possible.

rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"WT"},{'MediumChangeTimepoint',1}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_untreated = RecordingGroup(rec_array, rg_params);

%% Set parameters
%Save path
save_path = "/home/phornauer/Git/DeePhys/Data/Figure5"; %Optional // change

% Dimensionality reduction
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
act_unit_features = ["ActivityFeatures","RegularityFeatures","Catch22"];%ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
wf_unit_features = "ReferenceWaveform"; %"ReferenceWaveform","WaveformFeatures","
comb_unit_features = ["ReferenceWaveform", "ActivityFeatures","RegularityFeatures","Catch22"];
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

% Clustering 
clust_method = "louvain";
resolution_parameter = 0.4;

% Classification
N_hyper = 0;
K_fold = 5;

%% The following code describes how you can do the cell type classification based on activity features
% The same code works for the waveform or the combination of both, you just
% need to switch out act_unit_features for the corresponding unit_features
% variable

%% Activity dimensionality reduction
% We reduce the dimensionality based on the selected activity features

act.umap = wt_untreated.reduceDimensionality(dr_level, dr_method, n_dims, act_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
% And then use the louvain community detection algorithm to generate our
% clusters

act.cluster_idx = wt_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Plot waveforms
% We can then check if we see corresponding differences in the waveforms

[act.colors, act.mean_wfs] = wt_untreated.plot_cluster_waveforms(clust_method);

%% Classification
% Next, we use these cluster labels to train a classifier 

act.prediction = wt_untreated.classifyClusteredUnits(act.cluster_idx, act_unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);
act.prediction = rmfield(act.prediction,'Mdl'); %Optional // Remove the model if you want to save the results (pretty large file size)

%% Concatenate results from CVs
% We aggregate the results and store them in one struct variable. We also
% visualize the performance of our classification and visualize the feature
% "fingerprints" of the clusters.

act.y = vertcat(act.prediction.Y_test);
act.y_hat = vertcat(act.prediction.Y_pred);
figure('Color','w'); cm = confusionchart(act.y,act.y_hat,'RowSummary','row-normalized','ColumnSummary','column-normalized');

act.confusion_mat = cm.NormalizedValues;
act.features = act.prediction(1).Features;
act.pI.mean = mean(vertcat(act.prediction.predImp));
act.pI.sd = std(vertcat(act.prediction.predImp));
[act.pI.sorted,act.pI.sorted_idx] = sort(act.pI.mean,'descend');
% wt_taylor.plot_unit_cluster_features(act.cluster_idx, act.features(act.pI.sorted_idx(1:5)), colors)
act.feature_heatmap = wt_untreated.plot_unit_cluster_heatmap(act.cluster_idx, act.features(act.pI.sorted_idx(1:5)));

%% Save activity results
save(fullfile(save_path,'activity_prediction.mat'),'act') %Optional

%% Let's apply our classifier to some other lines
% Again, we will demonstrate the workflow using the activity classifier,
% but it works similarly for all other features

%% Full recording group
rg_params.Selection.Inclusion = {{'DIV',34,45}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%Cell array of cell arrays with fieldname + value
comb_rg = RecordingGroup(comb_array, rg_params);

%% Classification parameters
pooling_vals = {{"A53T"    "Control"    "GBA"    "SNCA"    "WT"},{["Untreated","ntLNA"],["LNA"]}};
[group_idx, group_labels] = comb_rg.combineMetadataIndices("Unit",["ShortName","Treatment"], pooling_vals);

train_idx = group_idx == 3; % Here you specify which of the lines you want to train the classifier on
% Check the group_labels variable to adjust this for your own dataset

test_idx = ~train_idx;  %Make all other lines the test set
clf = []; %Here we could also provide the classifier that we trained beforehand.
% But since we had CV enabled, we would not have trained on the full
% dataset, so we do that now 
act_pred = comb_rg.applyClusteredUnitClassifier(act.cluster_idx, train_idx, test_idx, act_unit_features, clf,grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Plot the results
% This is still a bit messy and needs some user input to adjust the order
% of conditions to be grouped appropriately
metadata_idx = findgroups(group_idx(test_idx));
label_order = [3,4,7,8,5,6,9,10,1,2]; %Change this
act_data = nan(max(metadata_idx),max(ctrl_act.cluster_idx));
for i = 1:max(metadata_idx)
    act_data(i,:) = histcounts(act_pred.Y_pred(metadata_idx == i));
end

act_data = [act_data([1,2],:); histcounts(vertcat(ctrl_act.prediction.Y_pred)); act_data([3:end],:)]; %Change this
% [~,sorted_idx] = sort(act_data(3,:),'descend');
act_data = act_data(label_order,:);
figure('Color','w');
bar(act_data./sum(act_data,2),'stacked');xticklabels(group_labels(label_order))