addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
% path_logic = {'SCR*','*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
path_logic = {'230123','How_*','M*','w*','sorted','day_1'}; %Each cell corresponds to one subdirectory
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
taylor_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'DeePhysS*','2*','*','w*','sorted'};
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

load('/home/phornauer/Git/DeePhys/Data/Figure4/top_features.mat')

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
fcdi_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Remove low unit recordings
comb_array = [taylor_array fcdi_array];
for r = 1:length(comb_array)
   keep = arrayfun(@(x) length(x.SpikeTimes) > 2,comb_array(r).Units);
   comb_array(r).Units(~keep) = [];
end

%% Filter recordings for Taylor cells
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"WT"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
ctrl_untreated = RecordingGroup(comb_array, rg_params);

%% Full recording group
rg_params.Selection.Inclusion = {{'DIV',34,45}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
comb_rg = RecordingGroup(comb_array, rg_params);

%% Set parameters
%Save path
save_path = "/home/phornauer/Git/DeePhys/Data/Figure6";

% Dimensionality reduction
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
act_unit_features = ["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
wf_unit_features = "ReferenceWaveform";
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

%% Activity dimensionality reduction
ctrl_act.reduction = ctrl_untreated.reduceDimensionality(dr_level, dr_method, n_dims, act_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
ctrl_act.cluster_idx = ctrl_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Train on untreated WT
ctrl_act.prediction = ctrl_untreated.classifyClusteredUnits(ctrl_act.cluster_idx, act_unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);

%% Apply classifiers to other conditions
pooling_vals = {{"A53T"    "Control"    "GBA"    "SNCA"    "WT"},{["Untreated","ntLNA"],["LNA"]}};
[group_idx, group_labels] = comb_rg.combineMetadataIndices("Unit",["ShortName","Treatment"], pooling_vals);

train_idx = group_idx == 3; test_idx = ~train_idx; clf = []; 
act_pred = comb_rg.applyClusteredUnitClassifier(ctrl_act.cluster_idx, train_idx, test_idx, act_unit_features, clf,grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Cluster proportions from predictions
metadata_idx = findgroups(group_idx(test_idx));
label_order = [3,4,7,8,5,6,9,10,1,2];
act_data = nan(max(metadata_idx),max(ctrl_act.cluster_idx));
for i = 1:max(metadata_idx)
    act_data(i,:) = histcounts(act_pred.Y_pred(metadata_idx == i));
end

act_data = [act_data([1,2],:); histcounts(vertcat(ctrl_act.prediction.Y_pred)); act_data([3:end],:)];
% [~,sorted_idx] = sort(act_data(3,:),'descend');
act_data = act_data(label_order,:);
figure('Color','w');
bar(act_data./sum(act_data,2),'stacked');xticklabels(group_labels(label_order))
act_table = array2table(act_data,'RowNames',group_labels(label_order),'VariableNames',"Cluster " + (1:size(act_data,2)));

%% Waveform dimensionality reduction
ctrl_wf.reduction = ctrl_untreated.reduceDimensionality(dr_level, dr_method, n_dims, wf_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
ctrl_wf.cluster_idx = ctrl_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Train on untreated WT
ctrl_wf.prediction = ctrl_untreated.classifyClusteredUnits(ctrl_wf.cluster_idx, wf_unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);

%% Classify untreated WT
wf_pred = comb_rg.applyClusteredUnitClassifier(ctrl_wf.cluster_idx, train_idx, test_idx, wf_unit_features, clf,grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Cluster proportions from predictions
metadata_idx = findgroups(group_idx(test_idx));
label_order = [3,4,7,8,5,6,9,10,1,2];
wf_data = nan(max(metadata_idx),max(ctrl_wf.cluster_idx));
for i = 1:max(metadata_idx)
    wf_data(i,:) = histcounts(wf_pred.Y_pred(metadata_idx == i));
end

wf_data = [wf_data([1,2],:); histcounts(ctrl_wf.cluster_idx); wf_data([3:end],:)];
% [~,sorted_idx] = sort(wf_data(3,:),'descend');
wf_data = wf_data(label_order,:);
figure('Color','w');
bar(wf_data./sum(wf_data,2),'stacked');xticklabels(group_labels(label_order)); box off
wf_table = array2table(wf_data,'RowNames',group_labels(label_order),'VariableNames',"Cluster " + (1:size(wf_data,2)));

%% Combined dimensionality reduction
ctrl_comb.reduction = ctrl_untreated.reduceDimensionality(dr_level, dr_method, n_dims, comb_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
ctrl_comb.cluster_idx = ctrl_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Train on untreated WT
ctrl_comb.prediction = ctrl_untreated.classifyClusteredUnits(ctrl_comb.cluster_idx, comb_unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);

%% Classify untreated WT
comb_pred = comb_rg.applyClusteredUnitClassifier(ctrl_comb.cluster_idx, train_idx, test_idx, comb_unit_features, clf,grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Cluster proportions from predictions
metadata_idx = findgroups(group_idx(test_idx));
label_order = [3,4,7,8,5,6,9,10,1,2];
comb_data = nan(max(metadata_idx),max(ctrl_comb.cluster_idx));
for i = 1:max(metadata_idx)
    comb_data(i,:) = histcounts(comb_pred.Y_pred(metadata_idx == i));
end

comb_data = [comb_data([1,2],:); histcounts(ctrl_comb.cluster_idx); comb_data([3:end],:)];
% [~,sorted_idx] = sort(comb_data(3,:),'descend');
comb_data = comb_data(label_order,:);
figure('Color','w');
bar(comb_data./sum(comb_data,2),'stacked');xticklabels(group_labels(label_order)); box off
comb_table = array2table(comb_data,'RowNames',group_labels(label_order),'VariableNames',"Cluster " + (1:size(comb_data,2)));

%% Save data
save('cluster_proportions.mat','wf_table','act_table','comb_table')

%% Activity dimensionality reduction
act.reduction = comb_rg.reduceDimensionality(dr_level, dr_method, n_dims, act_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
act.cluster_idx = comb_rg.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Plot waveforms
[act.colors, act.mean_wfs] = comb_rg.plot_cluster_waveforms(clust_method);

%%
separation_var = ["ShortName","Treatment"];
pooling_vals = {{"A53T"    "Control"    "GBA"    "SNCA"    "WT"},{["Untreated","ntLNA"],["LNA"]}};
comb_rg.plot_cluster_proportions(separation_var, pooling_vals, clust_method,act.colors);

