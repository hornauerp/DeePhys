addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate baseline sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
% path_logic = {'SCR*','*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
path_logic = {'230123','How_*','M*','w*','sorted','day_1'}; %Each cell corresponds to one subdirectory
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
mc_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Generate QP sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'2302*','How_*','M*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "full_min_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
qp_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'DeePhysS*','2*','*','w*','sorted'};
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
fcdi_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Set parameters
%Save path
save_path = "/home/phornauer/Git/DeePhys/Data/Figure5";

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
resolution_parameter = 0.5;

% Classification
N_hyper = 0;
K_fold = 5;

%% Filter recordings for ctrl Taylor cells
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"WT"},{'MediumChangeTimepoint',1}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Treatment',"LNA"}}; %Cell array of cell arrays with fieldname + value
ctrl_untreated = RecordingGroup(mc_array, rg_params);

rg_params.Selection.Inclusion = {{'Source','Taylor'},{'MediumChangeTimepoint',1}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Mutation',"GBA"}}; %Cell array of cell arrays with fieldname + value
full_mc = RecordingGroup(mc_array, rg_params);

%% Test dataset
rg_params.Selection.Inclusion = {{'DIV',56},{'Timepoint',0}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Mutation',"GBA"}}; %Cell array of cell arrays with fieldname + value
qp_bl_taylor = RecordingGroup(qp_array, rg_params);

rg_params.Selection.Inclusion = {{'DIV',56},{'Timepoint',10},{'ShortName',"Control"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings 
rg_params.Selection.Exclusion = {{'Treatment',"LNA"}}; %Cell array of cell arrays with fieldname + value 
qp_post_taylor = RecordingGroup(qp_array, rg_params);

rg_params.Selection.Inclusion = {{'DIV',34}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
test_fcdi = RecordingGroup(fcdi_array, rg_params);

rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
test_rg = RecordingGroup([qp_bl_taylor.Recordings, test_fcdi.Recordings], rg_params);

full_rg = RecordingGroup([test_rg.Recordings, ctrl_untreated.Recordings], rg_params);

%% Combined dimensionality reduction
tic;
ctrl_comb.reduction = ctrl_untreated.reduceDimensionality(dr_level, dr_method, n_dims, wf_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);
toc;
%% Generate clusters
tic;
ctrl_comb.cluster_idx = ctrl_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);
toc;
%% Plot waveforms
[ctrl_comb.colors, ctrl_comb.mean_wfs] = ctrl_untreated.plot_cluster_waveforms(clust_method);

%% Train on untreated WT
tic;
ctrl_comb.prediction = ctrl_untreated.classifyClusteredUnits(ctrl_comb.cluster_idx, wf_unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);
toc;
%% Classify untreated WT
[group_idx, group_labels] = full_rg.combineMetadataIndices("Unit","DIV");
train_idx = group_idx == 2; test_idx = ~train_idx; clf = []; 

comb_pred = full_rg.applyClusteredUnitClassifier(ctrl_comb.cluster_idx, train_idx, test_idx, wf_unit_features, clf, grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Apply classifiers to other conditions
pooling_vals = {{"A53T"    "Control"    "SNCA"    "WT"},{["Untreated","ntLNA"],["LNA"]},{34,45,56}};
[metadata_idx, group_labels] = full_rg.combineMetadataIndices("Unit",["ShortName","Treatment","DIV"], pooling_vals);

%% Cluster compositions
cl_id = [1:2,4:max(metadata_idx)];
comb_data = nan(max(metadata_idx),max(ctrl_comb.cluster_idx));
for i = cl_id
    comb_data(i,:) = histcounts(comb_pred.Y_pred(metadata_idx == i), 0.5:1:max(ctrl_comb.cluster_idx)+0.5);
end
arranged_data = comb_data([4,5,6,7,3,8,9,1,2],:);
figure('Color','w');
ax = axes;
bar(ax, arranged_data./sum(arranged_data,2),'stacked');
set(ax,'Box','off','ColorOrder',ctrl_comb.colors)

%% Firing rate change
metadata_id = [4]; %Corresponding to heterogeneous untreated controls
post_timerange = [1200, 1500];
bl_units = full_rg.Units(ismember(metadata_idx,metadata_id));
for i = 1:length(post_units)
   bl_fr(i) = length(bl_units(i).SpikeTimes)./ 600;
end
% bl_af = vertcat(bl_units.ActivityFeatures);
% bl_fr = 1./[bl_af.MeanInterSpikeInterval];
post_units= [qp_post_taylor.Units];
for i = 1:length(post_units)
    post_fr(i) = sum(qp_post_taylor.Units(i).SpikeTimes > post_timerange(1) & qp_post_taylor.Units(i).SpikeTimes < post_timerange(2))./ (post_timerange(2) - post_timerange(1));
end
% post_af = vertcat(post_units.ActivityFeatures);
% post_fr = 1./[post_af.MeanInterSpikeInterval];

fr_change = post_fr./bl_fr';
fr_idx = fr_change<1;
bl_cluster_idx = comb_pred.Y_pred(ismember(metadata_idx,metadata_id));

clust_counts = histcounts(bl_cluster_idx);
fr_counts = histcounts(bl_cluster_idx(fr_idx==1));

figure('Color','w');
ax = axes;
bar(ax, fr_counts./clust_counts,'FaceColor','flat')
set(ax.Children, 'CData', ctrl_comb.colors)