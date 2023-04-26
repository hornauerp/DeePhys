addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'2302*','How_*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "min_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for Taylor cells
rg_params.Selection.Inclusion = {{'Mutation',"WT"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_untreated = RecordingGroup(rec_array, rg_params);

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
grouping_var = "Timepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = "baseline"; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

% Clustering 
clust_method = "louvain";
resolution_parameter = 0.4;

% Classification
N_hyper = 0;
K_fold = 5;

% Plot trajectories
comp_var = "Mutation";
pooling_vals = {};
plot_normalization = [];
plot_feature_names = [];%"MeanInterSpikeInterval";

%% Activity dimensionality reduction
act.umap = wt_untreated.reduceDimensionality(dr_level, dr_method, n_dims, act_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
act.cluster_idx = wt_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Activity dimensionality reduction
comb.umap = wt_untreated.reduceDimensionality(dr_level, dr_method, n_dims, comb_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
comb.cluster_idx = wt_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%%
wt_untreated.assignUnitClusterIdx("louvain",true,true);

%% Plot waveforms
[comb.colors, comb.mean_wfs] = wt_untreated.plot_cluster_waveforms(clust_method);

%%
wt_untreated.plot_feature_trajectories(dr_level, grouping_var, grouping_values, network_features, act_unit_features, plot_normalization, plot_feature_names, comp_var, pooling_vals, useClustered, tolerance);
