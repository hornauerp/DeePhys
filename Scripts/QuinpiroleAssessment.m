addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'2301*','Quinpirole','*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "min_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
qp_array = recording_array_from_single_files(sorting_path_list, rm_con);
qp_array = remove_low_unit_recordings(qp_array, qp_array(1).Parameters.QC.N_Units);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Mutation',"LRRK2"},{'Treatment',"AS0"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
qp_fdci = RecordingGroup(qp_array, rg_params);


%% Check unsupervised recordings
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = [];%["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "Timepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = "baseline"; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
qp_fdci.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for mutation clusters
grouping_var = ["Mutation"];["Treatment"];["Treatment","Mutation"];
pooling_vals = {};{{["Untreated","ntLNA"], ["LNA","ASO"]}};{{["Untreated","ntLNA"], ["LNA","ASO"]},{"WT","A53T"}};
plot_centroid = false;
nodeSz = 20;
mapSz = 200;
figure("Color","w");
[cluster_idx, group_labels_comb] = qp_fdci.plot_true_clusters(dr_level, dr_method, grouping_var, pooling_vals, plot_centroid, nodeSz, mapSz);

%% Check units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "Timepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = "baseline"; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
qp_fdci.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for mutation clusters
grouping_var = ["Mutation"];
nodeSz = 5;
mapSz = 300;
figure("Color","w");
[cluster_idx, group_labels_comb] = qp_fdci.plot_true_clusters(dr_level, dr_method, grouping_var, nodeSz, mapSz);

%% Plot cluster densities
grouping_var = ["Mutation"];
n_bins = 100;
value_array = qp_fdci.plot_cluster_densities(dr_level,dr_method,grouping_var,n_bins);

%% Generate clusters
clust_method = "kmeans";
qp_fdci.clusterByFeatures(dr_method,dr_level,clust_method);

%% Plot waveforms
% figure("Color","w");
qp_fdci.plot_cluster_waveforms(clust_method);

%% Plot cluster shifts
figure('Color','w');
qp_fdci.plot_cluster_shifts(value_array, group_labels_comb);

%% Assign cluster IDs 
qp_fdci.assignUnitClusterIdx("kmeans",false);
level = "Unit";
grouping_var = "Timepoint";
grouping_values = nan;
network_features = [];
unit_features = ["ActivityFeatures","RegularityFeatures","Catch22"];
normalization = "baseline";
feature_names = [];
comp_var = "Mutation"; %Compare groups by this metadata field
useClustered = false;
tolerance = 0;
qp_fdci.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance);

%% Feature trajectories
level = "Recording";
grouping_var = "Timepoint";
grouping_values = nan;
network_features = ["Catch22"];
unit_features = ["ActivityFeatures"];
normalization = "baseline";
feature_names = [];
comp_var = "Mutation"; %Compare groups by this metadata field
useClustered = false;
tolerance = 1;

qp_fdci.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, useClustered, tolerance)

%% Feature heatmaps
level = "Recording";
grouping_var = "Timepoint";
grouping_values = nan;
network_features = ["Burst","Regularity","Catch22"];
unit_features = ["ActivityFeatures"];
normalization = "baseline";
comp_var = "Mutation"; %Compare groups by this metadata field
comp_val = "WT";
useClustered = false;
tolerance = 1;
color_lim = 3; %set maximum fold change to be visualized (clim capped at this value)

qp_fdci.plot_feature_heatmap(level, grouping_var, grouping_values, network_features, unit_features, normalization, comp_var, comp_val, useClustered, tolerance, color_lim)

%% Unit feature trajectories
level = "Unit";
grouping_var = "Timepoint";
grouping_values = nan;
network_features = ["Catch22"];
unit_features = ["ActivityFeatures"];
normalization = "baseline";
feature_names = [];
comp_var = "Mutation"; %Compare groups by this metadata field
useClustered = false;
tolerance = 1;

qp_fdci.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, useClustered, tolerance)