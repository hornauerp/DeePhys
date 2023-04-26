addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'230123','How_*','M*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "day_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
mc_array = recording_array_from_single_files(sorting_path_list, rm_con);
mc_array_filtered = remove_low_unit_recordings(mc_array, mc_array(1).Parameters.QC.N_Units);

%% Cluster purities over time
cluster_purities_sc_clustered = [];
tps = {-1,0,1,2,3,-1:3};

dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = true; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "MediumChangeTimepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

%% Loop over time points
for i = 1:length(tps)
    rg_params.Selection.Inclusion = {{'MediumChangeTimepoint',tps{i}}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
    rg_params.Selection.Exclusion = {{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
    single_mc = RecordingGroup(mc_array_filtered, rg_params);
    single_mc.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
        normalization_var, grouping_var, grouping_values, normalization, tolerance);
    cluster_purities_sc_clustered(i) = single_mc.calculateClusterPurity("UMAP", "kmeans");
end
%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname (char) + value (string)
mc_fdci = RecordingGroup(mc_array_filtered, rg_params);

%% Check unsupervised recordings
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = true; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "MediumChangeTimepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = -1; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
mc_fdci.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for mutation clusters
plot_groups = ["Mutation"];
pooling_vals = {};
plot_centroid = false;
nodeSz = 20;
mapSz = 300;
sigma = 5;
% cmap = [othercolor('Reds9',7); othercolor('Greens9',7); othercolor('Blues9',7)];
% cmap = cmap([3:7, 10:14, 17:end],:);
figure("Color","w");
[cluster_idx, group_labels_comb] = single_mc.plot_true_clusters(dr_level, dr_method, plot_groups, pooling_vals, plot_centroid, nodeSz, mapSz, sigma);%, cmap);

%% Check units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform","ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
single_mc.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Clusters in one plot
grouping_var = ["ShortName"];
figure("Color","w");
[cluster_idx, group_labels_comb] = single_mc.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Generate clusters
clust_method = "louvain";
single_mc.clusterByFeatures(dr_method,dr_level,clust_method,0.4);

%% Plot waveforms
single_mc.plot_cluster_waveforms(clust_method);

%% Assign cluster IDs 
single_mc.assignUnitClusterIdx("louvain",true,true);
level = "Unit";
grouping_var = "MediumChangeTimepoint";
grouping_values = nan;
network_features = [];
unit_features = ["ActivityFeatures","RegularityFeatures"];
normalization = [];
feature_names = [];
comp_var = "Mutation"; %Compare groups by this metadata field
pooling_vals = {};
useClustered = false;
tolerance = 0;
mc_fdci.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance);

%%
level = "Recording";
grouping_var = "MediumChangeTimepoint";
grouping_values = nan;
network_features = ["all"];
unit_features = [];
normalization = "baseline";
comp_var = "Mutation"; %Compare groups by this metadata field
comp_val = "WT";
useClustered = false;
tolerance = 0;
color_lim = 3; %set maximum fold change to be visualized (clim capped at this value)

mc_fdci.plot_feature_heatmap(level, grouping_var, grouping_values, network_features, unit_features, normalization, comp_var, comp_val, useClustered, tolerance, color_lim)