%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Mutation',"GBA","LRRK2"},{'Treatment',"AS0"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
qp_full = RecordingGroup(all_qp, rg_params);

%% Check unsupervised recordings
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "Timepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = ["baseline"]; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
qp_full.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Check for mutation clusters
grouping_var = ["ShortName"];["Treatment", "ShortName"];  ["Treatment"];
pooling_vals = {};{{["Untreated","ntLNA"], ["LNA","ASO"]},{["WT","Control"], ["SNCA","A53T"]}};{{["Untreated","ntLNA"], ["LNA","ASO"]},{"WT","A53T"}};
plot_centroid = false;
nodeSz = 20;
mapSz = 200;
figure("Color","w");
[cluster_idx, group_labels_comb] = qp_full.plot_true_clusters(dr_level, dr_method, grouping_var, pooling_vals, plot_centroid, nodeSz, mapSz);

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
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
qp_full.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
clust_method = "louvain";
qp_full.clusterByFeatures(dr_method,dr_level,clust_method,0.2);

%%
figure('Color','w');
qp_full.plot_cluster_outlines(qp_full.DimensionalityReduction.Unit.UMAP.Reduction, qp_full.Clustering.Unit.louvain.Index)

%% Plot waveforms
colors = qp_full.plot_cluster_waveforms(clust_method);

%% Assign cluster IDs 
qp_full.assignUnitClusterIdx("louvain",true,true);

%% Cluster proportions
separation_var = "ShortName";
clust_method = "louvain";
pooling_vals = {};

qp_full.plot_cluster_proportions(separation_var, pooling_vals, clust_method, colors);

%% Trajectories
level = "Unit";
grouping_var = "Timepoint";
grouping_values = nan;
network_features = [];
unit_features = ["ActivityFeatures","RegularityFeatures"];%,"Catch22"
normalization = "baseline";
feature_names = [];
comp_var = "Mutation"; %Compare groups by this metadata field
useClustered = false;
tolerance = 0;
qp_full.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance, colors);