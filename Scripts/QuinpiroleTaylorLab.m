addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'2302*','How_*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "min_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
qp_taylor_array = recording_array_from_single_files(sorting_path_list, rm_con);
qp_taylor_array = remove_low_unit_recordings(qp_taylor_array, qp_taylor_array(1).Parameters.QC.N_Units);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Mutation',"GBA"},{'Treatment',"AS0"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
qp_taylor = RecordingGroup(qp_taylor_array, rg_params);

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
qp_taylor.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Check for mutation clusters
grouping_var = ["Mutation","Treatment"];  
pooling_vals = {};{{["Untreated","ntLNA"], ["LNA","ASO"]}};{{["Untreated","ntLNA"], ["LNA","ASO"]},{"WT","A53T"}};
plot_centroid = false;
nodeSz = 20;
mapSz = 300;
figure("Color","w");
[cluster_idx, group_labels_comb] = qp_taylor.plot_true_clusters(dr_level, dr_method, grouping_var, pooling_vals, plot_centroid, nodeSz, mapSz, 4, qp_results.response.cmap([4,3,2,1],:));%, stacked_color);

%% Assign cluster IDs 
level = "Recording";
grouping_var = "Timepoint";
grouping_values = nan;
network_features = ["Regularity","Burst","Catch22"];
unit_features = [];["ActivityFeatures","RegularityFeatures","Catch22"];
normalization = [];
feature_names = [];
comp_var = ["Mutation","Treatment"]; %Compare groups by this metadata field
useClustered = false;
tolerance = 0;
qp_taylor.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance);%, colors);

%% Check baseline units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform","ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "Timepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = 0; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
[bl_reduction, bl_input_table] = qp_taylor.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
clust_method = "louvain";
bl_cluster_idx = qp_taylor.clusterByFeatures(dr_method,dr_level,clust_method,0.4);

%% Check unit responses
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = "baseline"; %"baseline" (divided by first data point) or "scaled" [0 1]
[norm_reduction, norm_input_table] = qp_taylor.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Check for firing rate changes
% fr_idx = ((norm_input_table.MeanInterSpikeInterval_35 > 1.2) & (norm_input_table.MeanInterSpikeInterval_10 > 1.2)) * 1 + 1;
fr_idx = (norm_input_table.MeanInterSpikeInterval_10 < 1) * 1 + 1;

%%
nodeSz = 2;
figure('Color','w');
umap_reduction = RecordingGroup.plot_cluster_outlines(bl_reduction, bl_cluster_idx);
ax = get(gca);
arrayfun(@(x) set(ax.Children(x),'Visible','off'),1:length(ax.Children) - 1)

scatter(gca, umap_reduction(fr_idx == 1,1),umap_reduction(fr_idx == 1,2),nodeSz,repmat([0 0 0],sum(fr_idx == 1),1),'o','filled', 'MarkerEdgeColor','none')
scatter(gca, umap_reduction(fr_idx == 2,1),umap_reduction(fr_idx == 2,2),nodeSz,repmat([1 1 1],sum(fr_idx == 2),1),'o','filled', 'MarkerEdgeColor','none')

%%
nodeSz = 5;
figure('Color','w');
ax = axes;
RecordingGroup.plot_cluster_outlines(bl_reduction, fr_idx, ax, plot_centroid, nodeSz, mapSz);%, sigma, cmap)

%% Generate clusters
clust_method = "louvain";
norm_cluster_idx = qp_taylor.clusterByFeatures(dr_method,dr_level,clust_method,0.4);

%% FR changes per clusters
sel_idx = bl_cluster_idx;
% fr_indicator = norm_fr; %norm_input_table.MeanInterSpikeInterval_35
clust_counts = histcounts(sel_idx);
fr_counts = histcounts(sel_idx(fr_idx == 1));
ratios = fr_counts./clust_counts;
figure('Color','w');
bar(ratios)
box off
xlabel('Cluster ID')
ylabel('Ratio of units with decreased FR')

%% Check for mutation clusters
grouping_var = ["Mutation"];
pooling_vals = {};
plot_centroid = false;
nodeSz = 5;
mapSz = 300;
figure("Color","w");
[cluster_idx, group_labels_comb] = qp_taylor.plot_true_clusters(dr_level, dr_method, grouping_var, pooling_vals, plot_centroid,nodeSz, mapSz);

%% Clusters split across subplots
qp_taylor.plot_dimensionality_reduction(qp_taylor.DimensionalityReduction.Unit.(dr_method).Reduction, cluster_idx,group_labels_comb);

%% Plot cluster densities
grouping_var = ["Mutation"];
n_bins = 100;
value_array = qp_taylor.plot_cluster_densities(dr_level,dr_method,grouping_var);

%% Plot waveforms
% figure("Color","w");
colors = qp_taylor.plot_cluster_waveforms(clust_method);

%% Plot cluster shifts
figure('Color','w');
qp_taylor.plot_cluster_shifts(value_array, group_labels_comb);

%% Assign cluster IDs 
qp_taylor.assignUnitClusterIdx("louvain",true,true);
level = "Unit";
grouping_var = "Timepoint";
grouping_values = nan;
network_features = [];
unit_features = ["ActivityFeatures"];%,"RegularityFeatures","Catch22"
normalization = "baseline";
feature_names = "MeanInterSpikeInterval";
comp_var = [];%"Mutation"; %Compare groups by this metadata field
useClustered = false;
tolerance = 0;
[feature_table, mean_mat, sd_mat, group_labels_comb] = qp_taylor.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, ...
    feature_names, comp_var, pooling_vals, useClustered, tolerance, colors);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Timepoint',0}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Mutation',"GBA"},{'Treatment',"AS0"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
qp_baseline = RecordingGroup(qp_taylor_array, rg_params);

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
normalization = "baseline"; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
qp_baseline.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
clust_method = "louvain";
qp_baseline.clusterByFeatures(dr_method,dr_level,clust_method,0.4);

%%
qp_baseline.assignUnitClusterIdx("louvain",true,false);

%%
figure;
RecordingGroup.plot_cluster_outlines(qp_baseline.DimensionalityReduction.Unit.UMAP.Reduction, qp.Unit.cluster_idx, gca, plot_centroid, 10, mapSz)

%% Clusters split across subplots
qp_baseline.plot_dimensionality_reduction(qp_baseline.DimensionalityReduction.Unit.UMAP.Reduction, qp_taylor.Clustering.Unit.louvain.Index,group_labels_comb);

%% Plot waveforms
% figure("Color","w");
colors = qp_baseline.plot_cluster_waveforms(clust_method);

%% Predict QP response cluster from baseline activity
N_hyper = 0;
K_fold = 5;
qp_response.prediction = qp_baseline.classifyClusteredUnits(qp.Unit.cluster_idx, unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);

%% Evaluate
qp_response.y = vertcat(qp_response.prediction.Y_test);
qp_response.y_hat = vertcat(qp_response.prediction.Y_pred);
figure('Color','w'); 
cm = confusionchart(qp_response.y,qp_response.y_hat,'RowSummary','row-normalized','ColumnSummary','column-normalized');