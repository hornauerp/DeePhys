% This tutorial shows how to use the extracted features on the single-cell
% level to find clusters that might correspond to different cell types etc.

% Ideally, this would be used together with some pharmacological
% intervention that indicates the correspondence of a certain cluster with
% a specific cell type (see original DeePhys paper).

% I would recommend to copy these scripts into a separate folder, where you
% can change them without being overwritten if you pull a newer version of
% DeePhys.

%% SKIP THIS PART IF YOU DID THE NETWORK TUTORIAL BEFORE

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
fprintf("Generated %i sorting paths\n",length(path_list))

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

%% SKIP UNTIL HERE IF YOU DID THE NETWORK TUTORIAL BEFORE

%% Check units (unsupervised)
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions

% Here, we specify the input features. This can be a list of the feature
% names or simply "all". Specifying "ReferenceWaveform" basically runs the
% WaveMAP algorithm.
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
feature_names = [];
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = ["PlatingDate"]; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. drug response)
grouping_values = []; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
pattern_rg.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

% Since the UMAP algorithm was developed by another group, there will be
% pop-ups that prompt you to download accelerants. If it is possible for
% you, I would recommend downloading them as it drastically speeds up
% computations. However, this is not mandatory. Also, you can ignore any
% error thrown by this function (typically caused by an attempted
% connection to the download servers for the accelerants/example datasets).

%% Clusters in one plot
grouping_var = ["Patterning"];
figure("Color","w");
[cluster_idx, group_labels_comb] = pattern_rg.plot_true_clusters(dr_level, dr_method, grouping_var);

%% Clusters split across subplots
% For the single-cell plots this is typically too overwhelming/crowded, so
% we can split the plot into different panels.
pattern_rg.plot_dimensionality_reduction(pattern_rg.DimensionalityReduction.Unit.UMAP.Reduction, cluster_idx, group_labels_comb);

% This might give us a first impression about some of the more obvious
% differences, but often they will be more subtle. 

%% Generating clusters
% To quantify these differences, we generate clusters from the
% dimensionality reduction results. This is basically an extension of the
% WaveMAP approach:

% For single-cell clustering, louvain works best out of all available
% methods, but doesnt allow specifying the number of clusters. However, we
% can tune that indirectly by changing the resolution parameter.

clust_method = "louvain"; %"kmeans", "hierarchical", "spectral", "gmm", "louvain"

% Either N clusters or resolution param for louvain. The default resolution
% parameter is 1, with smaller values yielding fewer clusters.
clust_param = 1; 
pattern_rg.clusterByFeatures(dr_method,dr_level,clust_method,clust_param);

% Only relevant for louvain clustering
n_clust = max(pattern_rg.Clustering.Unit.(clust_method).Index);
fprintf('Found %i clusters \n', n_clust)

%% Visualize clustering
% Coordinates of dimensionality reduction
reduction = pattern_rg.DimensionalityReduction.Unit.(dr_method).Reduction;
cluster_idx = pattern_rg.Clustering.Unit.(clust_method).Index;

plot_centroid = false; % Plots the centroid of each cluster on top
cluster_var = ["Patterning"]; % Can also be a list of several metadata vars

nodeSz = 5; % Scatter point size
mapSz = 500; % Size of the plot
sigma = 10; % Spread of the background coloring
cmap = othercolor('Set16',n_clust); % Adjust for the number of conditions you want to plot simultaneously

figure("Color","w")
ax = gca;
[reduction, cmap] = RecordingGroup.plot_cluster_outlines(reduction, cluster_idx, ax, plot_centroid, ...
    nodeSz, mapSz, sigma, cmap);

% Play around with the clust_param value to find a reasonable clustering, 
% e.g.if you want to isolate a particular cell cluster.

%% Alternative: Plot cluster densities
pooling_vals = {};
% cluster_idx or grouping_var depending on whether you want to see inferred
% or true clusters
cluster_idx = pattern_rg.Clustering.Unit.(clust_method).Index;

n_bins = 50; % Granularity of the 2D binning
smoothing_factor = 0.5; % Smoothing to obtain cleaner contours

[value_array, group_labels_comb] = pattern_rg.plot_cluster_densities(dr_level, dr_method, cluster_idx, ...
    pooling_vals, n_bins, smoothing_factor);

%% Changes in cluster distribution
% Not really applicable here, but if used over development it allows the
% visualization of the single-cell developmental trajectory. Ideally
% combined with different conditions to highlight distinct changes over
% time.
metadata_var = "Patterning";
[value_array, group_labels_comb] = pattern_rg.plot_cluster_densities(dr_level, dr_method, metadata_var, ...
    pooling_vals, n_bins, smoothing_factor);

figure("Color","w")
pattern_rg.plot_cluster_shifts(value_array, group_labels_comb)
% Here we see which patterning condition result in the biggest single-cell
% distribution changes.

%% Visualize features driving cluster differences
% This is a very intuitive way of identifying distinctive features.
feature_names = []; % We plot all features
pattern_rg.plot_unit_cluster_heatmap(cluster_idx,feature_names);

% Similarly, we can visualize this data as box plot. However, this is
% arguably better suited to visualize a pre-selection of features.
feature_names = ["FiringRate","CVInterSpikeInterval"]; % As an example, leave empty to see all

% We use the cmap from the clustering plot to match clusters
pattern_rg.plot_unit_cluster_features(cluster_idx,feature_names, cmap);
% This might cut off the upper part of the error bars (intentionally), to
% focus on the mean differences. Use 'axis tight' to adjust this.

%% Visualize cluster distribution across conditions
% This chart shows which percentage of units belong to each cluster for
% each condition specified.
metadata_var = "Patterning";
clust_method = "louvain";
colors = othercolor('Spectral9',n_clust);
pooling_vals = {};
figure('Color','w');
pattern_ratios = pattern_rg.plot_cluster_proportions(metadata_var, pooling_vals, clust_method, colors);
xticklabels(["CHIR","RA","BMP4","CHIR+BMP4","RA+BMP4","CTRL"])

%% And across indivual recordings
metadata_var = "ChipID";
figure('Color','w');
pattern_ratios = pattern_rg.plot_cluster_proportions(metadata_var, pooling_vals, clust_method, colors);
xticks(1:length(pattern_ratios))

%% Runing WaveMAP
% We expanded on the original WaveMAP approach by additionally
% incorporating activity-based features. To run (basically) the original:

wavemap_rg = RecordingGroup(mearec_array, rg_params);

dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = "ReferenceWaveform";
feature_names = [];
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = ["PlatingDate"]; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. drug response)
grouping_values = []; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
wavemap_rg.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

clust_param = 1;
wavemap_rg.clusterByFeatures(dr_method,dr_level,clust_method,clust_param);

%% Visualization
% We can visualize the different cluster waveforms:
wavemap_rg.plot_cluster_waveforms(clust_method);

%% Removing a cluster !USE WITH CARE!
% If you see a cluster/waveform that seems odd to you (too little activity,
% weird waveform shape), we can remove it. 

%%%%%%%%%%%%%%%%%%%
% THIS WILL REMOVE ALL INSTANCES OF THOSE UNITS FROM MEMORY, ALSO FROM ALL 
% OTHER VARIABLES. THE DATA NEEDS TO BE LOADED IN AGAIN FOR THE UNITS TO 
% REAPPEAR.
%%%%%%%%%%%%%%%%%%%

calc_feat = true; % This allows you to set useClustered to true when using UnitFeatures (see dimensionality reduction)
concatenated = false; % If the units result from a concatenated recording that was later split, we can remove this unit from every MEArecording
wavemap_rg.assignUnitClusterIdx(clust_method,calc_feat, concatenated);

ID = 2; % Cluster ID of units to remove
wavemap_rg.removeUnitsByCluster(clust_method, ID, concatenated);

% We can see that the cluster is now gone
wavemap_rg.plot_cluster_waveforms(clust_method);

% After you remove units, you should update spike times/features
for r = 1:length(wavemap_rg.Recordings)
    wavemap_rg.Recordings(r).updateSpikeTimes();
    wavemap_rg.Recordings(r).aggregateSingleCellFeatures();
end

%% Classification of cell clusters
% If - through a drug experiment or through dedicated cell type ratios in
% your cultures - you would like to generate the phenotype of your cells in
% a supervised fashion, this is also possible. This might be important if
% you want to identify a specific cell population in a different
% culture/cell line.

cluster_idx = pattern_rg.Clustering.Unit.louvain.Index; % Choose your clustering of choice
feature_groups = "all";
grouping_var = [];      % Only relevant when the recording was concatenated and units can be tracked (e.g. drug response)
grouping_values = [];   % Values corresponding to grouping_var (use nan if you want to use all available values)
feature_names = [];     % If you already know which features drive the difference between clusters
normalization = [];     % Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
normalization_var = []; 

% Indicates the number of iterations for the hyperparameter optimization, 0
% skips it. This increases the runtime A LOT, so I would only run this if
% you are not satisfied with the default hyperparameter performance.
N_hyper = 0;

% Indicates the k-fold split for cross-validation. -1 performs
% leave-one-out cross-validation, which we DO NOT want to use here, since 
% we typically work with 1000s of units. This value directly correlates 
% with the runtime, so increasing it will drastically increase performance. 
K_fold = 4;

tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))


result = pattern_rg.classifyClusteredUnits(cluster_idx, feature_groups, grouping_var, grouping_values, ...
    feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);

%% Evaluate classification results
% For a basic assessment we can run 
[metrics_special, train_acc, avg_pred_imp, conf_mat] = pattern_rg.assessClassifier(result);
% All metrics_special are informative, but the F1score typically gives a
% good balanced view (0 to 1).

figure('Color','w');
confusionchart(conf_mat)

%% Find important features
% If the performance of the classifier is good enough, we can further look
% for the important features that underlie it. 

% We sort them for visual clarity
[sorted_pI, sort_idx] = sort(avg_pred_imp,'descend');

figure('Color','w');
bar(sorted_pI); xticks(1:length(avg_pred_imp)); xticklabels(result(1).Features(sort_idx)); box off
ylabel('Feature importance'); setAllFontSizes(gcf,7)

% For a table/explanations of the different features, please take a look at
% the original publication (https://doi.org/10.1016/j.stemcr.2023.12.008)
% or the wiki on Github.

%% Apply a classifier
% While we don't really have a use case for this workflow here, we will
% showcase it anyways. First, we need to specify the cultures we want to
% train on (which might be some sort of ground truth), and the dataset we
% want to apply the classifier to (maybe some mixed cultures).
pooling_vals = {};
metadata_var = "Patterning";
[group_idx, group_labels] = pattern_rg.combineMetadataIndices("Unit", metadata_var, pooling_vals);

% Since this is just for display, we pick a random condition that we want
% to predict for
test_idx = group_idx == 3;

% And assign the rest to the training set.
train_idx = ~test_idx; 
train_cluster_idx = cluster_idx(train_idx); %Select only the cluster 
% Here we could also provide the classifier that we trained beforehand.
clf = [];
% But since we had CV enabled, we would not have trained on the full
% dataset, so we do that now 

result = pattern_rg.applyClusteredUnitClassifier(train_cluster_idx, train_idx, test_idx, feature_groups, ...
    clf, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

figure("Color","w")
conf_mat = confusionchart(cluster_idx(test_idx),result.Y_pred);

% That is it for the core analysis pipeline for single-cell analyses, if
% you have any questions, please do not hesitate to open an issue on the
% github page!