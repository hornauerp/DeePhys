addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list manually
% sorting_path_list = ["/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/SCR_rebuttal_week_4/230106/M05562/well000/sorted",
%     "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/SCR_rebuttal_week_4/230105/M04163/well005/sorted"];

%% Generate sorting_path_list automatically
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'SCR*','*','*','w*'}; %Each cell corresponds to one subdirectory
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Set parameters values
emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();
params.QC.Amplitude = [];
params.QC.FiringRate = [0.01 100];
params.QC.Axon = 1;
params.QC.Noise = 8;
params.Analyses.Bursts = 0;
params.Analyses.Connectivity.DDC = 0;
params.Outlier.Method = []; %No outlier removal

%%
rec_array = {};

parfor (iPath = 1:length(sorting_path_list),12)
    metadata = struct();
    metadata.LookupPath = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx";
    metadata.InputPath = sorting_path_list(iPath);
    temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
    if exist(temp_file,"file")
        spk_temp = readNPY(temp_file);
        if length(unique(spk_temp)) > 10
            rec_array{iPath} = MEArecording(metadata, params);
        end
    end
end
rec_array = [rec_array{:}];

%% Remove MEArecordings without units 
min_N_units = 20;
N_units = arrayfun(@(x) length(rec_array(x).Units),1:length(rec_array));
rec_array_filtered = rec_array(N_units >= min_N_units);
fprintf("Kept %i out of %i units\n", length(rec_array_filtered),length(rec_array))

%%
rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
%{'Source','Taylor'},{'Mutation','WT'}
rg_params.Selection.Exclusion = {{'DIV',7,12,14}}; %Cell array of cell arrays with fieldname + value
rec_group = RecordingGroup(rec_array_filtered, rg_params);

%%
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform"];%,"ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["Regularity","Burst"];%"all"; %Only relevant for level = Recording 
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
rec_group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%%
clust_method = "spectral";
rec_group.clusterByFeatures(dr_method,dr_level,clust_method);

%%
% figure("Color","w");
rec_group.plot_cluster_waveforms(clust_method);

%%
grouping_var = ["Source","Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = rec_group.plot_true_clusters(dr_level, dr_method, grouping_var);

%%
rec_group.plot_dimensionality_reduction(rec_group.DimensionalityReduction.Unit.UMAP.Reduction, cluster_idx,group_labels_comb)

%%
method = "spectral";
calc_feat = true;
assignUnitClusterIdx(rec_group,method,calc_feat);

%% Check unsupervised recordings
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = [];["ActivityFeatures","RegularityFeatures","WaveformFeatures"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Burst","Regularity","GraphFeatures"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = [];%"RecordingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = [];%Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = []; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
rec_group.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%%
grouping_var = ["EI_Ratio"];
figure("Color","w");
[cluster_idx, group_labels_comb] = rec_group.plot_true_clusters(dr_level, dr_method, grouping_var,50);

%%
rec_group.plot_dimensionality_reduction(rec_group.DimensionalityReduction.Recording.UMAP.Reduction, cluster_idx,group_labels_comb)

%%
dr_method = "tSNE";
dr_level = "Culture";
n_dims = 2; %Number of output dimensions
grouping_var = []; %Has to correspond to a MEArecording metadata field
unit_features = "all";
network_features = "all";
useClustered = false;
cluster_age = [21,28,35];
rec_group.reduceDimensionality(dr_level, dr_method, n_dims, grouping_var, unit_features, network_features, useClustered, cluster_age);

%%
grouping_var = ["DIV"];
figure("Color","w");
[cluster_idx, group_labels_comb] = rec_group.plot_true_clusters(dr_level, dr_method, grouping_var,50);


%% Check for correlations with unit number
object_group = [rec_group.Recordings];
input_table = object_group.getRecordingFeatures(network_features, unit_features, useClustered);
N_units = arrayfun(@(x) length(object_group(x).Units),1:length(object_group));
figure("Color","w");
for i = 1:size(input_table,2)
%     mdl = fitlm(N_units,input_table(:,i).Variables);
    nexttile
%     plot(mdl)
    scatter(N_units,input_table(:,i).Variables,10,cluster_idx,'filled')
    title(input_table.Properties.VariableNames(i))
    if i ==1
        ylim([0 0.5])
        cb = colorbar; cb.TickLabels = group_labels_comb;
    end
end