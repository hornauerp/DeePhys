addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list manually
sorting_path_list = ["/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/SCR_rebuttal_week_4/230106/M05562/well000/sorted",
    "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/SCR_rebuttal_week_4/230105/M04163/well005/sorted"];

%% Generate sorting_path_list automatically
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/SCR_rebuttal_week_4/";
path_logic = fullfile("*","*","w*"); %Each string corresponds to one subdirectory
sorting_path_list = generate_sorting_path_list(root_path, path_logic);

%% Set parameters values
emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();
params.QC.Amplitude = [];
params.QC.FiringRate = [0.1 100];
params.QC.Axon = 0.8;
params.QC.Noise = 10;
params.Analyses.Bursts = 0;
params.Analyses.Connectivity.DDC = 0;
params.Outlier.Method = []; %No outlier removal

%%
rec_array = {};

parfor iPath = 1:length(sorting_path_list)
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
empty_idx = arrayfun(@(x) isempty(rec_array(x).Units),1:length(rec_array));
rec_array(empty_idx) = [];

%%
rg_params.Selection.Inclusion = {{'Source','FCDI'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
rec_group = RecordingGroup(rec_array, rg_params);

%%
dr_method = "UMAP";
dr_level = "Unit";
n_dims = 2; %Number of output dimensions
grouping_var = []; %Has to correspond to a MEArecording metadata field // leave empty to normalize the whole dataset together
unit_features = ["ReferenceWaveform","ActivityFeatures"];
rec_group.reduceDimensionality(dr_level, dr_method, n_dims, grouping_var, unit_features);

%%
clust_method = "spectral";
rec_group.clusterByFeatures(dr_method,dr_level,clust_method);
%%
figure("Color","w");
rec_group.plot_cluster_waveforms(clust_method);

%%
grouping_var = "Treatment";
figure("Color","w");
rec_group.plot_true_clusters(dr_level, dr_method, grouping_var)

%%
method = "spectral";
calc_feat = true;
assignUnitClusterIdx(rec_group,method,calc_feat);

%%
dr_method = "PCA";
dr_level = "Recording";
n_dims = 2; %Number of output dimensions
grouping_var = []; %Has to correspond to a MEArecording metadata field
unit_features = "all";
network_features = "all";
useClustered = true;
rec_group.reduceDimensionality(dr_level, dr_method, n_dims, grouping_var, unit_features, network_features, useClustered);

%%
grouping_var = "Mutation";
figure("Color","w");
rec_group.plot_true_clusters(dr_level, dr_method, grouping_var)