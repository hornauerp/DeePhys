addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list automatically
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'23011*','Quinpirole','M*','w*'}; %Each cell corresponds to one subdirectory
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Set parameters values
emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();
params.QC.Amplitude = [];
params.QC.FiringRate = [0.01 100];
params.QC.Axon = 0.9;
params.QC.Noise = 8;
params.Analyses.Bursts = 0;
params.Analyses.Connectivity.DDC = 0;
params.Analyses.Catch22 = 1;
params.Outlier.Method = []; %No outlier removal

split_times = [0 10; 25 35; 50 60]; %Cutouts of the recording
output_folder = "min_" + split_times(:,1); %Name folders according to the minute the cutout starts

min_N_units = 20; %Minimum number of units

%%
rec_array = {};

parfor iPath = 1:length(sorting_path_list)
    sorting_path = sorting_path_list(iPath);
    
    if startsWith(sorting_path, "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/230118/Quinpirole/M05867/")
        split_times = [0 10; 10 20; 35 45];
    else
        split_times = [0 10; 25 35; 50 60]; %Cutouts of the recording
    end
    output_folder = "min_" + split_times(:,1); %Name folders according to the minute the cutout starts
    output_path_list = split_sortings(sorting_path, split_times, output_folder);
    
    for iSplit = numel(output_path_list):-1:1
        
        split_params = params;
        if iSplit == 3 %Perform QC on first part, keep the same units for the rest
            split_params.QC.GoodUnits = [];
        else
            split_params.QC.GoodUnits = mearec.Parameters.QC.GoodUnits; %Currently only supported in boolean/logic indexing
        end
        
        metadata = struct();
        metadata.LookupPath = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx";
        metadata.InputPath = output_path_list(iSplit);
        
        metadata.PathPartIndex = [10,12,13]; %[RecordingDate, PlateID, WellID]
        metadata.Timepoint = 25*(iSplit - 1);%split_times(iSplit,1);
        metadata.AcuteTreatment = "Quinpirole";
        temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
        if exist(temp_file,"file")
            spk_temp = readNPY(temp_file);
            if length(unique(spk_temp)) > min_N_units
                mearec = MEArecording(metadata, split_params);
                rec_array = [rec_array mearec];
            else
                break
            end
        end
    end
end

%%
min_N_units = 20;
filtered_array = remove_low_unit_recordings(rec_array, min_N_units);
rg_wf = RecordingGroup(filtered_array);

%%
dr_level = "Unit"; %Unit or Recording
dr_method = "UMAP";
n_dims = 2; %Number of output dimensions
unit_features = ["ReferenceWaveform","ActivityFeatures"];%"ReferenceWaveform",
network_features = "all"; %Only relevant for level = Recording 
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "Timepoint"; %Only relevant when the recording was concatenated and units can be tracked (e.g. Timepoint or DIV)
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = "baseline"; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
rg.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, useClustered,...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%%
clust_method = "kmeans";
rg.clusterByFeatures(dr_method,dr_level,clust_method,4);
rg.plot_cluster_waveforms(clust_method);

%%
grouping_var = ["Mutation"];
figure("Color","w");
[cluster_idx, group_labels_comb] = rg.plot_true_clusters(dr_level, dr_method, grouping_var);

%%
rg.plot_dimensionality_reduction(rg.DimensionalityReduction.Unit.UMAP.Reduction, cluster_idx,group_labels_comb)

%%
n_bins = 100;
value_array = rg.plot_cluster_densities(dr_level,dr_method,grouping_var,n_bins);

%%
figure;plot_cluster_shifts(rg, value_array, group_labels_comb)

%%
dr_method = "UMAP";
dr_level = "Recording";
n_dims = 2; %Number of output dimensions
grouping_var = []; %Has to correspond to a MEArecording metadata field
unit_features = [];
network_features = "all";
useClustered = false;
rg.reduceDimensionality(dr_level, dr_method, n_dims, grouping_var, unit_features, network_features, useClustered);

%%
grouping_var = ["Mutation","Timepoint"];
figure("Color","w");
[cluster_idx, group_labels_comb] = rg.plot_true_clusters(dr_level, dr_method, grouping_var,50);