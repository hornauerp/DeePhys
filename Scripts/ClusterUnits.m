addpath(genpath("/home/phornauer/Git/DeePhys"))

%%
sorting_path_list = ["/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/SCR_rebuttal_week_3/221230/M05562/well000/sorted"];
emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();
params.QC.Amplitude = [];
params.Analyses.Bursts = 0;
params.Analyses.Connectivity.DDC = 0;

%%
rec_array = {};
for iPath = 1:length(sorting_path_list)
    metadata = struct();
    metadata.InputPath = sorting_path_list(iPath);
    path_parts = cellstr(strsplit(sorting_path_list(iPath), filesep));
    metadata.ChipID = string(path_parts{end-2}) + "_" + string(str2double(path_parts{end-1}(end)) + 1);
    metadata.PlatingDate = "221209";
    rec_array{iPath} = MEArecording(metadata, params);
end

%%
rec_group = RecordingGroup([rec_array{:}]);

%%
dr_method = "UMAP";
dr_level = "Unit";

clust_method = "kmeans";

%%
rec_group.reduceDimensionality(dr_level, dr_method);
rec_group.clusterByFeatures(dr_method,dr_level,clust_method);
%%
figure("Color","w");
rec_group.plot_cluster_waveforms(clust_method);