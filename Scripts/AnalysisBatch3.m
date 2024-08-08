addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'DeePhysS*','*','*','w*','sorted'};
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

load('/home/phornauer/Git/DeePhys/Data/Figure4/top_features.mat')

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'PlatingDate',200121},{'DIV',27},{'Mutation',"WT"},{'Treatment',"Untreated"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_untreated = RecordingGroup(rec_array, rg_params);

rg_params.Selection.Inclusion = {{'PlatingDate',200121},{'DIV',27},{'Mutation',"WT"},{'Treatment',"LNA"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_lna = RecordingGroup(rec_array, rg_params);

rg_params.Selection.Inclusion = {{'PlatingDate',200121},{'DIV',27},{'Mutation',"A53T"},{'Treatment',"Untreated"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
a53t_untreated = RecordingGroup(rec_array, rg_params);

rg_params.Selection.Inclusion = {{'PlatingDate',200121},{'DIV',27},{'Mutation',"A53T"},{'Treatment',"LNA"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
a53t_lna = RecordingGroup(rec_array, rg_params);

%% Extract spiketimes and bin
bin_limits = [0 120];
spikes{1} = histcounts(wt_untreated.Recordings(end).Spikes.Times,'BinWidth',0.1,'BinLimits',bin_limits);
spikes{2} = histcounts(wt_lna.Recordings(end).Spikes.Times,'BinWidth',0.1,'BinLimits',bin_limits);
spikes{3} = histcounts(a53t_untreated.Recordings(1).Spikes.Times,'BinWidth',0.1,'BinLimits',bin_limits);
spikes{4} = histcounts(a53t_lna.Recordings(end).Spikes.Times,'BinWidth',0.1,'BinLimits',bin_limits);

spike_mat = vertcat(spikes{:})';
spike_tbl = array2table(spike_mat,'VariableNames',["WT Untreated", "WT LNA", "A53T Untreated","A53T LNA"]);
stacked_color = [blues(7,:); blues(4,:); reds(7,:); reds(4,:)];
y_labels = ["A53T" + newline + "LNA", "A53T" + newline + "Untreated" , "WT" + newline + "LNA", "WT" + newline + "Untreated"];
%% Plot coactivity
figure('Color','w');
sp = stackedplot(spike_tbl,'XLabel',"Time [s]");
arrayfun(@(x) set(sp.LineProperties(x),'Color',stacked_color(x,:)),1:length(sp.LineProperties));
arrayfun(@(x) set(sp.AxesProperties(x),'YLimits',[-50 550]),1:length(sp.LineProperties));
ax = findobj(sp.NodeChildren, 'Type','Axes');
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','FontSize',fontsz)
arrayfun(@(x) set(ax(x).YLabel, 'String', y_labels(x)),1:length(ax))
arrayfun(@(x) set(x, 'XTick',0:300:1200, 'FontSize', fontsz, 'XTickLabels',0:30:120,'YTick',[0 500],'YTickLabels',[0 500]),ax)

%%
rg_params.Selection.Inclusion = {{'PlatingDate',200121},{'Mutation',"WT"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt = RecordingGroup(rec_array, rg_params);

rg_params.Selection.Inclusion = {{'PlatingDate',200121},{'Mutation',"A53T"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
a53t = RecordingGroup(rec_array, rg_params);


%% Feature trajectories WT
level = "Recording";
grouping_var = "DIV";
grouping_values = [6:7:27];
network_features = ["all"];
unit_features = ["all"];
normalization = [];
feature_names = [sc_features_sorted(1:10) nw_features_sorted(1:10)];
comp_var = "Treatment"; %Compare groups by this metadata field
pooling_vals = {{["Untreated","ntLNA"],["LNA"]}};
useClustered = false;
tolerance = 0;
color_lim = 3;
[wt_feature_table, wt_mean_mat, wt_sd_mat, group_labels_comb] = wt.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance);

[wt_color_mat, features] = wt.plot_feature_heatmap(level, grouping_var, grouping_values, network_features, unit_features, feature_names, normalization, comp_var, pooling_vals, useClustered, tolerance, color_lim);
%% A53T
[a53t_feature_table, a53t_mean_mat, a53t_sd_mat, group_labels_comb] = a53t.plot_feature_trajectories(level, grouping_var, grouping_values, network_features, unit_features, normalization, feature_names, comp_var, pooling_vals, useClustered, tolerance);

[a53t_color_mat, features] = a53t.plot_feature_heatmap(level, grouping_var, grouping_values, network_features, unit_features, feature_names, normalization, comp_var, pooling_vals, useClustered, tolerance, color_lim);

%% MLM WT
for r = 1:length(wt.Recordings)
    if  wt.Recordings(r).Metadata.Treatment == "ntLNA"
        wt.Recordings(r).Metadata.Treatment = "Untreated";
    end
end
comparison_var = "Treatment";
[wt_mlm_result_array, wt_p_G, wt_p_GxT, features] = wt.runMLM(network_features, unit_features, feature_names, useClustered, grouping_var, comparison_var);

%% MLM A53T
for r = 1:length(a53t.Recordings)
    if  a53t.Recordings(r).Metadata.Treatment == "ntLNA"
        a53t.Recordings(r).Metadata.Treatment = "Untreated";
    end
end
[a53t_mlm_result_array, a53t_p_G, a53t_p_GxT, features] = a53t.runMLM(network_features, unit_features, feature_names, useClustered, grouping_var, comparison_var);


%% 
rg_params.Selection.Inclusion = {{'PlatingDate',200121}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
lna = RecordingGroup(rec_array, rg_params);

%% UMAP single-cell
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = [];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = sc_features_sorted;
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [13:7:27]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
sc_reduction = lna.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);
plot_vars = ["Mutation","Treatment"];
figure;
[sc_true_idx, sc_group_labels_comb] = lna.plot_true_clusters(dr_level, dr_method, plot_vars);

%% UMAP network
unit_features = [];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = nw_features_sorted;
nw_reduction = lna.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);
plot_vars = ["Mutation","Treatment"];
figure;
[nw_true_idx, nw_group_labels_comb] = lna.plot_true_clusters(dr_level, dr_method, plot_vars);

%% UMAP combined
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [sc_features_sorted nw_features_sorted];
comb_reduction = lna.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);
plot_vars = ["Mutation","Treatment"];
figure;
[comb_true_idx, comb_group_labels_comb] = lna.plot_true_clusters(dr_level, dr_method, plot_vars);

%% Batch 1-3
rg_params.Selection.Inclusion = {}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
batch123 = RecordingGroup(rec_array, rg_params);
for r = 1:length(batch123.Recordings)
    if  batch123.Recordings(r).Metadata.Treatment == "ntLNA"
        batch123.Recordings(r).Metadata.Treatment = "Untreated";
    end
end

%% UMAP across batches
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [sc_features_sorted nw_features_sorted];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "PlatingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [13:7:27]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))
full_lna_reduction = batch123.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);
plot_vars = ["Mutation","Treatment"];
figure('Color','w');
[full_lna_idx, full_lna_labels] = batch123.plot_true_clusters(dr_level, dr_method, plot_vars);

%% Age prediction
rg_params.Selection.Inclusion = {{'Mutation',"WT"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV', 34}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt = RecordingGroup(rec_array, rg_params);

rg_params.Selection.Inclusion = {{'Mutation',"A53T"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV', 34}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
a53t = RecordingGroup(rec_array, rg_params);
%% Age prediction WT 
level = "Recording"; %Unit or Recording
alg = "rf"; %rf, svm, cnb, knn
stratification_var = "Treatment"; %Specify the variable by which to split training and test dataset 
stratification_values = "Untreated"; %Corresponding to stratification_var, if a value is specified then this will be used as training data (e.g. train on untreated, test on treated)
pooling_vals = [];
network_features = ["Regularity","Burst","Catch22"];
unit_features = ["ActivityFeatures","RegularityFeatures","WaveformFeatures","Catch22"];
useClustered = false;
normalization_var = "PlatingDate";
N_hyper = 0; %If >0 do hyperparameter optimization
K_fold = 5; % number of K-fold CV

wt_age_result = wt.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);

%% Age regression plot
regression_var = "DIV";
color_var = "Mutation";
color_order = [];
wt.plot_regression_results(regression_var, color_var, color_order)

%% Age prediction A53T
a53t_age_result = a53t.predictAge(level, alg, stratification_var, stratification_values, pooling_vals, network_features, unit_features, useClustered, normalization_var, N_hyper, K_fold);
a53t.plot_regression_results(regression_var, color_var, color_order)

%% Prepare for saving
wt_age = rmfield(wt_age_result, 'Mdl');
wt_age = rmfield(wt_age, 'objects');
a53t_age = rmfield(a53t_age_result, 'Mdl');
a53t_age = rmfield(a53t_age, 'objects');