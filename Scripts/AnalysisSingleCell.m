addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
% path_logic = {'SCR*','*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
path_logic = {'230123*','How_*','M*','w*','sorted'}; %Each cell corresponds to one subdirectory
split_prefix = "day_*"; %Only use when previously split with split_sortings()
sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings
rm_con = true; %Remove connectivity data to save memory
rec_array = recording_array_from_single_files(sorting_path_list, rm_con);

%% Filter recordings for Taylor cells
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"WT"},{'MediumChangeTimepoint',1}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'Treatment',"LNA"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt_untreated = RecordingGroup(rec_array, rg_params);

%% Set parameters
%Save path
save_path = "/home/phornauer/Git/DeePhys/Data/Figure5";

% Dimensionality reduction
dr_level = "Unit"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
act_unit_features = ["ActivityFeatures","RegularityFeatures","Catch22"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
wf_unit_features = "ReferenceWaveform";
comb_unit_features = ["ReferenceWaveform", "ActivityFeatures","RegularityFeatures","Catch22"];
network_features = []; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
feature_names = [];
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = []; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = []; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = nan; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 0; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

% Clustering 
clust_method = "louvain";
resolution_parameter = 0.4;

% Classification
N_hyper = 0;
K_fold = 5;

%% Activity dimensionality reduction
act.umap = wt_untreated.reduceDimensionality(dr_level, dr_method, n_dims, act_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
act.cluster_idx = wt_untreated.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Plot waveforms
[act.colors, act.mean_wfs] = wt_untreated.plot_cluster_waveforms(clust_method);

%% Classification
act.prediction = wt_untreated.classifyClusteredUnits(act.cluster_idx, act_unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);
act.prediction = rmfield(act.prediction,'Mdl');

%% Concatenate results from CVs
act.y = vertcat(act.prediction.Y_test);
act.y_hat = vertcat(act.prediction.Y_pred);
figure('Color','w'); cm = confusionchart(act.y,act.y_hat,'RowSummary','row-normalized','ColumnSummary','column-normalized');

act.confusion_mat = cm.NormalizedValues;
act.features = act.prediction(1).Features;
act.pI.mean = mean(vertcat(act.prediction.predImp));
act.pI.sd = std(vertcat(act.prediction.predImp));
[act.pI.sorted,act.pI.sorted_idx] = sort(act.pI.mean,'descend');
% wt_taylor.plot_unit_cluster_features(act.cluster_idx, act.features(act.pI.sorted_idx(1:5)), colors)
act.feature_heatmap = wt_untreated.plot_unit_cluster_heatmap(act.cluster_idx, act.features(act.pI.sorted_idx(1:5)));

%% Save activity results
save(fullfile(save_path,'qp_activity_prediction.mat'),'act')

%% Waveform dimensionality reduction
wf.umap = wt_untreated.reduceDimensionality(dr_level, dr_method, n_dims, wf_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);
wf.Graph = wt_untreated.DimensionalityReduction.Unit.UMAP.Graph;

%% Generate clusters
wf.cluster_idx = wt_untreated.clusterByFeatures(dr_method, dr_level,clust_method, resolution_parameter);

%% Plot waveforms
[wf.colors, wf.mean_wfs] = wt_untreated.plot_cluster_waveforms(clust_method);

%% Classification
wf.prediction = wt_untreated.classifyClusteredUnits(wf.cluster_idx, wf_unit_features, grouping_var, grouping_values, feature_names, normalization, normalization_var, N_hyper, K_fold, tolerance);
wf.prediction = rmfield(wf.prediction,'Mdl');

%% Concatenate results from CVs
wf.y = vertcat(wf.prediction.Y_test);
wf.y_hat = vertcat(wf.prediction.Y_pred);
figure('Color','w'); cm = confusionchart(wf.y,wf.y_hat,'RowSummary','row-normalized','ColumnSummary','column-normalized');

wf.confusion_mat = cm.NormalizedValues;
wf.features = wf.prediction(1).Features;
wf.pI.mean = mean(vertcat(wf.prediction.predImp));
wf.pI.sd = std(vertcat(wf.prediction.predImp));
[wf.pI.sorted,wf.pI.sorted_idx] = sort(wf.pI.mean,'descend');

%% Save waveform results
save(fullfile(save_path,'qp_waveform_prediction.mat'),'wf')

%% Select top features for combined prediction
comb_features = [];%[wf.features(wf.pI.sorted_idx(1:10)), act.features(act.pI.sorted_idx(1:10))];

%% Combined dimensionality reduction
comb.umap = wt_untreated.reduceDimensionality(dr_level, dr_method, n_dims, comb_unit_features, network_features, comb_features, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);
comb.Graph = wt_untreated.DimensionalityReduction.Unit.UMAP.Graph;

%% Generate clusters
comb.cluster_idx = wt_untreated.clusterByFeatures(dr_method, dr_level, clust_method, resolution_parameter);

%% Plot waveforms
[comb.colors, comb.mean_wfs] = wt_untreated.plot_cluster_waveforms(clust_method);

%% Classification
comb.prediction = wt_untreated.classifyClusteredUnits(comb.cluster_idx, comb_unit_features, grouping_var, grouping_values, comb_features, normalization, normalization_var, N_hyper, K_fold, tolerance);
comb.prediction = rmfield(comb.prediction,'Mdl');

%% Concatenate results from CVs
comb.y = vertcat(comb.prediction.Y_test);
comb.y_hat = vertcat(comb.prediction.Y_pred);
figure('Color','w'); cm = confusionchart(comb.y,comb.y_hat,'RowSummary','row-normalized','ColumnSummary','column-normalized');

comb.confusion_mat = cm.NormalizedValues;
comb.features = comb.prediction(1).Features;
comb.pI.mean = mean(vertcat(comb.prediction.predImp));
comb.pI.sd = std(vertcat(comb.prediction.predImp));
[comb.pI.sorted,comb.pI.sorted_idx] = sort(comb.pI.mean,'descend');

%% Save combined results
save(fullfile(save_path,'qp_combined_prediction.mat'),'comb')

%%
wf_idx = startsWith(pred_results(1).Features,"Waveform");
act_feats = pred_results(1).Features(~wf_idx);
pI_act = pI(~wf_idx);
[top_pI_act,top_idx]= maxk(pI_act,30);
feature_names = act_feats(top_idx);
% feature_names = ["Asymmetry","AUC_peak_2","AUC_trough"];
wt_untreated.plot_unit_cluster_features(true_cluster_idx, feature_names, colors)

%%
data = wt_untreated.plot_unit_cluster_heatmap(true_cluster_idx);%, feature_names);

%%
pI_wf = pI(wf_idx);
x = 1:sum(wf_idx) + 1;
xx = 1:0.01:length(x); % Adjust for missing peak (0 variance)
y = [pI_wf(1:19) 0 pI_wf(20:60)];
yy = pchip(x,y,xx);
wf_table = wt_untreated.Recordings.getUnitFeatures("ReferenceWaveform");
mean_wf = mean(wf_table.Variables);
interp_wf = pchip(x,mean_wf,xx);

mean_wf_matrix = [mean_wfs{:}]';
interp_wfs = pchip(x,mean_wf_matrix,xx);
%%
figure('Color','w');plot(normalize(interp_wfs','range',[-1 1]),'LineWidth',1.5);axis off
figure('Color','w');plot(interp_wfs','LineWidth',1.5);axis off

%% Filter recordings for Taylor cells with LNA treatment
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'Mutation',"WT"},{'MediumChangeTimepoint',1}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
wt_all = RecordingGroup(rec_array, rg_params);

%% Activity dimensionality reduction
wt_all.reduceDimensionality(dr_level, dr_method, n_dims, act_unit_features, network_features, feature_names, useClustered,...
    normalization_var, grouping_var, grouping_values, normalization, tolerance);

%% Generate clusters
wt_all.clusterByFeatures(dr_method,dr_level,clust_method, resolution_parameter);

%% Remove units with 0 activity
wt_all.assignUnitClusterIdx("louvain",false,false);
wt_all.removeUnitsByCluster("louvain",5);

%% Apply classifiers to other 
[group_idx, group_labels] = wt_all.combineMetadataIndices("Unit","Treatment");
train_idx = group_idx == 2; test_idx = group_idx == 1; clf = []; 
clf_wt_lna_wf = wt_all.applyClusteredUnitClassifier(wf.cluster_idx, train_idx, test_idx, wf_unit_features, clf,grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Compare cluster proportions
[conc_idx,conc_label] = wt_all.combineMetadataIndices(wt_all.Units(group_idx == 1),"Concentration");
clust_wt_33_wf = histcounts(clf_wt_lna_wf.Y_pred(conc_idx==1));
clust_wt_100_wf = histcounts(clf_wt_lna_wf.Y_pred(conc_idx==2));
clust_wt_untreated_wf = histcounts(wf.y_hat);

%%
figure('Color','w');
subplot(1,2,1)
plot(horzcat(wf.mean_wfs{:}))
box off; xlim([0 50]); xticks(0:20:60); xticklabels(0:3); xlabel('Time [ms]')
subplot(1,2,2)
bar([clust_wt_untreated_wf./sum(clust_wt_untreated_wf); clust_wt_33_wf./sum(clust_wt_33_wf); clust_wt_100_wf./sum(clust_wt_100_wf)],'stacked')
box off
xticklabels([0,33,100])
xlabel('LNA concentration')
title('Cluster proportions')

%% Apply classifiers to other 
[group_idx, group_labels] = wt_all.combineMetadataIndices("Unit","Treatment");
train_idx = group_idx == 2; test_idx = group_idx == 1; clf = []; 
clf_wt_lna_act = wt_all.applyClusteredUnitClassifier(act.cluster_idx, train_idx, test_idx, act_unit_features, clf,grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Compare cluster proportions
clust_wt_33_act = histcounts(clf_wt_lna_act.Y_pred(conc_idx==1));
clust_wt_100_act = histcounts(clf_wt_lna_act.Y_pred(conc_idx==2));
clust_wt_untreated_act = histcounts(act.y_hat);
figure('Color','w'); bar([clust_wt_untreated_act./sum(clust_wt_untreated_act); clust_wt_33_act./sum(clust_wt_33_act); clust_wt_100_act./sum(clust_wt_100_act)],'stacked')
box off
xticklabels([0,33,100])
xlabel('LNA concentration')
title('Cluster proportions','FontWeight','normal')
setAllFontSizes(gcf, 6.5)

%% Filter recordings for Taylor cells with LNA treatment
rg_params.Selection.Inclusion = {{'Source','Taylor'},{'MediumChangeTimepoint',1}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {}; %Cell array of cell arrays with fieldname + value
all_taylor = RecordingGroup(rec_array, rg_params);

N_spikes = cellfun(@length, {all_taylor.Units.SpikeTimes});
all_taylor.DimensionalityReduction.Unit.ObjectGroup = all_taylor.Units;
fake_ids = ones(1,length(N_spikes)); fake_ids(N_spikes==0) = 99;
all_taylor.Clustering.Unit.Fake.Index = fake_ids;
all_taylor.removeUnitsByCluster("Fake",99);

%% Waveform prediction of other cultures
[group_idx, group_labels] = all_taylor.combineMetadataIndices("Unit",["Mutation","Concentration"]);
train_group_idx = 5;
train_idx = group_idx == train_group_idx;
test_idx = ~train_idx;
clf = [];
feature_names = [];
prediction_from_wf = all_taylor.applyClusteredUnitClassifier(wf.cluster_idx, train_idx, test_idx, wf_unit_features, clf, grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% 
metadata_idx = findgroups(group_idx(test_idx));

wf_data = nan(max(metadata_idx),max(wf.cluster_idx));
for i = 1:max(metadata_idx)
    wf_data(i,:) = histcounts(prediction_from_wf.Y_pred(metadata_idx == i));
end

wf_data = [histcounts(wf.y_hat); wf_data([5,6,3,4,1,2],:)];
wf_table = array2table(wf_data,'RowNames',group_labels([5,6,7,3,4,1,2]),'VariableNames',"Cluster " + (1:size(wf_data,2)));

%% Activity prediction
prediction_from_act = all_taylor.applyClusteredUnitClassifier(act.cluster_idx, train_idx, test_idx, act_unit_features, clf, grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% Plot
metadata_idx = findgroups(group_idx(test_idx));

act_data = nan(max(metadata_idx),max(act.cluster_idx));
for i = 1:max(metadata_idx)
    act_data(i,:) = histcounts(prediction_from_act.Y_pred(metadata_idx == i));
end

act_data = [histcounts(act.y_hat); act_data([5,6,3,4,1,2],:)];
act_table = array2table(act_data,'RowNames',group_labels([5,6,7,3,4,1,2]),'VariableNames',"Cluster " + (1:size(act_data,2)));

%% Combined prediction
prediction_from_comb = all_taylor.applyClusteredUnitClassifier(comb.cluster_idx, train_idx, test_idx, comb_unit_features, clf, grouping_var, ...
    grouping_values, feature_names, normalization, normalization_var, N_hyper, tolerance);

%% 
metadata_idx = findgroups(group_idx(test_idx));

comb_data = nan(max(metadata_idx),max(comb.cluster_idx));
for i = 1:max(metadata_idx)
    comb_data(i,:) = histcounts(prediction_from_comb.Y_pred(metadata_idx == i));
end

comb_data = [histcounts(comb.y_hat); comb_data];
comb_table = array2table(comb_data,'RowNames',group_labels([5,6,7,3,4,1,2]),'VariableNames',"Cluster " + (1:size(comb_data,2)));