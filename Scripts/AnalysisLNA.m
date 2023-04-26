addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'SCR*','*','*','w*','sorted'}; %Each cell corresponds to one subdirectory
sorting_path_list_new = generate_sorting_path_list(root_path, path_logic);
path_logic = {'DeePhysS*','*','*','w*','sorted'};
sorting_path_list_old = generate_sorting_path_list(root_path, path_logic);
sorting_path_list = [sorting_path_list_new sorting_path_list_old];
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Load processed MEArecordings  
full_rec_array = recording_array_from_single_files(sorting_path_list);
unit_threshold = 20;
rec_array = remove_low_unit_recordings(full_rec_array, unit_threshold);

%% Load top predictor names
data_path = "/home/phornauer/Git/DeePhys/Data/Figure4";
load(fullfile(data_path, "top_features.mat"))

%% Find relevant metadata information
metadata = [rec_array.Metadata];
plating_dates = unique([metadata.PlatingDate]);

%% Filter recordings for relevant LNA dataset
rg_params.Selection.Inclusion = {{'Source','FCDI'},{'PlatingDate',200121}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',9,12},{'Treatment',"ASO"},{'Mutation',"LRRK2"},{'Concentration',16,66}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
batch3 = RecordingGroup(rec_array, rg_params);
for i = 1:length(batch3.Recordings)
    if batch3.Recordings(i).Metadata.Treatment == "ntLNA"
        batch3.Recordings(i).Metadata.Concentration = 0;
    end
end

%% Check recordings unsupervised
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "RecordingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [14:7:28]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

%% Use top predictive features
feature_names = [sc_features_sorted, nw_features_sorted];
sc_reduction = batch3.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Plot UMAP
plot_var = ["Mutation","Concentration"];
plot_centroid = false;
nodeSz = 20;
mapSz = 300;
sigma = mapSz/60;
cmap = othercolor('RdBu9',8);
% cmap = cmap([6,5,1,2],:);
cmap = cmap([1,2,3,8,7,6],:);
% pooling_vals = {{"WT","A53T"},{["Untreated","ntLNA"], "LNA"}};
pooling_vals = {};
figure("Color","w");
[sc_true_idx, sc_group_labels_comb] = batch3.plot_true_clusters(dr_level, dr_method, plot_var, pooling_vals, plot_centroid, nodeSz, mapSz, sigma, cmap);


%% Filter recordings for relevant dataset
rg_params.Selection.Inclusion = {{'Source','FCDI'},{'PlatingDate',190308, 190903, 200121}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',9,12},{'Treatment',"ASO"}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
batch123 = RecordingGroup(rec_array, rg_params);

%% Use top predictive features
feature_names = [sc_features_sorted, nw_features_sorted];
sc_reduction = batch123.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Plot UMAP
plot_var = ["Mutation","Treatment"];
plot_centroid = false;
nodeSz = 20;
mapSz = 300;
sigma = mapSz/60;
cmap = othercolor('RdBu9',6);
cmap = cmap([6,5,1,2],:);
pooling_vals = {{"WT","A53T"},{["Untreated","ntLNA"], "LNA"}};
figure("Color","w");
[sc_true_idx, sc_group_labels_comb] = batch123.plot_true_clusters(dr_level, dr_method, plot_var, pooling_vals, plot_centroid, nodeSz, mapSz, sigma, cmap);

%% Apply classifier to 3. batch (LNA application)
level = "Recording"; %Unit or Recording
alg = "rf"; %"rf", "svm", "cnb"
classification_var = "Mutation"; %Metadata field that determines y_test/y_train
pooling_vals = {}; %Values corresponding to classification_var to pool together (e.g. ntLNA and Untreated)
test_var = "PlatingDate";
test_pooling_vals = {{[190308, 190903], 200121}};
network_features = ["all"];["Regularity","Burst"];%,"GraphFeatures"]
unit_features = ["all"];%["ActivityFeatures","RegularityFeatures","WaveformFeatures"]; %"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
useClustered = false;
grouping_var = "DIV"; %Metadata field that groups recordings to cultures
grouping_values = [14:7:28]; %Selected values corresponding to grouping_var
normalization = []; %Normalization within grouped units/recordings, "baseline" (divided by first datapoint) or "scaled" [0 1]
normalization_var = "PlatingDate"; % normalization by each value of normalization_var
N_hyper = 0; %If >0 do hyperparameter optimization
tolerance = 1;

feature_names = [];%[sc_features_sorted, nw_features_sorted];
result = batch123.applyClassifier(level, alg, classification_var, pooling_vals, test_var, test_pooling_vals, network_features, unit_features, feature_names, useClustered,...
    grouping_var, grouping_values, normalization, normalization_var, N_hyper, tolerance);

%% Check batch 3 results
assessment_var = "Treatment";
assess_pooling_vals = {{["Untreated","ntLNA"], "LNA"}};
accuracy_mat = batch123.assessAppliedClassifier(result, assessment_var, assess_pooling_vals);

%% Filter recordings for Batch 3 + 4
rg_params.Selection.Inclusion = {{'Source','FCDI'}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV',9,12},{'Treatment',"ASO"},{'Mutation',"LRRK2"},{'Concentration',16,66}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
batch1234 = RecordingGroup(rec_array, rg_params);
for i = 1:length(batch1234.Recordings)
    if batch1234.Recordings(i).Metadata.Treatment == "ntLNA"
        batch1234.Recordings(i).Metadata.Concentration = 0;
    end
    
    if batch1234.Recordings(i).Metadata.Treatment == "ASO"
        batch1234.Recordings(i).Metadata.Concentration = 16;
    end
end

%% 
test_pooling_vals = {{[190308, 190903], [200121, 221207, 221208]}};
feature_names = [];%[sc_features_sorted, nw_features_sorted];
result = batch1234.applyClassifier(level, alg, classification_var, pooling_vals, test_var, test_pooling_vals, network_features, unit_features, feature_names, useClustered,...
    grouping_var, grouping_values, normalization, normalization_var, N_hyper, tolerance);

%% Check batch 4 results
assessment_var = "Treatment";
assess_pooling_vals = {{["Untreated","ntLNA"], ["LNA","ASO"]}};
accuracy_mat = batch1234.assessAppliedClassifier(result, assessment_var, assess_pooling_vals);

%%
[conc_idx, conc_val] = batch1234.combineMetadataIndices(result.objects, "Concentration");
concentrations = double(conc_val(conc_idx));
[mut_idx, mut_val] = batch1234.combineMetadataIndices(result.objects, "Mutation");
colors = othercolor('RdBu9',2);

figure('Color','w');
for i = 1:length(mut_val)
    x = (concentrations(mut_idx == i)/16) + (i-1)*0.4;
    y = result.scores(mut_idx == i,i);
   boxchart(x,y, 'BoxFaceColor',colors(i,:),'BoxWidth',0.3)
   hold on
end
yline(0.5,'k:','Label',"Decision" + newline +  "threshold",'LabelVerticalAlignment','middle','FontSize',6.5,'LabelHorizontalAlignment','right')
leg = legend(mut_val,'Box','off','Location','southwest'); leg.ItemTokenSize = [5 5];
xticks([0.2 2.2 6.2])
xticklabels([0 33 100])
xlabel('LNA concentration')
xlim([-0.5 9])
ylim([0 1])
ylabel('Classifier confidence')

%% Concentration regression
regression_var = "Concentration";
K_fold = -1;
stratification_var = "Mutation";
stratification_values = [];
reg_pooling_vals = {};


result = batch1234.regressionByFeatureGroups(level, regression_var, stratification_var, stratification_values, reg_pooling_vals, grouping_var, grouping_values, ...
    network_features, unit_features, useClustered, normalization_var, normalization, tolerance, N_hyper, K_fold);

batch1234.plot_regression_results(regression_var);
%% Check recordings unsupervised
dr_level = "Recording"; %"Unit" or "Recording"
dr_method = "UMAP"; %"UMAP", "tSNE", "PCA"
n_dims = 2; %Number of output dimensions
unit_features = ["all"];%"ReferenceWaveform","WaveformFeatures","ActivityFeatures","RegularityFeatures","GraphFeatures","Catch22"
network_features = ["all"];%["Regularity","Burst","Catch22"]; %Only relevant for level = Recording ["Regularity","Burst","Catch22","GraphFeatures"]
useClustered = false; %Only usable when clustering and assignClusterIdx was run beforehand
normalization_var = "RecordingDate"; %Use to normalize by groups of the normalization_var (e.g. PlatingDate would normalize by individual batches)
grouping_var = "DIV"; %Only relevant when the recording was concatenated and units can be tracked (e.g. "Timepoint" or "DIV")
grouping_values = [14:7:28]; %Values corresponding to grouping_var (use nan if you want to use all available values)
normalization = []; %"baseline" (divided by first data point) or "scaled" [0 1]
tolerance = 1; %Tolerance for the grouping_values (e.g. if you want to group recordings with a DIV that is off by one (14 +/- 1))

%% Use top predictive features
feature_names = [];%[sc_features_sorted, nw_features_sorted];
sc_reduction = batch1234.reduceDimensionality(dr_level, dr_method, n_dims, unit_features, network_features, feature_names, useClustered, ...
                    normalization_var, grouping_var, grouping_values, normalization, tolerance);
                
%% Plot UMAP
plot_var = ["Mutation","Concentration"];
plot_centroid = true;
nodeSz = 10;
mapSz = 300;
sigma = mapSz/60;
cmap = othercolor('RdBu9',8);
cmap = cmap([1:3,8:-1:6],:);
pooling_vals = {{"WT","A53T"},{0, 33, 100}};
figure("Color","w");
[sc_true_idx, sc_group_labels_comb] = batch1234.plot_true_clusters(dr_level, dr_method, plot_var, pooling_vals, plot_centroid, nodeSz, mapSz, sigma, cmap);

