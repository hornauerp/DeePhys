

%% Age prediction
% Demonstrated on WT cultures, can be expanded to other cell lines
rg_params.Selection.Inclusion = {{'Mutation',"WT"}}; %Cell array of cell arrays with fieldname + value // empty defaults to including all recordings
rg_params.Selection.Exclusion = {{'DIV', 34}};%,{'ChipID',"4135_0","4043_0"}}; %Cell array of cell arrays with fieldname + value
wt = RecordingGroup(rec_array, rg_params);

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