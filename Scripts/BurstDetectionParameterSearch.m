addpath(genpath("/home/phornauer/Git/DeePhys"))
%% Generate sorting_path_list 
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'SCR*','*','*','w*'}; %Each cell corresponds to one subdirectory
% path_logic = {'DeePhysS*','*','*','w*'};
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Set parameters values
emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();
params.QC.Amplitude             = [];
params.QC.FiringRate            = [0.01 100];
params.QC.Axon                  = 0.8;
params.QC.Noise                 = 8;
params.QC.N_Units               = 10;
params.Analyses.SingleCell      = 1;
params.Analyses.Regularity      = 0;
params.Analyses.Catch22         = 0;
params.Analyses.Bursts          = 1;
params.Analyses.Connectivity    = [];
params.Outlier.Method           = []; %No outlier removal

lookup_path = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx";
path_part_idx = [11,12,13];
min_N_units = params.QC.N_Units;
parallel = false;

%% Randomly pick one sorting and check the burst detection (rerun to check more sortings/test out different parameter sets)
path_list = randsample(sorting_path_list,1);
params.Bursts.ISI_N         = 1;  %Change if needed
params.Bursts.Skewness_TH   = 2;  %Change if needed
params.Bursts.Merge_t       = 2;  %Change if needed

test_rec = generate_MEArecordings_from_sorting_list(path_list, lookup_path, path_part_idx, min_N_units, params, parallel);
test_rec.PlotBurstCheck();