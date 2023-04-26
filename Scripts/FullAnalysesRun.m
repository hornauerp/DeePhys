addpath(genpath("/home/phornauer/Git/DeePhys"))
%% Generate sorting_path_list 
% root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/mpriouret/iNeurons/";
% path_logic = {'SCR*','*','*','w*'}; %Each cell corresponds to one subdirectory
% path_logic = {'DeePhysS*','*','*','w*'};
% path_logic = {'M*','w*'};
path_logic = {'230302','Amplitude','*','w*','sorted'};

sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Set parameters values
emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();
params.QC.Amplitude             = [];
params.QC.FiringRate            = [0.01 100];
params.QC.Axon                  = 0.8;
params.QC.Noise                 = 8;
params.QC.N_Units               = 20;
params.Analyses.SingleCell      = 1;
params.Analyses.Regularity      = 1;
params.Analyses.Catch22         = 1;
params.Analyses.Bursts          = 1;
params.Analyses.Connectivity    = ["STTC","CCG"];
params.Outlier.Method           = []; %No outlier removal
params.Save.Flag                = 1; %Save individual MEArecordings to prevent data loss if the execution is interrupted

lookup_path = "/home/phornauer/Data/iNeurons/iNeurons_lookup.xlsx";
path_part_idx = [11,13,14]; %[recording_date, plate_id, well_id]
min_N_units = params.QC.N_Units;
parallel = true;
if parallel
    parpool(6);
end

%% Run full loop
rec_array = generate_MEArecordings_from_sorting_list(sorting_path_list, lookup_path, path_part_idx, min_N_units, params, parallel);
save("/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Maxtwo/mpriouret/iNeurons/230302/amplitude_based.mat", "rec_array")