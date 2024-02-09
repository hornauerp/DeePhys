%% General feature extraction
% This code serves as a template to extract features from a list of
% sortings. We have provided two sortings and one lookup table with metadata to showcase
% the general approach. 
% 

%% Add the DeePhys code to the search path
addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list automatically
% Assuming that you stored the sortings according to some tree logic, we
% can automatically generate a list with all sorting paths
% Alternatively, you can also supply a manually curated sorting_path_list

root_path = "/home/phornauer/Git/DeePhys/Data/SortingExamples/";        %Common path of all sorting paths
path_logic = {'200*','52*','well000','sorted'};                         %Each cell corresponds to one subdirectory // * represents a wildcard character
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))

%% Set inference parameters (can be left at default)
%Leaving any of these empty will skip that inference/QC

emptyRec = MEArecording();
params = emptyRec.returnDefaultParams();            %For a full list of available parameters check the params variable

% List of commonly used/changed paramters
params.QC.Amplitude             = [];               %Amplitude minimum and maximum
params.QC.FiringRate            = [0.01 100];       %Firing rate minimum and maximum
params.QC.Axon                  = 0.8;              %Maximum ratio between positive and negative peaks allowed
params.QC.Noise                 = 8;                %Maximum number of changes in the first 
params.QC.N_Units               = 20;
params.Analyses.SingleCell      = 1;                %Run single-cell analysis (always recommended)
params.Analyses.Regularity      = 1;                %Infer regularity features
params.Analyses.Catch22         = 1;                %Infer Catch22 features
params.Analyses.Bursts          = 1;                %Infer burst features
params.Analyses.Connectivity    = ["STTC","CCG"];   %Connectivity inferece algorithms to be run (string array of method names)
params.Outlier.Method           = [];               %No outlier removal if no method is specified
params.Save.Flag                = 1;                %Save individual MEArecordings to prevent data loss if the execution is interrupted (recommended)
params.Save.Overwrite           = 1;                %Overwrite previous analysis if you changed parameters


%% Parameters that most likely NEED TO BE CHANGED 
path_part_idx                   = [8,9,10];         %These numbers need to correspond to the index in the path after splitting at / [RecordingDate, PlateID, WellID]
lookup_path                     = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx"; %Path where the lookup tabel containing metainformation is stored
N_cores                         = 6;                %Parallelize feature extraction (recommended)

if N_cores > 1
    parallel = true;
    parpool(N_cores);
else 
    parallel = false;
end

%% Run full loop
rec_array = generate_MEArecordings_from_sorting_list(sorting_path_list, lookup_path, path_part_idx, params.QC.N_Units, params, parallel);

save("rec_array.mat", "rec_array")
% And we have the recording group we can use to perform all subsequent
% analyses!
