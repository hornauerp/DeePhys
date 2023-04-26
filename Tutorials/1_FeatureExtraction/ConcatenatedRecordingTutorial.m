%% Feature extraction of concatenated recordings
% This code serves as a template to extract features from a list of
% sortings that had been concatenated prior to spike sorting. We have provided two sortings and one lookup table with metadata to showcase
% the general approach. 
% 
addpath(genpath("/home/phornauer/Git/DeePhys")) %Add your personal path to the DeePhys folder

%% Generate sorting_path_list automatically
% Assuming that you stored the sortings according to some tree logic, we
% can automatically generate a list with all sorting paths
% Alternatively, you can also supply a manually curated sorting_path_list

root_path = "/home/phornauer/Git/DeePhys/Data/SortingExamples/";        %Common path of all sorting paths
path_logic = {'200*','52*','well000','sorted'};                         %Each cell corresponds to one subdirectory // * represents a wildcard character
sorting_path_list = generate_sorting_path_list(root_path, path_logic);
fprintf("Generated %i sorting paths\n",length(sorting_path_list))


%% Set inference parameters 
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


%% Set metadata information MOST LIKELY NEEDS TO BE CHANGED

split_times = [0 2; 4 6; 8 10]; %Cutouts of the recording // specify in pairs of starting time and end time (:,2) 
output_folder = "min_" + split_times(:,1); %Name the folders according to some logic, here we use the minute the cutout starts
lookup_path = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx"; %Path where the lookup tabel containing metainformation is stored
path_part_index = [8,9,10]; %These numbers need to correspond to the index in the path after splitting at / [RecordingDate, PlateID, WellID]
min_N_units = 20; %Minimum number of spike-sorted units for the feature inference to be run for a recording
acute_treatment = "Quinpirole"; %Change to your drug of choice

%% Now we run the main loop
% You should not need to change anything here

rec_array = {};
parfor iPath = 1:length(sorting_path_list)
    sorting_path = sorting_path_list(iPath);
    spkt_file = fullfile(sorting_path,'spike_times.npy'); %Make sure that we have the sorted data 
    if exist(spkt_file,'file')
        spk_t = readNPY(spkt_file);
        output_folder = "min_" + split_times(:,1); %Name folders according to the minute the cutout starts
        output_path_list = split_sortings(sorting_path, split_times, output_folder); %We split the recording into subfolders and distribute the sorted data accordingly
        
        for iSplit = 1:numel(output_path_list)
            
            split_params = params;
            if iSplit == 1 %Perform QC on first part, keep the same units for the rest
                split_params.QC.GoodUnits = [];
            else
                split_params.QC.GoodUnits = mearec.Parameters.QC.GoodUnits; %Currently only supported in boolean/logic indexing
            end
            
            metadata = struct();
            metadata.LookupPath = lookup_path;
            metadata.InputPath = output_path_list(iSplit);
            
            metadata.PathPartIndex = path_part_index;
            metadata.Timepoint = split_times(iSplit,1);
            metadata.AcuteTreatment = acute_treatment;
            temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
            if exist(temp_file,"file") %Make sure that the splitting worked and we have the sorted data
                spk_temp = readNPY(temp_file);
                if length(unique(spk_temp)) > min_N_units %Here, we check for the number of units before the quality control
                    mearec = MEArecording(metadata, split_params); %Main algorithm
                else
                    break
                end
            end
            rec_array = [rec_array mearec];
        end
    else
        continue
    end
end

%% Optional: Remove low unit recordings after the QC
filtered_array = remove_low_unit_recordings(rec_array, min_N_units);
rg = RecordingGroup(filtered_array);

% And we have the recording group we can use to perform all subsequent
% analyses!
