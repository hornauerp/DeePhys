addpath(genpath("/home/phornauer/Git/DeePhys"))

%% Generate sorting_path_list automatically
root_path = "/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/phornauer/";
path_logic = {'230123','How_*','M*','w*','sorted'}; %Each cell corresponds to one subdirectory
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
params.Save.Overwrite           = 1; %Overwrite previous analysis if you changed parameters

split_times = [0 15; 15 30; 30 45; 45 60; 60 75]; %Cutouts of the recording 
timepoints = -1:3;
output_folder = "day_" + timepoints; %Name folders according to cutout starts

min_N_units = 20; %Minimum number of units

%%

parfor iPath = 1:length(sorting_path_list)
    sorting_path = sorting_path_list(iPath);
    spkt_file = fullfile(sorting_path,'spike_templates.npy');
    if exist(spkt_file,'file') && length(readNPY(unique(spkt_file))) > 20
        output_path_list = split_sortings(sorting_path, split_times, output_folder);
        output_path_list = output_path_list([2,1,3,4,5]);
        for iSplit = 1:numel(output_path_list)
            
            split_params = params;
            if iSplit == 1 %Perform QC on first part, keep the same units for the rest
                split_params.QC.GoodUnits = [];
            else
                split_params.QC.GoodUnits = mearec.Parameters.QC.GoodUnits; %Currently only supported in boolean/logic indexing
            end
            
            metadata = struct();
            metadata.LookupPath = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx";
            metadata.InputPath = output_path_list(iSplit);
            
            metadata.PathPartIndex = [10,12,13]; 
            metadata.MediumChangeTimepoint = timepoints(iSplit);
            temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
            mearec_file = fullfile(metadata.InputPath,"MEArecording.mat");
            if exist(temp_file,"file") %&& ~exist(fullfile(metadata.InputPath, 'MEArecording.mat'),'file') 
                mearec = MEArecording(metadata, split_params);
%             elseif exist(fullfile(metadata.InputPath, 'MEArecording.mat'),'file')
%                 loaded_obj = load(fullfile(metadata.InputPath, 'MEArecording.mat'));
%                 mearec = loaded_obj.obj;
            else
                break
            end
        end
    else
        continue
    end
end