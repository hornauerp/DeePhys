function failed_sortings = generate_MEArecordings_from_sorting_list(sorting_path_list, metadata, params, parallel)
arguments
    sorting_path_list
    metadata struct
    params = struct()
    parallel logical = true
end
warning('off')
failed_sortings = [];

if parallel
    parfor iPath = 1:length(sorting_path_list)
        save_file = fullfile(sorting_path_list(iPath),'MEArecording.mat');
        md = metadata;
        if ~exist(save_file,'file') || params.Save.Overwrite
            md.InputPath = sorting_path_list(iPath);
            temp_file = fullfile(md.InputPath,"spike_templates.npy");
            if exist(temp_file,"file")
                spk_temp = readNPY(temp_file);
                if length(unique(spk_temp)) > params.QC.N_Units
                    try
                        MEArecording(md, params);
                    catch
                        disp("Failed for sorting: " + sorting_path_list(iPath))
                        failed_sortings = [failed_sortings, sorting_path_list(iPath)];
                    end
                else
                    disp("Not enough units for: " + sorting_path_list(iPath))
                end
            else
                disp("No sorting file found for: " + sorting_path_list(iPath))
            end
            disp("Finished recording " + num2str(iPath))
        end
    end
else
    for iPath = 1:length(sorting_path_list)
        save_file = fullfile(sorting_path_list(iPath),'MEArecording.mat');
        if ~exist(save_file,'file') || params.Save.Overwrite
            metadata.InputPath = sorting_path_list(iPath);
            temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
            if exist(temp_file,"file")
                spk_temp = readNPY(temp_file);
                if length(unique(spk_temp)) > params.QC.N_Units
                    try
                        MEArecording(metadata, params);
                    catch
                        disp("Failed for sorting: " + sorting_path_list(iPath))
                        failed_sortings = [failed_sortings, sorting_path_list(iPath)];
                    end
                else
                    disp("Not enough units for: " + sorting_path_list(iPath))
                end
            else
                disp("No sorting file found for: " + sorting_path_list(iPath))
            end
            disp("Finished recording " + num2str(iPath))
        end
    end
end
end