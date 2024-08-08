function rec_array = generate_MEArecordings_from_sorting_list(sorting_path_list, lookup_path, path_part_idx, min_N_units, params, parallel)
arguments
    sorting_path_list
    lookup_path %(1,1) string = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx"
    path_part_idx (1,:) double % = [11,12,13]
    min_N_units (1,1) double = 10
    params = struct()
    parallel logical = true
end

rec_array = {};
if parallel
    parfor (iPath = 1:length(sorting_path_list),24)
        save_file = fullfile(sorting_path_list(iPath),'MEArecording.mat');
        if ~exist(save_file,'file') || params.Save.Overwrite
            metadata = struct();
            metadata.LookupPath = lookup_path;
            metadata.InputPath = sorting_path_list(iPath);
            metadata.PathPartIndex = path_part_idx;
            temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
            if exist(temp_file,"file")
                spk_temp = readNPY(temp_file);
                if length(unique(spk_temp)) > min_N_units
                    rec = MEArecording(metadata, params);
                    rec_array{iPath} = rec;
                else
                    disp("Not enough units for: " + sorting_path_list(iPath))
                end
            else
                disp("No sorting file found for: " + sorting_path_list(iPath))
            end
            disp("Finished recording " + num2str(iPath))
        else
            rec_array{iPath} = load(save_file,"obj");
            disp("Loading existing object")
        end
    end
else
    for iPath = 1:length(sorting_path_list)
        save_file = fullfile(sorting_path_list(iPath),'MEArecording.mat');
        if ~exist(save_file,'file') || params.Save.Overwrite
            metadata = struct();
            metadata.LookupPath = lookup_path;
            metadata.InputPath = sorting_path_list(iPath);
            metadata.PathPartIndex = path_part_idx;
            temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
            if exist(temp_file,"file")
                spk_temp = readNPY(temp_file);
                if length(unique(spk_temp)) > min_N_units
                    rec_array{iPath} = MEArecording(metadata, params);
                else
                    disp("Not enough units for: " + sorting_path_list(iPath))
                end
            else
                disp("No sorting file found for: " + sorting_path_list(iPath))
            end
            disp("Finished recording " + num2str(iPath))
        else
            rec_array{iPath} = load(save_file,"obj");
            disp("Loading existing object")
        end
    end
end
rec_array = [rec_array{:}];
end