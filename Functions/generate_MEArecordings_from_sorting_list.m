function rec_array = generate_MEArecordings_from_sorting_list(sorting_path_list, lookup_path, path_part_idx, min_N_units, params, parallel)
arguments
    sorting_path_list
    lookup_path %(1,1) string = "/home/phornauer/Git/DeePhys/Data/cellline_lookup.xlsx"
    path_part_idx (1,:) double = [11,12,13]
    min_N_units (1,1) double = 10
    params = struct()
    parallel logical = true
end

rec_array = {};
if parallel
    parfor (iPath = 1:length(sorting_path_list),24)
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
            end
        end
    end
else
    for iPath = 1:length(sorting_path_list)
        metadata = struct();
        metadata.LookupPath = lookup_path;
        metadata.InputPath = sorting_path_list(iPath);
        metadata.PathPartIndex = path_part_idx;
        temp_file = fullfile(metadata.InputPath,"spike_templates.npy");
        if exist(temp_file,"file")
            spk_temp = readNPY(temp_file);
            if length(unique(spk_temp)) > min_N_units
                rec_array{iPath} = MEArecording(metadata, params);
            end
        end
    end
end
rec_array = [rec_array{:}];
end