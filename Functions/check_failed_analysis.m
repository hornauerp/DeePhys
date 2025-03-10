function check_failed_analysis(sorting_path_list, metadata, params)
for s = 1:length(sorting_path_list)
    metadata.InputPath = sorting_path_list(s);
    MEArecording(metadata, params);
end