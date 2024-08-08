function sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix)
arguments
    root_path string %Common denominator of all sorting paths
    path_logic cell = {'*','*','w*'} %Expression describing the subfolder structure
    split_prefix string = []% only use if a concatenated sorting was previously split by the split_sortings() function
end
sorted_dirs = dir(fullfile(root_path, path_logic{:}));
sorted_dirs = string(unique({sorted_dirs.folder}));
sorting_path_list = [];
if isempty(split_prefix)
    sorting_path_list = sorted_dirs;
else
    for i = 1:length(sorted_dirs)
        dirs = dir(fullfile(sorted_dirs(i), split_prefix));
        sorting_path_list = [sorting_path_list; arrayfun(@(x) fullfile(string(x.folder), x.name),dirs)];
    end
end

end