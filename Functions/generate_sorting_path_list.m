function sorting_path_list = generate_sorting_path_list(root_path, path_logic, split_prefix)
arguments
    root_path string %Common denominator of all sorting paths
    path_logic cell = {'*','*','w*'} %Expression describing the subfolder structure
    split_prefix string = []% only use if a concatenated sorting was previously split by the split_sortings() function
end
% dir() returns the matched entries for a wildcard last element, but the CONTENTS of the
% matched folders for a literal one (e.g. 'sorter_output'); the "." entry marks the latter.
% A bare '*' also lists "." and "..", so the last element itself decides which case applies.
entries = dir(fullfile(root_path, path_logic{:}));
is_dot = ismember({entries.name}, {'.', '..'});
last_is_wildcard = ~isempty(path_logic) && contains(string(path_logic{end}), ["*", "?"]);
if last_is_wildcard
    listed = strings(0, 1);
    rest = entries(~is_dot);                                                  % the matches
else
    listed = unique(string({entries(strcmp({entries.name}, '.')).folder}));  % literal folder
    rest = entries(~is_dot & ~ismember(string({entries.folder}), listed));   % literal file
end
matched = string(arrayfun(@(e) fullfile(e.folder, e.name), rest(:), 'UniformOutput', false));
matched_dirs = matched([rest.isdir]);
matched_files = string({rest(~[rest.isdir]).folder});                     % a file pattern: its folder
sorted_dirs = unique([listed(:); matched_dirs(:); matched_files(:)])';
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
