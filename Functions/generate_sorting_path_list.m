function sorting_path_list = generate_sorting_path_list(root_path, path_logic)
arguments
   root_path string %Common denominator of all sorting paths
   path_logic string = ["*","*","w*"] %Expression describing the subfolder structure
end

sorted_dirs = dir(fullfile(root_path, path_logic(:)));

sorting_path_list = arrayfun(@(x) fullfile(string(sorted_dirs(x).folder), sorted_dirs(x).name, "sorted"), 1:length(sorted_dirs));

end