function check_failed_analysis(sorting_path_list, metadata, params)
% CHECK_FAILED_ANALYSIS  Re-run MEArecording construction for a list of sorting paths.
%
% Iterates over sorting_path_list and constructs an MEArecording for each path,
% useful for retrying previously failed analyses.
%
% INPUTS:
%   sorting_path_list - String array of sorting output folder paths
%   metadata          - Struct with recording metadata (InputPath will be set per path)
%   params            - Parameters struct passed to MEArecording

for s = 1:length(sorting_path_list)
    metadata.InputPath = sorting_path_list(s);
    MEArecording(metadata, params);
end
end
