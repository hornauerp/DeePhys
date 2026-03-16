function failed_sortings = generate_MEArecordings_from_sorting_list(sorting_path_list, metadata, params, parallel)
% GENERATE_MEARECORDINGS_FROM_SORTING_LIST  Construct MEArecording objects from a list of sorting paths.
%
% Skips paths where a MEArecording.mat already exists (unless params.Save.Overwrite
% is true). Optionally processes paths in parallel via parfor.
%
% INPUTS:
%   sorting_path_list - String array of paths to sorting output folders
%   metadata          - Struct with recording metadata (InputPath will be set per path)
%   params            - (optional) Parameters struct passed to MEArecording; default struct()
%   parallel          - (optional) If true, uses parfor; default true
%
% OUTPUTS:
%   failed_sortings   - String array of paths that raised an error during construction

arguments
    sorting_path_list
    metadata struct
    params = struct()
    parallel logical = true
end
prev_warn_state = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
cleanup_warnings = onCleanup(@() warning(prev_warn_state));
failed_sortings = [];

% Pre-read the lookup table once to avoid N redundant xlsx reads in parfor
if isfield(metadata, 'LookupPath') && ~isempty(metadata.LookupPath) && ...
   ~isfield(metadata, 'MetadataTable')
    metadata.MetadataTable = readtable(metadata.LookupPath, ...
        "ReadVariableNames", true, "TextType", "string");
end

if parallel
    parfor iPath = 1:length(sorting_path_list)
        md = metadata;
        failed = processOneSorting(sorting_path_list(iPath), md, params, iPath);
        failed_sortings = [failed_sortings, failed]; %#ok<AGROW>
    end
else
    for iPath = 1:length(sorting_path_list)
        failed = processOneSorting(sorting_path_list(iPath), metadata, params, iPath);
        failed_sortings = [failed_sortings, failed]; %#ok<AGROW>
    end
end
end


function failed = processOneSorting(path, metadata, params, iPath)
% PROCESSONSORTING  Attempt to construct an MEArecording for a single sorting path.
%   Returns the path as a string if construction fails, or string.empty on success.
failed = string.empty;
save_file = fullfile(path, 'MEArecording.mat');
if ~exist(save_file, 'file') || params.Save.Overwrite
    metadata.InputPath = path;
    temp_file = fullfile(metadata.InputPath, "spike_templates.npy");
    if exist(temp_file, "file")
        spk_temp = readNPY(temp_file);
        if length(unique(spk_temp)) > params.QC.N_Units
            try
                MEArecording(metadata, params);
            catch
                disp("Failed for sorting: " + path)
                failed = path;
            end
        else
            disp("Not enough units for: " + path)
        end
    else
        disp("No sorting file found for: " + path)
    end
    disp("Finished recording " + num2str(iPath))
end
end
