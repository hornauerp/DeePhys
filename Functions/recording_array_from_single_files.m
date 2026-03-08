function rec_array = recording_array_from_single_files(sorting_path_list, rm_con)
% RECORDING_ARRAY_FROM_SINGLE_FILES  Load MEArecording objects from a list of sorting paths.
%
% Uses parfor to load MEArecording.mat files in parallel. Skips paths where
% no .mat file exists and removes recordings that have no units.
%
% INPUTS:
%   sorting_path_list - String array of paths to sorting folders
%   rm_con            - (optional) If true, clears Connectivity data after loading to save
%                       memory; default true
%
% OUTPUTS:
%   rec_array - Array of MEArecording objects

arguments
    sorting_path_list string
    rm_con logical = true %If true, removes connectivity data to save memory
end
rec_array = cell(1,length(sorting_path_list));
parfor iPath = 1:length(sorting_path_list)
    mearec = fullfile(sorting_path_list(iPath), "MEArecording.mat");
    if exist(mearec,'file')
        opened = load(mearec,"obj");
        if rm_con
            opened.obj.Connectivity = [];
        end
        rec_array{iPath} = opened.obj;
    end
end
rec_array = [rec_array{:}];
empty_idx = arrayfun(@(x) isempty(x.Units),rec_array); %Recordings with no units
rec_array(empty_idx) = []; %Remove them
fprintf("%i/%i objects retrieved\n", length(rec_array),length(sorting_path_list))
end
