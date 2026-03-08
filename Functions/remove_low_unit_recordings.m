function filtered_array = remove_low_unit_recordings(recording_array, unit_threshold)
% REMOVE_LOW_UNIT_RECORDINGS  Filter out recordings with too few units.
%
% INPUTS:
%   recording_array - Array of MEArecording objects
%   unit_threshold  - Minimum number of units required (recordings with <= this value are removed)
%
% OUTPUTS:
%   filtered_array  - Subset of recording_array with unit count > unit_threshold

N_units = arrayfun(@(x) length(recording_array(x).Units),1:length(recording_array));
filtered_array = recording_array(N_units > unit_threshold);
fprintf("Kept %i out of %i recordings\n", length(filtered_array),length(recording_array))
end
