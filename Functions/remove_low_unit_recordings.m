function filtered_array = remove_low_unit_recordings(recording_array, unit_threshold)
N_units = arrayfun(@(x) length(recording_array(x).Units),1:length(recording_array));
filtered_array = recording_array(N_units > unit_threshold);
fprintf("Kept %i out of %i units\n", length(filtered_array),length(recording_array))