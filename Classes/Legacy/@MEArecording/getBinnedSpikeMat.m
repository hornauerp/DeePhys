function binned_mat = getBinnedSpikeMat(obj, bin_size, time_window)
% GETBINNEDSPIKMAT  Bin spike times for all units into a units × time-bins matrix.
%
% INPUTS:
%   obj         - MEArecording object
%   bin_size    - Bin width in seconds
%   time_window - (optional) [start, end] in seconds; default [0, recording duration]
%
% OUTPUTS:
%   binned_mat  - (N_units x N_bins) spike count matrix

arguments
    obj         MEArecording
    bin_size    (1,1) double
    time_window (1,2) double = [0, obj.RecordingInfo.Duration]
end

n_units    = length(obj.Units);
bins       = time_window(1) : bin_size : time_window(2);
binned_mat = zeros(n_units, length(bins) - 1);

for u = 1:n_units
    spike_times = obj.Units(u).SpikeTimes;
    in_window   = spike_times >= time_window(1) & spike_times < time_window(2);
    binned_mat(u,:) = histcounts(spike_times(in_window), bins);
end
end
