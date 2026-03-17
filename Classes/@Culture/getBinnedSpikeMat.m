function binned_mat = getBinnedSpikeMat(culture, bin_size, time_window)
% GETBINNEDSPIKMAT  Bin spike times across all recordings into a units × time-bins matrix.
%
% Recordings are concatenated in order; each recording's spike times are offset
% by the cumulative duration of prior recordings. Bin size and analysis window
% are specified over the concatenated timeline.
%
% INPUTS:
%   culture     - Culture object
%   bin_size    - Bin width in seconds
%   time_window - (optional) [start, end] in seconds over the concatenated timeline;
%                 default [0, total duration]
%
% OUTPUTS:
%   binned_mat  - (N_units x N_bins) spike count matrix

arguments
    culture    Culture
    bin_size   (1,1) double
    time_window (1,2) double = [0, sum(arrayfun(@(r) r.RecordingInfo.Duration, culture.Recordings))]
end

baseline_units = culture.Units;           % canonical unit order from Recordings(1)
n_units    = length(baseline_units);
if n_units == 0
    binned_mat = zeros(0, max(0, floor((diff(time_window))/bin_size)));
    return
end
% Use TemplateID for matching — unitID is fragile after save/load
baseline_ids = [baseline_units.TemplateID];
bins       = time_window(1) : bin_size : time_window(2);
binned_mat = zeros(n_units, length(bins) - 1);
offset     = 0;

for r = 1:culture.N
    rec        = culture.Recordings(r);
    rec_ids    = [rec.Units.TemplateID];
    [~, row_idx] = ismember(baseline_ids, rec_ids);  % map baseline order → rec order
    for u = 1:n_units
        if row_idx(u) == 0; continue; end            % unit absent in this recording
        spike_times = rec.Units(row_idx(u)).SpikeTimes + offset;
        in_window   = spike_times >= time_window(1) & spike_times < time_window(2);
        binned_mat(u,:) = binned_mat(u,:) + histcounts(spike_times(in_window), bins);
    end
    offset = offset + rec.RecordingInfo.Duration;
end
end
