function [clean_spikes, clean_units] = removeTonic(spike_times, spike_units, firing_rates, cv_isi)
% REMOVETONIC  Remove tonic-firing units from spike data.
%
%   [clean_spikes, clean_units] = removeTonic(spike_times, spike_units, firing_rates, cv_isi)
%
% Tonic units are those with outlier-high firing rates AND below-average CV(ISI).
% This replicates the logic from MEArecording.removeTonic but accepts pre-computed
% activity features instead of requiring a MEArecording object.
%
% INPUTS:
%   spike_times   - (N x 1) double, sorted spike times in seconds
%   spike_units   - (N x 1) double, unit assignment per spike
%   firing_rates  - (M x 1) double, firing rate per unit
%   cv_isi        - (M x 1) double, CV of ISI per unit
%
% OUTPUTS:
%   clean_spikes  - spike times with tonic units removed
%   clean_units   - re-grouped unit assignments for remaining spikes

    arguments
        spike_times   (:,1) double
        spike_units   (:,1) double
        firing_rates  (:,1) double
        cv_isi        (:,1) double
    end

    tonic_ids = find((isoutlier(firing_rates, "ThresholdFactor", 1) & ...
        firing_rates > median(firing_rates)) & ...
        cv_isi < mean(cv_isi));

    clean_idx = ~ismember(spike_units, tonic_ids);
    clean_units = findgroups(spike_units(clean_idx));
    clean_spikes = spike_times(clean_idx);
end
