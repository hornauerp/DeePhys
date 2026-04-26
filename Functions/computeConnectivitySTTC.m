function connectivity = computeConnectivitySTTC(unit_spike_times_cell, spike_times, ...
    spike_units, n_units, duration, sampling_rate, sttc_params)
% COMPUTECONNECTIVITYSTTC  Compute STTC-based functional connectivity.
%
%   connectivity = computeConnectivitySTTC(unit_spike_times_cell, spike_times, ...
%       spike_units, n_units, duration, sampling_rate, sttc_params)
%
% INPUTS:
%   unit_spike_times_cell - {1 x M} cell array, spike time vectors per unit
%   spike_times           - (N x 1) double, all spike times
%   spike_units           - (N x 1) double, unit assignment per spike
%   n_units               - (1 x 1) double, number of units
%   duration              - (1 x 1) double, recording duration in seconds
%   sampling_rate         - (1 x 1) double, sampling rate in Hz
%   sttc_params           - struct: .MaxLag, .SurrogateMethod, .N_Surrogates, .Percentile
%
% OUTPUTS:
%   connectivity - struct with .wu (weighted undirected), .bu (binary undirected), .threshold

    arguments
        unit_spike_times_cell cell
        spike_times           (:,1) double
        spike_units           (:,1) double
        n_units               (1,1) double
        duration              (1,1) double
        sampling_rate         (1,1) double
        sttc_params           struct
    end

    if isempty(spike_times)
        connectivity.wu = zeros(n_units);
        connectivity.bu = zeros(n_units);
        return
    end

    empirical_sttc = zeros(n_units);
    for i = 1:n_units-1
        N1v = length(unit_spike_times_cell{i});
        for j = (i+1):n_units
            N2v = length(unit_spike_times_cell{j});
            empirical_sttc(i,j) = sttc(N1v, N2v, sttc_params.MaxLag, ...
                [0 duration], unit_spike_times_cell{i}, unit_spike_times_cell{j});
        end
    end
    empirical_sttc = empirical_sttc + triu(empirical_sttc, -1).';

    surrogate_mat = nan([size(empirical_sttc), sttc_params.N_Surrogates]);
    silent_units = find(cellfun(@isempty, unit_spike_times_cell));
    spks_data = convert2asdf2(spike_times, spike_units, sampling_rate);

    for su = fliplr(silent_units)
        spks_data.raster = [spks_data.raster(1:su-1); cell(1,1); spks_data.raster(su:end)];
    end
    spks_data.nchannels = length(spks_data.raster);

    for s = 1:sttc_params.N_Surrogates
        spks_rand = randomizeasdf2(spks_data, sttc_params.SurrogateMethod);
        surrogate_sttc = zeros(spks_data.nchannels);
        for i = 1:(spks_data.nchannels - 1)
            spike_times_1 = double(spks_rand.raster{i}) / sampling_rate;
            N1v = length(spike_times_1);
            for j = (i+1):spks_data.nchannels
                spike_times_2 = double(spks_rand.raster{j}) / sampling_rate;
                N2v = length(spike_times_2);
                surrogate_sttc(i,j) = sttc(N1v, N2v, sttc_params.MaxLag, ...
                    [0 duration], spike_times_1, spike_times_2);
            end
        end
        surrogate_sttc = surrogate_sttc + triu(surrogate_sttc, -1).';
        surrogate_mat(:,:,s) = surrogate_sttc;
    end

    connectivity.threshold = prctile(surrogate_mat, sttc_params.Percentile, 3);
    connectivity.wu = empirical_sttc;
    connectivity.bu = connectivity.wu > connectivity.threshold;
end
