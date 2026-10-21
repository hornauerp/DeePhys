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
%   sttc_params           - struct: .MaxLag, .SurrogateMethod, .N_Surrogates, .EdgeTest
%                           ("fdr", default, with .FDR; or "percentile" with .Percentile),
%                           .JitterSD (s, 'jitter' surrogates; default 2 x MaxLag)
%
% OUTPUTS:
%   connectivity - struct with .wu (weighted undirected), .bu (binary undirected), .threshold
%                  ("fdr" also: .p one-sided surrogate p-values, .q BH-adjusted)

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

    % Ascending order: each insertion assumes the units before su are already in place
    for su = silent_units(:)'
        spks_data.raster = [spks_data.raster(1:su-1); cell(1,1); spks_data.raster(su:end)];
    end
    spks_data.nchannels = length(spks_data.raster);

    % Jitter width of the surrogates: it must exceed the coincidence window MaxLag, or the
    % jittered spikes keep their +-MaxLag coincidences and the null equals the data (the former
    % fixed 10 ms jitter left almost no pair significant at MaxLag = 50 ms).
    jitter_sd = 2 * sttc_params.MaxLag;
    if isfield(sttc_params, 'JitterSD') && ~isempty(sttc_params.JitterSD)
        jitter_sd = sttc_params.JitterSD;
        if strcmp(sttc_params.SurrogateMethod, 'jitter') && jitter_sd <= sttc_params.MaxLag
            warning('computeConnectivitySTTC:jitterTooSmall', ...
                ['STTC.JitterSD (%g s) does not exceed MaxLag (%g s): the surrogates keep most ' ...
                 'coincidences and few pairs can be significant. Leave JitterSD empty (2 x MaxLag).'], ...
                jitter_sd, sttc_params.MaxLag);
        end
    end
    rand_args = {};
    if strcmp(sttc_params.SurrogateMethod, 'jitter')
        rand_args = {'stdtime', jitter_sd * 1000};   % randomizeasdf2 takes ms
    end
    for s = 1:sttc_params.N_Surrogates
        spks_rand = randomizeasdf2(spks_data, sttc_params.SurrogateMethod, rand_args{:});
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

    % Edge test. Each pair has its own surrogate null; testing every pair at a fixed
    % percentile leaves (100 - Percentile)% false edges among independent pairs, which then
    % drive density, clustering and modularity. EdgeTest "fdr" (default) z-scores each pair
    % against its surrogates (one-sided: positive synchrony) and keeps the pairs that pass
    % Benjamini-Hochberg control at STTC.FDR across all tested pairs; "percentile" keeps the
    % former per-pair percentile threshold.
    edge_test = "fdr";
    if isfield(sttc_params, 'EdgeTest') && ~isempty(sttc_params.EdgeTest)
        edge_test = lower(string(sttc_params.EdgeTest));
    end
    surr_mu = mean(surrogate_mat, 3, 'omitnan');
    surr_sd = std(surrogate_mat, 0, 3, 'omitnan');
    connectivity.wu = empirical_sttc;
    switch edge_test
        case "fdr"
            q_level = 0.05;
            if isfield(sttc_params, 'FDR') && ~isempty(sttc_params.FDR), q_level = sttc_params.FDR; end
            z = (empirical_sttc - surr_mu) ./ surr_sd;
            z(surr_sd == 0 & empirical_sttc > surr_mu) = Inf;    % constant null, data above it
            z(surr_sd == 0 & empirical_sttc <= surr_mu) = -Inf;
            p = 1 - normcdf(z);
            ut = triu(true(n_units), 1) & ~isnan(z);
            q = NaN(n_units);
            q(ut) = StatisticalTest.fdrCorrection(p(ut));
            q = min(q, q.');                                     % symmetric (one value per pair)
            connectivity.p  = min(p, p.');
            connectivity.q  = q;
            connectivity.bu = q <= q_level;
            % Smallest STTC that would pass for each pair (diagnostic; Inf when no pair passes)
            p_cut = max([p(ut & q <= q_level); 0]);
            connectivity.threshold = surr_mu + norminv(1 - p_cut) .* surr_sd;
            connectivity.threshold(isnan(z)) = NaN;
        case "percentile"
            connectivity.threshold = prctile(surrogate_mat, sttc_params.Percentile, 3);
            connectivity.bu = connectivity.wu > connectivity.threshold;
        otherwise
            error('computeConnectivitySTTC:unknownEdgeTest', ...
                'Unknown STTC.EdgeTest "%s" (use "fdr" or "percentile").', edge_test);
    end
    connectivity.bu(1:n_units+1:end) = false;
end
