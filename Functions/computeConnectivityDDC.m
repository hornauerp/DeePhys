function connectivity = computeConnectivityDDC(spike_times, spike_units, n_units, duration, ddc_params)
% COMPUTECONNECTIVITYDDC  Compute DDC (derivative cross-correlation) connectivity.
%
%   connectivity = computeConnectivityDDC(spike_times, spike_units, n_units, duration, ddc_params)
%
% INPUTS:
%   spike_times  - (N x 1) double, sorted spike times in seconds
%   spike_units  - (N x 1) double, unit assignment per spike
%   n_units      - (1 x 1) double, number of units
%   duration     - (1 x 1) double, recording duration in seconds
%   ddc_params   - struct: .BinSize, .Threshold, .N_Surrogates (default 100),
%                  .Percentile (default 95)
%
% OUTPUTS:
%   connectivity - struct with .wu (weighted directed matrix),
%                  .bd (binary directed, surrogate-thresholded), .threshold

    arguments
        spike_times  (:,1) double
        spike_units  (:,1) double
        n_units      (1,1) double
        duration     (1,1) double
        ddc_params   struct
    end

    if ~isfield(ddc_params, 'N_Surrogates')
        ddc_params.N_Surrogates = 100;
    end
    if ~isfield(ddc_params, 'Percentile')
        ddc_params.Percentile = 95;
    end

    bin_size = ddc_params.BinSize;
    n_bins = ceil(duration / bin_size);
    V_obs = zeros(n_bins, n_units);
    spiketimes_binned = max(1, min(n_bins, round(spike_times / bin_size) + 1));
    V_obs(sub2ind(size(V_obs), spiketimes_binned, spike_units)) = 1;

    connectivity.wu = ddcCore(V_obs, ddc_params.Threshold);

    if all(isnan(connectivity.wu(:)))
        connectivity.bd = zeros(n_units);
        connectivity.threshold = NaN(n_units);
        return
    end

    % Surrogate test: circularly shift each unit's binned activity independently
    surrogate_mat = zeros(n_units, n_units, ddc_params.N_Surrogates);
    for s = 1:ddc_params.N_Surrogates
        V_shuf = V_obs;
        for u = 1:n_units
            V_shuf(:, u) = circshift(V_obs(:, u), randi(n_bins));
        end
        surrogate_mat(:,:,s) = ddcCore(V_shuf, ddc_params.Threshold);
    end

    connectivity.threshold = prctile(abs(surrogate_mat), ddc_params.Percentile, 3);
    bd = double(abs(connectivity.wu) > connectivity.threshold);
    bd(1:n_units+1:end) = 0;
    connectivity.bd = bd;
end


function wu = ddcCore(V_obs, relu_threshold)
    TR = 1;
    [~, N] = size(V_obs);
    Fx = V_obs - relu_threshold;
    Fx(Fx < 0) = 0;
    tmp = cov([Fx V_obs]);
    B = tmp(1:N, N+1:end);
    dV = (-1/2*V_obs(1:end-2,:) + 1/2*V_obs(3:end,:)) / TR;
    dV = [mean(dV); dV; mean(dV)];
    tmp = cov([dV V_obs]);
    dCov = tmp(1:N, N+1:end);
    if rcond(B) < 1e-10
        wu = NaN(N);
    else
        wu = dCov / B;
    end
end
