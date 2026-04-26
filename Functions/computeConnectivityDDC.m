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
%   ddc_params   - struct: .BinSize, .Threshold
%
% OUTPUTS:
%   connectivity - struct with .wu (weighted undirected connectivity matrix)

    arguments
        spike_times  (:,1) double
        spike_units  (:,1) double
        n_units      (1,1) double
        duration     (1,1) double
        ddc_params   struct
    end

    TR = 1;
    bin_size = ddc_params.BinSize;
    n_bins = ceil(duration / bin_size);
    V_obs = zeros(n_bins, n_units);
    spiketimes_binned = max(1, min(n_bins, round(spike_times / bin_size) + 1));
    V_obs(sub2ind(size(V_obs), spiketimes_binned, spike_units)) = 1;

    [~, N] = size(V_obs);
    Fx = V_obs - ddc_params.Threshold;
    Fx(Fx < 0) = 0;
    tmp = cov([Fx V_obs]);
    B = tmp(1:N, N+1:end);
    dV = (-1/2*V_obs(1:end-2,:) + 1/2*V_obs(3:end,:)) / TR;
    dV = [mean(dV); dV; mean(dV)];
    tmp = cov([dV V_obs]);
    dCov = tmp(1:N, N+1:end);
    if rcond(B) < eps
        warning('computeConnectivityDDC:singularB', ...
            'B matrix is singular (rcond=%.2e). Setting connectivity to NaN.', rcond(B));
        connectivity.wu = NaN(size(dCov));
    else
        connectivity.wu = dCov / B;
    end
end
