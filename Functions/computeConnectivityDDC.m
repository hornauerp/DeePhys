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
%   connectivity - struct with
%     .wd        (n x n) weighted directed graph: the signed DDC coefficients of all
%                pairs, diagonal 0. Positive ~ excitatory, negative ~ inhibitory influence.
%     .wd_sig    the same, restricted to the surrogate-significant links (0 elsewhere)
%     .bd        unweighted (binary) directed graph of the significant links
%     .threshold (n x n) surrogate |DDC| percentile per link
%   computeGraphFeatures uses .bd for the binary graph features and |.wd_sig| for the
%   weighted ones (Parameters.GraphFeatures.Types). Legacy MEArecording stored the
%   unthresholded matrix as .wu.

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

    % Binned activity as a sparse binary matrix (n_bins x n_units). The dense version
    % took n_bins x n_units doubles per copy (600 s at 1 ms, 200 units: ~1 GB).
    bin_size = ddc_params.BinSize;
    n_bins = ceil(duration / bin_size);
    bins   = max(1, min(n_bins, round(spike_times / bin_size) + 1));
    V_obs  = spones(sparse(bins, spike_units, 1, n_bins, n_units));

    wd = ddcCore(V_obs, ddc_params.Threshold);
    wd(1:n_units+1:end) = 0;
    connectivity.wd = wd;

    if all(isnan(wd(:)))
        connectivity.wd_sig = zeros(n_units);
        connectivity.bd = zeros(n_units);
        connectivity.threshold = NaN(n_units);
        return
    end

    % Surrogate test: circularly shift each unit's binned activity independently
    [rows, cols] = find(V_obs);
    surrogate_mat = zeros(n_units, n_units, ddc_params.N_Surrogates);
    for s = 1:ddc_params.N_Surrogates
        shift  = zeros(n_units, 1);
        for u = 1:n_units
            shift(u) = randi(n_bins);   % same draws as circshift(V(:, u), randi(n_bins))
        end
        rows_s = mod(rows - 1 + shift(cols), n_bins) + 1;
        V_shuf = sparse(rows_s, cols, 1, n_bins, n_units);
        surrogate_mat(:,:,s) = ddcCore(V_shuf, ddc_params.Threshold);
    end

    connectivity.threshold = prctile(abs(surrogate_mat), ddc_params.Percentile, 3);
    bd = double(abs(wd) > connectivity.threshold);
    bd(1:n_units+1:end) = 0;
    connectivity.wd_sig = wd .* bd;
    connectivity.bd = bd;
end


function W = ddcCore(V, relu_threshold)
% DDC = cov(dV, V) / cov(ReLU(V - threshold), V) for a binary (sparse) activity matrix V,
% with dV the central difference of V padded by its column means. Covariances are
% computed from sparse products: cov(A, B) = (A'B - n mean(A)' mean(B)) / (n - 1).
    [n, N] = size(V);
    mV  = full(mean(V, 1));
    cVV = (full(V' * V) - n * (mV' * mV)) / (n - 1);
    % V is binary, so ReLU(V - t) = (1 - t) V for 0 <= t < 1, and 0 for t >= 1
    B   = max(1 - relu_threshold, 0) * cVV;
    if relu_threshold < 0
        % ReLU(V - t) = V - t for t < 0: a constant shift leaves the covariance unchanged
        B = cVV;
    end
    dV  = (V(3:end, :) - V(1:end-2, :)) / 2;              % central difference (TR = 1)
    mdV = full(mean(dV, 1));                               % padding rows (first, last)
    S   = full(dV' * V(2:end-1, :)) + mdV' * full(V(1, :) + V(end, :));
    dCov = (S - n * (mdV' * mV)) / (n - 1);                % the padded dV has mean mdV
    if rcond(B) < 1e-10
        W = NaN(N);
    else
        W = dCov / B;
    end
end
