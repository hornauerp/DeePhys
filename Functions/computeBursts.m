function [bursts, burst_features] = computeBursts(spike_times, spike_units, ...
    firing_rates, cv_isi, n_units, duration, params)
% COMPUTEBURSTS  Detect and characterize network bursts from spike data.
%
%   [bursts, burst_features] = computeBursts(spike_times, spike_units, ...
%       firing_rates, cv_isi, n_units, duration, params)
%
% INPUTS:
%   spike_times   - (N x 1) double, sorted spike times in seconds
%   spike_units   - (N x 1) double, unit assignment per spike
%   firing_rates  - (M x 1) double, firing rate per unit (for tonic removal)
%   cv_isi        - (M x 1) double, CV of ISI per unit (for tonic removal)
%   n_units       - (1 x 1) double, total number of units
%   duration      - (1 x 1) double, recording duration in seconds
%   params        - struct with fields:
%       .Method         "ISIN" or "logISI" (default "ISIN")
%       .MergeFactor    (default 0.5)
%       .Outlier        struct with .Method, .ThresholdFactor
%       .Binning        (default 0.01, for burst statistics)
%
% OUTPUTS:
%   bursts         - struct: T_start, T_end, best_N, best_ISI_N, Method
%   burst_features - (1 x 10) table from computeBurstStatistics

    arguments
        spike_times   (:,1) double
        spike_units   (:,1) double
        firing_rates  (:,1) double
        cv_isi        (:,1) double
        n_units       (1,1) double
        duration      (1,1) double
        params        struct = struct()
    end

    if ~isfield(params, 'Method'),      params.Method = "ISIN"; end
    if ~isfield(params, 'MergeFactor'), params.MergeFactor = 0.5; end
    if ~isfield(params, 'Outlier'),     params.Outlier = struct('Method', '', 'ThresholdFactor', 3); end
    if ~isfield(params, 'Binning'),     params.Binning = 0.01; end

    bursts = struct('T_start', [], 'T_end', [], 'Method', params.Method);

    if string(params.Method) == "logISI"
        B = detectBurstsLogISI_local(spike_times);
        bursts.T_start = B.T_start;
        bursts.T_end = B.T_end;
        if isfield(B, 'LogISI_Threshold')
            bursts.LogISI_Threshold = B.LogISI_Threshold;
        end
        bursts.best_ISI_N = NaN;
        bursts.best_N = NaN;
    else
        [best_N, best_ISI_N, N_array, ISI_N_array, dunns_coeffs] = ...
            inferOptimalBurstParameters_local(spike_times, spike_units, firing_rates, cv_isi);
        bursts.best_N = best_N;
        bursts.best_ISI_N = best_ISI_N;
        bursts.N = N_array;
        bursts.ISI_N = ISI_N_array;
        bursts.dunns = dunns_coeffs;

        if isnan(best_N)
            burst_features = computeBurstStatistics(bursts, spike_times, spike_units, n_units, params);
            return
        end

        [clean_spikes, ~] = removeTonic(spike_times, spike_units, firing_rates, cv_isi);
        B = detectBurstsISIN_local(clean_spikes, best_N, best_ISI_N);
        bursts.T_start = B.T_start;
        bursts.T_end = B.T_end;
    end

    bursts = pruneBursts_local(bursts, spike_times);
    bursts = mergeBursts_local(bursts, params.MergeFactor);

    burst_features = computeBurstStatistics(bursts, spike_times, spike_units, n_units, params);
end


%% ========================================================================
%  Local helper functions
%  ========================================================================

function Burst = detectBurstsLogISI_local(spike_times, min_spikes)
    if nargin < 2, min_spikes = 3; end

    Burst.T_start = [];
    Burst.T_end = [];

    if length(spike_times) < min_spikes * 2
        return
    end

    isi = diff(spike_times);
    log_isi = log10(isi);
    log_isi = log_isi(isfinite(log_isi));
    if length(log_isi) < 20
        return
    end

    try
        gm = fitgmdist(log_isi(:), 2, 'RegularizationValue', 1e-5, ...
            'Options', statset('MaxIter', 500, 'TolFun', 1e-6), 'Replicates', 3);
    catch
        return
    end

    [~, sort_idx] = sort(gm.mu);
    intra_comp = sort_idx(1);
    inter_comp = sort_idx(2);

    separation = abs(diff(gm.mu)) / mean(sqrt(gm.Sigma(:)));
    if separation < 1.0
        return
    end

    mu1 = gm.mu(intra_comp); sigma1 = sqrt(gm.Sigma(intra_comp));
    mu2 = gm.mu(inter_comp); sigma2 = sqrt(gm.Sigma(inter_comp));
    w1 = gm.ComponentProportion(intra_comp);
    w2 = gm.ComponentProportion(inter_comp);

    x_range = linspace(mu1, mu2, 1000);
    pdf1 = w1 * normpdf(x_range, mu1, sigma1);
    pdf2 = w2 * normpdf(x_range, mu2, sigma2);
    [~, cross_idx] = min(abs(pdf1 - pdf2));
    log_threshold = x_range(cross_idx);
    threshold = 10^log_threshold;

    is_intra = isi <= threshold;

    burst_starts = [];
    burst_ends = [];
    in_burst = false;
    b_start = 0;

    for i = 1:length(is_intra)
        if is_intra(i) && ~in_burst
            b_start = i;
            in_burst = true;
        elseif ~is_intra(i) && in_burst
            n_spikes_in_burst = i - b_start + 1;
            if n_spikes_in_burst >= min_spikes
                burst_starts = [burst_starts; spike_times(b_start)]; %#ok<AGROW>
                burst_ends = [burst_ends; spike_times(i)]; %#ok<AGROW>
            end
            in_burst = false;
        end
    end

    if in_burst
        n_spikes_in_burst = length(is_intra) - b_start + 2;
        if n_spikes_in_burst >= min_spikes
            burst_starts = [burst_starts; spike_times(b_start)];
            burst_ends = [burst_ends; spike_times(end)];
        end
    end

    Burst.T_start = burst_starts;
    Burst.T_end = burst_ends;
    Burst.LogISI_Threshold = threshold;
    Burst.GMModel = gm;
end


function Burst = detectBurstsISIN_local(spike_times, N, ISI_N)
    if N == 0
        Burst.T_start = [];
        Burst.T_end = [];
        return
    end

    dT = zeros(N, length(spike_times)) + inf;
    for j = 0:N-1
        dT(j+1, N:length(spike_times)-(N-1)) = spike_times((N:end-(N-1))+j) - ...
            spike_times((1:end-(N-1)*2)+j);
    end
    Criteria = zeros(size(spike_times));
    Criteria(min(dT) <= ISI_N) = 1;

    SpikeBurstNumber = zeros(size(spike_times)) - 1;
    INBURST = 0;
    NUM_ = 0;
    NUMBER = -1;
    BL = 0;
    for i = N:length(spike_times)
        if INBURST == 0
            if Criteria(i)
                INBURST = 1;
                NUM_ = NUM_ + 1;
                NUMBER = NUM_;
                BL = 1;
            else
                continue
            end
        else
            if ~Criteria(i)
                INBURST = 0;
                if BL < N
                    SpikeBurstNumber(SpikeBurstNumber == NUMBER) = -1;
                    NUM_ = NUM_ - 1;
                end
                NUMBER = -1;
            elseif diff(spike_times([i-(N-1) i])) > ISI_N && BL >= N
                NUM_ = NUM_ + 1;
                NUMBER = NUM_;
                BL = 1;
            else
                BL = BL + 1;
            end
        end
        SpikeBurstNumber(i) = NUMBER;
    end

    MaxBurstNumber = max(SpikeBurstNumber);
    if MaxBurstNumber < 1
        Burst.T_start = [];
        Burst.T_end = [];
        return
    end
    Burst.T_start = zeros(1, MaxBurstNumber);
    Burst.T_end = zeros(1, MaxBurstNumber);
    for i = 1:MaxBurstNumber
        ID = find(SpikeBurstNumber == i);
        Burst.T_start(i) = spike_times(ID(1));
        Burst.T_end(i) = spike_times(ID(end));
    end
end


function [best_N, best_ISI_N, N_array, ISI_N_array, dunns_coeffs] = ...
    inferOptimalBurstParameters_local(spike_times, spike_units, firing_rates, cv_isi)

    [min_isis, N_array, ~] = calculateMinISIN_local(spike_times, spike_units, firing_rates, cv_isi);
    ISI_N_array = min_isis / 1000;

    [clean_spikes, clean_units] = removeTonic(spike_times, spike_units, firing_rates, cv_isi);

    dunns_coeffs = zeros(size(N_array));
    N_spikes = 10000;
    valid = ISI_N_array > 0;
    N_array = N_array(valid);
    ISI_N_array = ISI_N_array(valid);
    dunns_coeffs = zeros(size(N_array));

    for ni = 1:length(N_array)
        i = 0;
        max_i = max([1, floor(length(clean_spikes) / N_spikes)]);
        ISI_N_ms = ISI_N_array(ni) * 1000;
        ISI_N_local = -1;

        while (min(ISI_N_local) > ISI_N_ms || max(ISI_N_local) < ISI_N_ms) && i < max_i
            spike_cutout = clean_spikes((i*N_spikes)+1:min([(i+1)*N_spikes, length(clean_spikes)]));
            ISI_N_local = 1e3 * (spike_cutout(N_array(ni):end) - spike_cutout(1:end-(N_array(ni)-1)));
            i = i + 1;
        end

        if min(ISI_N_local) > ISI_N_ms || max(ISI_N_local) < ISI_N_ms
            dunns_coeffs(ni) = 0;
        else
            idx = ones(size(ISI_N_local));
            idx(ISI_N_local > ISI_N_ms) = 2;
            distM = squareform(pdist(ISI_N_local(:)));
            dunns_coeffs(ni) = dunns(2, distM, idx);
        end
    end

    dunns_coeffs(ISI_N_array < 0.01) = 0;
    [max_vals, max_idx] = maxk(dunns_coeffs, 3);

    if ~isempty(max_vals)
        vals = sum(max_vals > max_vals(1) * 0.9);
        candidate_idx = max_idx(1:vals);
    else
        candidate_idx = [];
    end

    if isempty(candidate_idx)
        best_N = NaN;
        best_ISI_N = NaN;
    else
        best_N = max(N_array(candidate_idx));
        best_ISI_N = ISI_N_array(N_array == best_N);
    end
end


function [min_isis, N_array, prom] = calculateMinISIN_local(spike_times, spike_units, firing_rates, cv_isi)
    Steps = 10.^(-3:.025:2);

    [SpikeTimes, SpikeUnits] = removeTonic(spike_times, spike_units, firing_rates, cv_isi);

    unique_units = unique(SpikeUnits);
    n_unique = length(unique_units);
    min_N = max([4, ceil(n_unique / 30)]);
    max_N = ceil(n_unique / 2);
    if (max_N - min_N) < 20
        N_array = min_N:max_N;
    else
        N_array = round(linspace(min_N, max_N, 20));
    end

    prom = zeros(1, length(N_array));
    min_isis = zeros(1, length(N_array));
    for i = 1:length(N_array)
        FRnum = max([2, N_array(i)]);
        if FRnum > length(SpikeTimes)
            min_isis(i) = 0;
            continue
        end
        ISI_N = SpikeTimes(FRnum:end) - SpikeTimes(1:end-(FRnum-1));
        n = histc(ISI_N * 1000, Steps * 1000);
        n = smooth(n, 'lowess');
        [~, b, ~, p] = findpeaks(n);
        if length(b) < 2
            min_isis(i) = 0;
        else
            [sorted_p, sorted_idx] = sort(p, 'descend');
            b = b(sorted_idx(1:2));
            [~, m_idx] = min(n(min(b):max(b)));
            isi_idx = (m_idx + min(b) - 1);
            min_isis(i) = Steps(isi_idx) * 1000;
            prom(i) = sorted_p(2);
        end
    end
end


function bursts = pruneBursts_local(bursts, spike_times, ops)
    if nargin < 3
        ops.bins = 200;
        ops.min_width = 0.005;
    end

    if isempty(bursts.T_start)
        return
    end

    best_ISI_N = bursts.best_ISI_N;
    if isnan(best_ISI_N)
        best_ISI_N = 0.01;
    end

    pruned.T_start = nan(size(bursts.T_start));
    pruned.T_end = nan(size(bursts.T_end));
    for iburst = 1:length(bursts.T_start)
        padding = max([best_ISI_N, (bursts.T_end(iburst) - bursts.T_start(iburst)) * 0.5]);
        padded_start = bursts.T_start(iburst) - padding;
        padded_end = bursts.T_end(iburst) + padding;
        if ops.bins > 1
            bin_width = max([ops.min_width, (padded_end - padded_start) / ops.bins]);
        else
            bin_width = ops.bins;
        end
        burst_activity = spike_times(spike_times <= padded_end & spike_times >= padded_start);
        binned_activity = histcounts(burst_activity, padded_start:bin_width:padded_end);
        smoothed_activity = smoothdata(binned_activity, 'lowess', 'SmoothingFactor', 0.5) + 1;
        cum_activity_l = cumsum(smoothed_activity);
        cum_activity_r = cumsum(smoothed_activity(end:-1:1));
        threshold_l = triangle_threshold(cum_activity_l, 'L', false);
        threshold_r = triangle_threshold(cum_activity_r(end:-1:1), 'R', false);

        checks(1) = mean(binned_activity(threshold_l:threshold_r)) > mean(binned_activity);
        checks(2) = threshold_l > 1 && threshold_r < length(binned_activity);
        checks(3) = threshold_l < threshold_r;

        if all(checks)
            pruned.T_start(iburst) = padded_start + threshold_l * bin_width;
            pruned.T_end(iburst) = padded_end - (length(smoothed_activity) - threshold_r) * bin_width;
        else
            pruned.T_start(iburst) = bursts.T_start(iburst);
            pruned.T_end(iburst) = bursts.T_end(iburst);
        end
    end
    bursts.T_start = pruned.T_start;
    bursts.T_end = pruned.T_end;
end


function bursts = mergeBursts_local(bursts, merge_factor)
    if isempty(bursts.T_start)
        return
    end

    N_bursts = length(bursts.T_start);
    if isnan(bursts.best_ISI_N)
        return
    end
    IBI_merge = merge_factor * bursts.best_ISI_N;
    if N_bursts > 1
        keep_separate_idx = find(bursts.T_start(2:end) - bursts.T_end(1:end-1) > IBI_merge);
        bursts.T_start = bursts.T_start([1, keep_separate_idx + 1]);
        bursts.T_end = bursts.T_end([keep_separate_idx, N_bursts]);
        fprintf('Merged %i bursts\n', N_bursts - length(bursts.T_start))
    end
end
