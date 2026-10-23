function burst_features = computeBurstStatistics(bursts, spike_times, spike_units, n_units, params)
% COMPUTEBURSTSTATISTICS  Compute network-level burst statistics.
%
%   burst_features = computeBurstStatistics(bursts, spike_times, spike_units, n_units, params)
%
% INPUTS:
%   bursts       - struct with .T_start, .T_end vectors
%   spike_times  - (N x 1) double, sorted spike times
%   spike_units  - (N x 1) double, unit assignment per spike
%   n_units      - (1 x 1) double, total number of units
%   params       - struct with fields:
%       .Binning  (default 0.01)
%       .Outlier.Method (string, e.g. "median", or "" to skip)
%       .Outlier.ThresholdFactor (double)
%       .ShapeLevel (default 0.25, see burstShape)
%
% OUTPUTS:
%   burst_features - (1 x 14) table. BurstDuration / RiseTime / FallTime (legacy, unchanged):
%       detector span and the maximum bin within it, so they depend on how the detector sets
%       its edges (unrelated between detectors for the same bursts) and on Binning.
%       BurstWidth / BurstRise / BurstFall / BurstFWHM: medians of burstShape, measured at
%       ShapeLevel of each burst's own peak above baseline, independent of the detector edges.

    arguments
        bursts      struct
        spike_times (:,1) double
        spike_units (:,1) double
        n_units     (1,1) double
        params      struct = struct()
    end

    params = parseStructParameters(struct( ...
        'Binning', 0.01, ...
        'Outlier', struct('Method', '', 'ThresholdFactor', 3), ...
        'ShapeLevel', 0.25), params);

    if length(bursts.T_start) < 3
        Burst.MeanInterBurstInterval = NaN;
        Burst.BurstDuration = NaN;
        Burst.RiseTime = NaN;
        Burst.FallTime = NaN;
        Burst.PeakFR = NaN;
        Burst.StdFR = NaN;
        Burst.StdBurstDuration = NaN;
        Burst.StdInterBurstInterval = NaN;
        Burst.IntraFiringRate = NaN;
        Burst.InterFiringRate = NaN;
        Burst.BurstWidth = NaN;
        Burst.BurstRise = NaN;
        Burst.BurstFall = NaN;
        Burst.BurstFWHM = NaN;
        burst_features = struct2table(Burst);
        return
    end

    IBIs = bursts.T_start(2:end) - bursts.T_end(1:end-1);
    BDs = bursts.T_end - bursts.T_start;

    [RiseTime, FallTime, PeakFR, StdFR] = deal(nan(1, length(BDs)));
    Fs = round(1 / params.Binning);
    for iburst = 1:length(BDs)
        try
            burst_activity = spike_times(spike_times <= bursts.T_end(iburst) & spike_times >= bursts.T_start(iburst));
            binned_activity = histcounts(burst_activity, ...
                'BinLimits', [bursts.T_start(iburst), bursts.T_end(iburst)], ...
                'BinWidth', params.Binning);
            [max_binned, imax_binned] = max(binned_activity);
            norm_activity = smoothdata(binned_activity / max_binned, 'lowess');
            RiseTime(iburst) = imax_binned * params.Binning;
            FallTime(iburst) = (length(binned_activity) - imax_binned) * params.Binning;
            PeakFR(iburst) = max_binned / n_units * Fs;
            StdFR(iburst) = std(norm_activity);
        catch ME
            warning(ME.message)
        end
    end

    [inburst, outburst] = getBurstSpikes_local(bursts, spike_times, spike_units);
    inter_spikes = cellfun(@length, outburst.spikes);
    inter_fr = inter_spikes(2:end-1) ./ (IBIs * n_units);
    intra_spikes = cellfun(@length, inburst.spikes);
    intra_fr = intra_spikes ./ (BDs * n_units);

    cellfun_input = {IBIs, BDs, RiseTime, FallTime, PeakFR, StdFR};

    if strlength(string(params.Outlier.Method)) > 0   % "" (string) and '' both skip
        cellfun_input = cellfun(@(x) rmoutliers(x, params.Outlier.Method, ...
            'ThresholdFactor', params.Outlier.ThresholdFactor), cellfun_input, 'un', 0);
    end
    feature_means = cellfun(@(x) median(x, 'omitnan'), cellfun_input, 'un', 0);

    [Burst.MeanInterBurstInterval, ...
        Burst.BurstDuration, ...
        Burst.RiseTime, ...
        Burst.FallTime, ...
        Burst.PeakFR, ...
        Burst.StdFR] = feature_means{:};
    % Spread of the same (outlier-filtered) values the central values above summarise;
    % on the raw values one merged or very long burst dominated the SD
    Burst.StdBurstDuration = std(cellfun_input{2}, 'omitnan');
    Burst.StdInterBurstInterval = std(cellfun_input{1}, 'omitnan');
    Burst.IntraFiringRate = median(intra_fr);
    Burst.InterFiringRate = median(inter_fr);
    % detector-independent timing and shape (burstShape at ShapeLevel)
    S = burstShape(spike_times, bursts.T_start(:)', bursts.T_end(:)', ...
        max([spike_times; bursts.T_end(:)]) + 1, Level=params.ShapeLevel);
    Burst.BurstWidth = median(S.Duration, 'omitnan');
    Burst.BurstRise  = median(S.Rise, 'omitnan');
    Burst.BurstFall  = median(S.Fall, 'omitnan');
    Burst.BurstFWHM  = median(S.FWHM, 'omitnan');
    burst_features = struct2table(Burst);
end


function [inburst, outburst] = getBurstSpikes_local(bursts, spike_times, spike_units)
    n_bursts = length(bursts.T_start);
    inburst.spikes = cell(1, n_bursts);
    inburst.units = cell(1, n_bursts);
    outburst.spikes = cell(1, n_bursts + 1);
    outburst.units = cell(1, n_bursts + 1);

    for iburst = 1:n_bursts
        burst_idx = spike_times > bursts.T_start(iburst) & spike_times < bursts.T_end(iburst);
        inburst.spikes{iburst} = spike_times(burst_idx);
        inburst.units{iburst} = spike_units(burst_idx);
        if iburst == 1
            outburst_idx = spike_times > 0 & spike_times < bursts.T_start(iburst);
        else
            outburst_idx = spike_times > bursts.T_end(iburst-1) & spike_times < bursts.T_start(iburst);
        end
        outburst.spikes{iburst} = spike_times(outburst_idx);
        outburst.units{iburst} = spike_units(outburst_idx);
    end
    outburst_idx = spike_times > bursts.T_end(end);
    outburst.spikes{n_bursts + 1} = spike_times(outburst_idx);
    outburst.units{n_bursts + 1} = spike_units(outburst_idx);
end
