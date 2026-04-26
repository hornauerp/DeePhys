function [activity_table, regularity_table, catch22_table] = computeActivityFeatures(spike_times, duration, params)
% COMPUTEACTIVITYFEATURES  Compute per-unit activity, regularity, and catch22 features.
%
%   [act, reg, c22] = computeActivityFeatures(spike_times, duration, params)
%
% INPUTS:
%   spike_times - (N x 1) double, spike times in seconds
%   duration    - (1 x 1) double, recording duration in seconds
%   params      - struct with fields:
%       .FanoBinWidth      double (default 0.1)
%       .RefractoryPeriod  double (default 0.002)
%       .Regularity        struct with .Binning, .N_peaks (omit field to skip)
%       .Catch22           struct with .BinSize (omit field to skip)
%
% OUTPUTS:
%   activity_table   - (1 x 9) table
%   regularity_table - (1 x 3) table, or empty table if regularity not requested
%   catch22_table    - (1 x 22) table, or empty table if catch22 not requested

    arguments
        spike_times (:,1) double
        duration    (1,1) double
        params      struct = struct()
    end

    if isfield(params, 'FanoBinWidth')
        fano_bin_width = params.FanoBinWidth;
    else
        fano_bin_width = 0.1;
    end

    if isfield(params, 'RefractoryPeriod')
        refract = params.RefractoryPeriod;
    else
        refract = 0.002;
    end

    do_regularity = isfield(params, 'Regularity') && ~isempty(params.Regularity);
    do_catch22    = isfield(params, 'Catch22') && ~isempty(params.Catch22);

    act_names = ["FiringRate","MeanInterSpikeInterval","VarianceInterSpikeInterval","CVInterSpikeInterval", ...
        "CV2InterSpikeInterval","LocalVariation","RevisedLocalVariation","FanoFactor","PartialAutocorrelation"];

    if length(spike_times) <= 2
        activity_table = array2table(NaN(1, length(act_names)), 'VariableNames', act_names);
        if do_regularity
            reg_names = ["RegularityFrequency","RegularityMagnitude","RegularityFit"];
            regularity_table = array2table(NaN(1, length(reg_names)), 'VariableNames', reg_names);
        else
            regularity_table = table();
        end
        if do_catch22
            c22_names = "SC_" + string(GetAllFeatureNames());
            catch22_table = array2table(NaN(1, length(c22_names)), 'VariableNames', c22_names);
        else
            catch22_table = table();
        end
        return
    end

    isi = diff(spike_times);
    act_feat.FiringRate = length(spike_times) / duration;
    act_feat.MeanInterSpikeInterval = mean(isi);
    act_feat.VarianceInterSpikeInterval = var(isi);
    mean_isi = mean(isi);
    if mean_isi > 0
        act_feat.CVInterSpikeInterval = std(isi) / mean_isi;
    else
        act_feat.CVInterSpikeInterval = NaN;
    end

    if length(isi) > 1
        cv2_vals = 2 * abs(isi(2:end) - isi(1:end-1)) ./ (isi(2:end) + isi(1:end-1));
        act_feat.CV2InterSpikeInterval = mean(cv2_vals);
    else
        act_feat.CV2InterSpikeInterval = NaN;
    end

    if length(isi) > 1
        lv_vals = 3 * (isi(2:end) - isi(1:end-1)).^2 ./ (isi(2:end) + isi(1:end-1)).^2;
        act_feat.LocalVariation = mean(lv_vals);
    else
        act_feat.LocalVariation = NaN;
    end

    if length(isi) > 1
        lvr_vals = (1 - 4*isi(1:end-1).*isi(2:end) ./ (isi(1:end-1)+isi(2:end)).^2) ...
                   .* (1 + 4*refract ./ (isi(1:end-1)+isi(2:end)));
        act_feat.RevisedLocalVariation = 3 * mean(lvr_vals);
    else
        act_feat.RevisedLocalVariation = NaN;
    end

    bin_width = fano_bin_width;
    edges = 0:bin_width:duration;
    if length(edges) > 2
        spike_counts = histcounts(spike_times, edges);
        mean_count = mean(spike_counts);
        if mean_count > 0
            act_feat.FanoFactor = var(spike_counts) / mean_count;
        else
            act_feat.FanoFactor = NaN;
        end
    else
        act_feat.FanoFactor = NaN;
    end

    if length(isi) > 3
        pacf = parcorr(isi, 1);
        act_feat.PartialAutocorrelation = pacf(2);
    else
        act_feat.PartialAutocorrelation = NaN;
    end
    activity_table = struct2table(act_feat);

    if do_regularity
        binned_activity = histcounts(spike_times, ...
            'BinLimits', [0 duration], ...
            'BinWidth', params.Regularity.Binning);
        norm_activity = binned_activity / max(binned_activity);
        [freq, mag, fit_coeff] = computeRegularity(norm_activity, params.Regularity.Binning, params.Regularity.N_peaks);
        reg_feat.RegularityFrequency = freq;
        reg_feat.RegularityMagnitude = mag;
        reg_feat.RegularityFit = fit_coeff;
        regularity_table = struct2table(reg_feat);
    else
        regularity_table = table();
    end

    if do_catch22
        catch22_table = computeCatch22(spike_times, duration, params.Catch22.BinSize);
        catch22_table.Properties.VariableNames = "SC_" + string(catch22_table.Properties.VariableNames);
    else
        catch22_table = table();
    end
end
