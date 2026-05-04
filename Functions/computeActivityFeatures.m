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
%   activity_table   - (1 x 10) table
%   regularity_table - (1 x 3) table, or empty table if regularity not requested
%   catch22_table    - (1 x 22) table, or empty table if catch22 not requested

    arguments
        spike_times (:,1) double
        duration    (1,1) double
        params      struct = struct()
    end

    p = parseStructParameters(struct( ...
        'FanoBinWidth',    0.1, ...
        'RefractoryPeriod', 0.002, ...
        'Regularity',      [], ...
        'Catch22',         []), params);

    fano_bin_width = p.FanoBinWidth;
    refract        = p.RefractoryPeriod;
    do_regularity  = ~isempty(p.Regularity);
    do_catch22     = ~isempty(p.Catch22);

    act_names = ["FiringRate","MeanInterSpikeInterval","VarianceInterSpikeInterval","CVInterSpikeInterval", ...
        "CV2InterSpikeInterval","LocalVariation","RevisedLocalVariation","FanoFactor","PartialAutocorrelation", ...
        "ISIBimodality"];

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
        denom_cv2 = isi(2:end) + isi(1:end-1);
        cv2_vals = 2 * abs(isi(2:end) - isi(1:end-1)) ./ denom_cv2;
        cv2_vals(denom_cv2 == 0) = NaN;
        act_feat.CV2InterSpikeInterval = mean(cv2_vals, 'omitnan');
    else
        act_feat.CV2InterSpikeInterval = NaN;
    end

    if length(isi) > 1
        denom_lv = (isi(2:end) + isi(1:end-1)).^2;
        lv_vals = 3 * (isi(2:end) - isi(1:end-1)).^2 ./ denom_lv;
        lv_vals(denom_lv == 0) = NaN;
        act_feat.LocalVariation = mean(lv_vals, 'omitnan');
    else
        act_feat.LocalVariation = NaN;
    end

    if length(isi) > 1
        isi_a = isi(1:end-1);
        isi_b = isi(2:end);
        denom_lvr = isi_a + isi_b;
        valid_lvr = denom_lvr > 0;
        lvr_vals = nan(size(denom_lvr));
        lvr_vals(valid_lvr) = (1 - 4*isi_a(valid_lvr).*isi_b(valid_lvr) ./ denom_lvr(valid_lvr).^2) ...
                               .* (1 + 4*refract ./ denom_lvr(valid_lvr));
        act_feat.RevisedLocalVariation = 3 * mean(lvr_vals, 'omitnan');
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

    % Bimodality coefficient on log(ISI): BC = (skewness^2 + 1) / kurtosis
    % BC > 5/9 (~0.556) suggests bimodality (e.g. burst vs tonic ISI modes)
    if length(isi) >= 10
        log_isi = log(isi(isi > 0));
        if length(log_isi) >= 10
            k = kurtosis(log_isi);
            if isnan(k) || k < 0.1
                act_feat.ISIBimodality = NaN;
            else
                act_feat.ISIBimodality = (skewness(log_isi)^2 + 1) / k;
            end
        else
            act_feat.ISIBimodality = NaN;
        end
    else
        act_feat.ISIBimodality = NaN;
    end
    activity_table = struct2table(act_feat);

    if do_regularity
        binned_activity = histcounts(spike_times, ...
            'BinLimits', [0 duration], ...
            'BinWidth', params.Regularity.Binning);
        max_activity = max(binned_activity);
        if max_activity == 0
            norm_activity = zeros(size(binned_activity));
        else
            norm_activity = binned_activity / max_activity;
        end
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
