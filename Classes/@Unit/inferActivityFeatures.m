function inferActivityFeatures(unit, force)
% INFERACTIVITYFEATURES  Compute firing rate, ISI statistics, and optional regularity/catch22.
%
%   unit.inferActivityFeatures()          — use cache if available
%   unit.inferActivityFeatures(true)      — force recomputation
%
% Caches results in unit.FeatureCache to avoid redundant computation.

    arguments
        unit Unit
        force logical = false
    end

    unit.ensureCache();

    % Fano factor bin width — configurable via UnitParams.Activity.FanoBinWidth.
    % Fallback to 0.1 s for objects saved before this parameter was introduced.
    if isfield(unit.UnitParams, 'Activity') && isfield(unit.UnitParams.Activity, 'FanoBinWidth')
        fano_bin_width = unit.UnitParams.Activity.FanoBinWidth;
    else
        fano_bin_width = 0.1;  % legacy default: 100 ms
    end

    cache_key = "ActivityFeatures_" + paramHash(struct('Duration', unit.RecordingDuration, 'FanoBinWidth', fano_bin_width));

    if ~force && unit.FeatureCache.isKey(cache_key) && ~isempty(unit.ActivityFeatures)
        return
    end

    if length(unit.SpikeTimes) <= 2
        % Insufficient spikes — fill with NaN tables to prevent downstream vertcat errors
        act_names = ["FiringRate","MeanInterSpikeInterval","VarianceInterSpikeInterval","CVInterSpikeInterval", ...
            "CV2InterSpikeInterval","LocalVariation","RevisedLocalVariation","FanoFactor","PartialAutocorrelation"];
        unit.ActivityFeatures = array2table(NaN(1, length(act_names)), 'VariableNames', act_names);
        if unit.UnitParams.Analyses.Regularity
            reg_names = ["RegularityFrequency","RegularityMagnitude","RegularityFit"];
            unit.RegularityFeatures = array2table(NaN(1, length(reg_names)), 'VariableNames', reg_names);
        end
        if unit.UnitParams.Analyses.Catch22
            c22_names = "SC_" + string(GetAllFeatureNames());
            unit.Catch22 = array2table(NaN(1, length(c22_names)), 'VariableNames', c22_names);
        end
        unit.FeatureCache(cache_key) = true;
        return
    end

    isi = diff(unit.SpikeTimes);
    act_feat.FiringRate = length(unit.SpikeTimes) / unit.RecordingDuration;
    act_feat.MeanInterSpikeInterval = mean(isi);
    act_feat.VarianceInterSpikeInterval = var(isi);
    mean_isi = mean(isi);
    if mean_isi > 0
        act_feat.CVInterSpikeInterval = std(isi) / mean_isi;
    else
        act_feat.CVInterSpikeInterval = NaN;
    end
    % CV2: local ISI variability (Holt et al., 1996)
    % Insensitive to firing rate changes, captures local irregularity
    if length(isi) > 1
        cv2_vals = 2 * abs(isi(2:end) - isi(1:end-1)) ./ (isi(2:end) + isi(1:end-1));
        act_feat.CV2InterSpikeInterval = mean(cv2_vals);
    else
        act_feat.CV2InterSpikeInterval = NaN;
    end

    % LV: local variation (Shinomoto et al., 2003)
    if length(isi) > 1
        lv_vals = 3 * (isi(2:end) - isi(1:end-1)).^2 ./ (isi(2:end) + isi(1:end-1)).^2;
        act_feat.LocalVariation = mean(lv_vals);
    else
        act_feat.LocalVariation = NaN;
    end

    % LvR: revised local variation with refractory period correction (Shinomoto et al., 2009)
    refract = unit.UnitParams.QC.RefractoryPeriod;  % refractory period in seconds
    if length(isi) > 1
        lvr_vals = (1 - 4*isi(1:end-1).*isi(2:end) ./ (isi(1:end-1)+isi(2:end)).^2) ...
                   .* (1 + 4*refract ./ (isi(1:end-1)+isi(2:end)));
        act_feat.RevisedLocalVariation = 3 * mean(lvr_vals);
    else
        act_feat.RevisedLocalVariation = NaN;
    end

    % Fano factor: var(count)/mean(count) in time bins
    bin_width = fano_bin_width;
    edges = 0:bin_width:unit.RecordingDuration;
    if length(edges) > 2
        spike_counts = histcounts(unit.SpikeTimes, edges);
        mean_count = mean(spike_counts);
        if mean_count > 0
            act_feat.FanoFactor = var(spike_counts) / mean_count;
        else
            act_feat.FanoFactor = NaN;
        end
    else
        act_feat.FanoFactor = NaN;
    end

    if length(isi) > 3  % parcorr requires sufficient samples for lag=1
        pacf = parcorr(isi, 1);
        act_feat.PartialAutocorrelation = pacf(2);
    else
        act_feat.PartialAutocorrelation = NaN;
    end
    unit.ActivityFeatures = struct2table(act_feat);
    unit.FeatureCache(cache_key) = true;

    if unit.UnitParams.Analyses.Regularity
        unit.RegularityFeatures = unit.getRegularity();
    end
    if unit.UnitParams.Analyses.Catch22
        catch_22_table = unit.MEArecording.runCatch22(unit.SpikeTimes);
        catch_22_table.Properties.VariableNames = "SC_" + string(catch_22_table.Properties.VariableNames);
        unit.Catch22 = catch_22_table;
    end
end
