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
    cache_key = "ActivityFeatures_" + paramHash(struct('Duration', unit.RecordingDuration));

    if ~force && unit.FeatureCache.isKey(cache_key) && ~isempty(unit.ActivityFeatures)
        return
    end

    if length(unit.SpikeTimes) > 2
        isi = diff(unit.SpikeTimes);
        act_feat.FiringRate = length(unit.SpikeTimes) / unit.RecordingDuration;
        act_feat.MeanInterSpikeInterval = mean(isi);
        act_feat.VarianceInterSpikeInterval = var(isi);
        act_feat.CVInterSpikeInterval = std(isi) / mean(isi);
        pacf = parcorr(isi, 1);
        act_feat.PartialAutocorrelation = pacf(2);
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
end
