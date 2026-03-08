function inferActivityFeatures(unit)
    if length(unit.SpikeTimes) > 2
        isi = diff(unit.SpikeTimes);
        act_feat.FiringRate = length(unit.SpikeTimes) / unit.RecordingDuration;
        act_feat.MeanInterSpikeInterval = mean(isi);
        act_feat.VarianceInterSpikeInterval = var(isi);
        act_feat.CVInterSpikeInterval = std(isi) / mean(isi);
        pacf = parcorr(isi, 1);
        act_feat.PartialAutocorrelation = pacf(2);
        unit.ActivityFeatures = struct2table(act_feat);
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
