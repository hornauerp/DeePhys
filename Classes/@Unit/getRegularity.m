function regularity_table = getRegularity(unit, force)
% GETREGULARITY  Compute regularity features from binned spike activity.
%
%   regularity_table = unit.getRegularity()        — use cache if available
%   regularity_table = unit.getRegularity(true)    — force recomputation

    arguments
        unit Unit
        force logical = false
    end

    unit.ensureCache();
    cache_key = "Regularity_" + paramHash(unit.UnitParams.Regularity);

    if ~force && unit.FeatureCache.isKey(cache_key) && ~isempty(unit.RegularityFeatures)
        regularity_table = unit.RegularityFeatures;
        return
    end

    binned_activity = histcounts(unit.SpikeTimes, ...
        'BinLimits', [0 unit.RecordingDuration], ...
        'BinWidth', unit.UnitParams.Regularity.Binning);
    norm_activity = binned_activity / max(binned_activity);
    [freq, mag, fit_coeff] = computeRegularity(norm_activity, unit.UnitParams.Regularity.Binning, unit.UnitParams.Regularity.N_peaks);
    reg_feat.RegularityFrequency = freq;
    reg_feat.RegularityMagnitude = mag;
    reg_feat.RegularityFit = fit_coeff;
    regularity_table = struct2table(reg_feat);
    unit.FeatureCache(cache_key) = true;
end
