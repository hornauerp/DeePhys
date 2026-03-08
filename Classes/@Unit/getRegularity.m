function regularity_table = getRegularity(unit)
    binned_activity = histcounts(unit.SpikeTimes, ...
        'BinLimits', [0 unit.RecordingDuration], ...
        'BinWidth', unit.UnitParams.Regularity.Binning);
    norm_activity = binned_activity / max(binned_activity);
    [freq, mag, fit_coeff] = computeRegularity(norm_activity, unit.UnitParams.Regularity.Binning, unit.UnitParams.Regularity.N_peaks);
    reg_feat.RegularityFrequency = freq;
    reg_feat.RegularityMagnitude = mag;
    reg_feat.RegularityFit = fit_coeff;
    regularity_table = struct2table(reg_feat);
end
