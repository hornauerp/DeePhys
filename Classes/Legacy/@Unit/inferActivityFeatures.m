function inferActivityFeatures(unit, force)
% INFERACTIVITYFEATURES  Compute firing rate, ISI statistics, and optional regularity/catch22.
%
%   unit.inferActivityFeatures()          — use cache if available
%   unit.inferActivityFeatures(true)      — force recomputation
%
% Delegates core computation to computeActivityFeatures() and caches results.

    arguments
        unit Unit
        force logical = false
    end

    unit.ensureCache();

    if isfield(unit.UnitParams, 'Activity') && isfield(unit.UnitParams.Activity, 'FanoBinWidth')
        fano_bin_width = unit.UnitParams.Activity.FanoBinWidth;
    else
        fano_bin_width = 0.1;
    end

    cache_key = "ActivityFeatures_" + paramHash(struct('Duration', unit.RecordingDuration, 'FanoBinWidth', fano_bin_width));

    if ~force && unit.FeatureCache.isKey(cache_key) && ~isempty(unit.ActivityFeatures)
        return
    end

    params.FanoBinWidth     = fano_bin_width;
    params.RefractoryPeriod = unit.UnitParams.QC.RefractoryPeriod;
    if unit.UnitParams.Analyses.Regularity
        params.Regularity = unit.UnitParams.Regularity;
    end
    if unit.UnitParams.Analyses.Catch22
        params.Catch22 = unit.UnitParams.Catch22;
    end

    [act_tbl, reg_tbl, c22_tbl] = computeActivityFeatures(unit.SpikeTimes, unit.RecordingDuration, params);

    unit.ActivityFeatures = act_tbl;
    if unit.UnitParams.Analyses.Regularity && ~isempty(reg_tbl)
        unit.RegularityFeatures = reg_tbl;
    end
    if unit.UnitParams.Analyses.Catch22 && ~isempty(c22_tbl)
        unit.Catch22 = c22_tbl;
    end
    unit.FeatureCache(cache_key) = true;
end
