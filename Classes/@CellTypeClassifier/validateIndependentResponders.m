function res = validateIndependentResponders(ctc, metadata_filter, opts)
% VALIDATEINDEPENDENTRESPONDERS  Check the labels against responders of an independent condition.
%
%   res = ctc.validateIndependentResponders({'Construct', 2})
%
% Uses cultures that were NOT used for training labels, e.g. cultures with an inhibitory
% DREADD (hM4Di in interneurons): units that DEcrease under the actuator are silenced
% interneurons (expected inhibitory), units that INcrease are disinhibited excitatory
% neurons (expected excitatory). Calls use the control-calibrated method
% (calibratedResponseTable / calibrateCalls with Bootstrap.ControlFilter, ControlQuantile,
% CalibratedAlpha, HighDoses), so no thresholds are tuned. Reports the AUC of P(inhibitory)
% separating the two groups (95% CI by resampling cultures) and the fraction labelled
% inhibitory in each group. Run after classifyUnits().
%
% NAME-VALUE OPTIONS:
%   InhibitoryDirection - response direction of expected interneurons: "decrease" (default,
%                         inhibitory DREADD) or "increase" (activating DREADD)
%   NBoot               - bootstrap resamples for the CI (default 1000)
%
% OUTPUT: struct with AUC, AUC_CI, n_expected_inh, n_expected_exc, frac_inh_expected_inh,
%         frac_inh_expected_exc, base_rate, calibration (thresholds / estimated FDR), units
%         (table: UnitID, Culture, ExpectedInhibitory, PInhibitory of the validation units)
    arguments
        ctc CellTypeClassifier
        metadata_filter cell
        opts.InhibitoryDirection (1,1) string = "decrease"
        opts.NBoot (1,1) double = 1000
    end
    assert(~isempty(ctc.UnitLabels), 'Run classifyUnits() before validateIndependentResponders().');
    p    = ctc.Parameters.Bootstrap;
    fs   = ctc.FeatureStore;
    meta = fs.MetadataTable;
    keys = ctc.Parameters.CultureKeys;
    ud_all = ctc.UnitDataArray;
    if numel(ctc.OriginalUnitDataArray) > numel(ud_all), ud_all = ctc.OriginalUnitDataArray; end

    rec_cult = FeatureStore.buildCultureIDs(meta, keys);
    matches = @(field, val) arrayfun(@(i) isequal(string(meta.(field)(i)), string(val)), (1:height(meta))');
    assert(iscell(p.ControlFilter) && numel(p.ControlFilter) == 2, 'CellTypeClassifier:noControlFilter', ...
        'Set Bootstrap.ControlFilter = {field, value} for the control cultures.');
    controls = unique(rec_cult(matches(p.ControlFilter{1}, p.ControlFilter{2})));
    treated  = setdiff(unique(rec_cult(matches(metadata_filter{1}, metadata_filter{2}))), controls);
    T = calibratedResponseTable(ud_all, meta, keys, string(ctc.Parameters.UMAP.GroupingVar), ...
        ctc.Parameters.UMAP.GroupingValues, [treated(:); controls(:)], p);
    is_ctrl = ismember(T.Culture, controls);
    [inc, dec, info] = calibrateCalls(T, is_ctrl, p.ControlQuantile, p.CalibratedAlpha);
    if opts.InhibitoryDirection == "decrease"
        exp_inh = dec & ~is_ctrl; exp_exc = inc & ~is_ctrl;
    else
        exp_inh = inc & ~is_ctrl; exp_exc = dec & ~is_ctrl;
    end

    % classifier output per (UnitID, culture)
    ut = fs.UnitTable;
    cult_rows = FeatureStore.getCultureIDsForUnits(ut.RecordingID, meta, keys);
    lab = ctc.UnitLabels(:); conf = ctc.UnitConfidence(:);
    p_inh = conf; p_inh(lab == 1) = 1 - conf(lab == 1);
    [hit, loc] = ismember(T.UnitID + "|" + T.Culture, string(ut.UnitID) + "|" + cult_rows);
    P = nan(height(T), 1); L = P;
    P(hit) = p_inh(loc(hit)); L(hit) = lab(loc(hit));
    sel = (exp_inh | exp_exc) & ~isnan(P);
    y = exp_inh(sel); s = P(sel); cg = T.Culture(sel);
    res = struct();
    res.n_expected_inh = sum(y); res.n_expected_exc = sum(~y);
    res.frac_inh_expected_inh = mean(L(exp_inh & ~isnan(L)) == 2);
    res.frac_inh_expected_exc = mean(L(exp_exc & ~isnan(L)) == 2);
    in_treated = ismember(cult_rows, treated) & ~isnan(lab);
    res.base_rate = mean(lab(in_treated) == 2);
    res.AUC = aucScore(s, y);
    % CI: resample cultures with replacement
    uc = unique(cg); B = opts.NBoot; a = nan(B, 1);
    rs = RandStream('twister', 'Seed', 1);
    for b = 1:B
        pick = uc(randi(rs, numel(uc), numel(uc), 1));
        idx = cell2mat(arrayfun(@(c) find(cg == c), pick, 'UniformOutput', false));
        if any(y(idx)) && any(~y(idx)), a(b) = aucScore(s(idx), y(idx)); end
    end
    res.AUC_CI = prctile(a, [2.5 97.5]);
    res.calibration = info;
    % the validation units themselves, e.g. to test other features against them
    res.units = table(T.UnitID(sel), T.Culture(sel), y, s, ...
        'VariableNames', {'UnitID', 'Culture', 'ExpectedInhibitory', 'PInhibitory'});
    fprintf(['Independent validation (%s = %s, expected interneurons %s): AUC %.3f [%.3f %.3f]\n' ...
             '   expected inhibitory n=%d: %.1f%% labelled inhibitory | expected excitatory n=%d: %.1f%% | base rate %.1f%%\n'], ...
        string(metadata_filter{1}), string(metadata_filter{2}), opts.InhibitoryDirection, res.AUC, res.AUC_CI, ...
        res.n_expected_inh, 100 * res.frac_inh_expected_inh, res.n_expected_exc, 100 * res.frac_inh_expected_exc, ...
        100 * res.base_rate);
end

function a = aucScore(s, y)
% Mann-Whitney AUC with tied ranks
    y = logical(y(:)); r = tiedrank(s(:));
    n1 = sum(y); n0 = sum(~y);
    a = (sum(r(y)) - n1 * (n1 + 1) / 2) / (n1 * n0);
end
