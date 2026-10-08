function report = qualityReport(ctc, opts)
% QUALITYREPORT  Acceptance checks for a classification run (instead of parameter tuning).
%
%   report = ctc.qualityReport()
%   report = ctc.qualityReport(ConditionVar="Construct", ValidationFilter={'Construct', 2})
%
% Checks (PASS / WARN), all computed without ground truth:
%   1. Training seeds     — seed counts (>= 30 per class), enriched communities, and the
%                           agreement of the community-detection restarts with the best
%                           partition (reported; the "enrichment" seeds are a consensus
%                           over restarts).
%   2. Responder calls    — estimated FDR from control cultures (control_calibrated method).
%   3. Condition balance  — inhibitory fraction per value of ConditionVar. Cultures of one
%                           preparation should agree; a fraction that tracks the treatment
%                           means the features encode the drug response (leakage).
%   4. Independent check  — validateIndependentResponders(ValidationFilter), if given.
%   5. Inhibitory fraction with a 95% CI from resampling cultures.
% Run after classifyUnits().
%
% NAME-VALUE OPTIONS:
%   ConditionVar      - metadata column with the treatment condition (default "AAV" if present)
%   ValidationFilter  - {field, value} of an independent condition (default: none)
%   MaxConditionRatio - WARN above this max/min ratio of condition fractions (default 1.25)
%   MaxFDR            - WARN above this estimated responder FDR (default 0.2)
%   MinValidationAUC  - WARN if the AUC CI lower bound is below this (default 0.6)
    arguments
        ctc CellTypeClassifier
        opts.ConditionVar (1,1) string = ""
        opts.ValidationFilter cell = {}
        opts.MaxConditionRatio (1,1) double = 1.25
        opts.MaxFDR (1,1) double = 0.2
        opts.MinValidationAUC (1,1) double = 0.6
    end
    assert(~isempty(ctc.UnitLabels), 'Run classifyUnits() before qualityReport().');
    fs = ctc.FeatureStore; meta = fs.MetadataTable; ut = fs.UnitTable;
    lab = ctc.UnitLabels(:);
    report = struct('checks', table(strings(0, 1), strings(0, 1), strings(0, 1), ...
        'VariableNames', {'Check', 'Status', 'Detail'}));
    add = @(r, name, ok, detail) setfield(r, 'checks', [r.checks; {string(name), string(ok), string(detail)}]); %#ok<SFLD>
    status = @(ok) string(ifelse(ok, "PASS", "WARN"));

    % 1. training seeds
    tl = ctc.TrainLabels;
    n_e = sum(tl.sorted_y_train == 1); n_i = sum(tl.sorted_y_train == 2);
    detail = sprintf('%d inhibitory / %d excitatory seeds', n_i, n_e);
    ok = min(n_e, n_i) >= 30;
    if isfield(tl, 'enrichment')
        detail = sprintf('%s; %d enriched communities', detail, numel(tl.inh_comm_ids));
    end
    if isfield(tl, 'louvain') && isfield(tl.louvain, 'ari_median')
        % Seeds are a consensus over all restarts, so a degenerate partition landscape is
        % reported but does not fail the check.
        detail = sprintf('%s; community restarts vs best partition: median ARI %.2f (seeds by consensus)', ...
            detail, tl.louvain.ari_median);
    end
    report = add(report, "Training seeds", status(ok), detail);

    % 2. responder calls
    if isfield(ctc.ResponsivenessDetail, 'calibration')
        ci = ctc.ResponsivenessDetail.calibration;
        report = add(report, "Responder FDR", status(ci.fdr_inc <= opts.MaxFDR), ...
            sprintf('increase: est. FDR %.2f (%.1f%% treated vs %.1f%% control)', ci.fdr_inc, ...
            100 * ci.rate_inc_treated, 100 * ci.rate_inc_control));
    else
        report = add(report, "Responder FDR", "INFO", "not estimated (use Bootstrap.GroundTruthMethod = ""control_calibrated"")");
    end

    % 3. condition balance
    cv = opts.ConditionVar;
    if cv == "" && ismember("AAV", string(meta.Properties.VariableNames)), cv = "AAV"; end
    if cv ~= "" && ismember(cv, string(meta.Properties.VariableNames))
        [~, mi] = ismember(string(ut.RecordingID), string(meta.RecordingID));
        cond = string(meta.(cv)(mi));
        vals = unique(cond(~isnan(lab)));
        fr = arrayfun(@(v) mean(lab(cond == v & ~isnan(lab)) == 2), vals);
        ratio = max(fr) / max(min(fr), eps);
        report.condition = table(vals, 100 * fr, 'VariableNames', {char(cv), 'PercentInhibitory'});
        report = add(report, "Condition balance", status(ratio <= opts.MaxConditionRatio), ...
            sprintf('%s: %s (max/min %.2f)', cv, strjoin(compose('%s %.1f%%', vals, 100 * fr), ', '), ratio));
    end

    % 4. independent validation
    if ~isempty(opts.ValidationFilter)
        v = ctc.validateIndependentResponders(opts.ValidationFilter);
        report.validation = v;
        report = add(report, "Independent validation", status(v.AUC_CI(1) >= opts.MinValidationAUC), ...
            sprintf('AUC %.3f [%.3f %.3f]; expected I %.0f%% vs expected E %.0f%% labelled inhibitory', ...
            v.AUC, v.AUC_CI, 100 * v.frac_inh_expected_inh, 100 * v.frac_inh_expected_exc));
    end

    % 5. inhibitory fraction with a culture-resampling CI
    cult = FeatureStore.getCultureIDsForUnits(ut.RecordingID, meta, ctc.Parameters.CultureKeys);
    ok_l = ~isnan(lab);
    uc = unique(cult(ok_l));
    per = arrayfun(@(c) [sum(lab(cult == c & ok_l) == 2), sum(cult == c & ok_l)], uc, 'UniformOutput', false);
    per = vertcat(per{:});
    rs = RandStream('twister', 'Seed', 1); B = 1000; f = zeros(B, 1);
    for b = 1:B
        k = randi(rs, numel(uc), numel(uc), 1);
        f(b) = sum(per(k, 1)) / sum(per(k, 2));
    end
    report.fraction = sum(per(:, 1)) / sum(per(:, 2));
    report.fraction_CI = prctile(f, [2.5 97.5]);
    report.per_culture = table(uc, 100 * per(:, 1) ./ per(:, 2), per(:, 2), ...
        'VariableNames', {'Culture', 'PercentInhibitory', 'NUnits'});
    report = add(report, "Inhibitory fraction", "INFO", sprintf('%.1f%% [%.1f %.1f] over %d cultures', ...
        100 * report.fraction, 100 * report.fraction_CI, numel(uc)));

    fprintf('\nClassification quality report\n');
    for i = 1:height(report.checks)
        fprintf('  [%-4s] %-23s %s\n', report.checks.Status(i), report.checks.Check(i), report.checks.Detail(i));
    end
end

function out = ifelse(c, a, b)
    if c, out = a; else, out = b; end
end
