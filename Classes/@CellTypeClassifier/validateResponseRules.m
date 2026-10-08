function T = validateResponseRules(ctc, opts)
% VALIDATERESPONSERULES  Train on each response rule alone, validate on the others.
%
%   T = ctc.validateResponseRules()
%
% For setups with several perturbations of the same cell class (Bootstrap.ResponseRules, e.g.
% an activating and a silencing DREADD in interneurons): for each rule k the training seeds
% are recomputed from rule k's responders only (generateTrainLabels + classifyUnits on the same
% features and graph), and the labels are validated against the responders of every other
% rule j with validateIndependentResponders (expected interneurons respond in rule j's
% direction, control-calibrated calls). Agreement in both directions means the two
% perturbations identify the same class and pooling them (the default with ResponseRules) is
% justified; a one-sided failure points to a rule whose responders are not the target class.
% Run after identifyResponsiveUnits (with ResponseRules) and discardSegmentUnits. The pooled
% training state (TrainLabels, UnitLabels, ...) is restored afterwards.
%
% OPTIONS:
%   NBoot - bootstrap resamples for the AUC CI (default 200)
%
% OUTPUT: table with TrainRule, ValidationRule, NSeedsInh, NSeedsExc, AUC, AUCLow, AUCHigh,
%         FracInhExpectedInh, FracInhExpectedExc, NExpectedInh, NExpectedExc
arguments
    ctc CellTypeClassifier
    opts.NBoot (1,1) double = 200
end
det = ctc.ResponsivenessDetail;
assert(isstruct(det) && isfield(det, 'rules') && numel(det.rules) >= 2, ...
    'CellTypeClassifier:noResponseRules', ...
    'Run identifyResponsiveUnits() with at least two Bootstrap.ResponseRules first.');
rules = det.rules;
rule_of = det.rule_of_unit(:)';
assert(numel(rule_of) == numel(ctc.ResponsiveUnitIdx), 'CellTypeClassifier:ruleIndex', ...
    'ResponsivenessDetail.rule_of_unit does not match ResponsiveUnitIdx.');

% the state that generateTrainLabels / classifyUnits overwrite
keep = ["ResponsiveUnitIdx", "TrainLabels", "UnitLabels", "UnitConfidence", ...
    "UnitGraphConnectivity", "NormalizationParams", "Graph", "Reduction", "UMAP"];
saved = struct();
for f = keep, saved.(f) = ctc.(f); end
name = @(r) sprintf('%s=%s %s', string(r.Filter{1}), string(r.Filter{2}), r.Direction);

rows = {};
try
    for k = 1:numel(rules)
        ctc.ResponsiveUnitIdx = saved.ResponsiveUnitIdx & rule_of == k;
        evalc('ctc.generateTrainLabels(); ctc.classifyUnits();');
        tl = ctc.TrainLabels;
        for j = setdiff(1:numel(rules), k)
            inh_dir = rules(j).Direction;
            assert(inh_dir ~= "both", 'CellTypeClassifier:ruleDirection', ...
                'Rule %d has direction "both"; validation needs the direction of the target class.', j);
            [~, v] = evalc('ctc.validateIndependentResponders(rules(j).Filter, InhibitoryDirection=inh_dir, NBoot=opts.NBoot)');
            rows(end+1, :) = {string(name(rules(k))), string(name(rules(j))), sum(tl.sorted_y_train == 2), ...
                sum(tl.sorted_y_train == 1), v.AUC, v.AUC_CI(1), v.AUC_CI(2), v.frac_inh_expected_inh, ...
                v.frac_inh_expected_exc, v.n_expected_inh, v.n_expected_exc}; %#ok<AGROW>
            fprintf('train %-22s -> validate %-22s AUC %.3f [%.3f %.3f] (%d / %d validation units)\n', ...
                name(rules(k)), name(rules(j)), v.AUC, v.AUC_CI, v.n_expected_inh, v.n_expected_exc);
        end
    end
catch err
    for f = keep, ctc.(f) = saved.(f); end
    rethrow(err);
end
for f = keep, ctc.(f) = saved.(f); end
T = cell2table(rows, 'VariableNames', {'TrainRule', 'ValidationRule', 'NSeedsInh', 'NSeedsExc', ...
    'AUC', 'AUCLow', 'AUCHigh', 'FracInhExpectedInh', 'FracInhExpectedExc', 'NExpectedInh', 'NExpectedExc'});
end
