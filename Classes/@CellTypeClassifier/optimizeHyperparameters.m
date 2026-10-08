function results = optimizeHyperparameters(ctc, opts) %#ok<STOUT,INUSD>
% OPTIMIZEHYPERPARAMETERS  Removed: tuned the supervised UMAP of the former classifier.
%
% Labels now come from label propagation on the unit graph (classifyUnits), which
% has no supervised embedding, and its search settings (SupervisedNDimsRange,
% TargetWeightRange, ObjectiveMetric, ...) are no longer in the parameters, so this
% method could not run. The default pipeline needs no tuning: check a run with
% ctc.qualityReport(). optimizeUnsupervisedUMAP() remains for the optional graph
% and community settings. The former implementation is in the git history.
arguments
    ctc CellTypeClassifier %#ok<INUSA>
    opts.Verbose (1,1) logical = true
end
error('CellTypeClassifier:removed', ['optimizeHyperparameters was removed: the supervised ' ...
    'UMAP it tuned is no longer part of the classifier. Use ctc.qualityReport() to check ' ...
    'a run, or optimizeUnsupervisedUMAP() for the optional graph settings.']);
end
