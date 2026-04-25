function cv_results = evaluatePipeline(ctc, opts)
% EVALUATEPIPELINE  Leave-one-culture-out cross-validation of the CTC pipeline.
%
% For each fold, one culture is withheld from the UMAP training label set while
% all units remain embedded in the shared global UMAP (transductive setting).
% Classification performance is evaluated against drug-response ground truth:
%   responsive units (ResponsiveUnitIdx) → proxy inhibitory ground truth.
%
% The UMAP embedding is reused across folds by default (ReuseUMAP=true): since
% global normalization is fitted on all units and the unit set is identical in
% every fold, the embedding does not change — only the training label subset does.
%
% WORKFLOW:
%   ctc.identifyResponsiveUnits();   % must run first
%   ctc.generateTrainLabels();       % must run once to build initial UMAP
%   cv = ctc.evaluatePipeline();
%   % Labels/confidence are left from the last fold — re-run pipeline if needed.
%
% INPUTS:
%   ctc                    - CellTypeClassifier handle
%   opts.HeldOutCultureIdx - (1 x K) culture indices to hold out; default = all (LOO)
%   opts.ReuseUMAP         - (logical) reuse UMAP across folds (default: true)
%   opts.Verbose           - (logical) print per-fold summary (default: true)
%
% OUTPUTS:
%   cv_results - struct array (one element per fold) with fields:
%     fold_idx       - culture index held out
%     n_responsive   - # drug-confirmed inh units in held-out culture
%     n_labeled_inh  - # units labeled inh in held-out culture
%     true_positives - responsive units correctly labeled inh
%     precision_inh  - TP / n_labeled_inh  (0 when n_labeled_inh=0)
%     recall_inh     - TP / n_responsive   (NaN when n_responsive=0)
%     labels         - (1 x N) full UnitLabels for this fold
%     confidence     - (1 x N) full UnitConfidence for this fold
%
% NOTE: Ensemble classification (Parameters.Ensemble.Enabled) multiplies runtime
% by the number of seeds per fold. Disable for fast iteration:
%   ctc.Parameters.Ensemble.Enabled = false;

arguments
    ctc   CellTypeClassifier
    opts.HeldOutCultureIdx (1,:) double = []
    opts.ReuseUMAP         (1,1) logical = true
    opts.Verbose           (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before evaluatePipeline().');
assert(~isempty(ctc.UMAP) && ~isempty(ctc.UMAP.head), ...
    ['Run generateTrainLabels() at least once before evaluatePipeline() ' ...
     'so the initial UMAP is available for reuse across folds.']);

% -- Unique cultures from UnitTable (same ordering as TrainingCultureIdx) -----
unit_table      = ctc.FeatureStore.UnitTable;
meta            = ctc.FeatureStore.MetadataTable;
cult_ids        = FeatureStore.getCultureIDsForUnits( ...
    unit_table.RecordingID, meta, ctc.Parameters.CultureKeys);
unique_cultures = unique(cult_ids, 'stable');
n_cultures      = numel(unique_cultures);

if isempty(opts.HeldOutCultureIdx)
    held_out_indices = 1:n_cultures;
else
    assert(all(opts.HeldOutCultureIdx >= 1 & opts.HeldOutCultureIdx <= n_cultures), ...
        'HeldOutCultureIdx contains values outside [1, %d].', n_cultures);
    held_out_indices = opts.HeldOutCultureIdx;
end

resp_class_label   = ctc.Parameters.TrainLabels.ResponsiveClassLabel;
orig_training_idx  = ctc.Parameters.UMAP.TrainingCultureIdx;

cv_results = struct('fold_idx', {}, 'n_responsive', {}, 'n_labeled_inh', {}, ...
    'true_positives', {}, 'precision_inh', {}, 'recall_inh', {}, ...
    'labels', {}, 'confidence', {});

n_folds = numel(held_out_indices);

for fi = 1:n_folds
    c = held_out_indices(fi);

    % Hold out culture c: train on all other cultures
    ctc.Parameters.UMAP.TrainingCultureIdx = setdiff(1:n_cultures, c);

    if opts.Verbose
        fprintf('--- Fold %d/%d: holding out culture %d/%d ---\n', fi, n_folds, c, n_cultures);
    end

    % Recompute training labels (UMAP reused when ReuseUMAP=true).
    % clearCache invalidates NormalizedFeatures so buildNormalizedFeatures
    % recomputes subset_mask with the new TrainingCultureIdx; it does NOT
    % clear ctc.UMAP or ctc.Reduction.
    ctc.clearCache();
    ctc.generateTrainLabels(ReuseExistingUMAP=opts.ReuseUMAP);
    ctc.classifyUnits();

    % -- Held-out culture mask over UnitDataArray / UnitTable rows -------------
    % UnitDataArray and UnitTable are guaranteed same length and order.
    held_out_mask = (cult_ids == unique_cultures(c))';  % (1 x N) logical

    % -- Ground truth: responsive (drug-confirmed inh) in held-out culture -----
    resp_in_held  = ctc.ResponsiveUnitIdx & held_out_mask;
    n_responsive  = sum(resp_in_held);

    % -- Predicted: units labeled inhibitory in held-out culture ---------------
    pred_inh      = ~isnan(ctc.UnitLabels) & ctc.UnitLabels == resp_class_label;
    inh_in_held   = pred_inh & held_out_mask;
    n_labeled_inh = sum(inh_in_held);

    % -- Performance -----------------------------------------------------------
    tp        = sum(resp_in_held & pred_inh);
    precision = tp / max(n_labeled_inh, 1);
    if n_responsive > 0
        recall = tp / n_responsive;
    else
        recall = NaN;
    end

    if opts.Verbose
        fprintf('  Culture %d | responsive=%d | labeled_inh=%d | TP=%d | prec=%.2f | recall=%.2f\n', ...
            c, n_responsive, n_labeled_inh, tp, precision, recall);
    end

    cv_results(fi).fold_idx       = c;
    cv_results(fi).n_responsive   = n_responsive;
    cv_results(fi).n_labeled_inh  = n_labeled_inh;
    cv_results(fi).true_positives = tp;
    cv_results(fi).precision_inh  = precision;
    cv_results(fi).recall_inh     = recall;
    cv_results(fi).labels         = ctc.UnitLabels;
    cv_results(fi).confidence     = ctc.UnitConfidence;
end

% Restore original TrainingCultureIdx
ctc.Parameters.UMAP.TrainingCultureIdx = orig_training_idx;

if opts.Verbose
    valid_recall = [cv_results.recall_inh];
    valid_recall = valid_recall(~isnan(valid_recall));
    fprintf('--- Cross-validation complete (%d folds) ---\n', n_folds);
    fprintf('  Mean precision: %.3f | Mean recall: %.3f\n', ...
        mean([cv_results.precision_inh]), mean(valid_recall));
    fprintf('  NOTE: UnitLabels/Confidence reflect the last fold. Re-run pipeline to restore full-data results.\n');
end
end
