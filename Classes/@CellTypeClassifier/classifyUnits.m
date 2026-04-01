function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Feature extraction uses ctc.Parameters.Harmonization so DeePhys ACGs are
% recomputed at the configured bin size and lag (default: 0.5 ms / 100 ms).
%
% Normalisation pipeline:
%   1. Z-score within each NormalizationVar group (e.g. ChipID) — fit on all units
%   2. Global z-score: fit on train, apply train statistics to test
%   3. Drop columns with NaN or Inf in either train or test
%   4. Scale by max(abs(train))
%
% Requires ctc.TrainLabels to be set (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyUnits()');

labels = ctc.TrainLabels;
p_umap = ctc.Parameters.UMAP;

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr]                        = ctc.getOrExtract(ctc.UnitDataArray);
[X_raw, feat_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% ── Step 1: per-group normalisation using FeatureStore table column ───────────
X_raw(isnan(X_raw)) = 0;

norm_var = p_umap.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col = ctc.FeatureStore.UnitTable.(norm_var);
    [~, ~, iG] = unique(string(group_col), 'stable');
    for g = 1:max(iG)
        mask = (iG == g);
        X_raw(mask, :) = normalize(X_raw(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end

train_idx = logical(labels.umap_train_idx);
test_idx  = ~train_idx;

% ── Steps 2–4: global z-score, drop bad columns, scale ───────────────────────
[X_train, X_test, feat_names] = CellTypeClassifier.normalizeForClassification( ...
    X_raw(train_idx, :), X_raw(test_idx, :), feat_names);

% ── Supervised UMAP classification ───────────────────────────────────────────
[Y_pred, ~, ~, test_reduction, train_reduction] = supervisedUMAP(ctc, ...
    X_train, labels.sorted_y_train, feat_names, X_test);

ctc.Reduction.Train = train_reduction;
ctc.Reduction.Test  = test_reduction;

% Ground-truth labels for training units; predicted labels for all others
full_labels = nan(1, numel(ctc.UnitDataArray));
full_labels(labels.sorted_train_ids) = labels.sorted_y_train;
full_labels(test_idx)                = Y_pred;
ctc.UnitLabels = full_labels;

% ── kNN confidence in UMAP embedding space ─────────────────────────────────
% For each test unit, fraction of its k nearest training neighbors sharing
% the predicted label. Training units get confidence = 1.0 (ground truth).
conf_k = p_umap.ConfidenceK;
full_confidence = nan(1, numel(ctc.UnitDataArray));
full_confidence(labels.sorted_train_ids) = 1.0;

if ~isempty(test_reduction) && ~isempty(train_reduction)
    k_actual = min(conf_k, size(train_reduction, 1));
    [nn_idx, ~] = knnsearch(train_reduction, test_reduction, 'K', k_actual);
    test_global_idx = find(test_idx);
    for i = 1:numel(test_global_idx)
        neighbor_labels = labels.sorted_y_train(nn_idx(i,:));
        full_confidence(test_global_idx(i)) = sum(neighbor_labels == Y_pred(i)) / k_actual;
    end
end
ctc.UnitConfidence = full_confidence;

n_exc = sum(full_labels == 1, 'omitnan');
n_inh = sum(full_labels == 2, 'omitnan');
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, sum(isnan(full_labels)));
end
