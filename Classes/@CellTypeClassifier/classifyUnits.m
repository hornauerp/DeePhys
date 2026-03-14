function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Feature extraction uses ctc.Parameters.Harmonization so DeePhys ACGs are
% recomputed at the configured bin size and lag (default: 0.5 ms / 100 ms).
%
% Requires ctc.TrainLabels to be set (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyUnits()');

labels = ctc.TrainLabels;

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr]       = extractUnitWaveformsAndACGs(ctc, ctc.UnitList);
[X_raw, feat_names] = buildFeatureMatrix(ctc, wf, acg, sr);

X_train_raw = X_raw(labels.umap_train_idx, :);
X_test_raw  = X_raw(labels.umap_test_idx,  :);

[Y_pred, train_reduction, test_reduction] = classifyWithFeatureMatrices(ctc, ...
    X_train_raw, labels.sorted_y_train, X_test_raw, feat_names);

ctc.Reduction.Train = train_reduction;
ctc.Reduction.Test  = test_reduction;

full_labels = nan(1, length(ctc.UnitList));
full_labels(labels.sorted_train_ids) = labels.sorted_y_train;
full_labels(labels.umap_test_idx)    = Y_pred;
ctc.UnitLabels = full_labels;

n_exc = sum(full_labels == 1, 'omitnan');
n_inh = sum(full_labels == 2, 'omitnan');
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, sum(isnan(full_labels)));
end
