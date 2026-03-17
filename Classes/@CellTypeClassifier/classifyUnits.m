function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Feature extraction uses ctc.Parameters.Harmonization so DeePhys ACGs are
% recomputed at the configured bin size and lag (default: 0.5 ms / 100 ms).
%
% Normalisation matches the main-branch pipeline:
%   1. Z-score within each NormalizationVar group (e.g. ChipID)
%   2. Global z-score: fit on train, apply train stats to test
%   3. Remove NaN columns, scale by max(abs(train))
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
rg     = ctc.RecordingGroup;

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr]                        = ctc.getOrExtract(ctc.UnitList);
[X_raw, feat_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% ── Chip-level normalisation (matching main branch prepareInputMatrix) ───────
X_raw(isnan(X_raw)) = 0;

% Step 1: z-score within each NormalizationVar group (all units per group)
[iG, G] = rg.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx = unique(iG);
for g = 1:length(g_idx)
    mask = (iG == g_idx(g));
    X_raw(mask, :) = normalize(X_raw(mask, :));
end

train_idx = logical(labels.umap_train_idx);
test_idx  = logical(labels.umap_test_idx);

% Step 2: global z-score — fit on train, apply train stats to test
if length(G) > 1
    [X_train, mu, sigma] = normalize(X_raw(train_idx, :));
    X_test = normalize(X_raw(test_idx, :), 'center', mu, 'scale', sigma);
else
    X_train = X_raw(train_idx, :);
    X_test  = X_raw(test_idx, :);
end

% Step 3: remove NaN columns (from train or test)
nan_cols = any(isnan(X_train), 1) | any(isnan(X_test), 1);
X_train(:, nan_cols)   = [];
X_test(:, nan_cols)    = [];
feat_names(nan_cols)   = [];

% Step 4: scale by max(abs(train))
scale = max(abs(X_train), [], 1);
scale(scale == 0) = 1;
X_train = X_train ./ scale;
X_test  = X_test  ./ scale;

% ── Supervised UMAP classification ───────────────────────────────────────────
[Y_pred, ~, ~, test_reduction, train_reduction] = supervisedUMAP(ctc, ...
    X_train, labels.sorted_y_train, feat_names, X_test);

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
