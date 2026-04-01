function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Feature extraction uses ctc.Parameters.Harmonization so DeePhys ACGs are
% recomputed at the configured bin size and lag (default: 0.5 ms / 100 ms).
%
% Normalisation uses parameters stored by generateTrainLabels for consistency:
%   1. Z-score within each NormalizationVar group (e.g. ChipID)
%   2. Apply stored global z-score (mu, sigma from generateTrainLabels subset)
%   3. Remove same NaN columns as generateTrainLabels
%   4. Apply same max(abs) scaling
%
% This ensures the supervised UMAP sees the same feature space as the
% unsupervised UMAP used for label generation.
%
% Requires ctc.TrainLabels and ctc.NormalizationParams to be set
% (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyUnits()');
assert(~isempty(ctc.NormalizationParams), ...
    'NormalizationParams not set — re-run generateTrainLabels() to store them.');

labels = ctc.TrainLabels;
p_umap = ctc.Parameters.UMAP;
rg     = ctc.RecordingGroup;
np     = ctc.NormalizationParams;

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr]                        = ctc.getOrExtract(ctc.UnitList);
[X_raw, feat_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% ── Chip-level normalisation (same as generateTrainLabels) ────────────────────
[iG, ~] = rg.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx = unique(iG);
for g = 1:length(g_idx)
    mask = (iG == g_idx(g));
    X_raw(mask, :) = normalize(X_raw(mask, :));
end

% ── Apply stored normalisation from generateTrainLabels ───────────────────────
% This ensures the supervised UMAP operates in the same feature space as the
% unsupervised UMAP used for label generation.
X_all = normalize(X_raw, 'center', np.mu_global, 'scale', np.sigma_global);
X_all(:, np.nan_cols)  = [];
feat_names(np.nan_cols) = [];
X_all = X_all ./ np.scale;

train_idx = logical(labels.umap_train_idx);
test_idx  = logical(labels.umap_test_idx);

X_train = X_all(train_idx, :);
X_test  = X_all(test_idx, :);

% ── Supervised UMAP classification ───────────────────────────────────────────
[Y_pred, ~, ~, test_reduction, train_reduction, extras] = supervisedUMAP(ctc, ...
    X_train, labels.sorted_y_train, feat_names, X_test);

ctc.Reduction.Train  = train_reduction;
ctc.Reduction.Test   = test_reduction;
ctc.Reduction.Extras = extras;

% ── Confidence estimation (kNN distance to same-class training points) ───────
% For each test unit, compute distance to the k nearest training points of the
% predicted class. Units far from their assigned class in embedding space are
% unreliable. Confidence = 1 - normalised distance (0 = no confidence, 1 = max).
k_conf = min(5, size(train_reduction, 1));
confidence = nan(1, length(Y_pred));
for cls = unique(Y_pred(~isnan(Y_pred)))'
    cls_train_mask = (labels.sorted_y_train == cls);
    cls_test_mask  = (Y_pred == cls);
    if sum(cls_train_mask) == 0 || sum(cls_test_mask) == 0; continue; end
    D = pdist2(test_reduction(cls_test_mask, :), ...
               train_reduction(cls_train_mask, :), 'euclidean');
    k_use = min(k_conf, size(D, 2));
    D_sorted = sort(D, 2, 'ascend');
    confidence(cls_test_mask) = mean(D_sorted(:, 1:k_use), 2)';
end
% Normalise to [0, 1] where 1 = most confident (closest)
if any(~isnan(confidence))
    max_d = max(confidence(~isnan(confidence)));
    if max_d > 0
        confidence = 1 - confidence / max_d;
    end
end

% ── Assemble labels ──────────────────────────────────────────────────────────
full_labels     = nan(1, length(ctc.UnitList));
full_confidence = zeros(1, length(ctc.UnitList));

full_labels(labels.sorted_train_ids) = labels.sorted_y_train;
full_confidence(labels.sorted_train_ids) = 1;  % training labels are certain
full_labels(labels.umap_test_idx)    = Y_pred;
full_confidence(labels.umap_test_idx) = confidence;

% Apply optional confidence threshold — low-confidence predictions become NaN
p_conf = ctc.Parameters.Classification;
if p_conf.UseConfidenceThreshold
    low_conf = full_confidence < p_conf.ConfidenceThreshold & ...
               ~labels.umap_train_idx;  % never filter training labels
    full_labels(low_conf) = NaN;
end

ctc.UnitLabels  = full_labels;
ctc.Confidence  = full_confidence;

n_exc = sum(full_labels == 1, 'omitnan');
n_inh = sum(full_labels == 2, 'omitnan');
n_unc = sum(isnan(full_labels));
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, n_unc);
end
