function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Normalisation uses parameters stored by generateTrainLabels for consistency:
%   1. Z-score within each NormalizationVar group (e.g. ChipID)
%   2. Apply stored global z-score (mu, sigma from generateTrainLabels subset)
%   3. Remove same NaN columns as generateTrainLabels
%   4. Apply same max(abs) scaling
%   5. Apply feature_selection_mask if FeatureSelection was enabled
%
% This ensures the supervised UMAP sees the same feature space as the
% unsupervised UMAP used for label generation.
%
% Requires ctc.TrainLabels and ctc.NormalizationParams to be set
% (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)
%       ctc.UnitConfidence  kNN vote fraction in [0,1]; 1.0 for training units

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyUnits()');
assert(~isempty(ctc.NormalizationParams), ...
    'NormalizationParams not set — re-run generateTrainLabels() to store them.');

labels = ctc.TrainLabels;
p_umap = ctc.Parameters.UMAP;
np     = ctc.NormalizationParams;

% -- Use same unique-unit deduplication as generateTrainLabels ----------------
% labels.all_to_unique maps each full-array row to its unique-unit index.
% We rebuild unique_ud via uniqueUnitMap (same logic, same baseline selection).
[unique_ud, all_to_unique, ~] = ctc.uniqueUnitMap();
n_unique     = numel(unique_ud);
n_units_full = numel(ctc.UnitDataArray);

% Recover unique_to_rep (representative row index for each unique unit)
unique_to_rep = zeros(1, n_unique);
for u = 1:n_unique
    rows = find(all_to_unique == u);
    unique_to_rep(u) = rows(1);
end

% -- Extract features via harmonized path (unique units only) -----------------
[wf, acg, sr]                        = ctc.getOrExtract(unique_ud);
[X_raw, feat_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% -- Step 1: per-group normalisation using representative FeatureStore rows ---
X_raw(isnan(X_raw)) = 0;

norm_var = p_umap.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col        = ctc.FeatureStore.UnitTable.(norm_var);
    group_col_unique = group_col(unique_to_rep);
    [~, ~, iG] = unique(string(group_col_unique), 'stable');
    for g = 1:max(iG)
        mask = (iG == g);
        X_raw(mask, :) = normalize(X_raw(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end

% -- Steps 2-5: apply stored normalisation from generateTrainLabels -----------
% Ensures the supervised UMAP operates in the same feature space as the
% unsupervised UMAP used for label generation.
X_all = normalize(X_raw, 'center', np.mu_global, 'scale', np.sigma_global);
X_all(:, np.nan_cols)  = [];
feat_names(np.nan_cols) = [];
X_all = X_all ./ np.scale;

% Apply feature selection mask if generateTrainLabels ran FeatureSelection
if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X_all      = X_all(:, np.feature_selection_mask);
    feat_names = feat_names(np.feature_selection_mask);
end

train_idx = logical(labels.umap_train_idx);
test_idx  = ~train_idx;
X_train   = X_all(train_idx, :);
X_test    = X_all(test_idx, :);

% -- AutoNNeighbors + per-scenario TargetWeight for supervised UMAP -----------
orig_sup_nneigh    = p_umap.SupervisedNNeighbors;
orig_target_weight = p_umap.TargetWeight;

if p_umap.AutoNNeighbors
    N_train = size(X_train, 1);
    auto_sup_n = max(p_umap.MinNNeighbors, round(sqrt(N_train)));
    fprintf('Auto SupervisedNNeighbors: %d (N_train=%d)\n', auto_sup_n, N_train);
    ctc.Parameters.UMAP.SupervisedNNeighbors = auto_sup_n;
end

% Metadata labels are clean ground truth — higher TargetWeight is appropriate
% (more label-driven topology). Inferred drug-response labels use the default
% (0.5) because label noise warrants a more data-driven embedding.
if isfield(labels, 'has_explicit_ce') && labels.has_explicit_ce
    ctc.Parameters.UMAP.TargetWeight = p_umap.MetadataTargetWeight;
    fprintf('Metadata labels: using TargetWeight=%.2f\n', p_umap.MetadataTargetWeight);
end

% -- Supervised UMAP classification ------------------------------------------
[Y_pred, ~, ~, test_reduction, train_reduction, extras] = supervisedUMAP(ctc, ...
    X_train, labels.sorted_y_train, feat_names, X_test);

% Restore overridden params
ctc.Parameters.UMAP.SupervisedNNeighbors = orig_sup_nneigh;
ctc.Parameters.UMAP.TargetWeight         = orig_target_weight;

ctc.Reduction.Train  = train_reduction;
ctc.Reduction.Test   = test_reduction;
ctc.Reduction.Extras = extras;

% -- Distance-weighted kNN confidence in UMAP embedding space ----------------
% For each test unit: sum of inverse-distance weights from k nearest training
% neighbors sharing the predicted label, divided by total weight. Closer
% neighbors contribute more, giving a softer signal near the decision boundary.
% Training units get confidence = 1.0 (ground truth).
if p_umap.AutoConfidenceK
    conf_k = max(5, round(sqrt(size(train_reduction, 1))));
    fprintf('Auto ConfidenceK: %d (N_train=%d)\n', conf_k, size(train_reduction, 1));
else
    conf_k = p_umap.ConfidenceK;
end

unique_confidence = nan(1, n_unique);
unique_confidence(labels.sorted_train_ids) = 1.0;

test_global_idx = find(test_idx);
if ~isempty(test_reduction) && ~isempty(train_reduction)
    k_actual = min(conf_k, size(train_reduction, 1));
    [nn_idx, nn_dists] = knnsearch(train_reduction, test_reduction, 'K', k_actual);
    for i = 1:numel(test_global_idx)
        neighbor_labels = labels.sorted_y_train(nn_idx(i, :));
        d = nn_dists(i, :);
        % Inverse-distance weights: w_j = 1 / (d_j + eps) where eps avoids
        % division by zero for co-located points.
        eps_d  = max(d) * 1e-6;
        if eps_d == 0; eps_d = 1e-10; end
        w      = 1 ./ (d + eps_d);
        w_total = sum(w);
        w_match = sum(w(neighbor_labels == Y_pred(i)));
        unique_confidence(test_global_idx(i)) = w_match / w_total;
    end
end

% -- Assemble labels in unique-unit space -------------------------------------
unique_labels = nan(1, n_unique);
unique_labels(labels.sorted_train_ids) = labels.sorted_y_train;
unique_labels(test_idx)                = Y_pred;

% Optional confidence threshold
p_conf = ctc.Parameters.Classification;
if p_conf.UseConfidenceThreshold
    low_conf = unique_confidence < p_conf.ConfidenceThreshold & ~train_idx;
    unique_labels(low_conf) = NaN;
    fprintf('Confidence threshold %.2f: %d units set to NaN\n', ...
        p_conf.ConfidenceThreshold, sum(low_conf));
end

% -- Broadcast unique-unit results back to all (unit x recording) rows --------
full_labels     = unique_labels(all_to_unique);
full_confidence = unique_confidence(all_to_unique);

ctc.UnitLabels     = full_labels;
ctc.UnitConfidence = full_confidence;

n_exc = sum(unique_labels == 1, 'omitnan');
n_inh = sum(unique_labels == 2, 'omitnan');
n_unc = sum(isnan(unique_labels));
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, n_unc);
end
