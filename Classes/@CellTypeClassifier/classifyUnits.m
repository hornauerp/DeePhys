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

% -- Extract features via harmonized path -------------------------------------
[wf, acg, sr]                        = ctc.getOrExtract(ctc.UnitDataArray);
[X_raw, feat_names, aligned_wf, norm_acgs] = buildFeatureMatrix(ctc, wf, acg, sr);

% Store harmonized data for downstream inspection / plotting
ctc.HarmonizedWaveforms = aligned_wf;
ctc.HarmonizedACGs      = norm_acgs;
ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% -- Step 1: per-group normalisation using FeatureStore table column ----------
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

% -- AutoNNeighbors for supervised UMAP (Phase 1) -----------------------------
orig_sup_nneigh = p_umap.SupervisedNNeighbors;
if p_umap.AutoNNeighbors
    N_train = size(X_train, 1);
    auto_sup_n = max(p_umap.MinNNeighbors, round(sqrt(N_train)));
    fprintf('Auto SupervisedNNeighbors: %d (N_train=%d)\n', auto_sup_n, N_train);
    ctc.Parameters.UMAP.SupervisedNNeighbors = auto_sup_n;
end

% -- AutoTargetWeight based on train/test distributional divergence (Phase 3) --
orig_target_weight = p_umap.TargetWeight;
if p_umap.AutoTargetWeight
    n_pcs = min(5, size(X_train, 2));
    [coeff, ~] = pca(X_train, 'NumComponents', n_pcs);
    pc_train   = X_train * coeff;
    pc_test    = X_test  * coeff;
    ks_stats   = zeros(1, n_pcs);
    for dim = 1:n_pcs
        [~, ~, ks_stats(dim)] = kstest2(pc_train(:, dim), pc_test(:, dim));
    end
    divergence     = mean(ks_stats);
    target_weight  = max(0.1, min(0.9, 0.5 * (1 - divergence)));
    fprintf('Auto TargetWeight: %.3f (mean KS divergence=%.3f)\n', target_weight, divergence);
    ctc.Parameters.UMAP.TargetWeight = target_weight;
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

% -- kNN confidence in UMAP embedding space ----------------------------------
% For each test unit: fraction of its k nearest training neighbors sharing
% the predicted label. Training units get confidence = 1.0 (ground truth).
if p_umap.AutoConfidenceK
    conf_k = max(5, round(sqrt(size(train_reduction, 1))));
    fprintf('Auto ConfidenceK: %d (N_train=%d)\n', conf_k, size(train_reduction, 1));
else
    conf_k = p_umap.ConfidenceK;
end

full_confidence = nan(1, numel(ctc.UnitDataArray));
full_confidence(labels.sorted_train_ids) = 1.0;

if ~isempty(test_reduction) && ~isempty(train_reduction)
    k_actual = min(conf_k, size(train_reduction, 1));
    [nn_idx, ~] = knnsearch(train_reduction, test_reduction, 'K', k_actual);
    test_global_idx = find(test_idx);
    for i = 1:numel(test_global_idx)
        neighbor_labels = labels.sorted_y_train(nn_idx(i, :));
        full_confidence(test_global_idx(i)) = sum(neighbor_labels == Y_pred(i)) / k_actual;
    end
end

% -- Assemble labels ----------------------------------------------------------
full_labels = nan(1, numel(ctc.UnitDataArray));
full_labels(labels.sorted_train_ids) = labels.sorted_y_train;
full_labels(test_idx)                = Y_pred;

% Optional confidence threshold: mark low-confidence test predictions as NaN
p_conf = ctc.Parameters.Classification;
if p_conf.UseConfidenceThreshold
    low_conf = full_confidence < p_conf.ConfidenceThreshold & ~train_idx;
    full_labels(low_conf) = NaN;
    fprintf('Confidence threshold %.2f: %d units set to NaN\n', ...
        p_conf.ConfidenceThreshold, sum(low_conf));
end

ctc.UnitLabels     = full_labels;
ctc.UnitConfidence = full_confidence;

n_exc = sum(full_labels == 1, 'omitnan');
n_inh = sum(full_labels == 2, 'omitnan');
n_unc = sum(isnan(full_labels));
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, n_unc);
end
