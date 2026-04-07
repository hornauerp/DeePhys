function [target_labels, target_ctc] = applyTo(ctc, target_fs, target_ud, opts)
% APPLYTO  Apply a trained CellTypeClassifier to a different dataset.
%
% Takes the training labels from this CTC and classifies all units in the
% target FeatureStore + UnitData array, producing a new CellTypeClassifier.
%
% Uses stored NormalizationParams from generateTrainLabels for the source data,
% ensuring the same feature space as the unsupervised and supervised UMAPs.
% Target data receives its own chip-level z-score, then the source's global
% z-score parameters are applied for cross-dataset consistency.
%
% USAGE:
%   ctc_source.identifyResponsiveUnits();
%   ctc_source.generateTrainLabels();
%
%   % Apply to a completely different dataset:
%   [labels, ctc_target] = ctc_source.applyTo(fs_target, ud_target);
%
%   % Apply with a different normalization variable on the target:
%   [labels, ctc_target] = ctc_source.applyTo(fs_target, ud_target, ...
%                              NormalizationVar="ChipID");
%
% INPUTS:
%   ctc       - trained CellTypeClassifier (source)
%   target_fs - FeatureStore for the target dataset
%   target_ud - UnitData array for the target dataset
%
% NAME-VALUE:
%   NormalizationVar - metadata field for within-group z-score on target
%                      (default: ctc.Parameters.UMAP.NormalizationVar)
%
% OUTPUTS:
%   target_labels - (1 x N_target) 1=excitatory, 2=inhibitory
%   target_ctc    - CellTypeClassifier for the target dataset with labels,
%                   reductions, and harmonized data populated

arguments
    ctc        CellTypeClassifier
    target_fs  FeatureStore
    target_ud  UnitData
    opts.NormalizationVar string = ctc.Parameters.UMAP.NormalizationVar
end

assert(~isempty(ctc.TrainLabels), ...
    'Source CTC has no TrainLabels — run generateTrainLabels() first.');
assert(~isempty(ctc.UnitDataArray), ...
    'Source CTC has no UnitDataArray — run identifyResponsiveUnits() first.');
assert(~isempty(ctc.NormalizationParams), ...
    'Source CTC has no NormalizationParams — re-run generateTrainLabels() to store them.');

labels = ctc.TrainLabels;
p_umap = ctc.Parameters.UMAP;
np     = ctc.NormalizationParams;

% -- Create target CTC (inherits all harmonization parameters) ---------------
target_ctc = CellTypeClassifier(target_fs, target_ud, ctc.Parameters);

fprintf('Applying trained classifier (%d train units) to %d target units\n', ...
    sum(labels.umap_train_idx), numel(target_ud));

% -- Extract training features (from source CTC, cached) ---------------------
[wf_train, acg_train, sr_train] = ctc.getOrExtract(ctc.UnitDataArray);
[X_train_raw, feat_names]       = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);

% -- Extract target features --------------------------------------------------
[wf_target, acg_target, sr_target] = target_ctc.getOrExtract(target_ud);
[X_target_raw, feat_names_target, aligned_wf, norm_acgs] = buildFeatureMatrix( ...
    target_ctc, wf_target, acg_target, sr_target);

assert(isequal(feat_names, feat_names_target), ...
    'Feature mismatch: source (%d) vs target (%d). Check Harmonization parameters.', ...
    numel(feat_names), numel(feat_names_target));

% Store harmonized data on target CTC
target_ctc.HarmonizedWaveforms = aligned_wf;
target_ctc.HarmonizedACGs      = norm_acgs;
target_ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% -- Step 1: per-group normalisation (source and target normalized independently) --
X_train_raw(isnan(X_train_raw))   = 0;
X_target_raw(isnan(X_target_raw)) = 0;

norm_var = opts.NormalizationVar;

if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col = ctc.FeatureStore.UnitTable.(norm_var);
    [~, ~, iG] = unique(string(group_col), 'stable');
    for g = 1:max(iG)
        mask = (iG == g);
        X_train_raw(mask, :) = normalize(X_train_raw(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in source UnitTable — per-group z-score skipped.', norm_var);
end

if ~isempty(norm_var) && ismember(norm_var, string(target_fs.UnitTable.Properties.VariableNames))
    group_col = target_fs.UnitTable.(norm_var);
    [~, ~, iG] = unique(string(group_col), 'stable');
    for g = 1:max(iG)
        mask = (iG == g);
        X_target_raw(mask, :) = normalize(X_target_raw(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in target UnitTable — per-group z-score skipped.', norm_var);
end

% -- Steps 2-5: apply stored normalisation from generateTrainLabels ----------
% Using source's stored NormalizationParams ensures train and test are in the
% same feature space as the unsupervised UMAP used for label generation.
X_train_all = normalize(X_train_raw, 'center', np.mu_global, 'scale', np.sigma_global);
X_target    = normalize(X_target_raw, 'center', np.mu_global, 'scale', np.sigma_global);

X_train_all(:, np.nan_cols) = [];
X_target(:, np.nan_cols)    = [];
feat_names(np.nan_cols)     = [];

X_train_all = X_train_all ./ np.scale;
X_target    = X_target    ./ np.scale;

% Apply feature selection mask if generateTrainLabels ran FeatureSelection
if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X_train_all = X_train_all(:, np.feature_selection_mask);
    X_target    = X_target(:, np.feature_selection_mask);
    feat_names  = feat_names(np.feature_selection_mask);
end

% Handle any residual NaN/Inf (edge cases from target data with different distribution)
bad_cols = any(isnan(X_train_all) | isinf(X_train_all), 1) | ...
           any(isnan(X_target)    | isinf(X_target),    1);
if any(bad_cols)
    X_train_all(:, bad_cols) = [];
    X_target(:, bad_cols)    = [];
    feat_names(bad_cols)     = [];
end

% Select training subset
X_train = X_train_all(labels.umap_train_idx, :);
y_train = labels.sorted_y_train;

% -- AutoConfidenceK for target (Phase 1) ------------------------------------
orig_conf_k = p_umap.ConfidenceK;
if p_umap.AutoConfidenceK
    conf_k = max(5, round(sqrt(size(X_train, 1))));
    fprintf('Auto ConfidenceK: %d (N_train=%d)\n', conf_k, size(X_train, 1));
else
    conf_k = orig_conf_k;
end

% -- Supervised UMAP classification ------------------------------------------
[target_labels, ~, ~, target_reduction, train_reduction, extras] = supervisedUMAP(ctc, ...
    X_train, y_train, feat_names, X_target);

target_ctc.UnitLabels          = target_labels;
target_ctc.TrainLabels         = labels;  % preserve source training info
target_ctc.NormalizationParams = np;      % preserve source normalization
target_ctc.Reduction.Train     = train_reduction;
target_ctc.Reduction.Test      = target_reduction;
target_ctc.Reduction.Extras    = extras;

% -- Distance-weighted kNN confidence for target units -----------------------
k_actual = min(conf_k, size(train_reduction, 1));
[nn_idx, nn_dists] = knnsearch(train_reduction, target_reduction, 'K', k_actual);
target_confidence = zeros(1, numel(target_labels));
for i = 1:numel(target_labels)
    neighbor_labels = y_train(nn_idx(i, :));
    d = nn_dists(i, :);
    eps_d  = max(d) * 1e-6;
    if eps_d == 0; eps_d = 1e-10; end
    w      = 1 ./ (d + eps_d);
    w_total = sum(w);
    w_match = sum(w(neighbor_labels == target_labels(i)));
    target_confidence(i) = w_match / w_total;
end
target_ctc.UnitConfidence = target_confidence;

n_exc = sum(target_labels == 1);
n_inh = sum(target_labels == 2);
fprintf('Target classification: %d excitatory, %d inhibitory (%d total)\n', ...
    n_exc, n_inh, numel(target_labels));
end
