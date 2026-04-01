function [target_labels, target_ctc] = applyTo(ctc, target_fs, target_ud, opts)
% APPLYTO  Apply a trained CellTypeClassifier to a different dataset.
%
% Takes the training labels from this CTC and classifies all units in the
% target FeatureStore + UnitData array, producing a new CellTypeClassifier.
%
<<<<<<< HEAD
% The trained CTC must have TrainLabels and UnitDataArray set (i.e., you must
% have run generateTrainLabels() or equivalent). All target units are treated
% as test — no responsive-unit identification is performed on the target.
=======
% Uses stored NormalizationParams from generateTrainLabels for the source data,
% ensuring the same feature space as the unsupervised and supervised UMAPs.
% Target data receives its own chip-level z-score, then the source's global
% z-score parameters are applied for cross-dataset consistency.
>>>>>>> claude/dreamy-varahamihira
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
%   target_fs - FeatureStore for the target dataset (provides metadata for normalisation)
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
<<<<<<< HEAD
assert(~isempty(ctc.UnitDataArray), ...
    'Source CTC has no UnitDataArray — run identifyResponsiveUnits() first.');

labels  = ctc.TrainLabels;

% ── Create target CTC (inherits all harmonization parameters) ──────────────
target_ctc = CellTypeClassifier(target_fs, target_ud, ctc.Parameters);
=======
assert(~isempty(ctc.UnitList), ...
    'Source CTC has no UnitList — run identifyResponsiveUnits() first.');
assert(~isempty(ctc.NormalizationParams), ...
    'Source CTC has no NormalizationParams — re-run generateTrainLabels() to store them.');

labels  = ctc.TrainLabels;
p_umap  = ctc.Parameters.UMAP;
np      = ctc.NormalizationParams;

% ── Create target CTC (inherits harmonization parameters) ──────────────
target_ctc = CellTypeClassifier(target_rg, ctc.Parameters);

% ── Build unit list for target RG ──────────────────────────────────────
target_units = [];
for r = 1:numel(target_rg.Recordings)
    rec = target_rg.Recordings(r);
    if ~isempty(rec.Units)
        target_units = [target_units, rec.Units]; %#ok<AGROW>
    end
end
assert(~isempty(target_units), 'Target RecordingGroup has no units.');
target_ctc.UnitList = target_units;
>>>>>>> claude/dreamy-varahamihira

fprintf('Applying trained classifier (%d train units) to %d target units\n', ...
    sum(labels.umap_train_idx), numel(target_ud));

% ── Extract training features (from source CTC, cached) ───────────────
[wf_train, acg_train, sr_train] = ctc.getOrExtract(ctc.UnitDataArray);
[X_train_raw, feat_names]       = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);

% ── Extract target features ────────────────────────────────────────────
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

<<<<<<< HEAD
% ── Step 1: per-group normalisation (source and target normalized independently) ─
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
=======
% ── Normalization ──────────────────────────────────────────────────────

% Step 1a: within-group z-score on training data (source RG)
rg_source = ctc.RecordingGroup;
[iG_train, ~] = rg_source.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx_train = unique(iG_train);
for g = 1:length(g_idx_train)
    mask = (iG_train == g_idx_train(g));
    X_train_raw(mask, :) = normalize(X_train_raw(mask, :));
end

% Step 1b: within-group z-score on target data (target RG)
[iG_target, ~] = target_rg.combineMetadataIndices(target_units, opts.NormalizationVar);
g_idx_target = unique(iG_target);
for g = 1:length(g_idx_target)
    mask = (iG_target == g_idx_target(g));
    X_target_raw(mask, :) = normalize(X_target_raw(mask, :));
>>>>>>> claude/dreamy-varahamihira
end

% Step 2: Apply stored global z-score from generateTrainLabels
% Source training data uses stored parameters for consistency
X_train_all = normalize(X_train_raw, 'center', np.mu_global, 'scale', np.sigma_global);
X_target = normalize(X_target_raw, 'center', np.mu_global, 'scale', np.sigma_global);

% Step 3: Remove same NaN columns as generateTrainLabels
X_train_all(:, np.nan_cols) = [];
X_target(:, np.nan_cols)    = [];
feat_names(np.nan_cols)     = [];

% Step 4: Apply same scaling
X_train_all = X_train_all ./ np.scale;
X_target    = X_target    ./ np.scale;

% Select training subset
X_train = X_train_all(labels.umap_train_idx, :);
y_train = labels.sorted_y_train;

<<<<<<< HEAD
% ── Steps 2–4: global z-score, drop bad columns, scale ────────────────────
[X_train, X_target, feat_names] = CellTypeClassifier.normalizeForClassification( ...
    X_train, X_target_raw, feat_names);
=======
% Handle any remaining NaN/Inf from target (different data may have edge cases)
bad_cols = any(isnan(X_train) | isinf(X_train), 1) | ...
           any(isnan(X_target) | isinf(X_target), 1);
if any(bad_cols)
    X_train(:, bad_cols)    = [];
    X_target(:, bad_cols)   = [];
    feat_names(bad_cols)    = [];
end
>>>>>>> claude/dreamy-varahamihira

% ── Supervised UMAP classification ─────────────────────────────────────
[target_labels, ~, ~, target_reduction, train_reduction, extras] = supervisedUMAP(ctc, ...
    X_train, y_train, feat_names, X_target);

target_ctc.UnitLabels      = target_labels;
target_ctc.TrainLabels     = labels;  % preserve source training info
target_ctc.NormalizationParams = np;  % preserve source normalization
target_ctc.Reduction.Train  = train_reduction;
target_ctc.Reduction.Test   = target_reduction;
target_ctc.Reduction.Extras = extras;

% ── kNN confidence for target units ─────────────────────────────────────
conf_k   = ctc.Parameters.UMAP.ConfidenceK;
k_actual = min(conf_k, size(train_reduction, 1));
[nn_idx, ~] = knnsearch(train_reduction, target_reduction, 'K', k_actual);
target_confidence = zeros(1, numel(target_labels));
for i = 1:numel(target_labels)
    neighbor_labels = y_train(nn_idx(i,:));
    target_confidence(i) = sum(neighbor_labels == target_labels(i)) / k_actual;
end
target_ctc.UnitConfidence = target_confidence;

n_exc = sum(target_labels == 1);
n_inh = sum(target_labels == 2);
fprintf('Target classification: %d excitatory, %d inhibitory (%d total)\n', ...
    n_exc, n_inh, numel(target_labels));
end
