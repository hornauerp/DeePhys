function [target_labels, target_ctc] = applyTo(ctc, target_rg, opts)
% APPLYTO  Apply a trained CellTypeClassifier to a different RecordingGroup.
%
% Takes the training labels from this CTC and classifies all units in the
% target RecordingGroup, producing a new CellTypeClassifier for the target.
%
% The trained CTC must have TrainLabels and UnitList set (i.e., you must
% have run generateTrainLabels() or equivalent). The target RG does not
% need responsive unit identification — all its units are treated as test.
%
% USAGE:
%   ctc_source = CellTypeClassifier(rg_source, params);
%   ctc_source.identifyResponsiveUnits();
%   ctc_source.generateTrainLabels();
%
%   % Apply to a completely different RecordingGroup:
%   [labels, ctc_target] = ctc_source.applyTo(rg_target);
%
%   % Apply with metadata-based normalization on the target:
%   [labels, ctc_target] = ctc_source.applyTo(rg_target, NormalizationVar="ChipID");
%
% INPUTS:
%   ctc       - trained CellTypeClassifier (source)
%   target_rg - RecordingGroup to classify
%
% NAME-VALUE:
%   NormalizationVar - metadata field for within-group z-score on target
%                      (default: ctc.Parameters.UMAP.NormalizationVar)
%
% OUTPUTS:
%   target_labels - (1 x N_target) 1=excitatory, 2=inhibitory
%   target_ctc    - CellTypeClassifier for the target RG with labels,
%                   reductions, and harmonized data populated

arguments
    ctc        CellTypeClassifier
    target_rg  RecordingGroup
    opts.NormalizationVar string = ctc.Parameters.UMAP.NormalizationVar
end

assert(~isempty(ctc.TrainLabels), ...
    'Source CTC has no TrainLabels — run generateTrainLabels() first.');
assert(~isempty(ctc.UnitList), ...
    'Source CTC has no UnitList — run identifyResponsiveUnits() first.');

labels  = ctc.TrainLabels;
p_umap  = ctc.Parameters.UMAP;

% ── Create target CTC (inherits harmonization parameters) ──────────────
target_ctc = CellTypeClassifier(target_rg, ctc.Parameters);

% ── Build unit list for target RG ──────────────────────────────────────
% Collect all units across all recordings (same as identifyResponsiveUnits)
target_units = [];
for r = 1:numel(target_rg.Recordings)
    rec = target_rg.Recordings(r);
    if ~isempty(rec.Units)
        target_units = [target_units, rec.Units]; %#ok<AGROW>
    end
end
assert(~isempty(target_units), 'Target RecordingGroup has no units.');
target_ctc.UnitList = target_units;

fprintf('Applying trained classifier (%d train units) to %d target units\n', ...
    sum(labels.umap_train_idx), numel(target_units));

% ── Extract training features (from source CTC, cached) ───────────────
[wf_train, acg_train, sr_train] = ctc.getOrExtract(ctc.UnitList);
[X_train_raw, feat_names]       = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);

% ── Extract target features ────────────────────────────────────────────
[wf_target, acg_target, sr_target] = extractUnitWaveformsAndACGs(target_ctc, target_units);
[X_target_raw, feat_names_target, aligned_wf, norm_acgs] = buildFeatureMatrix(target_ctc, wf_target, acg_target, sr_target);

assert(isequal(feat_names, feat_names_target), ...
    'Feature mismatch: source (%d) vs target (%d). Check Harmonization parameters.', ...
    numel(feat_names), numel(feat_names_target));

% Store harmonized data on target CTC
target_ctc.HarmonizedWaveforms = aligned_wf;
target_ctc.HarmonizedACGs      = norm_acgs;
target_ctc.HarmonizedSR        = ctc.Parameters.Harmonization.WaveformTargetSamplingRate;

% ── Normalization ──────────────────────────────────────────────────────
X_train_raw(isnan(X_train_raw)) = 0;
X_target_raw(isnan(X_target_raw)) = 0;

% Step 1a: within-group z-score on training data (source RG)
rg_source = ctc.RecordingGroup;
[iG_train, G_train] = rg_source.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx_train = unique(iG_train);
for g = 1:length(g_idx_train)
    mask = (iG_train == g_idx_train(g));
    X_train_raw(mask, :) = normalize(X_train_raw(mask, :));
end

% Step 1b: within-group z-score on target data (target RG)
[iG_target, G_target] = target_rg.combineMetadataIndices(target_units, opts.NormalizationVar);
g_idx_target = unique(iG_target);
for g = 1:length(g_idx_target)
    mask = (iG_target == g_idx_target(g));
    X_target_raw(mask, :) = normalize(X_target_raw(mask, :));
end

% Select training subset
X_train = X_train_raw(labels.umap_train_idx, :);
y_train = labels.sorted_y_train;

% Step 2: global z-score — fit on train, apply to target
[X_train, mu, sigma] = normalize(X_train);
X_target = normalize(X_target_raw, 'center', mu, 'scale', sigma);

% Step 3: remove NaN/Inf columns
bad_cols = any(isnan(X_train) | isinf(X_train), 1) | ...
           any(isnan(X_target) | isinf(X_target), 1);
X_train(:, bad_cols)    = [];
X_target(:, bad_cols)   = [];
feat_names(bad_cols)    = [];

% Step 4: scale by max(abs(train))
scale = max(abs(X_train), [], 1);
scale(scale == 0) = 1;
X_train  = X_train  ./ scale;
X_target = X_target ./ scale;

% ── Supervised UMAP classification ─────────────────────────────────────
[target_labels, ~, ~, target_reduction, train_reduction] = supervisedUMAP(ctc, ...
    X_train, y_train, feat_names, X_target);

target_ctc.UnitLabels      = target_labels;
target_ctc.TrainLabels     = labels;  % preserve source training info
target_ctc.Reduction.Train = train_reduction;
target_ctc.Reduction.Test  = target_reduction;

n_exc = sum(target_labels == 1);
n_inh = sum(target_labels == 2);
fprintf('Target classification: %d excitatory, %d inhibitory (%d total)\n', ...
    n_exc, n_inh, numel(target_labels));
end
