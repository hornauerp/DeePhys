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

warning('CellTypeClassifier:deprecated', ...
    ['applyTo is deprecated. Include all units in one CTC instead — ' ...
     'generateTrainLabels() and classifyUnits() now operate on the full NxN UMAP graph, ' ...
     'so cross-dataset transfer is no longer needed.']);

assert(~isempty(ctc.TrainLabels), ...
    'Source CTC has no TrainLabels — run generateTrainLabels() first.');
assert(~isempty(ctc.UnitDataArray), ...
    'Source CTC has no UnitDataArray — run identifyResponsiveUnits() first.');
assert(~isempty(ctc.NormalizationParams), ...
    'Source CTC has no NormalizationParams — re-run generateTrainLabels() to store them.');

labels = ctc.TrainLabels;
p_umap = ctc.Parameters.UMAP;
np     = ctc.NormalizationParams;

% -- Create target CTC with harmonization matched to source ------------------
% The source's actual waveform window is stored in NormalizationParams (set by
% generateTrainLabels). If the source used "trim" edge mode, the window may be
% narrower than the requested [−PreTrough, +PostTrough]. Force the target to
% use the same window with "zero" edge mode so short waveforms are padded
% rather than triggering a different trim that changes the feature count.
target_params = ctc.Parameters;
if isfield(np, 'wf_pre_trough_ms')
    target_params.Harmonization.WaveformPreTrough  = np.wf_pre_trough_ms;
    target_params.Harmonization.WaveformPostTrough = np.wf_post_trough_ms;
end
target_params.Harmonization.WaveformEdgeMode = "zero";
target_ctc = CellTypeClassifier(target_fs, target_ud, target_params);

pre_ms  = target_params.Harmonization.WaveformPreTrough;
post_ms = target_params.Harmonization.WaveformPostTrough;
fprintf('Applying trained classifier (%d train units) to %d target units\n', ...
    sum(labels.umap_train_idx), numel(target_ud));
fprintf('applyTo: using waveform window [-%g, +%g] ms with zero-padding\n', ...
    pre_ms, post_ms);

% -- Extract training features (from source CTC, cached) ---------------------
[wf_train, acg_train, sr_train] = ctc.getOrExtract(ctc.UnitDataArray);
[X_train_raw, feat_names]       = buildFeatureMatrix(ctc, wf_train, acg_train, sr_train);

% -- Extract target features --------------------------------------------------
[wf_target, acg_target, sr_target] = target_ctc.getOrExtract(target_ud);
[X_target_raw, feat_names_target, aligned_wf, norm_acgs] = buildFeatureMatrix( ...
    target_ctc, wf_target, acg_target, sr_target);

% -- Verify feature dimensions match ------------------------------------------
if ~isequal(feat_names, feat_names_target)
    n_src = numel(feat_names);
    n_tgt = numel(feat_names_target);
    n_wf_src = sum(startsWith(feat_names, "Waveform"));
    n_wf_tgt = sum(startsWith(feat_names_target, "Waveform"));
    n_acg_src = sum(startsWith(feat_names, "FullACG"));
    n_acg_tgt = sum(startsWith(feat_names_target, "FullACG"));
    error(['CellTypeClassifier:featureMismatch  ' ...
        'Feature mismatch: source has %d features (%d ACG + %d waveform), ' ...
        'target has %d features (%d ACG + %d waveform). ' ...
        'This usually means different ACG parameters. ' ...
        'Check Harmonization.ACGBinSize and ACGLag.'], ...
        n_src, n_acg_src, n_wf_src, n_tgt, n_acg_tgt, n_wf_tgt);
end

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

% -- Dimension-normalised feature group weighting (same as classifyUnits) ----
is_acg = startsWith(feat_names, "FullACG") | startsWith(feat_names, "ACG");
is_wf  = startsWith(feat_names, "Waveform") | startsWith(feat_names, "ReferenceWaveform");
n_acg  = sum(is_acg);
n_wf   = sum(is_wf);
if n_acg > 0
    X_train_all(:, is_acg) = X_train_all(:, is_acg) * sqrt(p_umap.ACGWeight / n_acg);
    X_target(:, is_acg)    = X_target(:, is_acg)    * sqrt(p_umap.ACGWeight / n_acg);
end
if n_wf > 0
    X_train_all(:, is_wf) = X_train_all(:, is_wf) * sqrt(p_umap.WaveformWeight / n_wf);
    X_target(:, is_wf)    = X_target(:, is_wf)    * sqrt(p_umap.WaveformWeight / n_wf);
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

% -- AutoNNeighbors + minority class ceiling for supervised UMAP -------------
orig_sup_nneigh = p_umap.SupervisedNNeighbors;

if p_umap.AutoNNeighbors
    N_train_loc = size(X_train, 1);
    auto_sup_n  = max(p_umap.MinNNeighbors, round(sqrt(N_train_loc)));
    fprintf('Auto SupervisedNNeighbors: %d (N_train=%d)\n', auto_sup_n, N_train_loc);
    ctc.Parameters.UMAP.SupervisedNNeighbors = auto_sup_n;
end

n_inh   = sum(y_train == 2);
safe_nn = floor(n_inh / 5);
if safe_nn > 0 && ctc.Parameters.UMAP.SupervisedNNeighbors > safe_nn
    fprintf('Minority class ceiling: SupervisedNNeighbors %d -> %d (n_inh=%d)\n', ...
        ctc.Parameters.UMAP.SupervisedNNeighbors, safe_nn, n_inh);
    ctc.Parameters.UMAP.SupervisedNNeighbors = safe_nn;
end

% -- Supervised UMAP (embedding only, labels from kNN/graph below) -----------
[~, ~, ~, target_reduction, train_reduction, extras] = supervisedUMAP(ctc, ...
    X_train, y_train, feat_names, X_target);

% Restore overridden params
ctc.Parameters.UMAP.SupervisedNNeighbors = orig_sup_nneigh;

p_conf = ctc.Parameters.Classification;
classify_method = lower(string(p_conf.Method));

N_target = size(X_target, 1);
target_labels     = nan(1, N_target);
target_confidence = nan(1, N_target);

if classify_method == "graph"
    % ── Graph-based label propagation on supervised UMAP embedding ──────
    all_reduction = [train_reduction; target_reduction];
    N_all = size(all_reduction, 1);
    N_train_loc = size(train_reduction, 1);

    if p_umap.AutoConfidenceK
        graph_k = max(5, round(sqrt(N_all)));
    else
        graph_k = p_umap.ConfidenceK;
    end
    graph_k = min(graph_k, N_all - 1);

    [knn_idx, knn_dists] = knnsearch(all_reduction, all_reduction, 'K', graph_k + 1);
    knn_idx   = knn_idx(:, 2:end);
    knn_dists = knn_dists(:, 2:end);

    sigma = median(knn_dists(:, end));
    if sigma == 0; sigma = 1; end
    rows_g = repmat((1:N_all)', 1, graph_k);
    G = sparse(rows_g(:), knn_idx(:), exp(-knn_dists(:).^2 / (2 * sigma^2)), N_all, N_all);
    G = (G + G') / 2;
    row_sums = sum(G, 2);
    row_sums(row_sums == 0) = 1;
    W = G ./ row_sums;

    fprintf('Graph classification (transfer): %d nodes, k=%d\n', N_all, graph_k);

    F_label = ones(N_all, 2) * 0.5;
    for ti = 1:numel(y_train)
        if y_train(ti) == 1
            F_label(ti, :) = [1, 0];
        else
            F_label(ti, :) = [0, 1];
        end
    end

    max_iter = p_conf.GraphMaxIter;
    tol      = p_conf.GraphConvergenceTol;
    for iter = 1:max_iter
        F_new = W * F_label;
        F_row_sum = sum(F_new, 2);
        F_row_sum(F_row_sum == 0) = 1;
        F_new = F_new ./ F_row_sum;
        for ti = 1:numel(y_train)
            if y_train(ti) == 1
                F_new(ti, :) = [1, 0];
            else
                F_new(ti, :) = [0, 1];
            end
        end
        if max(abs(F_new(:) - F_label(:))) < tol
            fprintf('Graph label propagation converged at iteration %d\n', iter);
            F_label = F_new;
            break
        end
        F_label = F_new;
    end

    for i = 1:N_target
        gi = N_train_loc + i;
        if F_label(gi, 1) >= F_label(gi, 2)
            target_labels(i) = 1;
            target_confidence(i) = F_label(gi, 1);
        else
            target_labels(i) = 2;
            target_confidence(i) = F_label(gi, 2);
        end
    end

else
    % ── Feature-space kNN classification (same as classifyUnits) ────────
    k_actual = min(conf_k, size(X_train, 1));
    [nn_idx, nn_dists] = knnsearch(X_train, X_target, 'K', k_actual);
    for i = 1:N_target
        neighbor_labels = y_train(nn_idx(i, :));
        d       = nn_dists(i, :);
        eps_d   = max(d) * 1e-6;
        if eps_d == 0; eps_d = 1e-10; end
        w       = 1 ./ (d + eps_d);
        w_total = sum(w);
        w_exc   = sum(w(neighbor_labels == 1)) / w_total;
        w_inh   = sum(w(neighbor_labels == 2)) / w_total;
        if w_exc >= w_inh
            target_labels(i) = 1;
            target_confidence(i) = w_exc;
        else
            target_labels(i) = 2;
            target_confidence(i) = w_inh;
        end
    end
end

% -- Store results on target CTC --------------------------------------------
target_ctc.UnitLabels          = target_labels;
target_ctc.UnitConfidence      = target_confidence;
target_labels_transfer         = labels;
target_labels_transfer.is_transfer = true;
target_ctc.TrainLabels         = target_labels_transfer;
target_ctc.NormalizationParams = np;
target_ctc.Reduction           = struct('Train', train_reduction, ...
                                        'Test',  target_reduction, ...
                                        'Extras', extras, ...
                                        'X_train', X_train, ...
                                        'X_test',  X_target, ...
                                        'feat_names', feat_names);

n_exc = sum(target_labels == 1);
n_inh = sum(target_labels == 2);
fprintf('Target classification: %d excitatory, %d inhibitory (%d total)\n', ...
    n_exc, n_inh, numel(target_labels));
end
