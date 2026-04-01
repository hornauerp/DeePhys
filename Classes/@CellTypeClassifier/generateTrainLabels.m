function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels using UMAP embedding and isolation forest.
%
% Projects units into a low-dimensional UMAP space, then uses an isolation
% forest to remove outliers from the responsive-unit candidates. Selects
% counterexamples (non-responsive units) via farthest-point sampling in
% normalised feature space to form a balanced training set
% (class 1 = excitatory, class 2 = inhibitory).
%
% Counterexample selection uses farthest-point sampling (k-center greedy)
% starting from the responsive cluster centroid. This ensures spatial
% diversity across the excitatory population without relying on UMAP
% distances (which distort global structure). The approach is deterministic
% and avoids arbitrary distance-percentile thresholds.
%
% When TrainingCultureIdx is set, only units from those cultures are used for
% the unsupervised UMAP and label selection. When empty (default), all units
% are used.
%
% Normalisation pipeline:
%   1. Z-score within each NormalizationVar group (e.g. ChipID) across all units
%   2. Global z-score on the UMAP subset
%   3. Remove NaN columns, scale by max(abs)
%   4. Store normalisation parameters for carryover to classifyUnits/applyTo
%
% Feature extraction uses ctc.Parameters.Harmonization so DeePhys ACGs are
% recomputed at the configured bin size and lag (default: 0.5 ms / 100 ms).
%
% Requires ctc.ResponsiveUnitIdx to be set (run identifyResponsiveUnits first).
% Sets: ctc.TrainLabels, ctc.UMAP, ctc.NormalizationParams

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
rg      = ctc.RecordingGroup;

% ── Seed RNG for reproducibility ──────────────────────────────────────────────
rng(ctc.Parameters.RNGSeed, 'twister');

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr] = ctc.getOrExtract(ctc.UnitList);
[X_raw, ~]    = buildFeatureMatrix(ctc, wf, acg, sr);

% ── Chip-level normalisation (all units) ──────────────────────────────────────
% NaN handling: leave as NaN so they propagate to column removal (step 3).
% This avoids bias from replacing NaN with 0 before z-scoring.
[iG, G] = rg.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx = unique(iG);
for g = 1:length(g_idx)
    mask = (iG == g_idx(g));
    X_raw(mask, :) = normalize(X_raw(mask, :));
end

% ── Determine training subset ────────────────────────────────────────────────
% When TrainingCultureIdx is set, restrict unsupervised UMAP to those cultures.
if ~isempty(p_umap.TrainingCultureIdx)
    subset_mask = buildCultureUnitMask(ctc, p_umap.TrainingCultureIdx);
else
    subset_mask = true(1, numel(ctc.UnitList));
end

X_subset = X_raw(subset_mask, :);

% Global z-score on subset — store parameters for carryover
if length(G) > 1
    [X_subset, mu_global, sigma_global] = normalize(X_subset);
else
    mu_global    = zeros(1, size(X_subset, 2));
    sigma_global = ones(1, size(X_subset, 2));
end

% Remove NaN columns and store which ones
nan_cols = any(isnan(X_subset), 1);
X_subset(:, nan_cols) = [];
scale = max(abs(X_subset), [], 1);
scale(scale == 0) = 1;
X_subset = X_subset ./ scale;

% ── Store normalisation parameters for carryover ─────────────────────────────
ctc.NormalizationParams = struct( ...
    'mu_global',    mu_global, ...
    'sigma_global', sigma_global, ...
    'nan_cols',     nan_cols, ...
    'scale',        scale);

% ── UMAP embedding (unsupervised, on subset only) ───────────────────────────
template_file = fullfile(p_umap.TemplateDir, 'ctc_umap_template.mat');
[reduction, umap_model, ~, ~] = run_umap(X_subset, ...
    'n_components',  p_umap.NDims, ...
    'n_neighbors',   p_umap.NNeighbors, ...
    'min_dist',      p_umap.MinDist, ...
    'spread',        p_umap.Spread, ...
    'sgd_tasks',     20, ...
    'method', 'Java', 'verbose', 'none', ...
    'save_template_file', template_file);
ctc.UMAP = umap_model;
ctc.Reduction.Unsupervised = reduction;

% ── Map responsive units to subset-local indices ─────────────────────────────
subset_responsive = ctc.ResponsiveUnitIdx(subset_mask);

% Isolation forest on responsive-unit candidates to remove outliers
candidate_reduction = reduction(subset_responsive, :);
[~, tf_forest, ~] = iforest(candidate_reduction, ...
    'ContaminationFraction', p_outlr.ContaminationFraction, ...
    'NumObservationsPerLearner', p_outlr.NObsPerLearner);

clean_candidates = candidate_reduction(~tf_forest, :);
if isempty(clean_candidates)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'All %i responsive candidates flagged as outliers by iforest — falling back to all candidates.', ...
        size(candidate_reduction, 1));
    clean_candidates = candidate_reduction;
    tf_forest = false(size(tf_forest));
end

% ── Counterexample selection via farthest-point sampling in feature space ─────
% Work in normalised feature space (X_subset) where distances are meaningful,
% rather than UMAP space where global distances are distorted.
n_candidates     = sum(~tf_forest);
n_counterexamples = p_outlr.CounterexampleRatio * n_candidates;

% Identify non-responsive units in subset (pool for counterexamples)
responsive_in_subset = find(subset_responsive);
clean_responsive_local = responsive_in_subset(~tf_forest(:)');
non_responsive_local = find(~subset_responsive);

% Centroid of clean responsive units in feature space
centroid = mean(X_subset(clean_responsive_local, :), 1);

% Farthest-point sampling (k-center greedy):
% Start from responsive centroid, iteratively pick the non-responsive unit
% farthest from all already-selected points. This maximises spatial diversity
% and naturally biases away from the inhibitory cluster.
n_pool = length(non_responsive_local);
n_pick = min(n_counterexamples, n_pool);
X_pool = X_subset(non_responsive_local, :);

% Initialise minimum distances from the centroid
min_dists = pdist2(centroid, X_pool, 'euclidean')';
selected = false(n_pool, 1);
pick_order = zeros(n_pick, 1);

for k = 1:n_pick
    % Pick the farthest point
    candidates_dists = min_dists;
    candidates_dists(selected) = -Inf;
    [~, idx] = max(candidates_dists);
    selected(idx) = true;
    pick_order(k) = idx;

    % Update minimum distances to include the newly selected point
    new_dists = pdist2(X_pool(idx, :), X_pool, 'euclidean')';
    min_dists = min(min_dists, new_dists);
end

counterexample_idx_local = non_responsive_local(pick_order);

% ── Map local indices back to global UnitList indices ────────────────────────
subset_global = find(subset_mask);  % global indices of subset units

in_train_id = subset_global(clean_responsive_local);
ex_train_id = subset_global(counterexample_idx_local);

train_ids = [ex_train_id, in_train_id];
y_train   = [ones(1, length(ex_train_id)), 2*ones(1, length(in_train_id))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids = sorted_train_ids;
labels.sorted_y_train   = y_train(sort_idx);
labels.umap_train_idx   = false(1, length(ctc.UnitList));
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx    = ~labels.umap_train_idx;

% Store outlier and counterexample info for diagnostic plotting
outlier_local = responsive_in_subset(tf_forest(:)');
labels.outlier_global_idx        = subset_global(outlier_local);
labels.excitatory_candidate_global_idx = subset_global(counterexample_idx_local);

% Store responsive strength for the training inhibitory candidates (for diagnostics)
if ~isempty(ctc.ResponsiveStrength)
    labels.inhibitory_strength = ctc.ResponsiveStrength(in_train_id);
else
    labels.inhibitory_strength = ones(1, length(in_train_id));
end

ctc.TrainLabels = labels;
fprintf('Training set: %i excitatory, %i inhibitory candidates\n', ...
    sum(y_train==1), sum(y_train==2));
end

function mask = buildCultureUnitMask(ctc, culture_indices)
% Build a logical mask over ctc.UnitList selecting units from specified cultures.
    rg = ctc.RecordingGroup;
    mask = false(1, numel(ctc.UnitList));
    offset = 0;
    for c = 1:length(rg.Cultures)
        culture = rg.Cultures(c);
        n_units_c = numel(culture.Units);
        if ismember(c, culture_indices)
            mask(offset+1 : offset+n_units_c) = true;
        end
        offset = offset + n_units_c;
    end
end
