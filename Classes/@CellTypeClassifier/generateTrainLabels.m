function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels using UMAP embedding and isolation forest.
%
% Projects units into a low-dimensional UMAP space, then uses an isolation
% forest to remove outliers from the responsive-unit candidates. Selects
% counterexamples (non-responsive units far from the responsive cluster centroid)
% to form a balanced training set (class 1 = excitatory, class 2 = inhibitory).
%
% When TrainingCultureIdx is set, only units from those cultures are used for
% the unsupervised UMAP and label selection (matching the main-branch behaviour
% where only specific recording dates feed label generation). When empty (default),
% all units are used.
%
% Normalisation matches the main-branch pipeline:
%   1. Z-score within each NormalizationVar group (e.g. ChipID) across all units
%   2. Global z-score on the UMAP subset
%   3. Remove NaN columns, scale by max(abs)
%
% Feature extraction uses ctc.Parameters.Harmonization so DeePhys ACGs are
% recomputed at the configured bin size and lag (default: 0.5 ms / 100 ms).
%
% Requires ctc.ResponsiveUnitIdx to be set (run identifyResponsiveUnits first).
% Sets: ctc.TrainLabels, ctc.UMAP

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
[wf, acg, sr] = extractUnitWaveformsAndACGs(ctc, ctc.UnitList);
[X_raw, ~]    = buildFeatureMatrix(ctc, wf, acg, sr);

% ── Chip-level normalisation (all units, matching main branch) ───────────────
X_raw(isnan(X_raw)) = 0;
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

% Global z-score on subset (matching main branch: train-only z-score + scale)
if length(G) > 1
    [X_subset, ~, ~] = normalize(X_subset);
end
nan_cols = any(isnan(X_subset), 1);
X_subset(:, nan_cols) = [];
scale = max(abs(X_subset), [], 1);
scale(scale == 0) = 1;
X_subset = X_subset ./ scale;

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
centroid = mean(clean_candidates, 1);

% Find counterexamples: units far from interneuron centroid (in subset space)
dists = pdist2(centroid, reduction);
candidate_dists = pdist2(centroid, clean_candidates);
dist_threshold = prctile(candidate_dists, p_outlr.DistancePercentile);

far_idx_local = find(dists > dist_threshold);
n_candidates = sum(~tf_forest);
n_counterexamples = p_outlr.CounterexampleRatio * n_candidates;
counterexample_idx_local = randsample(far_idx_local, min(n_counterexamples, length(far_idx_local)), false);

% ── Map local indices back to global UnitList indices ────────────────────────
subset_global = find(subset_mask);  % global indices of subset units
responsive_in_subset = find(subset_responsive);  % local indices of responsive units in subset
clean_responsive_local = responsive_in_subset(~tf_forest(:)');

in_train_id = subset_global(clean_responsive_local);
ex_train_id_global = subset_global(counterexample_idx_local);
ex_train_id = ex_train_id_global(~ismember(ex_train_id_global, in_train_id));

train_ids = [ex_train_id, in_train_id];
y_train   = [ones(1, length(ex_train_id)), 2*ones(1, length(in_train_id))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids = sorted_train_ids;
labels.sorted_y_train   = y_train(sort_idx);
labels.umap_train_idx   = false(1, length(ctc.UnitList));
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx    = ~labels.umap_train_idx;

ctc.TrainLabels = labels;
fprintf('Training set: %i excitatory, %i inhibitory candidates\n', ...
    sum(y_train==1), sum(y_train==2));
end

function mask = buildCultureUnitMask(ctc, culture_indices)
% Build a logical mask over ctc.UnitList selecting units from specified cultures.
% Culture indices refer to ctc.RecordingGroup.Cultures.
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
