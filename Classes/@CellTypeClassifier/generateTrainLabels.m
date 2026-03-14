function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels using UMAP embedding and isolation forest.
%
% Projects all units into a low-dimensional UMAP space, then uses an isolation
% forest to remove outliers from the responsive-unit candidates. Selects
% counterexamples (non-responsive units far from the responsive cluster centroid)
% to form a balanced training set (class 1 = excitatory, class 2 = inhibitory).
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

% ── Extract features via harmonized path ─────────────────────────────────────
[wf, acg, sr] = extractUnitWaveformsAndACGs(ctc, ctc.UnitList);
[X_raw, ~]    = buildFeatureMatrix(ctc, wf, acg, sr);

% ── Normalise: joint z-score all units, scale by column max ──────────────────
[X_norm, ~, ~] = normalize(X_raw, 1, 'zscore');
nan_cols = any(isnan(X_norm), 1);
X_norm(:, nan_cols) = [];
scale = max(abs(X_norm), [], 1);
scale(scale == 0) = 1;
X_all = X_norm ./ scale;

% ── UMAP embedding of all units (unsupervised) ───────────────────────────────
template_file = fullfile(p_umap.TemplateDir, 'ctc_umap_template.mat');
[reduction, umap_model, ~, ~] = run_umap(X_all, ...
    'n_components',  p_umap.NDims, ...
    'n_neighbors',   p_umap.NNeighbors, ...
    'min_dist',      p_umap.MinDist, ...
    'spread',        p_umap.Spread, ...
    'sgd_tasks',     20, ...
    'method', 'Java', 'verbose', 'none', ...
    'save_template_file', template_file);
ctc.UMAP = umap_model;
ctc.Reduction.Unsupervised = reduction;

% Isolation forest on responsive-unit candidates to remove outliers
candidate_reduction = reduction(ctc.ResponsiveUnitIdx, :);
[~, tf_forest, ~] = iforest(candidate_reduction, ...
    'ContaminationFraction', p_outlr.ContaminationFraction, ...
    'NumObservationsPerLearner', p_outlr.NObsPerLearner);

clean_candidates = candidate_reduction(~tf_forest, :);
if isempty(clean_candidates)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'All %i responsive candidates flagged as outliers by iforest — falling back to all candidates.', ...
        size(candidate_reduction, 1));
    clean_candidates = candidate_reduction;
end
centroid = mean(clean_candidates, 1);

% Find counterexamples: units far from interneuron centroid
dists = pdist2(centroid, reduction);
candidate_dists = pdist2(centroid, clean_candidates);
dist_threshold = prctile(candidate_dists, p_outlr.DistancePercentile);

far_idx = find(dists > dist_threshold);
n_candidates = sum(~tf_forest);
counterexample_idx = randsample(far_idx, min(n_candidates, length(far_idx)), false);

responsive_indices = find(ctc.ResponsiveUnitIdx);
in_train_id = responsive_indices(~tf_forest(:)');
ex_train_id = counterexample_idx(~ismember(counterexample_idx, in_train_id));

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
