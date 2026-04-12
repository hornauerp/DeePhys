function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels via outlier detection and feature-space operations.
%
% Pipeline:
%   1. Build normalized feature cache (shared with classifyUnits via buildNormalizedFeatures)
%   2. Global z-score on training subset; store NormalizationParams for classifyUnits
%   3. Optional feature selection
%   4. Unsupervised UMAP embedding (always computed; used for diagnostics and "umap" domain)
%   5. Outlier detection on responsive candidates (domain: "pca" or "umap")
%   6. Geometric consistency filter (feature space)
%   7. Counterexample selection: farthest-from-inhibitory-centroid (feature space)
%   8. Counterexample outlier detection (domain: "pca" or "umap")
%   9. Assemble training labels
%
% Outlier detection domain (OutlierDetection.Domain):
%   "pca"  — top NPCAComponents of candidate feature matrix (default).
%            Decouples labels from UMAP hyperparameters; PCA is deterministic.
%   "umap" — unsupervised UMAP embedding (legacy behavior).
%
% Two label-source strategies:
%
%   Drug response (inferred counterexamples):
%     Responsive units form one-class ground truth (e.g. inhibitory).
%     Counterexamples selected farthest from inhibitory centroid in normalized
%     feature space (correlation distance), ensuring unambiguously excitatory units.
%
%   Pure culture / metadata (explicit counterexamples):
%     Both classes have ground truth from ctc.CounterexampleUnitIdx.
%     Both classes are cleaned by outlier detection independently.
%
% Requires ctc.ResponsiveUnitIdx to be set (run identifyResponsiveUnits first).
% Sets: ctc.TrainLabels, ctc.NormalizationParams, ctc.Reduction.Unsupervised

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
p_train = ctc.Parameters.TrainLabels;

rng(ctc.Parameters.RNGSeed, 'twister');

% -- Build normalized feature cache (shared with classifyUnits) ---------------
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;

n_unique      = nf.n_unique;
n_units_full  = nf.n_units_full;
all_to_unique = nf.all_to_unique;
unique_to_rep = nf.unique_to_rep;
subset_mask   = nf.subset_mask;

% Map responsive flag to unique units by UnitID matching (order-independent).
% ResponsiveUnitIdx is in FeatureStore order; unique_to_rep is in UnitDataArray
% order. These arrays may not be aligned, so positional indexing would return
% wrong units. UnitID matching is immune to order differences.
unit_ids_table = string(ctc.FeatureStore.UnitTable.UnitID);
resp_uids      = unique(unit_ids_table(ctc.ResponsiveUnitIdx));
unique_ud_ids  = string({nf.unique_ud.UnitID});   % row vector
resp_unique    = ismember(unique_ud_ids, resp_uids);  % row vector (1 x n_unique)

% -- Global z-score on training subset ----------------------------------------
X_subset = nf.X_pergroup(subset_mask, :);
[X_subset, mu_global, sigma_global] = normalize(X_subset);

% -- Remove NaN columns + scale -----------------------------------------------
nan_cols = any(isnan(X_subset), 1);
X_subset(:, nan_cols) = [];
scale = max(abs(X_subset), [], 1);
scale(scale == 0) = 1;
X_subset = X_subset ./ scale;

% Track feature names through removal
feat_names_trimmed = nf.feat_names(~nan_cols);

% -- Optional feature selection -----------------------------------------------
if p_umap.FeatureSelection && size(X_subset, 2) > 1
    feat_var   = var(X_subset, 0, 1);
    var_thresh = prctile(feat_var, p_umap.MinVariancePercentile);
    low_var    = feat_var <= var_thresh;

    R = abs(corrcoef(X_subset));
    R(logical(eye(size(R)))) = 0;
    high_corr = false(1, size(X_subset, 2));
    for j = 2:size(R, 2)
        if any(R(1:j-1, j) > p_umap.MaxCorrelation & ~high_corr(1:j-1)')
            high_corr(j) = true;
        end
    end

    remove_feat = low_var | high_corr;
    feature_selection_mask = ~remove_feat;
    X_feat = X_subset(:, feature_selection_mask);
    feat_names_trimmed = feat_names_trimmed(feature_selection_mask);
    fprintf('Feature selection: removed %d low-var + %d redundant = %d total (kept %d)\n', ...
        sum(low_var), sum(high_corr & ~low_var), sum(remove_feat), size(X_feat, 2));
else
    feature_selection_mask = true(1, size(X_subset, 2));
    X_feat = X_subset;
end

% -- Store NormalizationParams -------------------------------------------------
ctc.NormalizationParams = struct( ...
    'mu_global',              mu_global, ...
    'sigma_global',           sigma_global, ...
    'nan_cols',               nan_cols, ...
    'scale',                  scale, ...
    'feature_selection_mask', feature_selection_mask, ...
    'feat_names_trimmed',     feat_names_trimmed);

% -- Normalize ALL unique units (for full-dataset UMAP) -----------------------
X_all = normalize(nf.X_pergroup, 'center', mu_global, 'scale', sigma_global);
X_all(:, nan_cols) = [];
X_all = X_all ./ scale;
if any(~feature_selection_mask)
    X_all = X_all(:, feature_selection_mask);
end

% -- Unsupervised UMAP embedding on ALL unique units --------------------------
% Run on the full dataset so the UMAP graph covers all units (not just the
% training subset). This enables graph-based label propagation in classifyUnits
% and provides a global embedding for diagnostics.
assert(~isempty(p_umap.TemplateDir), ...
    ['CellTypeClassifier:noTemplateDir  ', ...
    'Parameters.UMAP.TemplateDir must be set before calling generateTrainLabels. ' ...
    'Example: parameters.UMAP.TemplateDir = ''C:/path/to/output'';']);
assert(isfolder(p_umap.TemplateDir), ...
    'CellTypeClassifier:templateDirMissing  %s does not exist.', p_umap.TemplateDir);

N_all = size(X_all, 1);
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_all)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_all);
else
    n_neighbors = p_umap.NNeighbors;
end

[reduction, umap_model, ~, ~] = run_umap(X_all, ...
    'n_components',  p_umap.NDims, ...
    'n_neighbors',   n_neighbors, ...
    'min_dist',      p_umap.MinDist, ...
    'spread',        p_umap.Spread, ...
    'method', 'Java', 'verbose', 'none');
ctc.UMAP = umap_model;
if isempty(ctc.Reduction); ctc.Reduction = struct(); end
ctc.Reduction.Unsupervised = reduction;
fprintf('Unsupervised UMAP: %d units, %d features, %d dimensions\n', ...
    N_all, size(X_all, 2), p_umap.NDims);

% -- Isolation Forest on inhibitory candidates in PCA space -------------------
% Decouples training labels from UMAP hyperparameters (prerequisite for stable
% Bayesian optimization). Near-zero variance columns removed before PCA.
subset_responsive = resp_unique(subset_mask);
responsive_local  = find(subset_responsive);

if isempty(responsive_local)
    n_resp_total = sum(ctc.ResponsiveUnitIdx);
    error('CellTypeClassifier:noResponsiveUnits', ...
        ['No responsive units found in the training culture subset.\n' ...
         '  Total responsive units across all cultures: %d\n' ...
         '  Training subset size: %d unique units\n' ...
         'Possible causes:\n' ...
         '  1. Parameters.UMAP.TrainingCultureIdx excludes cultures with responsive units.\n' ...
         '  2. Parameters.UMAP.GroupingValues does not match the baseline dose in the data.\n' ...
         '  3. identifyResponsiveUnits() was not run, or produced 0 responsive units.\n' ...
         'Diagnostic: check sum(ctc.ResponsiveUnitIdx) and ctc.Parameters.UMAP.TrainingCultureIdx.'], ...
        n_resp_total, sum(subset_mask));
end

candidate_features     = X_feat(responsive_local, :);
col_var                = var(candidate_features, 0, 1);
candidate_features_pca = candidate_features(:, col_var > 1e-10);
n_pca = min([p_outlr.NPCAComponents, size(candidate_features_pca, 1) - 1, size(candidate_features_pca, 2)]);
if n_pca >= 1
    [~, inh_pca_scores] = pca(candidate_features_pca, 'NumComponents', n_pca, 'Algorithm', 'eig');
else
    inh_pca_scores = candidate_features_pca;
end

[~, tf_resp_outlier] = iforest(inh_pca_scores, ...
    'ContaminationFraction', p_outlr.ContaminationFraction);

n_rejected = sum(tf_resp_outlier);
n_total    = numel(tf_resp_outlier);
fprintf('Isolation forest (inh, PCA space): %d/%d flagged (%.0f%%)\n', ...
    n_rejected, n_total, n_rejected / n_total * 100);

clean_resp_local = responsive_local(~tf_resp_outlier(:)');

if isempty(clean_resp_local)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'Isolation forest flagged all responsive candidates — using all candidates (reduce ContaminationFraction).');
    clean_resp_local = responsive_local;
    tf_resp_outlier  = false(size(tf_resp_outlier));
end

% -- Geometric consistency filter on responsive candidates --------------------
% Remove candidates whose feature profile is closer to the excitatory centroid
% than the inhibitory centroid (correlation distance). These are likely false
% positives whose FR changed due to network effects, not direct DREADD activation.
inh_centroid = mean(X_feat(clean_resp_local, :), 1);
exc_centroid = mean(X_feat(~subset_responsive, :), 1);
d_to_inh = pdist2(inh_centroid, X_feat(clean_resp_local, :), 'correlation');
d_to_exc = pdist2(exc_centroid, X_feat(clean_resp_local, :), 'correlation');
geom_consistent = d_to_inh(:) < d_to_exc(:);
n_geom_removed = sum(~geom_consistent);
if n_geom_removed > 0
    fprintf('Geometric consistency: removed %d/%d responsive candidates\n', ...
        n_geom_removed, numel(clean_resp_local));
end
resp_geom_mask = ~geom_consistent(:)';   % true = removed by geometric filter
clean_resp_local = clean_resp_local(geom_consistent);

if isempty(clean_resp_local)
    warning('CellTypeClassifier:generateTrainLabels', ...
        ['All responsive candidates removed by geometric filter — skipping filter.\n' ...
         'This often means the responsive units are geometrically closer to the\n' ...
         'excitatory population centroid, suggesting label direction may be inverted.\n' ...
         'Check Parameters.TrainLabels.ResponsiveClassLabel and Parameters.Bootstrap.Direction.']);
    clean_resp_local = responsive_local;   % ultimate fallback: use all responsive units
    resp_geom_mask   = false(1, numel(clean_resp_local));
    inh_centroid     = mean(X_feat(clean_resp_local, :), 1);
end

% -- Counterexample selection --------------------------------------------------
n_responsive = numel(clean_resp_local);

has_explicit_ce = ~isempty(ctc.CounterexampleUnitIdx) && any(ctc.CounterexampleUnitIdx);

% Initialize diagnostic storage for CE
ce_distances_from_inh = [];
ce_outlier_mask_diag  = logical([]);

if has_explicit_ce
    % Metadata method: explicit ground-truth counterexamples.
    % UnitID-based matching (order-independent, like resp_unique above)
    ce_uid_table = unique(unit_ids_table(ctc.CounterexampleUnitIdx));
    ce_unique    = ismember(unique_ud_ids, ce_uid_table);
    ce_subset    = ce_unique(subset_mask);
    ce_local_all = find(ce_subset);

    % Isolation forest on counterexample class (PCA space).
    ce_feat     = X_feat(ce_local_all, :);
    col_var_ce  = var(ce_feat, 0, 1);
    ce_feat_pca = ce_feat(:, col_var_ce > 1e-10);
    n_pca_ce    = min([p_outlr.NPCAComponents, size(ce_feat_pca,1)-1, size(ce_feat_pca,2)]);
    if n_pca_ce >= 1
        [~, ce_pca_scores] = pca(ce_feat_pca, 'NumComponents', n_pca_ce, 'Algorithm', 'eig');
    else
        ce_pca_scores = ce_feat_pca;
    end
    [~, tf_ce_outlier] = iforest(ce_pca_scores, 'ContaminationFraction', p_outlr.ContaminationFraction);
    n_ce_rejected = sum(tf_ce_outlier);
    if n_ce_rejected > 0
        fprintf('Isolation forest (CE, PCA space): %d/%d flagged (%.0f%%)\n', ...
            n_ce_rejected, numel(tf_ce_outlier), n_ce_rejected / numel(tf_ce_outlier) * 100);
    end
    counterexample_idx_local = ce_local_all(~tf_ce_outlier(:)');
    ce_outlier_mask_diag = tf_ce_outlier(:)';

    % Distances for diagnostic
    ce_distances_from_inh = pdist2(inh_centroid, X_feat(ce_local_all, :), 'correlation');

    if isempty(counterexample_idx_local)
        warning('CellTypeClassifier:generateTrainLabels', ...
            'All counterexample candidates flagged as outliers — using all candidates.');
        counterexample_idx_local = ce_local_all;
    end
    fprintf('Metadata counterexamples: %d labeled, %d after outlier cleaning\n', ...
        numel(ce_local_all), numel(counterexample_idx_local));
else
    % Inferred method: farthest-from-inhibitory-centroid selection,
    % then outlier detection with top-up from reserve.
    non_responsive_local = find(~subset_responsive);
    ce_ratio             = p_outlr.CounterexampleRatio;
    n_target             = ce_ratio * n_responsive;
    n_pool               = numel(non_responsive_local);
    n_pick               = min(n_target, n_pool);

    if n_pick > 0
        % Distance-percentile counterexample selection.
        % Select only non-responsive units farther from the inhibitory centroid
        % than the CounterexampleDistancePercentile of inhibitory units' own
        % self-distances. This guarantees the excitatory pool is clearly separated
        % from the inhibitory cluster without using extreme outliers.
        inh_self_dists = pdist2(inh_centroid, X_feat(clean_resp_local, :), 'correlation');
        dist_threshold = prctile(inh_self_dists, p_outlr.CounterexampleDistancePercentile);

        pool_dists = pdist2(inh_centroid, X_feat(non_responsive_local, :), 'correlation');
        ce_distances_from_inh = pool_dists;   % keep for diagnostics
        far_mask  = pool_dists(:) >= dist_threshold;
        far_pool  = non_responsive_local(far_mask);

        fprintf('Distance-percentile pool: %d/%d units beyond %.0fth pctile (thresh=%.4f)\n', ...
            numel(far_pool), n_pool, p_outlr.CounterexampleDistancePercentile, dist_threshold);

        if numel(far_pool) < n_pick
            fprintf('Pool smaller than n_pick (%d < %d) — falling back to full non-responsive pool\n', ...
                numel(far_pool), n_pick);
            far_pool = non_responsive_local;
        end

        perm_order    = randperm(numel(far_pool));
        n_pick_actual = min(n_pick, numel(far_pool));
        sampled_local = far_pool(perm_order(1:n_pick_actual));

        % Isolation forest on selected counterexamples (PCA space).
        ce_feat     = X_feat(sampled_local, :);
        col_var_ce  = var(ce_feat, 0, 1);
        ce_feat_pca = ce_feat(:, col_var_ce > 1e-10);
        n_pca_ce    = min([p_outlr.NPCAComponents, size(ce_feat_pca,1)-1, size(ce_feat_pca,2)]);
        if n_pca_ce >= 1
            [~, ce_pca_scores] = pca(ce_feat_pca, 'NumComponents', n_pca_ce, 'Algorithm', 'eig');
        else
            ce_pca_scores = ce_feat_pca;
        end
        [~, tf_ce_outlier]   = iforest(ce_pca_scores, 'ContaminationFraction', p_outlr.ContaminationFraction);
        n_ce_rejected        = sum(tf_ce_outlier);
        ce_outlier_mask_diag = tf_ce_outlier(:)';

        if n_ce_rejected > 0
            fprintf('Isolation forest (CE, PCA space): %d/%d flagged\n', ...
                n_ce_rejected, n_pick_actual);
            counterexample_idx_local = sampled_local(~tf_ce_outlier(:)');
        else
            counterexample_idx_local = sampled_local;
        end

        if isempty(counterexample_idx_local)
            warning('CellTypeClassifier:generateTrainLabels', ...
                'All CE candidates flagged as outliers — using original sample.');
            counterexample_idx_local = sampled_local;
        end
    else
        counterexample_idx_local = zeros(1, 0);
        ce_outlier_mask_diag     = logical([]);
        ce_distances_from_inh    = [];
    end
end

% -- Map local indices back to global unique-unit indices ---------------------
subset_global = find(subset_mask);

resp_label    = p_train.ResponsiveClassLabel;
counter_label = 3 - resp_label;

in_train_id = subset_global(clean_resp_local);
ex_train_id = subset_global(counterexample_idx_local);

train_ids = [ex_train_id, in_train_id];
y_train   = [counter_label * ones(1, numel(ex_train_id)), ...
             resp_label    * ones(1, numel(in_train_id))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids  = sorted_train_ids;
labels.sorted_y_train    = y_train(sort_idx);
labels.umap_train_idx    = false(1, n_unique);
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx     = ~labels.umap_train_idx;

outlier_local = responsive_local(tf_resp_outlier(:)');
labels.outlier_global_idx              = subset_global(outlier_local);
labels.excitatory_candidate_global_idx = subset_global(counterexample_idx_local);

% Store mapping so classifyUnits can broadcast results back to all rows
labels.n_unique        = n_unique;
labels.all_to_unique   = all_to_unique;
labels.n_units_full    = n_units_full;
labels.has_explicit_ce = has_explicit_ce;

% -- Diagnostic storage -------------------------------------------------------
% resp_outlier_mask: (1 x numel(responsive_local)) — true = removed by outlier detection
labels.resp_outlier_mask = tf_resp_outlier(:)';

% resp_geom_mask: (1 x numel(clean_resp_after_outlier)) — true = removed by geometric filter
labels.resp_geom_mask = resp_geom_mask;

% ce_outlier_mask: (1 x n_ce_sampled) — true = removed by outlier detection
labels.ce_outlier_mask = ce_outlier_mask_diag;

% ce_distances_from_inh_centroid: correlation distances of all CE pool from inh centroid
labels.ce_distances_from_inh_centroid = ce_distances_from_inh;

% inh_centroid: (1 x F) centroid of clean inhibitory training set in feature space
labels.inh_centroid = inh_centroid;

% filter_counts for diagnostic funnel
labels.filter_counts = struct( ...
    'n_responsive_in',    numel(responsive_local), ...
    'n_after_outlier',    numel(responsive_local) - sum(tf_resp_outlier), ...
    'n_after_geometric',  numel(clean_resp_local), ...
    'n_ce_in',            n_responsive, ...
    'n_ce_after_outlier', numel(counterexample_idx_local));

% Responsive strength for training inhibitory candidates (UnitID-based lookup)
if ~isempty(ctc.ResponsiveStrength)
    in_train_ud_ids = string({nf.unique_ud(in_train_id).UnitID});
    [~, fs_locs]    = ismember(in_train_ud_ids, unit_ids_table);
    labels.inhibitory_strength = ctc.ResponsiveStrength(fs_locs);
else
    labels.inhibitory_strength = ones(1, numel(in_train_id));
end

ctc.TrainLabels = labels;
fprintf('Training set: %i %s, %i %s candidates\n', ...
    sum(y_train == counter_label), labelName(counter_label), ...
    sum(y_train == resp_label),    labelName(resp_label));

% -- Optional diagnostic ------------------------------------------------------
if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticTrainLabels();
end
end

%% -- Helper functions --------------------------------------------------------

function name = labelName(label)
    if label == 1; name = 'excitatory';
    elseif label == 2; name = 'inhibitory';
    else; name = sprintf('class_%d', label);
    end
end

