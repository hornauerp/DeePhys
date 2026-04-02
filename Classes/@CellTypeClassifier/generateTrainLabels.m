function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels via UMAP embedding and adaptive outlier detection.
%
% Projects all units into a low-dimensional UMAP space, then applies adaptive
% outlier detection on responsive-unit candidates:
%   1. Dip test on first PC of responsive units — determines uni/multimodality
%   2. Unimodal: robust Mahalanobis distance + chi-squared threshold
%   3. Multimodal: GMM with BIC-selected k + low-posterior outliers
%
% Counterexample selection uses farthest-point sampling (k-center greedy)
% in normalised feature space to maximise spatial diversity across the
% excitatory training set. This is deterministic and avoids arbitrary
% distance-percentile thresholds.
%
% ResponsiveClassLabel controls which class responsive units are assigned to:
%   2 (default) — responsive = inhibitory   (stimulation activates inhibitory neurons)
%   1           — responsive = excitatory   (direct excitatory stimulation paradigm)
%
% Normalisation note: the unsupervised UMAP stage fits z-score on the entire
% unit subset (no train/test split exists yet). Parameters are stored in
% ctc.NormalizationParams for carryover to classifyUnits and applyTo.
%
% Requires ctc.ResponsiveUnitIdx to be set (run identifyResponsiveUnits first).
% Sets: ctc.TrainLabels, ctc.UMAP, ctc.Reduction.Unsupervised, ctc.NormalizationParams

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
p_train = ctc.Parameters.TrainLabels;

rng(ctc.Parameters.RNGSeed, 'twister');

% -- Extract harmonized waveforms + ACGs from UnitDataArray ------------------
[wf, acg, sr] = ctc.getOrExtract(ctc.UnitDataArray);
[X_raw, ~]    = buildFeatureMatrix(ctc, wf, acg, sr);

% -- Step 1: per-group normalisation using FeatureStore table column ---------
X_raw(isnan(X_raw)) = 0;
norm_var = p_umap.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col = ctc.FeatureStore.UnitTable.(norm_var);
    [~, ~, iG] = unique(string(group_col), 'stable');
    for g = 1:max(iG)
        mask = iG == g;
        X_raw(mask, :) = normalize(X_raw(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end

% -- Determine training subset (by culture index or use all) -----------------
n_units = size(X_raw, 1);
if ~isempty(p_umap.TrainingCultureIdx)
    subset_mask = CellTypeClassifier.buildCultureSubsetMask( ...
        ctc.FeatureStore, p_umap.TrainingCultureIdx, ctc.Parameters.CultureKeys);
else
    subset_mask = true(1, n_units);
end

X_subset = X_raw(subset_mask, :);

% -- Step 2: Global z-score on subset — store parameters for carryover -------
[X_subset, mu_global, sigma_global] = normalize(X_subset);

% -- Step 3: Remove NaN columns and store which ones -------------------------
nan_cols = any(isnan(X_subset), 1);
X_subset(:, nan_cols) = [];
scale = max(abs(X_subset), [], 1);
scale(scale == 0) = 1;
X_subset = X_subset ./ scale;

% -- Step 4: Optional feature selection (Phase 3) ----------------------------
% Removes low-variance features and highly correlated feature pairs before UMAP.
% The selection mask is stored for consistent application in classifyUnits/applyTo.
if p_umap.FeatureSelection && size(X_subset, 2) > 1
    feat_var  = var(X_subset, 0, 1);
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
    X_umap = X_subset(:, feature_selection_mask);
    fprintf('Feature selection: removed %d low-var + %d redundant = %d total (kept %d)\n', ...
        sum(low_var), sum(high_corr & ~low_var), sum(remove_feat), size(X_umap, 2));
else
    feature_selection_mask = true(1, size(X_subset, 2));
    X_umap = X_subset;
end

% -- Store normalisation parameters for carryover ----------------------------
ctc.NormalizationParams = struct( ...
    'mu_global',             mu_global, ...
    'sigma_global',          sigma_global, ...
    'nan_cols',              nan_cols, ...
    'scale',                 scale, ...
    'feature_selection_mask', feature_selection_mask);

% -- AutoNDims: set from PCA variance (Phase 2) ------------------------------
if p_umap.AutoNDims
    [~, ~, latent] = pca(X_umap, 'Economy', true);
    cum_var = cumsum(latent) / sum(latent);
    n_dims  = find(cum_var >= p_umap.VarianceThreshold, 1, 'first');
    n_dims  = max(3, min(20, n_dims));
    fprintf('Auto NDims: %d (%.1f%% variance at threshold %.2f)\n', ...
        n_dims, cum_var(n_dims) * 100, p_umap.VarianceThreshold);
else
    n_dims = p_umap.NDims;
end

% -- AutoNNeighbors: scale with sqrt(N) (Phase 1) ----------------------------
N_umap = size(X_umap, 1);
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_umap)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_umap);
else
    n_neighbors = p_umap.NNeighbors;
end

% -- Unsupervised UMAP embedding ---------------------------------------------
assert(~isempty(p_umap.TemplateDir), ...
    'CellTypeClassifier:noTemplateDir', ...
    'Parameters.UMAP.TemplateDir must be set before calling generateTrainLabels. ' ...
    'Example: parameters.UMAP.TemplateDir = ''C:/path/to/output'';');
assert(isfolder(p_umap.TemplateDir), ...
    'CellTypeClassifier:templateDirMissing', ...
    'Parameters.UMAP.TemplateDir does not exist: %s', p_umap.TemplateDir);
template_file = fullfile(p_umap.TemplateDir, 'ctc_umap_template.mat');
[reduction, umap_model, ~, ~] = run_umap(X_umap, ...
    'n_components',  n_dims, ...
    'n_neighbors',   n_neighbors, ...
    'min_dist',      p_umap.MinDist, ...
    'spread',        p_umap.Spread, ...
    'sgd_tasks',     20, ...
    'method', 'Java', 'verbose', 'none', ...
    'save_template_file', template_file);
ctc.UMAP = umap_model;
ctc.Reduction.Unsupervised = reduction;

% -- Adaptive outlier detection on responsive candidates ---------------------
subset_responsive = ctc.ResponsiveUnitIdx(subset_mask);
responsive_local  = find(subset_responsive);
candidate_reduction = reduction(subset_responsive, :);

tf_outlier = detectResponsiveOutliers(candidate_reduction, p_outlr);

n_rejected = sum(tf_outlier);
n_total    = numel(tf_outlier);
if n_rejected > 0
    fprintf('Outlier detection: %d/%d responsive units flagged (%.0f%%)\n', ...
        n_rejected, n_total, n_rejected / n_total * 100);
end

clean_resp_local = responsive_local(~tf_outlier(:)');

if isempty(clean_resp_local)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'All responsive candidates flagged as outliers — using all candidates.');
    clean_resp_local = responsive_local;
    tf_outlier = false(size(tf_outlier));
end

% -- AutoCounterexampleRatio: match class balance to dataset (Phase 1) -------
non_responsive_local = find(~subset_responsive);
n_responsive = numel(clean_resp_local);

if p_outlr.AutoCounterexampleRatio
    n_pool = numel(non_responsive_local);
    ce_ratio = max(1, min(5, round(n_pool / max(n_responsive, 1))));
    fprintf('Auto CounterexampleRatio: %d (responsive=%d, pool=%d)\n', ...
        ce_ratio, n_responsive, n_pool);
else
    ce_ratio = p_outlr.CounterexampleRatio;
end
n_counterexamples = ce_ratio * n_responsive;

% -- AutoDistanceThreshold: robust pool filtering (Phase 1) ------------------
% Compute responsive centroid in the UMAP-input feature space (X_umap).
centroid = mean(X_umap(clean_resp_local, :), 1);
candidate_dists = pdist2(centroid, X_umap(clean_resp_local, :), 'euclidean')';

if p_outlr.AutoDistanceThreshold
    med_d = median(candidate_dists);
    mad_d = median(abs(candidate_dists - med_d));
    min_pool_dist = med_d + p_outlr.DistanceMADMultiplier * mad_d;

    X_pool_all  = X_umap(non_responsive_local, :);
    pool_dists  = pdist2(centroid, X_pool_all, 'euclidean')';
    pool_keep   = pool_dists > min_pool_dist;

    if any(pool_keep)
        X_pool_filt          = X_pool_all(pool_keep, :);
        non_responsive_filt  = non_responsive_local(pool_keep);
        fprintf('Auto pool filter: kept %d/%d units (dist > %.3f)\n', ...
            sum(pool_keep), numel(pool_keep), min_pool_dist);
    else
        X_pool_filt         = X_pool_all;
        non_responsive_filt = non_responsive_local;
    end
else
    X_pool_filt         = X_umap(non_responsive_local, :);
    non_responsive_filt = non_responsive_local;
end

% -- Farthest-point sampling: k-center greedy for diverse counterexamples ----
% Starting from the responsive centroid, iteratively pick the non-responsive
% unit farthest from all already-selected points. Guarantees spatial diversity
% across the excitatory population and biases away from the inhibitory cluster.
n_pool = size(X_pool_filt, 1);
n_pick = min(n_counterexamples, n_pool);

if n_pick > 0
    min_dists = pdist2(centroid, X_pool_filt, 'euclidean')';
    selected  = false(n_pool, 1);
    pick_order = zeros(n_pick, 1);

    for k = 1:n_pick
        cand_dists = min_dists;
        cand_dists(selected) = -Inf;
        [~, idx] = max(cand_dists);
        selected(idx) = true;
        pick_order(k) = idx;
        new_dists = pdist2(X_pool_filt(idx, :), X_pool_filt, 'euclidean')';
        min_dists = min(min_dists, new_dists);
    end

    counterexample_idx_local = non_responsive_filt(pick_order);
else
    counterexample_idx_local = zeros(1, 0);
end

% -- Map local indices back to global UnitDataArray indices ------------------
subset_global = find(subset_mask);

resp_label    = p_train.ResponsiveClassLabel;  % default 2 (inhibitory)
counter_label = 3 - resp_label;                % 1 if resp=2, 2 if resp=1

in_train_id = subset_global(clean_resp_local);
ex_train_id = subset_global(counterexample_idx_local);

train_ids = [ex_train_id, in_train_id];
y_train   = [counter_label * ones(1, numel(ex_train_id)), ...
             resp_label    * ones(1, numel(in_train_id))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids  = sorted_train_ids;
labels.sorted_y_train    = y_train(sort_idx);
labels.umap_train_idx    = false(1, n_units);
labels.umap_train_idx(sorted_train_ids) = true;
labels.umap_test_idx     = ~labels.umap_train_idx;

outlier_local = responsive_local(tf_outlier(:)');
labels.outlier_global_idx              = subset_global(outlier_local);
labels.excitatory_candidate_global_idx = subset_global(counterexample_idx_local);

% Store responsive strength for training inhibitory candidates (for diagnostics)
if ~isempty(ctc.ResponsiveStrength)
    labels.inhibitory_strength = ctc.ResponsiveStrength(in_train_id);
else
    labels.inhibitory_strength = ones(1, numel(in_train_id));
end

ctc.TrainLabels = labels;
fprintf('Training set: %i %s, %i %s candidates\n', ...
    sum(y_train == counter_label), labelName(counter_label), ...
    sum(y_train == resp_label),    labelName(resp_label));
end

%% -- Helper functions --------------------------------------------------------

function name = labelName(label)
% Map numeric label to readable name.
    if label == 1; name = 'excitatory';
    elseif label == 2; name = 'inhibitory';
    else; name = sprintf('class_%d', label);
    end
end

function tf_outlier = detectResponsiveOutliers(candidate_reduction, p)
% DETECTRESPONSIVEOUTLIERS  Adaptive outlier detection for responsive units.
%
% Two-stage approach:
%   1. Dip test on the first PC of candidate_reduction to test for multimodality
%   2a. Unimodal:   robust Mahalanobis distance + chi-squared threshold
%   2b. Multimodal: GMM with BIC-selected k + low-posterior outlier detection
%
% Uses DipTestAlpha for the dip test (conventional 0.05) and OutlierAlpha for
% the Mahalanobis/posterior threshold. When AutoOutlierAlpha is true, the
% outlier threshold is adapted from the compactness of the responsive cluster.

    n = size(candidate_reduction, 1);
    d = size(candidate_reduction, 2);

    if n < 4
        tf_outlier = false(n, 1);
        return
    end

    % -- Stage 1: dip test for multimodality (uses DipTestAlpha) -------------
    [~, scores] = pca(candidate_reduction, 'NumComponents', 1);
    pc1 = scores(:, 1);
    [~, p_dip] = diptest(pc1);

    dip_alpha = p.DipTestAlpha;  % default 0.05 (conventional)
    is_multimodal = p_dip < dip_alpha;

    % -- Adaptive outlier alpha (Phase 2) -------------------------------------
    if p.AutoOutlierAlpha && n >= 10
        try
            sil_vals = silhouette(candidate_reduction, ones(n, 1));
            mean_sil = mean(sil_vals);
            % Map silhouette [0.2, 0.5] -> alpha [0.05, 0.001] (log-linear)
            mean_sil_clamped = max(0.2, min(0.5, mean_sil));
            outlier_alpha = 10^interp1([0.2, 0.5], [log10(0.05), log10(0.001)], ...
                mean_sil_clamped, 'linear', 'extrap');
            fprintf('Auto OutlierAlpha: %.4f (mean silhouette=%.3f)\n', outlier_alpha, mean_sil);
        catch
            outlier_alpha = p.OutlierAlpha;
        end
    else
        outlier_alpha = p.OutlierAlpha;
    end

    if is_multimodal
        % -- Stage 2b: GMM + posterior-based outlier detection ----------------
        tf_outlier = gmmOutliers(candidate_reduction, p, outlier_alpha);
    else
        % -- Stage 2a: robust Mahalanobis + chi-squared -----------------------
        tf_outlier = mahalOutliers(candidate_reduction, d, outlier_alpha);
    end
end

function tf_outlier = mahalOutliers(data, d, outlier_alpha)
% Robust Mahalanobis distance outlier detection.
    n = size(data, 1);
    try
        [~, ~, mahal_d] = robustcov(data);
        thresh     = chi2inv(1 - outlier_alpha, d);
        tf_outlier = mahal_d > thresh;
    catch
        warning('CellTypeClassifier:robustcovFailed', ...
            'Robust covariance estimation failed — no outliers removed.');
        tf_outlier = false(n, 1);
    end
end

function tf_outlier = gmmOutliers(data, p, outlier_alpha)
% GMM-based outlier detection for multimodal responsive populations.
    n         = size(data, 1);
    max_k     = min(p.MaxResponsiveComponents, floor(n / 3));
    bic_vals  = inf(1, max_k);
    gmm_fits  = cell(1, max_k);

    for k = 1:max_k
        try
            gmm_fits{k} = fitgmdist(data, k, ...
                'RegularizationValue', 1e-5, ...
                'Options', statset('MaxIter', 500), ...
                'Replicates', 3);
            bic_vals(k) = gmm_fits{k}.BIC;
        catch
        end
    end

    [~, best_k] = min(bic_vals);

    if isinf(bic_vals(best_k))
        warning('CellTypeClassifier:gmmFailed', ...
            'All GMM fits failed — no outliers removed.');
        tf_outlier = false(n, 1);
        return
    end

    gmm = gmm_fits{best_k};
    posterior  = gmm.posterior(data);
    max_post   = max(posterior, [], 2);
    tf_outlier = max_post < outlier_alpha;

    fprintf('GMM selected k=%d components (BIC=%.1f); %d/%d outliers detected.\n', ...
        best_k, bic_vals(best_k), sum(tf_outlier), n);
end

function [dip, p_value] = diptest(x)
% DIPTEST  Hartigan's dip test for unimodality.
%
% Computes the dip statistic D_n and an approximate p-value for testing
% whether a 1D sample is drawn from a unimodal distribution.
% Small p_value -> reject unimodality (evidence of multimodality).
%
% Reference:
%   Hartigan, J.A. & Hartigan, P.M. (1985). The Dip Test of Unimodality.
%   Annals of Statistics, 13(1), 70-84.

    x = sort(x(:));
    n = numel(x);

    if n < 4
        dip = 0;
        p_value = 1;
        return
    end

    ecdf_y = (1:n)' / n;
    gcm = ecdfGCM(ecdf_y);
    lcm = ecdfLCM(ecdf_y);

    dip = max(lcm - gcm) / 2;

    dip_scaled = sqrt(n) * dip;

    alpha_table = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00];
    crit_table  = [1.20, 1.06, 0.89, 0.78, 0.66, 0.50, 0.35];

    if dip_scaled >= crit_table(1)
        p_value = alpha_table(1) / 2;
    elseif dip_scaled <= crit_table(end)
        p_value = 1.0;
    else
        p_value = interp1(crit_table, alpha_table, dip_scaled, 'linear', 'extrap');
        p_value = max(0, min(1, p_value));
    end
end

function gcm = ecdfGCM(y)
% Greatest convex minorant of a vector.
    n   = numel(y);
    gcm = y;
    i   = 1;
    while i < n
        j = i + 1;
        while j <= n
            slope_ij = (y(j) - gcm(i)) / (j - i);
            if j < n
                slope_next = (y(j+1) - gcm(i)) / (j + 1 - i);
                if slope_next < slope_ij
                    j = j + 1;
                    continue
                end
            end
            break
        end
        for k = i+1:j
            gcm(k) = gcm(i) + (y(j) - gcm(i)) * (k - i) / (j - i);
        end
        i = j;
    end
end

function lcm = ecdfLCM(y)
% Least concave majorant of a vector.
    n   = numel(y);
    lcm = y;
    i   = 1;
    while i < n
        j = i + 1;
        while j <= n
            slope_ij = (y(j) - lcm(i)) / (j - i);
            if j < n
                slope_next = (y(j+1) - lcm(i)) / (j + 1 - i);
                if slope_next > slope_ij
                    j = j + 1;
                    continue
                end
            end
            break
        end
        for k = i+1:j
            lcm(k) = lcm(i) + (y(j) - lcm(i)) * (k - i) / (j - i);
        end
        i = j;
    end
end
