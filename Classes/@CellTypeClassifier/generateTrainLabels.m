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

resp_unique = ctc.ResponsiveUnitIdx(unique_to_rep);

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

% -- Unsupervised UMAP embedding (always computed) ----------------------------
assert(~isempty(p_umap.TemplateDir), ...
    ['CellTypeClassifier:noTemplateDir  ', ...
    'Parameters.UMAP.TemplateDir must be set before calling generateTrainLabels. ' ...
    'Example: parameters.UMAP.TemplateDir = ''C:/path/to/output'';']);
assert(isfolder(p_umap.TemplateDir), ...
    'CellTypeClassifier:templateDirMissing  %s does not exist.', p_umap.TemplateDir);

N_feat = size(X_feat, 1);
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_feat)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_feat);
else
    n_neighbors = p_umap.NNeighbors;
end

template_file = fullfile(p_umap.TemplateDir, 'ctc_umap_template.mat');
[reduction, umap_model, ~, ~] = run_umap(X_feat, ...
    'n_components',  p_umap.NDims, ...
    'n_neighbors',   n_neighbors, ...
    'min_dist',      p_umap.MinDist, ...
    'spread',        p_umap.Spread, ...
    'sgd_tasks',     20, ...
    'method', 'Java', 'verbose', 'none', ...
    'save_template_file', template_file);
ctc.UMAP = umap_model;
if isempty(ctc.Reduction); ctc.Reduction = struct(); end
ctc.Reduction.Unsupervised = reduction;

% -- Determine outlier detection inputs ---------------------------------------
subset_responsive = resp_unique(subset_mask);
responsive_local  = find(subset_responsive);
use_pca = string(p_outlr.Domain) == "pca";

if use_pca
    % PCA on candidate feature matrix; near-zero variance columns removed first
    candidate_features = X_feat(responsive_local, :);
    col_var = var(candidate_features, 0, 1);
    candidate_features_pca = candidate_features(:, col_var > 1e-10);
    n_pca = min([p_outlr.NPCAComponents, size(candidate_features_pca, 1) - 1, size(candidate_features_pca, 2)]);
    if n_pca >= 1
        [~, inh_detection_data] = pca(candidate_features_pca, 'NumComponents', n_pca, 'Algorithm', 'eig');
    else
        inh_detection_data = candidate_features_pca;
    end
else
    inh_detection_data = reduction(subset_responsive, :);
end

% -- Outlier detection on responsive candidates -------------------------------
tf_resp_outlier = detectOutliers(inh_detection_data, p_outlr);

n_rejected = sum(tf_resp_outlier);
n_total    = numel(tf_resp_outlier);
if n_rejected > 0
    fprintf('Responsive outlier detection (%s space): %d/%d flagged (%.0f%%)\n', ...
        p_outlr.Domain, n_rejected, n_total, n_rejected / n_total * 100);
end

clean_resp_local = responsive_local(~tf_resp_outlier(:)');

if isempty(clean_resp_local)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'All responsive candidates flagged as outliers — using all candidates.');
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
        'All responsive candidates removed by geometric filter — skipping filter.');
    clean_resp_local = responsive_local(~tf_resp_outlier(:)');
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
    ce_unique    = ctc.CounterexampleUnitIdx(unique_to_rep);
    ce_subset    = ce_unique(subset_mask);
    ce_local_all = find(ce_subset);

    % Outlier detection on counterexample class.
    if use_pca
        ce_feat = X_feat(ce_local_all, :);
        col_var_ce = var(ce_feat, 0, 1);
        ce_feat_pca = ce_feat(:, col_var_ce > 1e-10);
        n_pca_ce = min([p_outlr.NPCAComponents, size(ce_feat_pca, 1) - 1, size(ce_feat_pca, 2)]);
        if n_pca_ce >= 1
            [~, ce_detection_data] = pca(ce_feat_pca, 'NumComponents', n_pca_ce, 'Algorithm', 'eig');
        else
            ce_detection_data = ce_feat_pca;
        end
    else
        ce_detection_data = reduction(ce_subset, :);
    end

    tf_ce_outlier = detectOutliers(ce_detection_data, p_outlr);
    n_ce_rejected = sum(tf_ce_outlier);
    if n_ce_rejected > 0
        fprintf('Counterexample outlier detection (%s space): %d/%d flagged (%.0f%%)\n', ...
            p_outlr.Domain, n_ce_rejected, numel(tf_ce_outlier), n_ce_rejected / numel(tf_ce_outlier) * 100);
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
        % Select farthest from inhibitory centroid (most unambiguously excitatory)
        pool_dists = pdist2(inh_centroid, X_feat(non_responsive_local, :), 'correlation');
        [~, dist_order] = sort(pool_dists, 'descend');
        sampled_pos = dist_order(1:n_pick);
        reserve_pos = dist_order(n_pick+1:end);

        % Store all distances for diagnostic
        ce_distances_from_inh = pool_dists;

        sampled_local = non_responsive_local(sampled_pos);

        % Outlier detection on initial sample
        if use_pca
            ce_feat = X_feat(sampled_local, :);
            col_var_ce = var(ce_feat, 0, 1);
            ce_feat_pca = ce_feat(:, col_var_ce > 1e-10);
            n_pca_ce = min([p_outlr.NPCAComponents, size(ce_feat_pca, 1) - 1, size(ce_feat_pca, 2)]);
            if n_pca_ce >= 1
                [~, ce_detection_data] = pca(ce_feat_pca, 'NumComponents', n_pca_ce, 'Algorithm', 'eig');
            else
                ce_detection_data = ce_feat_pca;
            end
        else
            ce_detection_data = reduction(sampled_local, :);
        end

        tf_ce_outlier = detectOutliers(ce_detection_data, p_outlr);
        n_ce_rejected = sum(tf_ce_outlier);
        ce_outlier_mask_diag = tf_ce_outlier(:)';

        if n_ce_rejected > 0
            fprintf('Counterexample outlier detection (%s space): %d/%d flagged (%.0f%%)\n', ...
                p_outlr.Domain, n_ce_rejected, n_pick, n_ce_rejected / n_pick * 100);

            clean_local = sampled_local(~tf_ce_outlier(:)');
            n_needed    = n_pick - numel(clean_local);

            if n_needed > 0 && ~isempty(reserve_pos)
                n_topup     = min(n_needed, numel(reserve_pos));
                topup_local = non_responsive_local(reserve_pos(1:n_topup));

                if use_pca
                    tu_feat = X_feat(topup_local, :);
                    col_var_tu = var(tu_feat, 0, 1);
                    tu_feat_pca = tu_feat(:, col_var_tu > 1e-10);
                    n_pca_tu = min([p_outlr.NPCAComponents, size(tu_feat_pca, 1) - 1, size(tu_feat_pca, 2)]);
                    if n_pca_tu >= 1
                        [~, topup_detection_data] = pca(tu_feat_pca, 'NumComponents', n_pca_tu, 'Algorithm', 'eig');
                    else
                        topup_detection_data = tu_feat_pca;
                    end
                else
                    topup_detection_data = reduction(topup_local, :);
                end

                tf_topup_outlier = detectOutliers(topup_detection_data, p_outlr);
                topup_clean      = topup_local(~tf_topup_outlier(:)');

                clean_local = [clean_local, topup_clean];
                fprintf('Top-up: added %d/%d reserve candidates (%d clean)\n', ...
                    n_topup, numel(reserve_pos), numel(topup_clean));
            end

            counterexample_idx_local = clean_local;
        else
            counterexample_idx_local = sampled_local;
        end

        if isempty(counterexample_idx_local)
            warning('CellTypeClassifier:generateTrainLabels', ...
                'All counterexample candidates flagged as outliers — using original sample.');
            counterexample_idx_local = sampled_local;
        end
    else
        counterexample_idx_local = zeros(1, 0);
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

% Responsive strength for training inhibitory candidates
if ~isempty(ctc.ResponsiveStrength)
    labels.inhibitory_strength = ctc.ResponsiveStrength(unique_to_rep(in_train_id));
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

function tf_outlier = detectOutliers(data, p)
% DETECTOUTLIERS  Adaptive outlier detection for a candidate population.
%
% Two-stage approach:
%   1. Dip test on the first PC of data to test for multimodality
%   2a. Unimodal:   robust Mahalanobis distance + chi-squared threshold
%   2b. Multimodal: GMM with BIC-selected k + low-posterior outlier detection

    n = size(data, 1);
    d = size(data, 2);

    if n < 4
        tf_outlier = false(n, 1);
        return
    end

    n_pcs = min(d, n - 1);
    if n_pcs >= 1
        [~, scores] = pca(data, 'NumComponents', min(n_pcs, 10));
        pc1 = scores(:, 1);
    else
        pc1 = data(:, 1);
    end
    [~, p_dip] = diptest(pc1);

    is_multimodal = p_dip < p.DipTestAlpha;
    outlier_alpha = p.OutlierAlpha;

    if is_multimodal
        tf_outlier = gmmOutliers(data, p, outlier_alpha);
    else
        tf_outlier = mahalOutliers(data, d, outlier_alpha);
    end
end

function tf_outlier = mahalOutliers(data, d, outlier_alpha)
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
    n        = size(data, 1);
    max_k    = min(p.MaxResponsiveComponents, floor(n / 3));
    bic_vals = inf(1, max_k);
    gmm_fits = cell(1, max_k);

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

    gmm       = gmm_fits{best_k};
    posterior  = gmm.posterior(data);
    max_post   = max(posterior, [], 2);
    tf_outlier = max_post < outlier_alpha;

    fprintf('GMM selected k=%d components (BIC=%.1f); %d/%d outliers detected.\n', ...
        best_k, bic_vals(best_k), sum(tf_outlier), n);
end

function [dip, p_value] = diptest(x)
% DIPTEST  Hartigan's dip test for unimodality.
    x = sort(x(:));
    n = numel(x);

    if n < 4
        dip     = 0;
        p_value = 1;
        return
    end

    ecdf_y = (1:n)' / n;
    gcm    = ecdfGCM(ecdf_y);
    lcm    = ecdfLCM(ecdf_y);

    dip        = max(lcm - gcm) / 2;
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
