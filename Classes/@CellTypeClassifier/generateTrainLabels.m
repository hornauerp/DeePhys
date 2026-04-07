function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels via UMAP embedding and adaptive outlier detection.
%
% Pipeline:
%   1. Normalise features (per-group z-score, global z-score, NaN removal, scaling)
%   2. Unsupervised UMAP embedding
%   3. Outlier detection on responsive candidates in UMAP space
%   4. Counterexample selection (uniform random for drug-response; explicit for metadata)
%   5. Outlier detection on counterexamples in UMAP space (metadata path only)
%   6. Assemble training labels
%
% Outlier detection in UMAP space (dip test + Mahalanobis/GMM):
%   Since classification happens through UMAP (supervised embedding), outlier
%   detection should operate in the same representation. This finds the
%   well-defined main population of responsive units, deliberately accepting
%   loss of smaller inhibitory subpopulations.
%
% Two label-source strategies:
%
%   Drug response (inferred counterexamples):
%     Responsive units form one-class ground truth (e.g. inhibitory).
%     Counterexamples are uniformly randomly sampled from the non-responsive
%     pool. Uniform sampling preserves the natural excitatory distribution;
%     contaminating inhibitory non-responders are caught by the subsequent
%     outlier detection step.
%
%   Pure culture / metadata (explicit counterexamples):
%     Both classes have ground truth from ctc.CounterexampleUnitIdx.
%     Both classes are cleaned by outlier detection independently.
%
% AutoSupervisedNDims: estimates the intrinsic dimensionality of the
% normalised feature manifold via the TWO-NN estimator (Facco et al. 2017),
% then sets SupervisedNDims accordingly.
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

% -- Deduplicate: one representative UnitData per unique UnitID --------------
[unique_ud, all_to_unique, ~] = ctc.uniqueUnitMap();
n_units_full = numel(ctc.UnitDataArray);

n_unique      = numel(unique_ud);
unique_to_rep = zeros(1, n_unique);
for u = 1:n_unique
    rows = find(all_to_unique == u);
    unique_to_rep(u) = rows(1);
end
resp_unique = ctc.ResponsiveUnitIdx(unique_to_rep);

% -- Extract harmonized waveforms + ACGs from unique units only --------------
[wf, acg, sr] = ctc.getOrExtract(unique_ud);
[X_raw, ~]    = buildFeatureMatrix(ctc, wf, acg, sr);

% -- Step 1: per-group normalisation -----------------------------------------
X_raw(isnan(X_raw)) = 0;
norm_var = p_umap.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col = ctc.FeatureStore.UnitTable.(norm_var);
    group_col_unique = group_col(unique_to_rep);
    [~, ~, iG] = unique(string(group_col_unique), 'stable');
    for g = 1:max(iG)
        mask = iG == g;
        X_raw(mask, :) = normalize(X_raw(mask, :));
    end
elseif ~isempty(norm_var)
    warning('CellTypeClassifier:normVarNotFound', ...
        'NormalizationVar "%s" not found in UnitTable — per-group z-score skipped.', norm_var);
end

% -- Determine training subset -----------------------------------------------
n_units = n_unique;
if ~isempty(p_umap.TrainingCultureIdx)
    subset_mask_full = CellTypeClassifier.buildCultureSubsetMask( ...
        ctc.FeatureStore, p_umap.TrainingCultureIdx, ctc.Parameters.CultureKeys);
    subset_mask = subset_mask_full(unique_to_rep);
else
    subset_mask = true(1, n_units);
end

X_subset = X_raw(subset_mask, :);

% -- Step 2: Global z-score --------------------------------------------------
[X_subset, mu_global, sigma_global] = normalize(X_subset);

% -- Step 3: Remove NaN columns + scale --------------------------------------
nan_cols = any(isnan(X_subset), 1);
X_subset(:, nan_cols) = [];
scale = max(abs(X_subset), [], 1);
scale(scale == 0) = 1;
X_subset = X_subset ./ scale;

% -- Step 4: Optional feature selection --------------------------------------
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
    X_feat = X_subset(:, feature_selection_mask);
    fprintf('Feature selection: removed %d low-var + %d redundant = %d total (kept %d)\n', ...
        sum(low_var), sum(high_corr & ~low_var), sum(remove_feat), size(X_feat, 2));
else
    feature_selection_mask = true(1, size(X_subset, 2));
    X_feat = X_subset;
end

% -- Store normalisation parameters for carryover ----------------------------
ctc.NormalizationParams = struct( ...
    'mu_global',             mu_global, ...
    'sigma_global',          sigma_global, ...
    'nan_cols',              nan_cols, ...
    'scale',                 scale, ...
    'feature_selection_mask', feature_selection_mask);

% -- AutoSupervisedNDims: TWO-NN intrinsic dimensionality estimation ----------
if p_umap.AutoSupervisedNDims
    id_est = estimateIntrinsicDimensionality(X_feat);
    sup_ndims = max(p_umap.MinSupervisedNDims, min(p_umap.MaxSupervisedNDims, round(id_est)));
    fprintf('TWO-NN intrinsic dimensionality estimate: %.1f -> SupervisedNDims = %d\n', id_est, sup_ndims);
    ctc.Parameters.UMAP.SupervisedNDims = sup_ndims;
end

% -- Unsupervised UMAP embedding (mandatory for outlier detection) -----------
assert(~isempty(p_umap.TemplateDir), ...
    ['CellTypeClassifier:noTemplateDir  ', ...
    'Parameters.UMAP.TemplateDir must be set before calling generateTrainLabels. ' ...
    'Example: parameters.UMAP.TemplateDir = ''C:/path/to/output'';']);
assert(isfolder(p_umap.TemplateDir), ...
    'CellTypeClassifier:templateDirMissing  %s does not exist.', p_umap.TemplateDir);

n_dims = p_umap.NDims;
N_feat = size(X_feat, 1);
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_feat)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_feat);
else
    n_neighbors = p_umap.NNeighbors;
end

template_file = fullfile(p_umap.TemplateDir, 'ctc_umap_template.mat');
[reduction, umap_model, ~, ~] = run_umap(X_feat, ...
    'n_components',  n_dims, ...
    'n_neighbors',   n_neighbors, ...
    'min_dist',      p_umap.MinDist, ...
    'spread',        p_umap.Spread, ...
    'sgd_tasks',     20, ...
    'method', 'Java', 'verbose', 'none', ...
    'save_template_file', template_file);
ctc.UMAP = umap_model;
if isempty(ctc.Reduction); ctc.Reduction = struct(); end
ctc.Reduction.Unsupervised = reduction;

% -- Adaptive outlier detection on responsive candidates (UMAP space) --------
% Outlier detection in UMAP space matches the representation where
% classification will happen. This finds the well-defined main population
% and accepts loss of smaller inhibitory subpopulations.
subset_responsive = resp_unique(subset_mask);
responsive_local  = find(subset_responsive);
candidate_reduction = reduction(subset_responsive, :);

tf_resp_outlier = detectOutliers(candidate_reduction, p_outlr);

n_rejected = sum(tf_resp_outlier);
n_total    = numel(tf_resp_outlier);
if n_rejected > 0
    fprintf('Responsive outlier detection (UMAP space): %d/%d flagged (%.0f%%)\n', ...
        n_rejected, n_total, n_rejected / n_total * 100);
end

clean_resp_local = responsive_local(~tf_resp_outlier(:)');

if isempty(clean_resp_local)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'All responsive candidates flagged as outliers — using all candidates.');
    clean_resp_local = responsive_local;
    tf_resp_outlier = false(size(tf_resp_outlier));
end

% -- Counterexample selection ---------------------------------------------------
n_responsive = numel(clean_resp_local);

has_explicit_ce = ~isempty(ctc.CounterexampleUnitIdx) && any(ctc.CounterexampleUnitIdx);

if has_explicit_ce
    % Metadata method: use explicit ground-truth counterexamples directly.
    ce_unique  = ctc.CounterexampleUnitIdx(unique_to_rep);
    ce_subset  = ce_unique(subset_mask);
    ce_local_all = find(ce_subset);

    % Outlier detection on counterexample class in UMAP space.
    ce_reduction = reduction(ce_subset, :);
    tf_ce_outlier = detectOutliers(ce_reduction, p_outlr);
    n_ce_rejected = sum(tf_ce_outlier);
    if n_ce_rejected > 0
        fprintf('Counterexample outlier detection (UMAP space): %d/%d flagged (%.0f%%)\n', ...
            n_ce_rejected, numel(tf_ce_outlier), n_ce_rejected / numel(tf_ce_outlier) * 100);
    end
    counterexample_idx_local = ce_local_all(~tf_ce_outlier(:)');

    if isempty(counterexample_idx_local)
        warning('CellTypeClassifier:generateTrainLabels', ...
            'All counterexample candidates flagged as outliers — using all candidates.');
        counterexample_idx_local = ce_local_all;
    end
    fprintf('Metadata counterexamples: %d labeled, %d after outlier cleaning\n', ...
        numel(ce_local_all), numel(counterexample_idx_local));
else
    % Inferred method: uniform random sampling from non-responsive pool,
    % followed by outlier detection in UMAP space. Removed units are replaced
    % by drawing additional candidates from the remaining pool so that the
    % target count (ce_ratio * n_responsive) is preserved where possible.
    % Contaminating inhibitory units that failed to respond cluster near the
    % inhibitory region in UMAP space and appear as a separate mode; the
    % dip test + GMM approach detects and removes them.
    non_responsive_local = find(~subset_responsive);
    ce_ratio          = p_outlr.CounterexampleRatio;
    n_target          = ce_ratio * n_responsive;
    n_pool            = numel(non_responsive_local);
    n_pick            = min(n_target, n_pool);

    if n_pick > 0
        % Draw initial sample
        pool_order    = randperm(n_pool);
        sampled_pos   = pool_order(1:n_pick);       % positions into non_responsive_local
        reserve_pos   = pool_order(n_pick+1:end);   % remaining pool for top-up

        sampled_local = non_responsive_local(sampled_pos);

        % Outlier detection on initial sample
        ce_reduction  = reduction(sampled_local, :);
        tf_ce_outlier = detectOutliers(ce_reduction, p_outlr);
        n_ce_rejected = sum(tf_ce_outlier);

        if n_ce_rejected > 0
            fprintf('Counterexample outlier detection (UMAP space): %d/%d flagged (%.0f%%)\n', ...
                n_ce_rejected, n_pick, n_ce_rejected / n_pick * 100);

            clean_local = sampled_local(~tf_ce_outlier(:)');
            n_needed    = n_pick - numel(clean_local);

            % Top-up from reserve pool if units were removed
            if n_needed > 0 && ~isempty(reserve_pos)
                n_topup    = min(n_needed, numel(reserve_pos));
                topup_local = non_responsive_local(reserve_pos(1:n_topup));

                % Outlier-check top-up candidates before accepting
                topup_reduction  = reduction(topup_local, :);
                tf_topup_outlier = detectOutliers(topup_reduction, p_outlr);
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

outlier_local = responsive_local(tf_resp_outlier(:)');
labels.outlier_global_idx              = subset_global(outlier_local);
labels.excitatory_candidate_global_idx = subset_global(counterexample_idx_local);

% Store mapping so classifyUnits can broadcast unique-unit results back to all rows
labels.n_unique         = n_unique;
labels.all_to_unique    = all_to_unique;
labels.n_units_full     = n_units_full;
labels.has_explicit_ce  = has_explicit_ce;

% Store responsive strength for training inhibitory candidates (for diagnostics)
if ~isempty(ctc.ResponsiveStrength)
    labels.inhibitory_strength = ctc.ResponsiveStrength(unique_to_rep(in_train_id));
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
    if label == 1; name = 'excitatory';
    elseif label == 2; name = 'inhibitory';
    else; name = sprintf('class_%d', label);
    end
end

function d_est = estimateIntrinsicDimensionality(X)
% ESTIMATEINTRINSCDIMENSIONALITY  TWO-NN estimator (Facco et al. 2017).
%
% Estimates intrinsic dimensionality from the ratio of distances to the
% second and first nearest neighbors. Fits a Pareto distribution to the
% mu = r2/r1 ratios via maximum likelihood.
%
% Reference:
%   Facco, E., d'Errico, M., Rodriguez, A. & Laio, A. (2017).
%   Estimating the intrinsic dimension of datasets by a minimal
%   neighborhood information. Scientific Reports, 7, 12140.

    [n, ~] = size(X);
    if n < 10
        d_est = 2;
        return
    end

    % Find 2 nearest neighbors (excluding self)
    [~, nn_dists] = knnsearch(X, X, 'K', 3);  % col 1 = self (dist 0)
    r1 = nn_dists(:, 2);  % distance to 1st neighbor
    r2 = nn_dists(:, 3);  % distance to 2nd neighbor

    % Remove degenerate cases (identical points)
    valid = r1 > 0;
    r1 = r1(valid);
    r2 = r2(valid);

    if numel(r1) < 5
        d_est = 2;
        return
    end

    mu = r2 ./ r1;

    % MLE for Pareto exponent: d = n / sum(log(mu))
    % The CDF of mu under d-dimensional uniformity is P(mu <= x) = 1 - x^(-d)
    d_est = numel(mu) / sum(log(mu));
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

    dip_alpha = p.DipTestAlpha;
    is_multimodal = p_dip < dip_alpha;
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
