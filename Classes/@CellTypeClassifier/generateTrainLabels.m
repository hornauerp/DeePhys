function generateTrainLabels(ctc)
% GENERATETRAINLABELS  Build training labels via UMAP embedding and adaptive outlier detection.
%
<<<<<<< HEAD
% Projects all units into a low-dimensional UMAP space, then applies adaptive
% outlier detection on responsive-unit candidates:
%   1. Dip test on first PC of responsive units — determines uni/multimodality
%   2. Unimodal: robust Mahalanobis distance + chi-squared threshold
%   3. Multimodal: GMM with BIC-selected k + low-posterior outliers
%
% Counterexample selection (non-responsive units far from the responsive centroid)
% is performed in the normalized feature space (deterministic, hyperparameter-free)
% rather than the UMAP embedding (which changes with UMAP hyperparameters).
%
% ResponsiveClassLabel controls which class responsive units are assigned to:
%   2 (default) — responsive = inhibitory   (stimulation activates inhibitory neurons)
%   1           — responsive = excitatory   (direct excitatory stimulation paradigm)
=======
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
>>>>>>> claude/dreamy-varahamihira
%
% Normalisation note: the unsupervised UMAP stage fits z-score on the entire
% unit subset (no train/test split exists yet). classifyUnits later fits only
% on the training partition — the difference is intentional.
%
% Requires ctc.ResponsiveUnitIdx to be set (run identifyResponsiveUnits first).
<<<<<<< HEAD
% Sets: ctc.TrainLabels, ctc.UMAP, ctc.Reduction.Unsupervised
=======
% Sets: ctc.TrainLabels, ctc.UMAP, ctc.NormalizationParams
>>>>>>> claude/dreamy-varahamihira

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
p_train = ctc.Parameters.TrainLabels;

rng(ctc.Parameters.RNGSeed, 'twister');

% ── Extract harmonized waveforms + ACGs from UnitDataArray ─────────────────
[wf, acg, sr] = ctc.getOrExtract(ctc.UnitDataArray);
[X_raw, ~]    = buildFeatureMatrix(ctc, wf, acg, sr);

<<<<<<< HEAD
% ── Step 1: per-group normalisation using FeatureStore table column ────────
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
=======
% ── Chip-level normalisation (all units) ──────────────────────────────────────
% NaN handling: leave as NaN so they propagate to column removal (step 3).
% This avoids bias from replacing NaN with 0 before z-scoring.
[iG, G] = rg.combineMetadataIndices(ctc.UnitList, p_umap.NormalizationVar);
g_idx = unique(iG);
for g = 1:length(g_idx)
    mask = (iG == g_idx(g));
    X_raw(mask, :) = normalize(X_raw(mask, :));
>>>>>>> claude/dreamy-varahamihira
end

% ── Determine training subset (by culture index or use all) ────────────────
n_units = size(X_raw, 1);
if ~isempty(p_umap.TrainingCultureIdx)
    subset_mask = CellTypeClassifier.buildCultureSubsetMask( ...
        ctc.FeatureStore, p_umap.TrainingCultureIdx, ctc.Parameters.CultureKeys);
else
    subset_mask = true(1, n_units);
end

X_subset = X_raw(subset_mask, :);

<<<<<<< HEAD
% Global z-score + scale on subset (whole subset — no train/test split at this stage)
[X_subset, ~, ~] = normalize(X_subset);
=======
% Global z-score on subset — store parameters for carryover
if length(G) > 1
    [X_subset, mu_global, sigma_global] = normalize(X_subset);
else
    mu_global    = zeros(1, size(X_subset, 2));
    sigma_global = ones(1, size(X_subset, 2));
end

% Remove NaN columns and store which ones
>>>>>>> claude/dreamy-varahamihira
nan_cols = any(isnan(X_subset), 1);
X_subset(:, nan_cols) = [];
scale = max(abs(X_subset), [], 1);
scale(scale == 0) = 1;
X_subset = X_subset ./ scale;

<<<<<<< HEAD
% ── Unsupervised UMAP embedding ─────────────────────────────────────────────
=======
% ── Store normalisation parameters for carryover ─────────────────────────────
ctc.NormalizationParams = struct( ...
    'mu_global',    mu_global, ...
    'sigma_global', sigma_global, ...
    'nan_cols',     nan_cols, ...
    'scale',        scale);

% ── UMAP embedding (unsupervised, on subset only) ───────────────────────────
>>>>>>> claude/dreamy-varahamihira
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

% ── Adaptive outlier detection on responsive candidates ───────────────────
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

<<<<<<< HEAD
clean_resp_local = responsive_local(~tf_outlier(:)');

if isempty(clean_resp_local)
    warning('CellTypeClassifier:generateTrainLabels', ...
        'All responsive candidates flagged as outliers — using all candidates.');
    clean_resp_local = responsive_local;
    tf_outlier = false(size(tf_outlier));
end

% ── Select counterexamples in feature space (not UMAP space) ──────────────
% Feature-space distance is deterministic and independent of UMAP hyperparameters.
centroid_feat    = mean(X_subset(clean_resp_local, :), 1);
dists            = pdist2(centroid_feat, X_subset);
candidate_dists  = pdist2(centroid_feat, X_subset(clean_resp_local, :));
dist_threshold   = prctile(candidate_dists, p_outlr.DistancePercentile);

far_idx_local       = find(dists > dist_threshold);
n_candidates        = numel(clean_resp_local);
n_counterexamples   = p_outlr.CounterexampleRatio * n_candidates;
counterexample_idx  = randsample(far_idx_local, min(n_counterexamples, numel(far_idx_local)), false);

% ── Map local indices back to global UnitDataArray indices ─────────────────
subset_global = find(subset_mask);

resp_label    = p_train.ResponsiveClassLabel;  % default 2 (inhibitory)
counter_label = 3 - resp_label;                % 1 if resp=2, 2 if resp=1

in_train_id  = subset_global(clean_resp_local);
ex_train_raw = subset_global(counterexample_idx);
ex_train_id  = ex_train_raw(~ismember(ex_train_raw, in_train_id));
=======
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
>>>>>>> claude/dreamy-varahamihira

train_ids = [ex_train_id, in_train_id];
y_train   = [counter_label * ones(1, numel(ex_train_id)), ...
             resp_label    * ones(1, numel(in_train_id))];
[sorted_train_ids, sort_idx] = sort(train_ids, 'ascend');

labels.sorted_train_ids  = sorted_train_ids;
labels.sorted_y_train    = y_train(sort_idx);
labels.umap_train_idx    = false(1, n_units);
labels.umap_train_idx(sorted_train_ids) = true;

<<<<<<< HEAD
outlier_local = responsive_local(tf_outlier(:)');
labels.outlier_global_idx              = subset_global(outlier_local);
labels.excitatory_candidate_global_idx = subset_global(counterexample_idx);
=======
% Store outlier and counterexample info for diagnostic plotting
outlier_local = responsive_in_subset(tf_forest(:)');
labels.outlier_global_idx        = subset_global(outlier_local);
labels.excitatory_candidate_global_idx = subset_global(counterexample_idx_local);
>>>>>>> claude/dreamy-varahamihira

% Store responsive strength for the training inhibitory candidates (for diagnostics)
if ~isempty(ctc.ResponsiveStrength)
    labels.inhibitory_strength = ctc.ResponsiveStrength(in_train_id);
else
    labels.inhibitory_strength = ones(1, length(in_train_id));
end

ctc.TrainLabels = labels;
fprintf('Training set: %i %s, %i %s candidates\n', ...
    sum(y_train == counter_label), labelName(counter_label), ...
    sum(y_train == resp_label),    labelName(resp_label));
end

<<<<<<< HEAD
%% ── Helper functions ────────────────────────────────────────────────────────

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

    n = size(candidate_reduction, 1);
    d = size(candidate_reduction, 2);

    if n < 4
        % Too few candidates for meaningful outlier detection
        tf_outlier = false(n, 1);
        return
    end

    % ── Stage 1: dip test for multimodality ──────────────────────────────
    % Project onto first principal component for 1D modality test
    [~, scores] = pca(candidate_reduction, 'NumComponents', 1);
    pc1 = scores(:, 1);
    [~, p_dip] = diptest(pc1);
    is_multimodal = p_dip < p.OutlierAlpha;

    if is_multimodal
        % ── Stage 2b: GMM + posterior-based outlier detection ────────────
        tf_outlier = gmmOutliers(candidate_reduction, p);
    else
        % ── Stage 2a: robust Mahalanobis + chi-squared ──────────────────
        tf_outlier = mahalOutliers(candidate_reduction, d, p);
    end
end

function tf_outlier = mahalOutliers(data, d, p)
% Robust Mahalanobis distance outlier detection.
% Uses robustcov (minimum covariance determinant) for robust center/covariance.
% robustcov returns [sigma, mu, mahal_distances, outlier_flags].
    n = size(data, 1);
    try
        [~, ~, mahal_d] = robustcov(data);
        thresh  = chi2inv(1 - p.OutlierAlpha, d);
        tf_outlier = mahal_d > thresh;
    catch
        % robustcov can fail for near-singular data (too few points or low rank)
        warning('CellTypeClassifier:robustcovFailed', ...
            'Robust covariance estimation failed — no outliers removed.');
        tf_outlier = false(n, 1);
    end
end

function tf_outlier = gmmOutliers(data, p)
% GMM-based outlier detection for multimodal responsive populations.
% Fits GMMs with k=1..MaxResponsiveComponents, selects by BIC.
% Outliers: units with max posterior < OutlierAlpha under all components.
    n         = size(data, 1);
    max_k     = min(p.MaxResponsiveComponents, floor(n / 3));  % need at least 3 points per component
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
            % fitgmdist may fail for certain k values
        end
    end

    [~, best_k] = min(bic_vals);

    if isinf(bic_vals(best_k))
        % All GMM fits failed
        warning('CellTypeClassifier:gmmFailed', ...
            'All GMM fits failed — no outliers removed.');
        tf_outlier = false(n, 1);
        return
    end

    gmm = gmm_fits{best_k};
    posterior = gmm.posterior(data);        % (n x k) posterior probabilities
    max_post  = max(posterior, [], 2);      % max posterior across all components
    tf_outlier = max_post < p.OutlierAlpha;

    fprintf('GMM selected k=%d components (BIC=%.1f); %d/%d outliers detected.\n', ...
        best_k, bic_vals(best_k), sum(tf_outlier), n);
end

function [dip, p_value] = diptest(x)
% DIPTEST  Hartigan's dip test for unimodality.
%
% Computes the dip statistic D_n and an approximate p-value for testing
% whether a 1D sample is drawn from a unimodal distribution.
%
% Small p_value → reject unimodality (evidence of multimodality).
%
% INPUTS:
%   x - (N x 1) or (1 x N) sample vector
%
% OUTPUTS:
%   dip     - Hartigan's dip statistic
%   p_value - approximate p-value (Hartigan & Hartigan 1985)
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

    % Empirical CDF
    ecdf_y = (1:n)' / n;

    % Greatest convex minorant (GCM) and least concave majorant (LCM)
    % of the empirical CDF
    gcm = ecdfGCM(ecdf_y);
    lcm = ecdfLCM(ecdf_y);

    % Dip = max distance between GCM and LCM, divided by 2
    dip = max(lcm - gcm) / 2;

    % Approximate p-value using the asymptotic distribution
    % (Hartigan & Hartigan 1985, Table 1: critical values for sqrt(n)*D)
    dip_scaled = sqrt(n) * dip;

    % Interpolation from Table 1 of Hartigan & Hartigan (1985)
    % Critical values of sqrt(n)*D for significance levels
    alpha_table = [0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00];
    crit_table  = [1.20, 1.06, 0.89, 0.78, 0.66, 0.50, 0.35];

    if dip_scaled >= crit_table(1)
        p_value = alpha_table(1) / 2;  % < 0.01
    elseif dip_scaled <= crit_table(end)
        p_value = 1.0;
    else
        p_value = interp1(crit_table, alpha_table, dip_scaled, 'linear', 'extrap');
        p_value = max(0, min(1, p_value));
    end
end

function gcm = ecdfGCM(y)
% Greatest convex minorant of a vector (isotonic regression from below).
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
        % Fill in the line segment from i to j
        for k = i+1:j
            gcm(k) = gcm(i) + (y(j) - gcm(i)) * (k - i) / (j - i);
        end
        i = j;
    end
end

function lcm = ecdfLCM(y)
% Least concave majorant of a vector (isotonic regression from above).
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
=======
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
>>>>>>> claude/dreamy-varahamihira
    end
end
