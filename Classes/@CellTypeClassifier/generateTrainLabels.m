function generateTrainLabels(ctc, opts)
% GENERATETRAINLABELS  Build training labels via outlier detection and feature-space operations.
%
% Pipeline:
%   1. Build normalized feature cache (shared with classifyUnits via buildNormalizedFeatures)
%   2. Global z-score on ALL unique units; store NormalizationParams for classifyUnits
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
    opts.ReuseExistingUMAP (1,1) logical = true
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before generateTrainLabels()');

p_umap  = ctc.Parameters.UMAP;
p_outlr = ctc.Parameters.OutlierDetection;
p_train = ctc.Parameters.TrainLabels;

rng(ctc.Parameters.RNGSeed, 'twister');

% -- Clear cached state so parameters are always honoured ---------------------
% buildNormalizedFeatures is idempotent within a single pipeline run
% (classifyUnits reuses the same cache). Clearing here ensures that calling
% generateTrainLabels again — e.g. after changing Parameters — starts fresh.
ctc.clearCache();

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

% -- Global z-score on ALL unique units ----------------------------------------
% Fit normalization statistics on the full population, not the training subset.
% In the transductive setting (one CTC for all units), subset-fit statistics
% would scale the embedding around the training cultures alone, biasing the
% global UMAP and the label propagation step that follows.
[X_all, mu_global, sigma_global] = normalize(nf.X_pergroup);

% -- Remove NaN columns + scale -----------------------------------------------
nan_cols = any(isnan(X_all), 1);
X_all(:, nan_cols) = [];
scale = max(abs(X_all), [], 1);
scale(scale == 0) = 1;
X_all = X_all ./ scale;

% Track feature names through removal
feat_names_trimmed = nf.feat_names(~nan_cols);

% -- Optional feature selection -----------------------------------------------
if p_umap.FeatureSelection && size(X_all, 2) > 1
    feat_var   = var(X_all, 0, 1);
    var_thresh = prctile(feat_var, p_umap.MinVariancePercentile);
    low_var    = feat_var <= var_thresh;

    % Graph-based correlation removal: iteratively remove the feature with the
    % most high-correlation neighbors. This is order-independent and avoids
    % the systematic bias of sequential scanning (which always kept early ACG
    % bins at the expense of late ACG / waveform features).
    R = abs(corrcoef(X_all));
    R(logical(eye(size(R)))) = 0;
    high_corr = false(1, size(X_all, 2));
    adj = R > p_umap.MaxCorrelation;
    while any(adj(:))
        n_edges = sum(adj & ~high_corr' & ~high_corr, 2);
        n_edges(high_corr) = 0;
        [~, worst] = max(n_edges);
        high_corr(worst) = true;
        adj(worst, :) = false;
        adj(:, worst) = false;
    end

    remove_feat = low_var | high_corr;
    feature_selection_mask = ~remove_feat;
    X_all = X_all(:, feature_selection_mask);
    feat_names_trimmed = feat_names_trimmed(feature_selection_mask);
    fprintf('Feature selection: removed %d low-var + %d redundant = %d total (kept %d)\n', ...
        sum(low_var), sum(high_corr & ~low_var), sum(remove_feat), size(X_all, 2));
else
    feature_selection_mask = true(1, size(X_all, 2));
end

% -- Store NormalizationParams -------------------------------------------------
% Record the actual waveform window produced by buildFeatureMatrix so that
% classifyExternalUnits can verify the external feature matrix has the same
% dimensions (avoiding feature count mismatches across platforms).
n_wf_features = sum(startsWith(nf.feat_names, "Waveform"));
ph = ctc.Parameters.Harmonization;
target_sr = ph.WaveformTargetSamplingRate;
dt_ms = 1000 / target_sr;
n_pre_requested  = round(ph.WaveformPreTrough  * target_sr / 1000);
n_post_requested = round(ph.WaveformPostTrough * target_sr / 1000);
n_requested = n_pre_requested + n_post_requested + 1;

if n_wf_features == n_requested
    wf_pre_ms  = ph.WaveformPreTrough;
    wf_post_ms = ph.WaveformPostTrough;
else
    % Trim mode narrowed the window; back-compute actual pre/post from features.
    % Trough position is at n_out_pre+1; total = n_out_pre + n_out_post + 1.
    % We keep the original ratio between pre and post.
    ratio = n_pre_requested / (n_pre_requested + n_post_requested);
    n_out_pre_actual  = round((n_wf_features - 1) * ratio);
    n_out_post_actual = n_wf_features - 1 - n_out_pre_actual;
    wf_pre_ms  = n_out_pre_actual  * dt_ms;
    wf_post_ms = n_out_post_actual * dt_ms;
end

ctc.NormalizationParams = struct( ...
    'mu_global',              mu_global, ...
    'sigma_global',           sigma_global, ...
    'nan_cols',               nan_cols, ...
    'scale',                  scale, ...
    'feature_selection_mask', feature_selection_mask, ...
    'feat_names_trimmed',     feat_names_trimmed, ...
    'wf_pre_trough_ms',       wf_pre_ms, ...
    'wf_post_trough_ms',      wf_post_ms);

% -- Unsupervised UMAP embedding on ALL unique units --------------------------
% Run on the full dataset so the UMAP graph covers all units. Enables graph-based
% label propagation in classifyUnits and provides a global embedding for diagnostics.
% Skipped when ReuseExistingUMAP=true and a valid UMAP is already stored (e.g.
% after optimizeUnsupervisedUMAP), saving several minutes on large datasets.

N_all = size(X_all, 1);
if p_umap.AutoNNeighbors
    n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(N_all)));
    fprintf('Auto NNeighbors: %d (N=%d)\n', n_neighbors, N_all);
else
    n_neighbors = p_umap.NNeighbors;
end

umap_ready = ~isempty(ctc.UMAP) && ~isempty(ctc.Reduction) ...
             && isfield(ctc.Reduction, 'Unsupervised') ...
             && ~isempty(ctc.Reduction.Unsupervised) ...
             && size(ctc.Reduction.Unsupervised, 1) == N_all;

if opts.ReuseExistingUMAP && umap_ready
    reduction = ctc.Reduction.Unsupervised;
    fprintf('generateTrainLabels: reusing existing UMAP (N=%d). Pass ReuseExistingUMAP=false to force rerun.\n', N_all);
else
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
end

% -- Extract training-subset view (outlier detection and CE selection) ---------
% Normalization was fit on all units; subset is a row-slice of the full matrix.
X_subset = X_all(subset_mask, :);
X_feat   = X_subset;   % feature selection already applied to X_all above

% -- Setup for outlier detection and CE selection -----------------------------
subset_global_idx = find(subset_mask);  % maps subset-local → global unique-unit index
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

p_comm   = ctc.Parameters.Community;
has_explicit_ce = ~isempty(ctc.CounterexampleUnitIdx) && any(ctc.CounterexampleUnitIdx);

% Initialize diagnostic fields (overwritten below)
tf_resp_outlier             = false(1, numel(responsive_local));
resp_geom_mask              = logical([]);
resp_community_outlier_mask = logical([]);
resp_purity_outlier_mask    = logical([]);
ce_distances_from_inh       = [];
ce_outlier_mask_diag        = logical([]);
community_ids               = [];
inh_comm_ids                = [];
Q_modularity                = NaN;
inh_centroid                = [];
use_community_path          = false;
n_ce_pool_diag              = 0;

if has_explicit_ce
    % ── Explicit counterexample path (metadata labels) ────────────────────────
    % Inhibitory candidates: use isolation forest (same as before)
    outlier_domain = lower(string(p_outlr.Domain));
    contamination  = p_outlr.ContaminationFraction;
    gmm_sep        = p_outlr.AutoGMMSeparation;

    inh_iforest_input = prepareIforestInput(X_feat, responsive_local, ...
        reduction, subset_global_idx, outlier_domain, p_outlr);
    [tf_resp_outlier, resp_outlier_info] = runOutlierDetection( ...
        inh_iforest_input, contamination, gmm_sep);
    fprintf('Isolation forest (inh, %s, %s): %d/%d flagged (%.0f%%)%s\n', ...
        outlier_domain, resp_outlier_info.mode, sum(tf_resp_outlier), ...
        numel(tf_resp_outlier), sum(tf_resp_outlier)/numel(tf_resp_outlier)*100, ...
        resp_outlier_info.detail);

    clean_resp_local = responsive_local(~tf_resp_outlier(:)');
    if isempty(clean_resp_local)
        clean_resp_local = responsive_local;
        tf_resp_outlier  = false(size(tf_resp_outlier));
    end

    % Geometric consistency filter
    inh_centroid = mean(X_feat(clean_resp_local, :), 1);
    exc_centroid = mean(X_feat(~subset_responsive, :), 1);
    d_to_inh = pdist2(inh_centroid, X_feat(clean_resp_local, :), 'correlation');
    d_to_exc = pdist2(exc_centroid, X_feat(clean_resp_local, :), 'correlation');
    geom_consistent  = d_to_inh(:) < d_to_exc(:);
    resp_geom_mask   = ~geom_consistent(:)';
    clean_resp_local = clean_resp_local(geom_consistent);
    if isempty(clean_resp_local)
        clean_resp_local = responsive_local;
        resp_geom_mask   = false(1, numel(clean_resp_local));
        inh_centroid     = mean(X_feat(clean_resp_local, :), 1);
    end

    % Explicit CE units — clean with isolation forest
    ce_uid_table = unique(unit_ids_table(ctc.CounterexampleUnitIdx));
    ce_unique    = ismember(unique_ud_ids, ce_uid_table);
    ce_subset    = ce_unique(subset_mask);
    ce_local_all = find(ce_subset);
    ce_iforest_input = prepareIforestInput(X_feat, ce_local_all, ...
        reduction, subset_global_idx, outlier_domain, p_outlr);
    [tf_ce_outlier, ~] = runOutlierDetection(ce_iforest_input, contamination, gmm_sep);
    ce_outlier_mask_diag     = tf_ce_outlier(:)';
    counterexample_idx_local = ce_local_all(~tf_ce_outlier(:)');
    if ~isempty(inh_centroid)
        ce_distances_from_inh = pdist2(inh_centroid, X_feat(ce_local_all,:), 'correlation');
    end
    if isempty(counterexample_idx_local)
        counterexample_idx_local = ce_local_all;
    end
    n_ce_pool_diag = numel(ce_local_all);
    fprintf('Metadata counterexamples: %d labeled, %d after outlier cleaning\n', ...
        numel(ce_local_all), numel(counterexample_idx_local));

else
    % ── Inferred counterexample path (no metadata labels) ────────────────────
    odm = lower(string(p_outlr.Method));   % "community" | "iforest" | "none"

    if odm == "community"
        % ── Community path: Louvain + graphPurityOutlier + k-medoids CE ──────
        assert(~isempty(ctc.UMAP) && ~isempty(ctc.UMAP.head), ...
            ['ctc.UMAP graph data is empty. ' ...
             'Run optimizeUnsupervisedUMAP() or generateTrainLabels() first, ' ...
             'or set Parameters.OutlierDetection.Method = "iforest" to skip community detection.']);
        G_umap  = sparse(ctc.UMAP.head, ctc.UMAP.tail, ctc.UMAP.weights, N_all, N_all);
        G_umap  = (G_umap + G_umap') / 2;
        k_graph = min(n_neighbors, N_all - 1);

        W_lou  = full(G_umap);
        best_Q = -Inf;
        best_M = ones(N_all, 1);
        for r_lou = 1:p_comm.LouvainRestarts
            rng(ctc.Parameters.RNGSeed + r_lou, 'twister');
            try
                [M_r, Q_r] = community_louvain(W_lou, p_comm.LouvainResolution);
                if Q_r > best_Q; best_Q = Q_r; best_M = M_r(:); end
            catch; end
        end
        community_ids = best_M;
        Q_modularity  = best_Q;
        n_comm        = max(best_M);
        fprintf('Louvain: %d communities, Q=%.3f (resolution=%.2f, %d restarts)\n', ...
            n_comm, best_Q, p_comm.LouvainResolution, p_comm.LouvainRestarts);

        resp_global_idx_all = subset_global_idx(responsive_local);
        resp_global_all     = false(1, N_all);
        resp_global_all(resp_global_idx_all) = true;

        comm_resp_frac = arrayfun(@(c) ...
            sum(resp_global_all(best_M == c)) / max(sum(best_M == c), 1), 1:n_comm);
        max_frac   = max(comm_resp_frac);
        p_resp_gen = sum(resp_unique) / N_all;
        inh_by_rel = comm_resp_frac >= max_frac * p_comm.InhibitoryCommunityRelThresh;
        inh_by_abs = comm_resp_frac >= p_comm.EnrichmentFactor * p_resp_gen;
        inh_comm_ids = find(inh_by_rel & inh_by_abs);
        rel_thr = max_frac * p_comm.InhibitoryCommunityRelThresh;
        abs_thr = p_comm.EnrichmentFactor * p_resp_gen;
        fprintf('Inhibitory communities: %d/%d (rel_thr=%.3f [%.2f×max=%.3f], abs_thr=%.3f [%.1f×p_resp=%.3f])\n', ...
            numel(inh_comm_ids), n_comm, rel_thr, p_comm.InhibitoryCommunityRelThresh, max_frac, ...
            abs_thr, p_comm.EnrichmentFactor, p_resp_gen);

        in_inh_comm_global = ismember(best_M', inh_comm_ids);
        frac_resp_in_inh   = sum(resp_global_all & in_inh_comm_global) / max(sum(resp_global_all), 1);
        fprintf('Responsive units in inhibitory communities: %.1f%%\n', 100*frac_resp_in_inh);

        if frac_resp_in_inh < p_comm.CommunityFallbackThreshold
            warning('CellTypeClassifier:communityFallback', ...
                ['Only %.0f%% of responsive units in inhibitory communities ' ...
                 '(threshold %.0f%%).\nFalling back to iforest path.\n' ...
                 'Consider re-running optimizeUnsupervisedUMAP(), adjusting ' ...
                 'Parameters.Community settings, or setting ' ...
                 'Parameters.OutlierDetection.Method = "iforest".'], ...
                100*frac_resp_in_inh, 100*p_comm.CommunityFallbackThreshold);
            odm = "iforest";   % fall through to iforest block below
        else
            % Community path succeeded
            use_community_path = true;

            in_inh_resp = ismember(best_M(resp_global_idx_all)', inh_comm_ids);
            resp_community_outlier_mask = ~in_inh_resp;
            resp_in_inh_global = resp_global_idx_all(in_inh_resp);
            fprintf('Community outlier removal: %d/%d responsive units excluded\n', ...
                sum(resp_community_outlier_mask), numel(resp_global_idx_all));

            [resp_clean_global, resp_purity_outlier_mask] = graphPurityOutlier( ...
                reduction, resp_in_inh_global, resp_global_all, k_graph, ...
                p_comm.PuritySigmaThreshold);
            fprintf('Graph purity outlier removal: %d/%d remaining responsive excluded\n', ...
                sum(resp_purity_outlier_mask), numel(resp_in_inh_global));

            if isempty(resp_clean_global)
                warning('CellTypeClassifier:generateTrainLabels', ...
                    'All responsive candidates removed by outlier detection — using community filter only.');
                resp_clean_global        = resp_in_inh_global;
                resp_purity_outlier_mask = false(1, numel(resp_in_inh_global));
            end

            [~, loc] = ismember(resp_clean_global, subset_global_idx);
            clean_resp_local = loc(loc > 0);

            tf_resp_outlier = resp_community_outlier_mask(:)';
            pass_comm_idx   = find(~resp_community_outlier_mask);
            tf_resp_outlier(pass_comm_idx(resp_purity_outlier_mask)) = true;

            inh_centroid = mean(X_feat(clean_resp_local, :), 1);

            non_resp_global  = subset_global_idx(find(~subset_responsive));
            in_ce_comm       = ~ismember(best_M(non_resp_global)', inh_comm_ids);
            ce_pool_global   = non_resp_global(in_ce_comm);
            ce_pool_comm_ids = best_M(ce_pool_global)';
            ce_communities_u = unique(ce_pool_comm_ids);
            n_ce_pool        = numel(ce_pool_global);
            n_ce_pool_diag   = n_ce_pool;
            fprintf('CE pool: %d non-responsive units in %d communities\n', ...
                n_ce_pool, numel(ce_communities_u));

            n_clean_resp = numel(clean_resp_local);
            k_alloc      = zeros(1, numel(ce_communities_u));
            for ci = 1:numel(ce_communities_u)
                k_alloc(ci) = max(p_comm.MinCEPerCommunity, ...
                    round(n_clean_resp * sum(ce_pool_comm_ids == ce_communities_u(ci)) / n_ce_pool));
            end
            k_diff = sum(k_alloc) - n_clean_resp;
            if k_diff ~= 0
                [~, largest] = max(k_alloc);
                k_alloc(largest) = max(p_comm.MinCEPerCommunity, k_alloc(largest) - k_diff);
            end

            counterexample_global = zeros(1, 0);
            for ci = 1:numel(ce_communities_u)
                c        = ce_communities_u(ci);
                c_mask   = ce_pool_comm_ids == c;
                c_global = ce_pool_global(c_mask);
                k_c      = min(k_alloc(ci), sum(c_mask));
                if k_c < 1; continue; end
                umap_c   = reduction(c_global, :);
                rng(ctc.Parameters.RNGSeed, 'twister');
                [~, med_coords] = kmedoids(umap_c, k_c, 'Distance', 'sqeuclidean', 'Replicates', 3);
                med_idx = knnsearch(umap_c, med_coords);
                counterexample_global = [counterexample_global, c_global(med_idx(:))]; %#ok<AGROW>
            end

            [~, ce_loc] = ismember(counterexample_global, subset_global_idx);
            counterexample_idx_local = ce_loc(ce_loc > 0);
            fprintf('CE selection: %d medoids (target: %d = n_clean_resp)\n', ...
                numel(counterexample_idx_local), n_clean_resp);
        end
    end

    if odm == "iforest"
        % ── Isolation forest path ─────────────────────────────────────────────
        outlier_domain = lower(string(p_outlr.Domain));
        contamination  = p_outlr.ContaminationFraction;
        gmm_sep        = p_outlr.AutoGMMSeparation;

        inh_iforest_input = prepareIforestInput(X_feat, responsive_local, ...
            reduction, subset_global_idx, outlier_domain, p_outlr);
        [tf_resp_outlier, resp_outlier_info] = runOutlierDetection( ...
            inh_iforest_input, contamination, gmm_sep);
        fprintf('Iforest (inh, %s): %d/%d flagged%s\n', outlier_domain, ...
            sum(tf_resp_outlier), numel(tf_resp_outlier), resp_outlier_info.detail);
        clean_resp_local = responsive_local(~tf_resp_outlier(:)');
        if isempty(clean_resp_local); clean_resp_local = responsive_local; end

        inh_centroid = mean(X_feat(clean_resp_local, :), 1);
        exc_centroid = mean(X_feat(~subset_responsive, :), 1);
        d_to_inh = pdist2(inh_centroid, X_feat(clean_resp_local, :), 'correlation');
        d_to_exc = pdist2(exc_centroid, X_feat(clean_resp_local, :), 'correlation');
        geom_consistent  = d_to_inh(:) < d_to_exc(:);
        resp_geom_mask   = ~geom_consistent(:)';
        clean_resp_local = clean_resp_local(geom_consistent);
        if isempty(clean_resp_local)
            clean_resp_local = responsive_local;
            resp_geom_mask   = false(1, numel(clean_resp_local));
            inh_centroid     = mean(X_feat(clean_resp_local, :), 1);
        end

        non_responsive_local = find(~subset_responsive);
        ce_ratio    = p_outlr.CounterexampleRatio;
        n_pick      = min(ce_ratio * numel(clean_resp_local), numel(non_responsive_local));
        pool_dists  = pdist2(inh_centroid, X_feat(non_responsive_local,:), 'correlation');
        ce_distances_from_inh = pool_dists;
        dist_threshold = prctile(pdist2(inh_centroid, X_feat(clean_resp_local,:), 'correlation'), ...
            p_outlr.CounterexampleDistancePercentile);
        far_pool = non_responsive_local(pool_dists(:) >= dist_threshold);
        if numel(far_pool) < n_pick; far_pool = non_responsive_local; end
        n_ce_pool_diag = numel(far_pool);
        perm_order    = randperm(numel(far_pool));
        sampled_local = far_pool(perm_order(1:min(n_pick, numel(far_pool))));
        ce_iforest_input = prepareIforestInput(X_feat, sampled_local, ...
            reduction, subset_global_idx, outlier_domain, p_outlr);
        [tf_ce_outlier, ~] = runOutlierDetection(ce_iforest_input, contamination, gmm_sep);
        ce_outlier_mask_diag     = tf_ce_outlier(:)';
        counterexample_idx_local = sampled_local(~tf_ce_outlier(:)');
        if isempty(counterexample_idx_local); counterexample_idx_local = sampled_local; end

    elseif odm == "none"
        % ── No outlier detection ──────────────────────────────────────────────
        clean_resp_local = responsive_local;
        tf_resp_outlier  = false(1, numel(responsive_local));
        inh_centroid     = mean(X_feat(clean_resp_local, :), 1);
        fprintf('Outlier detection disabled: using all %d responsive units\n', numel(clean_resp_local));

        non_responsive_local = find(~subset_responsive);
        ce_ratio    = p_outlr.CounterexampleRatio;
        n_pick      = min(ce_ratio * numel(clean_resp_local), numel(non_responsive_local));
        pool_dists  = pdist2(inh_centroid, X_feat(non_responsive_local,:), 'correlation');
        ce_distances_from_inh = pool_dists;
        dist_threshold = prctile(pdist2(inh_centroid, X_feat(clean_resp_local,:), 'correlation'), ...
            p_outlr.CounterexampleDistancePercentile);
        far_pool = non_responsive_local(pool_dists(:) >= dist_threshold);
        if numel(far_pool) < n_pick; far_pool = non_responsive_local; end
        n_ce_pool_diag = numel(far_pool);
        perm_order    = randperm(numel(far_pool));
        counterexample_idx_local = far_pool(perm_order(1:min(n_pick, numel(far_pool))));
        fprintf('CE selection (no filtering): %d units from distance pool\n', numel(counterexample_idx_local));
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
% (non-empty only in explicit CE and fallback paths)
labels.resp_geom_mask = resp_geom_mask;

% Community-based path diagnostics (empty in explicit CE and fallback paths)
labels.resp_community_outlier_mask = resp_community_outlier_mask;
labels.resp_purity_outlier_mask    = resp_purity_outlier_mask;
labels.community_ids               = community_ids;     % (N_all x 1) community assignment
labels.inh_comm_ids                = inh_comm_ids;      % community IDs identified as inhibitory
labels.Q_modularity                = Q_modularity;      % best Louvain modularity
labels.use_community_path          = use_community_path;

% ce_outlier_mask: (1 x n_ce_sampled) — true = removed by outlier detection (explicit CE / fallback)
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
    'n_ce_in',            n_ce_pool_diag, ...
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

function coords = prepareIforestInput(X_feat, local_indices, ...
        reduction, subset_global_idx, domain, p_outlr)
% PREPAREIFORESTINPUT  Get coordinates for iforest based on domain.
%   "umap" — index into the full UMAP embedding using global indices
%   "pca"  — PCA on candidate features (legacy behavior)

    if domain == "umap"
        global_idx = subset_global_idx(local_indices);
        coords = reduction(global_idx, :);
    else
        candidate_features = X_feat(local_indices, :);
        col_var = var(candidate_features, 0, 1);
        candidate_features_clean = candidate_features(:, col_var > 1e-10);
        n_pca = min([p_outlr.NPCAComponents, ...
                     size(candidate_features_clean, 1) - 1, ...
                     size(candidate_features_clean, 2)]);
        if n_pca >= 1
            [~, coords] = pca(candidate_features_clean, ...
                'NumComponents', n_pca, 'Algorithm', 'eig');
        else
            coords = candidate_features_clean;
        end
    end
end

function [clean_global, outlier_mask] = graphPurityOutlier( ...
        reduction, candidate_global, is_responsive_global, k, sigma_thresh)
% GRAPHPURITYOUTLIER  Remove responsive units whose UMAP neighbourhood is
%   dominated by non-responsive units (likely mislabelled or boundary units).
%
% For each unit in candidate_global, computes the fraction of its k nearest
% neighbours (in the full UMAP embedding) that are also responsive. Units
% whose purity is more than sigma_thresh robust-z-scores below the median
% purity are flagged as outliers.
%
% INPUTS:
%   reduction          — (N_all x D) UMAP embedding of all unique units
%   candidate_global   — (1 x M) global indices into reduction for responsive candidates
%   is_responsive_global — (1 x N_all) logical, true = responsive unit
%   k                  — number of neighbours to use (clipped to N_all-1)
%   sigma_thresh       — robust z-score cutoff (default: 2.5)
%
% OUTPUTS:
%   clean_global  — (1 x M') global indices of candidates passing the filter
%   outlier_mask  — (1 x M) logical, true = flagged as outlier

    if isempty(candidate_global)
        clean_global = candidate_global;
        outlier_mask = false(1, 0);
        return
    end

    N_all = size(reduction, 1);
    k_act = min(k, N_all - 1);
    M     = numel(candidate_global);

    if k_act < 1 || M == 0
        clean_global = candidate_global;
        outlier_mask = false(1, M);
        return
    end

    % For each candidate, find k nearest neighbours in the full embedding
    [nn_all, ~] = knnsearch(reduction, reduction(candidate_global, :), 'K', k_act + 1);
    nn_all = nn_all(:, 2:end);   % exclude self (first column)

    % Purity: fraction of k-NN that are also responsive
    purity = mean(is_responsive_global(nn_all), 2);   % (M x 1)

    % Robust z-score: (purity - median) / (1.4826 * MAD)
    med_p = median(purity);
    mad_p = median(abs(purity - med_p));
    if mad_p == 0
        % No variation — no outliers
        outlier_mask = false(1, M);
        clean_global = candidate_global;
        return
    end

    robust_z  = (purity - med_p) / (1.4826 * mad_p);
    outlier_mask = (robust_z < -sigma_thresh)';   % (1 x M)
    clean_global = candidate_global(~outlier_mask);
end

function [tf_outlier, info] = runOutlierDetection(candidate_coords, contamination, gmm_sep)
% RUNOUTLIERDETECTION  Isolation forest with fixed or auto contamination.
%   Fixed mode (contamination is numeric): iforest with that ContaminationFraction.
%   Auto mode (contamination == "auto"): iforest without ContaminationFraction,
%     then fit 2-component GMM on anomaly scores. Flag the high-score component
%     only if the two components are well-separated (> gmm_sep pooled stds apart).

    n = size(candidate_coords, 1);

    if isstring(contamination) || ischar(contamination)
        is_auto = lower(string(contamination)) == "auto";
    else
        is_auto = false;
    end

    if ~is_auto
        % Fixed contamination fraction
        [forest, tf_outlier] = iforest(candidate_coords, ...
            'ContaminationFraction', contamination);
        info.mode   = sprintf('fixed=%.2f', contamination);
        info.detail = '';
        return
    end

    % Auto mode: infer contamination from anomaly score distribution
    if n < 10
        % Too few samples for reliable GMM — skip outlier detection
        tf_outlier = false(n, 1);
        info.mode   = 'auto';
        info.detail = ' — skipped (n < 10)';
        return
    end

    % Train iforest and get anomaly scores
    [forest, ~] = iforest(candidate_coords);
    [~, scores] = isanomaly(forest, candidate_coords);

    % Fit 2-component GMM on anomaly scores
    try
        gmm = fitgmdist(scores(:), 2, 'RegularizationValue', 1e-6, ...
            'Options', statset('MaxIter', 200, 'TolFun', 1e-6));
    catch
        % GMM failed — no outliers
        tf_outlier = false(n, 1);
        info.mode   = 'auto';
        info.detail = ' — GMM fit failed, no outliers removed';
        return
    end

    % Identify the high-score (outlier) component
    [~, outlier_comp] = max(gmm.mu);
    inlier_comp       = 3 - outlier_comp;

    % Compute separation: distance between means / pooled std
    pooled_std = sqrt(mean(gmm.Sigma(:)));
    if pooled_std == 0; pooled_std = 1; end
    separation = abs(gmm.mu(1) - gmm.mu(2)) / pooled_std;

    if separation >= gmm_sep
        % Components are well-separated — flag outlier component
        posterior  = gmm.posterior(scores(:));
        tf_outlier = posterior(:, outlier_comp) > 0.5;
        info.mode   = 'auto';
        info.detail = sprintf(' — GMM separation=%.1f', separation);
    else
        % No clear bimodal structure — no outliers
        tf_outlier = false(n, 1);
        info.mode   = 'auto';
        info.detail = sprintf(' — no bimodal separation (%.1f < %.1f)', ...
            separation, gmm_sep);
    end
end

