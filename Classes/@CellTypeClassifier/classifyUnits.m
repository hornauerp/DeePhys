function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Two classification methods (Classification.Method):
%   "knn"   — distance-weighted kNN in feature space (default)
%   "graph" — label propagation on a kNN graph built in the supervised UMAP
%             embedding. The supervised UMAP reshapes the topology so that
%             same-class units cluster together; a kNN graph in this space
%             captures label-aware nonlinear neighborhoods. Training labels
%             are clamped and iteratively propagated through edge weights.
%
% Both methods run a supervised UMAP. For "knn" the embedding is visualization
% only; for "graph" it drives classification.
%
% Requires ctc.TrainLabels and ctc.NormalizationParams to be set
% (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)
%       ctc.UnitConfidence  class probability (graph) or vote fraction (kNN);
%                           1.0 for training units

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyUnits()');
assert(~isempty(ctc.NormalizationParams), ...
    'NormalizationParams not set — re-run generateTrainLabels() to store them.');

labels = ctc.TrainLabels;
p_umap = ctc.Parameters.UMAP;
p_conf = ctc.Parameters.Classification;
np     = ctc.NormalizationParams;

% -- Get normalized features (shared with generateTrainLabels) ----------------
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;

n_unique      = nf.n_unique;
n_units_full  = nf.n_units_full;  %#ok<NASGU>
all_to_unique = nf.all_to_unique;
unique_to_rep = nf.unique_to_rep; %#ok<NASGU>

% -- Apply stored NormalizationParams to ALL unique units ---------------------
X_all = normalize(nf.X_pergroup, 'center', np.mu_global, 'scale', np.sigma_global);
X_all(:, np.nan_cols)  = [];
X_all = X_all ./ np.scale;

if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X_all = X_all(:, np.feature_selection_mask);
end

% -- Dimension-normalised feature group weighting ----------------------------
% Scales ACG and Waveform groups so their total L2 contribution is
% proportional to ACGWeight / WaveformWeight regardless of group size.
% Applied here (not in generateTrainLabels) so training labels are unaffected
% and BayOpt can optimize weights independently.
if isfield(np, 'feat_names_trimmed') && ~isempty(np.feat_names_trimmed)
    feat_names_final = np.feat_names_trimmed;
    is_acg = startsWith(feat_names_final, "FullACG") | startsWith(feat_names_final, "ACG");
    is_wf  = startsWith(feat_names_final, "Waveform") | startsWith(feat_names_final, "ReferenceWaveform");
    n_acg  = sum(is_acg);
    n_wf   = sum(is_wf);
    if n_acg > 0
        X_all(:, is_acg) = X_all(:, is_acg) * sqrt(p_umap.ACGWeight / n_acg);
    end
    if n_wf > 0
        X_all(:, is_wf)  = X_all(:, is_wf)  * sqrt(p_umap.WaveformWeight / n_wf);
    end
end

train_idx = logical(labels.umap_train_idx);
test_idx  = ~train_idx;
X_train   = X_all(train_idx, :);
X_test    = X_all(test_idx, :);

% -- AutoNNeighbors for supervised UMAP ---------------------------------------
orig_sup_nneigh = p_umap.SupervisedNNeighbors;

if p_umap.AutoNNeighbors
    N_train    = size(X_train, 1);
    auto_sup_n = max(p_umap.MinNNeighbors, round(sqrt(N_train)));
    fprintf('Auto SupervisedNNeighbors: %d (N_train=%d)\n', auto_sup_n, N_train);
    ctc.Parameters.UMAP.SupervisedNNeighbors = auto_sup_n;
end

% -- Minority class ceiling on SupervisedNNeighbors ---------------------------
% Prevents n_neighbors from spanning a large fraction of the minority class,
% which causes UMAP to tear the manifold into disconnected components.
n_inh   = sum(labels.sorted_y_train == 2);
safe_nn = floor(n_inh / 5);
if safe_nn > 0 && ctc.Parameters.UMAP.SupervisedNNeighbors > safe_nn
    fprintf('Minority class ceiling: SupervisedNNeighbors %d -> %d (n_inh=%d)\n', ...
        ctc.Parameters.UMAP.SupervisedNNeighbors, safe_nn, n_inh);
    ctc.Parameters.UMAP.SupervisedNNeighbors = safe_nn;
end

% -- Supervised UMAP: get embedding (supervisorMatchedLabels not used) --------
% Classification is done via kNN in feature space below; UMAP embedding is
% stored for visualization only.
fprintf('supervisedUMAP input: X_train=%dx%d, X_test=%dx%d\n', ...
    size(X_train,1), size(X_train,2), size(X_test,1), size(X_test,2));
[~, ~, ~, test_reduction, train_reduction, extras] = supervisedUMAP(ctc, ...
    X_train, labels.sorted_y_train, nf.feat_names, X_test);

% Restore overridden params
ctc.Parameters.UMAP.SupervisedNNeighbors = orig_sup_nneigh;

ctc.Reduction.Train  = train_reduction;
ctc.Reduction.Test   = test_reduction;
ctc.Reduction.Extras = extras;

% -- Classification -----------------------------------------------------------
classify_method = lower(string(p_conf.Method));

unique_confidence = nan(1, n_unique);
unique_confidence(labels.sorted_train_ids) = 1.0;

if classify_method == "graph"
    % ── Graph-based label propagation on supervised UMAP embedding ──────
    % The supervised UMAP reshapes the topology so same-class units are
    % closer. Build a kNN graph in this label-aware space, then propagate
    % training labels through edge weights to classify test units.

    % Reconstruct full embedding in unique-unit order
    all_reduction = zeros(n_unique, size(train_reduction, 2));
    all_reduction(train_idx, :) = train_reduction;
    all_reduction(test_idx, :)  = test_reduction;

    if p_umap.AutoConfidenceK
        graph_k = max(5, round(sqrt(n_unique)));
    else
        graph_k = p_umap.ConfidenceK;
    end
    graph_k = min(graph_k, n_unique - 1);

    [knn_idx, knn_dists] = knnsearch(all_reduction, all_reduction, 'K', graph_k + 1);
    knn_idx   = knn_idx(:, 2:end);    % remove self
    knn_dists = knn_dists(:, 2:end);

    % Build sparse weighted adjacency (Gaussian kernel on distances)
    sigma = median(knn_dists(:, end));  % bandwidth from k-th neighbor distances
    if sigma == 0; sigma = 1; end
    rows = repmat((1:n_unique)', 1, graph_k);
    G = sparse(rows(:), knn_idx(:), exp(-knn_dists(:).^2 / (2 * sigma^2)), n_unique, n_unique);

    % Symmetrise and row-normalise
    G = (G + G') / 2;
    row_sums = sum(G, 2);
    row_sums(row_sums == 0) = 1;
    W = G ./ row_sums;  % row-stochastic transition matrix

    fprintf('Graph classification: %d nodes, k=%d, supervised UMAP %dD\n', ...
        n_unique, graph_k, size(all_reduction, 2));

    % Initialise label distributions: [P(exc), P(inh)] per unit
    F_label = ones(n_unique, 2) * 0.5;  % uninformed prior
    for ti = 1:numel(labels.sorted_train_ids)
        idx = labels.sorted_train_ids(ti);
        if labels.sorted_y_train(ti) == 1
            F_label(idx, :) = [1, 0];
        else
            F_label(idx, :) = [0, 1];
        end
    end

    % Iterative label propagation (clamp labeled nodes each iteration)
    max_iter = p_conf.GraphMaxIter;
    tol      = p_conf.GraphConvergenceTol;

    for iter = 1:max_iter
        F_new = W * F_label;

        % Re-normalise rows to valid distributions
        F_row_sum = sum(F_new, 2);
        F_row_sum(F_row_sum == 0) = 1;
        F_new = F_new ./ F_row_sum;

        % Clamp training labels
        for ti = 1:numel(labels.sorted_train_ids)
            idx = labels.sorted_train_ids(ti);
            if labels.sorted_y_train(ti) == 1
                F_new(idx, :) = [1, 0];
            else
                F_new(idx, :) = [0, 1];
            end
        end

        % Check convergence
        if max(abs(F_new(:) - F_label(:))) < tol
            fprintf('Graph label propagation converged at iteration %d\n', iter);
            F_label = F_new;
            break
        end
        F_label = F_new;
    end
    if iter == max_iter
        fprintf('Graph label propagation reached max iterations (%d)\n', max_iter);
    end

    % Extract predictions and confidence for test units
    test_global_idx = find(test_idx);
    Y_pred = nan(1, numel(test_global_idx));
    for i = 1:numel(test_global_idx)
        gi = test_global_idx(i);
        p_exc = F_label(gi, 1);
        p_inh = F_label(gi, 2);
        if p_exc >= p_inh
            Y_pred(i) = 1;
            unique_confidence(gi) = p_exc;
        else
            Y_pred(i) = 2;
            unique_confidence(gi) = p_inh;
        end
    end

else
    % ── Feature-space kNN classification ────────────────────────────────
    if p_umap.AutoConfidenceK
        conf_k = max(5, round(sqrt(size(train_reduction, 1))));
        fprintf('Auto ConfidenceK: %d (N_train=%d)\n', conf_k, size(train_reduction, 1));
    else
        conf_k = p_umap.ConfidenceK;
    end

    Y_pred          = nan(1, size(test_reduction, 1));
    test_global_idx = find(test_idx);

    if ~isempty(test_reduction) && ~isempty(train_reduction)
        k_actual = min(conf_k, size(X_train, 1));
        [nn_idx, nn_dists] = knnsearch(X_train, X_test, 'K', k_actual);
        for i = 1:numel(test_global_idx)
            neighbor_labels = labels.sorted_y_train(nn_idx(i, :));
            d       = nn_dists(i, :);
            eps_d   = max(d) * 1e-6;
            if eps_d == 0; eps_d = 1e-10; end
            w       = 1 ./ (d + eps_d);
            w_total = sum(w);
            w_exc   = sum(w(neighbor_labels == 1)) / w_total;
            w_inh   = sum(w(neighbor_labels == 2)) / w_total;
            if w_exc >= w_inh
                Y_pred(i) = 1;
                unique_confidence(test_global_idx(i)) = w_exc;
            else
                Y_pred(i) = 2;
                unique_confidence(test_global_idx(i)) = w_inh;
            end
        end
    end
end

% -- Assemble labels in unique-unit space -------------------------------------
unique_labels = nan(1, n_unique);
unique_labels(labels.sorted_train_ids) = labels.sorted_y_train;
unique_labels(test_idx)                = Y_pred;

% Optional confidence threshold
if p_conf.UseConfidenceThreshold
    low_conf = unique_confidence < p_conf.ConfidenceThreshold & ~train_idx;
    unique_labels(low_conf) = NaN;
    fprintf('Confidence threshold %.2f: %d units set to NaN\n', ...
        p_conf.ConfidenceThreshold, sum(low_conf));
end

% -- Broadcast unique-unit results back to all (unit x recording) rows --------
full_labels     = unique_labels(all_to_unique);
full_confidence = unique_confidence(all_to_unique);

ctc.UnitLabels     = full_labels;
ctc.UnitConfidence = full_confidence;

n_exc = sum(unique_labels == 1, 'omitnan');
n_inh_final = sum(unique_labels == 2, 'omitnan');
n_unc = sum(isnan(unique_labels));
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh_final, n_unc);

% -- Optional diagnostic ------------------------------------------------------
if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticClassification();
end
end
