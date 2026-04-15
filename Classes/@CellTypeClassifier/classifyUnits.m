function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all unique units using the UMAP graph or feature-space kNN.
%
% Uses the NxN fuzzy simplicial set reconstructed from ctc.UMAP.head/tail/weights
% (set by simplicial_set_embedding) for transductive label propagation.
% No supervised UMAP is run.
%
% Classification methods (Classification.Method):
%   "graph" — label propagation on the UMAP graph (transductive, recommended)
%   "knn"   — distance-weighted kNN in feature space
%
% Requires: ctc.TrainLabels, ctc.NormalizationParams, ctc.UMAP (with head/tail/weights)
% (run generateTrainLabels first).
%
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
all_to_unique = nf.all_to_unique;

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

% -- Classification -----------------------------------------------------------
classify_method = lower(string(p_conf.Method));

unique_confidence = nan(1, n_unique);
unique_confidence(labels.sorted_train_ids) = 1.0;

if classify_method == "graph"
    % ── Graph-based label propagation on UMAP's own NxN fuzzy simplicial set ──
    % Reconstructed from ctc.UMAP.head/tail/weights (COO edge list set by
    % simplicial_set_embedding). Covers all N unique units. Symmetrise →
    % row-normalise → propagate, clamping training labels each iteration.
    assert(~isempty(ctc.UMAP) && ~isempty(ctc.UMAP.head), ...
        ['ctc.UMAP graph data is empty — run generateTrainLabels() before classifyUnits(). ' ...
         'The graph is populated by the unsupervised UMAP run in generateTrainLabels.']);

    % Reconstruct weighted graph from COO edge list (head/tail/weights are always
    % populated by simplicial_set_embedding; graph field is not reliable with Java SGD).
    G = sparse(ctc.UMAP.head, ctc.UMAP.tail, ctc.UMAP.weights, n_unique, n_unique);
    G = (G + G') / 2;
    row_sums = sum(G, 2);
    row_sums(row_sums == 0) = 1;
    W = G ./ row_sums;  % row-stochastic transition matrix

    % Per-unit count of direct training-unit neighbors in the UMAP graph.
    % Stored in UnitGraphConnectivity so callers can distinguish "robustly
    % classified" units from ones whose label comes from 1-2 training neighbors.
    graph_train_neighbors = full(sum(G(:, train_idx) > 0, 2))';  % (1 x n_unique)

    low_conn = sum(graph_train_neighbors(test_idx) < 3);
    if low_conn > 0.1 * sum(test_idx)
        warning('CellTypeClassifier:sparseGraphNeighborhood', ...
            ['%d/%d test units have fewer than 3 training neighbors in the ' ...
             'UMAP graph. Their confidence scores may reflect sparse ' ...
             'connectivity rather than robust classification.'], ...
            low_conn, sum(test_idx));
    end

    fprintf('Graph classification: %d nodes, propagating on UMAP NxN graph\n', n_unique);

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
    % Connectivity: for kNN path, count training neighbors within ConfidenceK.
    graph_train_neighbors = zeros(1, n_unique);
    if p_umap.AutoConfidenceK
        conf_k = max(5, round(sqrt(size(X_train, 1))));
        fprintf('Auto ConfidenceK: %d (N_train=%d)\n', conf_k, size(X_train, 1));
    else
        conf_k = p_umap.ConfidenceK;
    end

    Y_pred          = nan(1, size(X_test, 1));
    test_global_idx = find(test_idx);

    if ~isempty(X_test) && ~isempty(X_train)
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
            % kNN connectivity: k_actual training neighbors by definition
            graph_train_neighbors(test_global_idx(i)) = k_actual;
        end
        % Training units themselves have full k_actual neighbors
        graph_train_neighbors(train_idx) = k_actual;
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
full_labels       = unique_labels(all_to_unique);
full_confidence   = unique_confidence(all_to_unique);
full_connectivity = graph_train_neighbors(all_to_unique);

ctc.UnitLabels            = full_labels;
ctc.UnitConfidence        = full_confidence;
ctc.UnitGraphConnectivity = full_connectivity;

n_exc       = sum(unique_labels == 1, 'omitnan');
n_inh_final = sum(unique_labels == 2, 'omitnan');
n_unc       = sum(isnan(unique_labels));
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh_final, n_unc);

% -- Optional diagnostic ------------------------------------------------------
if ctc.Parameters.Diagnostics.Enable
    ctc.diagnosticClassification();
end
end
