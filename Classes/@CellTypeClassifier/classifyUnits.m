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

% Delegate to ensemble path when enabled (majority vote across multiple seeds).
if ctc.Parameters.Ensemble.Enabled
    ctc.classifyUnitsEnsemble();
    return
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

train_idx = logical(labels.umap_train_idx);
test_idx  = ~train_idx;

% -- Classification -----------------------------------------------------------
classify_method = lower(string(p_conf.Method));

unique_confidence = nan(1, n_unique);
unique_confidence(labels.sorted_train_ids) = 1.0;

if classify_method == "graph"
    % ── Graph-based label propagation on the unit graph ──────────────────
    % ctc.graphMatrix(): the exact kNN fuzzy graph (Graph.Method "exact") or the UMAP
    % model's graph. ACGWeight and WaveformWeight were applied when it was built in
    % generateTrainLabels, so they are not re-applied here. Covers all N unique units;
    % row-normalise → propagate, clamping training labels each iteration.
    G = ctc.graphMatrix();
    assert(size(G, 1) == n_unique, 'CellTypeClassifier:graphSize', ...
        'Unit graph has %d nodes but there are %d unique units — rerun generateTrainLabels().', ...
        size(G, 1), n_unique);

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

    % Label propagation with the training labels clamped; unreached units (no graph path
    % to a training unit) are left unlabelled (NaN) instead of being called excitatory by
    % the tie rule. F is seed-count normalised (see CellTypeClassifier.propagateLabels).
    [F_label, unreached, converged, iter] = CellTypeClassifier.propagateLabels(G, ...
        labels.sorted_train_ids, labels.sorted_y_train, p_conf.GraphMaxIter, p_conf.GraphConvergenceTol);
    if converged
        fprintf('Graph label propagation converged at iteration %d\n', iter);
    else
        fprintf('Graph label propagation reached max iterations (%d)\n', p_conf.GraphMaxIter);
    end
    if any(unreached(test_idx))
        warning('CellTypeClassifier:unreached', ['%d units have no graph path to a ' ...
            'training unit and are left unlabelled (NaN).'], sum(unreached(test_idx)));
    end

    % Extract predictions and confidence for test units
    test_global_idx = find(test_idx);
    Y_pred = nan(1, numel(test_global_idx));
    for i = 1:numel(test_global_idx)
        gi = test_global_idx(i);
        if unreached(gi)
            unique_confidence(gi) = NaN;
            continue   % Y_pred(i) stays NaN
        end
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
    % Apply dimension-normalised feature group weighting before kNN search so
    % ACG and Waveform groups contribute proportionally to their configured weights.
    if isfield(np, 'feat_names_trimmed') && ~isempty(np.feat_names_trimmed)
        X_all = applyFeatureGroupWeights(X_all, np.feat_names_trimmed, p_umap);
    end
    X_train = X_all(train_idx, :);
    X_test  = X_all(test_idx, :);

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
        % Per-class seed weights: votes are divided by the class's training count so an
        % unbalanced training set does not tilt every vote (no-op when balanced).
        n_seed   = [sum(labels.sorted_y_train == 1), sum(labels.sorted_y_train == 2)];
        cls_w    = ones(1, 2);
        if all(n_seed > 0), cls_w = mean(n_seed) ./ n_seed; end
        for i = 1:numel(test_global_idx)
            neighbor_labels = labels.sorted_y_train(nn_idx(i, :));
            d       = nn_dists(i, :);
            eps_d   = max(d) * 1e-6;
            if eps_d == 0; eps_d = 1e-10; end
            w       = 1 ./ (d + eps_d);
            w_e     = cls_w(1) * sum(w(neighbor_labels == 1));
            w_i     = cls_w(2) * sum(w(neighbor_labels == 2));
            w_total = w_e + w_i;
            w_exc   = w_e / w_total;
            w_inh   = w_i / w_total;
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
