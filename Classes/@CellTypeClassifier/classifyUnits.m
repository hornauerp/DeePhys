function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Feature extraction and per-group z-score are shared with generateTrainLabels
% via the NormalizedFeatures cache (buildNormalizedFeatures). This function
% applies the NormalizationParams stored by generateTrainLabels to all units,
% ensuring an identical feature space across both UMAP runs:
%   1. Per-group z-score (done in buildNormalizedFeatures, shared)
%   2. Apply stored global z-score (mu, sigma from generateTrainLabels subset)
%   3. Remove same NaN columns as generateTrainLabels
%   4. Apply same max(abs) scaling
%   5. Apply feature_selection_mask if FeatureSelection was enabled
%
% Requires ctc.TrainLabels and ctc.NormalizationParams to be set
% (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)
%       ctc.UnitConfidence  distance-weighted kNN vote fraction in [0,1];
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
np     = ctc.NormalizationParams;

% -- Get normalized features (shared with generateTrainLabels) ----------------
% buildNormalizedFeatures is idempotent: returns immediately if already run.
% The cache holds X_pergroup (after per-group z-score) + deduplication info.
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;

n_unique      = nf.n_unique;
n_units_full  = nf.n_units_full;  %#ok<NASGU>
all_to_unique = nf.all_to_unique;
unique_to_rep = nf.unique_to_rep; %#ok<NASGU>

% -- Apply stored NormalizationParams to ALL unique units ---------------------
% The global z-score, NaN removal, and scaling were fit on the training subset
% in generateTrainLabels. Applying the same parameters here (rather than
% refitting) ensures the supervised UMAP sees the same feature space as the
% unsupervised UMAP used for label generation.
X_all = normalize(nf.X_pergroup, 'center', np.mu_global, 'scale', np.sigma_global);
X_all(:, np.nan_cols)  = [];
X_all = X_all ./ np.scale;

if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X_all = X_all(:, np.feature_selection_mask);
end

train_idx = logical(labels.umap_train_idx);
test_idx  = ~train_idx;
X_train   = X_all(train_idx, :);
X_test    = X_all(test_idx, :);

% -- AutoNNeighbors for supervised UMAP ---------------------------------------
orig_sup_nneigh    = p_umap.SupervisedNNeighbors;
orig_target_weight = p_umap.TargetWeight;

if p_umap.AutoNNeighbors
    N_train    = size(X_train, 1);
    auto_sup_n = max(p_umap.MinNNeighbors, round(sqrt(N_train)));
    fprintf('Auto SupervisedNNeighbors: %d (N_train=%d)\n', auto_sup_n, N_train);
    ctc.Parameters.UMAP.SupervisedNNeighbors = auto_sup_n;
end

% Metadata labels are clean ground truth — higher TargetWeight is appropriate.
% Inferred drug-response labels use the default (0.5) because label noise
% warrants a more data-driven embedding.
if isfield(labels, 'has_explicit_ce') && labels.has_explicit_ce
    ctc.Parameters.UMAP.TargetWeight = p_umap.MetadataTargetWeight;
    fprintf('Metadata labels: using TargetWeight=%.2f\n', p_umap.MetadataTargetWeight);
end

% -- Supervised UMAP classification ------------------------------------------
[Y_pred, ~, ~, test_reduction, train_reduction, extras] = supervisedUMAP(ctc, ...
    X_train, labels.sorted_y_train, nf.feat_names, X_test);

% Restore overridden params
ctc.Parameters.UMAP.SupervisedNNeighbors = orig_sup_nneigh;
ctc.Parameters.UMAP.TargetWeight         = orig_target_weight;

ctc.Reduction.Train  = train_reduction;
ctc.Reduction.Test   = test_reduction;
ctc.Reduction.Extras = extras;

% -- Distance-weighted kNN confidence in UMAP embedding space -----------------
% For each test unit: sum of inverse-distance weights from k nearest training
% neighbors sharing the predicted label, divided by total weight. Closer
% neighbors contribute more, giving a softer signal near the decision boundary.
% Training units get confidence = 1.0 (ground truth).
if p_umap.AutoConfidenceK
    conf_k = max(5, round(sqrt(size(train_reduction, 1))));
    fprintf('Auto ConfidenceK: %d (N_train=%d)\n', conf_k, size(train_reduction, 1));
else
    conf_k = p_umap.ConfidenceK;
end

unique_confidence = nan(1, n_unique);
unique_confidence(labels.sorted_train_ids) = 1.0;

test_global_idx = find(test_idx);
if ~isempty(test_reduction) && ~isempty(train_reduction)
    k_actual = min(conf_k, size(train_reduction, 1));
    [nn_idx, nn_dists] = knnsearch(train_reduction, test_reduction, 'K', k_actual);
    for i = 1:numel(test_global_idx)
        neighbor_labels = labels.sorted_y_train(nn_idx(i, :));
        d      = nn_dists(i, :);
        eps_d  = max(d) * 1e-6;
        if eps_d == 0; eps_d = 1e-10; end
        w      = 1 ./ (d + eps_d);
        w_total = sum(w);
        w_match = sum(w(neighbor_labels == Y_pred(i)));
        unique_confidence(test_global_idx(i)) = w_match / w_total;
    end
end

% -- Assemble labels in unique-unit space -------------------------------------
unique_labels = nan(1, n_unique);
unique_labels(labels.sorted_train_ids) = labels.sorted_y_train;
unique_labels(test_idx)                = Y_pred;

% Optional confidence threshold
p_conf = ctc.Parameters.Classification;
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
n_inh = sum(unique_labels == 2, 'omitnan');
n_unc = sum(isnan(unique_labels));
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, n_unc);
end
