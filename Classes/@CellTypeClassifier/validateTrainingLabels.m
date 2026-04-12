function results = validateTrainingLabels(ctc)
% VALIDATETRAININGLABELS  Leave-one-culture-out cross-validation on training labels.
%
% For each culture that contributes training units, trains a supervised UMAP
% on units from all other cultures and classifies the held-out culture's
% training units. Reports per-culture accuracy and overall agreement between
% predicted and ground-truth training labels.
%
% This gives an empirical estimate of label reliability without access to
% external ground truth. Low accuracy in a held-out culture indicates that
% the responsive units from that culture differ systematically from the
% others — a sign that the responsive-unit labels may be noisy or that
% the culture has an unusual cell-type composition.
%
% Requires generateTrainLabels() to have been run first.
%
% OUTPUTS:
%   results - struct with fields:
%     .overall_accuracy   - fraction of training units correctly classified
%     .per_culture        - table with columns: CultureID, N, Accuracy, ConfusionMatrix
%     .predicted_labels   - (1 x N_train) predicted labels in training order
%     .true_labels        - (1 x N_train) true labels from TrainLabels
%     .confusion          - 2x2 confusion matrix (rows=true, cols=predicted)

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before validateTrainingLabels()');
assert(~isempty(ctc.NormalizationParams), ...
    'NormalizationParams missing — re-run generateTrainLabels().');

rng(ctc.Parameters.RNGSeed, 'twister');

labels  = ctc.TrainLabels;
p_umap  = ctc.Parameters.UMAP;
np      = ctc.NormalizationParams;

% -- Rebuild unique unit space + full feature matrix --------------------------
[unique_ud, all_to_unique, ~] = ctc.uniqueUnitMap();
n_unique = numel(unique_ud);

unique_to_rep = zeros(1, n_unique);
for u = 1:n_unique
    rows = find(all_to_unique == u);
    unique_to_rep(u) = rows(1);
end

[wf, acg, sr]  = ctc.getOrExtract(unique_ud);
[X_raw, ~]     = buildFeatureMatrix(ctc, wf, acg, sr);
X_raw(isnan(X_raw)) = 0;

% Per-group z-score (same as generateTrainLabels / classifyUnits)
% Map unique units to FeatureStore rows via UnitID+RecordingID (order-independent)
ud = ctc.UnitDataArray;
rep_ud_keys_v = string({ud(unique_to_rep).UnitID}) + "|" + string({ud(unique_to_rep).RecordingID});
fs_keys_v     = string(ctc.FeatureStore.UnitTable.UnitID(:))' + "|" + ...
                string(ctc.FeatureStore.UnitTable.RecordingID(:))';
[~, unique_to_fs_v] = ismember(rep_ud_keys_v, fs_keys_v);

norm_var = p_umap.NormalizationVar;
if ~isempty(norm_var) && ismember(norm_var, string(ctc.FeatureStore.UnitTable.Properties.VariableNames))
    group_col        = ctc.FeatureStore.UnitTable.(norm_var);
    group_col_unique = group_col(unique_to_fs_v);
    [~, ~, iG] = unique(string(group_col_unique), 'stable');
    for g = 1:max(iG)
        mask = iG == g;
        X_raw(mask, :) = normalize(X_raw(mask, :));
    end
end

% Apply stored global normalization from generateTrainLabels
X_all = normalize(X_raw, 'center', np.mu_global, 'scale', np.sigma_global);
X_all(:, np.nan_cols)  = [];
X_all = X_all ./ np.scale;
if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X_all = X_all(:, np.feature_selection_mask);
end

% -- Identify training units and their cultures ------------------------------
train_idx   = logical(labels.umap_train_idx);   % in unique-unit space
train_ids   = labels.sorted_train_ids;           % indices into unique-unit space
y_train     = labels.sorted_y_train;
X_train_all = X_all(train_idx, :);

% Get culture ID for each training unit via FeatureStore metadata
meta         = ctc.FeatureStore.MetadataTable;
unit_rec_ids = ctc.FeatureStore.UnitTable.RecordingID;
culture_ids  = FeatureStore.getCultureIDsForUnits( ...
    unit_rec_ids, meta, ctc.Parameters.CultureKeys);

% Map training units to their cultures via UnitID matching (order-independent).
% culture_ids is in FeatureStore order; unique_to_rep is in UnitDataArray order.
unit_ids_table = string(ctc.FeatureStore.UnitTable.UnitID);
train_ud_ids   = string({unique_ud(train_ids).UnitID});
[~, table_loc] = ismember(train_ud_ids, unit_ids_table);
train_culture_ids = culture_ids(table_loc);
unique_cultures   = unique(train_culture_ids, 'stable');
n_cultures        = numel(unique_cultures);

if n_cultures < 2
    warning('CellTypeClassifier:validateTrainingLabels', ...
        'Only %d culture(s) in training set — leave-one-out CV requires at least 2. Skipping.', ...
        n_cultures);
    results = struct('overall_accuracy', NaN, 'per_culture', table(), ...
        'predicted_labels', nan(size(y_train)), 'true_labels', y_train, ...
        'confusion', nan(2));
    return
end

% -- AutoNNeighbors for supervised UMAP (matches classifyUnits logic) ---------
if p_umap.AutoNNeighbors
    cv_n_neighbors = max(p_umap.MinNNeighbors, round(sqrt(numel(train_ids))));
else
    cv_n_neighbors = p_umap.SupervisedNNeighbors;
end

cv_target_weight = p_umap.TargetWeight;

% -- Leave-one-culture-out loop -----------------------------------------------
predicted_labels = nan(1, numel(train_ids));
per_culture_acc  = nan(1, n_cultures);
per_culture_n    = zeros(1, n_cultures);
per_culture_cids = unique_cultures;

tdir = p_umap.TemplateDir;
assert(~isempty(tdir) && isfolder(tdir), ...
    'CellTypeClassifier:validateTrainingLabels', ...
    'Parameters.UMAP.TemplateDir must be set and exist for cross-validation.');

for ci = 1:n_cultures
    cid       = unique_cultures(ci);
    held_mask = train_culture_ids == cid;
    src_mask  = ~held_mask;

    if sum(src_mask) < 4 || sum(held_mask) < 1
        fprintf('Culture %s: too few source or held-out units — skipping.\n', cid);
        continue
    end

    % Check both classes represented in source set
    src_labels = y_train(src_mask);
    if numel(unique(src_labels)) < 2
        fprintf('Culture %s: source set lacks both classes — skipping.\n', cid);
        continue
    end

    X_src  = X_train_all(src_mask, :);
    y_src  = y_train(src_mask);
    X_held = X_train_all(held_mask, :);
    y_held = y_train(held_mask);

    % Brief supervised UMAP: source as training, held-out as test
    uid           = char(java.util.UUID.randomUUID);
    template_file = fullfile(tdir, ['ctc_cv_template_', uid, '.mat']);
    cleanup       = onCleanup(@() deleteIfExists(template_file));

    try
        train_args = { ...
            'label_column',       'end', ...
            'n_components',       p_umap.SupervisedNDims, ...
            'n_neighbors',        cv_n_neighbors, ...
            'min_dist',           p_umap.MinDist, ...
            'spread',             p_umap.Spread, ...
            'target_weight',      cv_target_weight, ...
            'method',            'java', ...
            'verbose',            'none', ...
            'save_template_file', template_file};
        run_umap([X_src, y_src'], train_args{:});

        X_held_with_label = [X_held, zeros(size(X_held, 1), 1)];
        test_args = { ...
            'label_column',      'end', ...
            'method',            'java', ...
            'verbose',           'none', ...
            'cluster_detail',    'adaptive', ...
            'match_supervisors', 1, ...
            'match_scenarios',   4, ...
            'template_file',     template_file};
        [~, ~, ~, extras_cv] = run_umap(X_held_with_label, test_args{:});

        y_pred_held = extras_cv.supervisorMatchedLabels;
        predicted_labels(held_mask) = y_pred_held;

        acc = mean(y_pred_held(:) == y_held(:));
        per_culture_acc(ci) = acc;
        per_culture_n(ci)   = sum(held_mask);

        fprintf('Culture %s: %d units, accuracy = %.1f%%\n', cid, sum(held_mask), 100 * acc);
    catch ME
        warning('CellTypeClassifier:validateTrainingLabels', ...
            'Culture %s CV failed: %s', cid, ME.message);
    end
end

% -- Aggregate results ---------------------------------------------------------
valid = ~isnan(per_culture_acc);
if any(valid)
    overall_accuracy = sum(per_culture_acc(valid) .* per_culture_n(valid)) / sum(per_culture_n(valid));
else
    overall_accuracy = NaN;
end

% 2x2 confusion matrix (true class 1 or 2 vs predicted class 1 or 2)
class_labels = [1, 2];
confusion    = zeros(2, 2);
valid_pred   = ~isnan(predicted_labels);
for ri = 1:2
    for ci = 1:2
        confusion(ri, ci) = sum( ...
            y_train(valid_pred) == class_labels(ri) & ...
            predicted_labels(valid_pred) == class_labels(ci));
    end
end

per_culture_table = table( ...
    per_culture_cids(:), per_culture_n(:), per_culture_acc(:), ...
    'VariableNames', {'CultureID', 'N', 'Accuracy'});

results = struct( ...
    'overall_accuracy', overall_accuracy, ...
    'per_culture',      per_culture_table, ...
    'predicted_labels', predicted_labels, ...
    'true_labels',      y_train, ...
    'confusion',        confusion);

fprintf('\nLeave-one-culture-out CV: overall accuracy = %.1f%% (%d/%d cultures evaluated)\n', ...
    100 * overall_accuracy, sum(valid), n_cultures);
fprintf('Confusion matrix (rows=true, cols=predicted):\n');
fprintf('              Pred exc  Pred inh\n');
fprintf('  True exc:   %6d    %6d\n', confusion(1, 1), confusion(1, 2));
fprintf('  True inh:   %6d    %6d\n', confusion(2, 1), confusion(2, 2));
end

function deleteIfExists(f)
    if isfile(f); delete(f); end
end
