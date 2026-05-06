# Classifier

Static classification methods operating on plain feature tables. **No instantiation needed** — all methods are static and have zero coupling to legacy recording objects.

## Example

```matlab
% Get feature matrix and labels from FeatureStore
[X, unit_ids] = fs.unitMatrix('all');
Y = fs.UnitTable.Mutation;

% Basic 5-fold random forest classification
result = Classifier.classify(X, Y, struct('Algorithm', 'rf', 'KFold', 5));

% Hierarchy-aware CV (units from same recording always in same fold)
opts.CVGroups = fs.UnitTable.RecordingID;
opts.NormalizationGroups = fs.UnitTable.PlatingDate;  % per-plate z-score
result = Classifier.classify(X, Y, opts);

% Summarise results
summary = ClassificationResult.summarizeFolds(result);
fprintf('Accuracy: %.2f ± %.2f\n', summary.mean_accuracy, summary.std_accuracy);

% Feature importance (RF only)
T = ClassificationResult.summarizeImportance(result);
disp(T(1:10, :));   % top 10 features

% Permutation test (null distribution of accuracy)
[p, null_accs, obs_acc] = Classifier.permutationTest(X, Y, opts, 500);

% Per-feature-group classification
group_defs = struct('Activity', ["FiringRate","MeanISI"], 'Waveform', ["T2Pdelay","HalfWidth"]);
results_by_group = Classifier.classifyByGroups(X, Y, group_defs, opts);

% Apply a trained model to new data
[Y_pred, scores] = Classifier.predict(result(1).Mdl, X_new);
```

## Static methods

### `result = Classifier.classify(X, Y, opts)`

K-fold cross-validated classification on a feature table.

| Argument | Type | Description |
|---|---|---|
| `X` | table | (N × F) feature table from `fs.unitMatrix()` etc. |
| `Y` | vector | (N × 1) categorical / string / numeric labels |
| `opts` | struct | Options (all fields optional; see below) |

Returns a `(1 × K) ClassificationResult` array, one object per fold.

### `[p, null_accs, obs_acc] = Classifier.permutationTest(X, Y, opts, n_permutations)`

Permutation test for classification significance. Shuffles Y labels `n_permutations` times and builds a null accuracy distribution.

| Return | Description |
|---|---|
| `p` | p-value: fraction of null accuracies ≥ observed |
| `null_accs` | (1 × n_permutations) null accuracy distribution |
| `obs_acc` | Observed mean accuracy |

### `results = Classifier.classifyByGroups(X, Y, feature_group_defs, opts)`

Runs separate classification for each named feature group. Returns a struct where each field is a `ClassificationResult` array.

```matlab
group_defs = struct('Activity', ["FiringRate","MeanISI","CVISI"], ...
                    'Waveform', ["T2Pdelay","HalfWidth","PeakAsymmetry"]);
results = Classifier.classifyByGroups(X, Y, group_defs, opts);
% results.Activity  — ClassificationResult array for activity features only
% results.Waveform  — ClassificationResult array for waveform features only
```

### `[Y_pred, scores] = Classifier.predict(trained_model, X)`

Apply a trained classifier to new data (no cross-validation).

## Options (`opts` struct fields)

All fields are optional.

### Algorithm

| Field | Default | Description |
|---|---|---|
| `Algorithm` | `"rf"` | `"rf"` (random forest), `"svm"`, `"cnb"` (complement naive Bayes), `"knn"` |
| `KFold` | `5` | Number of CV folds; `-1` = leave-one-group-out |
| `NHyper` | `0` | Bayesian hyperparameter optimisation evaluations (0 = skip) |
| `Prior` | `"empirical"` | `"empirical"` or `"uniform"` (equal class weighting for imbalanced data) |
| `Seed` | `[]` | RNG seed for reproducibility |

### Cross-validation grouping

| Field | Default | Description |
|---|---|---|
| `CVGroups` | `[]` | (N × 1) group IDs for hierarchy-aware CV (e.g. `fs.UnitTable.RecordingID`). When set, units from the same group are never split across folds. |
| `CVLevel` | `"recording"` | `"recording"` or `"culture"` — context label for serialization |

### Normalization

| Field | Default | Description |
|---|---|---|
| `NormalizationGroups` | `[]` | (N × 1) group IDs for per-group z-score (e.g. `fs.UnitTable.PlatingDate`). `[]` = global z-score. |
| `NormalizationPipeline` | `""` | `""` (z-score) or `"combat"` (ComBat batch correction) |
| `CovariateLabels` | `[]` | Biological covariates to protect from ComBat removal (defaults to `Y` when `NormalizationPipeline = "combat"`) |

### Label pooling

| Field | Default | Description |
|---|---|---|
| `PoolingValues` | `{}` | Cell array to merge label values, e.g. `{{'HET','HOM'}}` merges both into one class |

## Notes

- Requires the **Statistics and Machine Learning Toolbox** (`fitcensemble`, `fitcsvm`, `fitcnb`, `fitcknn`).
- OOB (out-of-bag) predictor importance is available for `"rf"` only; `result.predImp` is `[]` for other algorithms.
- Training accuracy uses OOB error for `"rf"` and resubstitution loss for `"svm"`/`"cnb"`/`"knn"`.
- `CVGroups` is passed to `GroupedCV.byGroups()` internally. See [GroupedCV](GroupedCV.md) for details.
- NaN values in `X` are imputed with the training-set column median before each fold (via `MLUtils`).
