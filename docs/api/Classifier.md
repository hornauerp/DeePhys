# Classifier

Static classification methods operating on plain feature tables. **No instantiation needed** — all methods are static and have zero coupling to legacy recording objects.

## Example

```matlab
% Get feature matrix and labels from FeatureStore
[X, unit_ids] = fs.unitMatrix('all');
Y = fs.UnitTable.Mutation;

% Basic 5-fold random forest classification
result = Classifier.classify(X, Y, struct('Algorithm', 'rf', 'KFold', 5));

% Hierarchy-aware CV (rows of the same group are always in the same fold)
opts.CVGroups = fs.UnitTable.RecordingID;   % or culture IDs, see Tutorial 3
opts.NormalizationGroups = fs.UnitTable.PlatingDate;  % per-plate z-score
result = Classifier.classify(X, Y, opts);

% Summarise results
summary = ClassificationResult.summarizeFolds(result);
fprintf('Balanced accuracy: %.2f\n', summary.balanced_accuracy);
fprintf('Fold accuracy: %.2f ± %.2f\n', summary.mean_accuracy, summary.std_accuracy);

% Feature importance (RF only)
T = ClassificationResult.summarizeImportance(result);
disp(T(1:10, :));   % top 10 features

% Permutation test (null distribution of balanced accuracy)
[p, null_accs, obs_acc] = Classifier.permutationTest(X, Y, opts, 500);

% Per-feature-group classification: struct array with .Name and .Columns
group_defs = struct('Name',    {'Activity', 'Waveform'}, ...
                    'Columns', {["FiringRate","MeanInterSpikeInterval"], ["T2Pdelay","HalfWidth"]});
results_by_group = Classifier.classifyByGroups(X, Y, group_defs, opts);
results_by_group(1).Result            % ClassificationResult array of the "Activity" group

% Apply a trained fold model to new raw data (pass the result, not result.Mdl)
[Y_pred, scores] = Classifier.predict(result(1), X_new);
```

`Classifier.predict` takes the `ClassificationResult` of one fold, not its raw model: the result carries the fitted preprocessing (imputation, normalisation, column mask), which `predict` replays on `X_new`, and it maps the predicted class indices back to the label values (`result(1).GroupLabels`). Calling `predict(result(1).Mdl, X_new)` skips the normalisation and returns class indices instead of labels.

## Static methods

### `result = Classifier.classify(X, Y, opts)`

K-fold cross-validated classification on a feature table.

| Argument | Type | Description |
|---|---|---|
| `X` | table | (N × F) feature table from `fs.unitMatrix()` etc. |
| `Y` | vector | (N × 1) categorical / string / numeric labels |
| `opts` | struct | Options (all fields optional; see below) |

Returns a `(1 × K) ClassificationResult` array, one object per fold. Rows with a missing label are left out (with a warning).

### `[p, null_accs, obs_acc] = Classifier.permutationTest(X, Y, opts, n_permutations)`

Permutation test for classification significance. Shuffles the labels `n_permutations` times (default 100), runs `classify` on each, and builds a null distribution of **balanced accuracy** (pooled out-of-fold mean per-class recall). Feature importance is not computed in the null runs.

| Return | Description |
|---|---|
| `p` | One-sided p-value with the +1 correction, `(1 + #{null ≥ observed}) / (1 + n_permutations)` (never exactly 0) |
| `null_accs` | (`n_permutations` × 1) column vector of null balanced accuracies |
| `obs_acc` | Observed balanced accuracy of the unpermuted `classify()` run |

With `opts.CVGroups` set and `Y` constant within every group (e.g. genotype per culture), labels are permuted between groups, not between observations; otherwise observations are shuffled within each group. With `opts.Seed` set, every permutation uses a derived seed (`Seed + p`).

### `results = Classifier.classifyByGroups(X, Y, feature_group_defs, opts)`

Runs a separate classification for each named feature group.

| Argument | Description |
|---|---|
| `feature_group_defs` | Struct array; each element has `.Name` (group label) and `.Columns` (string array of exact column names in `X`) |

Returns a struct array with fields `.Name` and `.Result` (the `ClassificationResult` array of that group). A group none of whose columns is found in `X` is skipped with a warning.

```matlab
group_defs = struct('Name',    {'Activity', 'Waveform'}, ...
                    'Columns', {["FiringRate","MeanInterSpikeInterval","CVInterSpikeInterval"], ...
                                ["T2Pdelay","HalfWidth"]});
results = Classifier.classifyByGroups(X, Y, group_defs, opts);
results(1).Name     % "Activity"
results(1).Result   % ClassificationResult array for the activity columns only
```

### `[Y_pred, scores] = Classifier.predict(trained_model, X, norm_groups)`

Apply a trained classifier to new data (no cross-validation).

| Argument | Description |
|---|---|
| `trained_model` | A single-fold `ClassificationResult`, e.g. `result(1)` (recommended). A raw MATLAB classifier is also accepted, but is then applied to `X.Variables` unchanged, so `X` must already be preprocessed like the training data |
| `X` | (N × F) table of **raw** feature values (columns matched by name) |
| `norm_groups` | (N × 1) normalisation groups of `X`; required when the model was trained with `NormalizationGroups` or ComBat |

`Y_pred` is returned in the original label space (`GroupLabels`); `scores` has one column per class of the training data (`NaN` for a class the fold's model never saw).

## Options (`opts` struct fields)

All fields are optional.

### Algorithm

| Field | Default | Description |
|---|---|---|
| `Algorithm` | `"rf"` | `"rf"` (random forest), `"svm"`, `"cnb"` (naive Bayes, `fitcnb`), `"knn"` |
| `KFold` | `5` | Number of CV folds. `-1` = leave-one-group-out when `CVGroups` is set; without `CVGroups`, `-1` means leave-one-out |
| `NHyper` | `0` | Bayesian hyperparameter optimisation evaluations (0 = skip). With `CVGroups`, the inner folds of the search also keep groups together |
| `Prior` | `"empirical"` | `"empirical"` or `"uniform"` (equal class weighting for imbalanced data) |
| `Seed` | `[]` | RNG seed for reproducibility |
| `RF` | `struct()` | Random forest overrides of `MLPipeline.returnDefaultParams().RF`: `NumCycles` (500), `MinLeafSize` (5), `MinParentSize` (`[]` = `templateTree` default 10), `NumVariablesToSample` (`"auto"`), `Surrogate`, `Reproducible`. Unknown fields are an error; the prior is set with `Prior`. Ignored by the hyperparameter search (`NHyper > 0`) |

The RF defaults suit unit-level data. A tree cannot split a node smaller than `MinParentSize` (10) or create a leaf smaller than `MinLeafSize` (5). When the training rows of a fold number about 10 or fewer, as with one row per well in a small experiment, the trees stay single leaves. A fold then predicts one class and its importances are zero. Lower both in that case:

```matlab
opts = struct('Algorithm', 'rf', 'KFold', 5, 'Prior', 'uniform', ...
              'RF', struct('MinLeafSize', 1, 'MinParentSize', 2));
result = Classifier.classify(X_wells, patterning, opts);
```

### Cross-validation grouping

| Field | Default | Description |
|---|---|---|
| `CVGroups` | `[]` | (N × 1) group IDs for hierarchy-aware CV (e.g. `fs.UnitTable.RecordingID`). When set, units from the same group are never split across folds. |

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

### Feature importance (RF only)

| Field | Default | Description |
|---|---|---|
| `ComputeImportance` | `true` | Compute permutation importance per fold. It is the slowest step; switch it off when only accuracy is needed (`permutationTest` does so for the null runs) |
| `KeepFullModel` | `false` | `result(k).Mdl` holds the compact model (it predicts, but has no training data or out-of-bag methods). `true` keeps the full model; full bagged ensembles store a copy of the training rows, which takes gigabytes in large unit-level runs |
| `ImportanceMethod` | `"auto"` | `"auto"`: `"test"` when `CVGroups` is set, otherwise `"oob"`. `"test"`: shuffle one feature at a time in the held-out fold and record the increase in the classification error. `"oob"`: out-of-bag permuted-predictor importance |

With grouped data the out-of-bag rows come from the same groups as the training rows, so `"oob"` favours features that merely identify a group; `"test"` measures what generalises to unseen groups. The test-fold importance uses a single shuffle per feature and fold, so it is noisy for small folds (values can be negative); interpret the mean and spread over folds from `ClassificationResult.summarizeImportance`, not a single fold. The two methods are on different scales and not comparable.

## Notes

- Requires the **Statistics and Machine Learning Toolbox** (`fitcensemble`, `fitcsvm`/`fitcecoc`, `fitcnb`, `fitcknn`).
- Predictor importance is available for `"rf"` only; `result.predImp` is `[]` for other algorithms or with `ComputeImportance = false`.
- `result(k).train_acc` uses OOB error for `"rf"` and resubstitution loss for `"svm"`/`"cnb"`/`"knn"`. It is optimistic whenever rows of a group are correlated; report the held-out `balanced_accuracy` from `ClassificationResult.summarizeFolds`.
- `CVGroups` is passed to `GroupedCV.byGroups()` internally. See [GroupedCV](GroupedCV.md) for details.
- NaN values in `X` are imputed with the training-set column median before each fold (via `MLUtils`).
- An `opts` field that no part of the ML layer reads (e.g. a misspelt `Kfold`) gives a `Classifier:unknownOption` warning instead of being ignored silently.
