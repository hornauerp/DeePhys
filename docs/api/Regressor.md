# Regressor

Static regression methods operating on plain feature tables. **No instantiation needed** — all methods are static.

## Example

```matlab
% Feature matrix and numeric targets
[X, rec_ids] = fs.recordingMatrix('all');
Y = fs.RecordingTable.DIV;

% 5-fold random forest regression with grouped CV
opts.Algorithm   = 'rf';
opts.KFold       = 5;
opts.CVGroups    = fs.RecordingTable.ChipID;  % never split same chip across folds
opts.NormalizationGroups = fs.RecordingTable.PlatingDate;
result = Regressor.regress(X, Y, opts);

% Summarise: pooled out-of-fold R² is the headline number
summary = RegressionResult.summarizeFolds(result);
fprintf('R² (pooled out-of-fold): %.3f\n', summary.R2_pooled);
fprintf('Fold R²: %.3f ± %.3f\n', summary.mean_R2, std(summary.per_fold.R2, 'omitnan'));

% Feature importance
T = RegressionResult.summarizeImportance(result);
disp(T(1:10, :));

% Permutation test
[p, null_r2s, obs_r2] = Regressor.permutationTest(X, Y, opts, 500);

% Per-feature-group regression: struct array with .Name and .Columns
Xa = fs.recordingMatrix("ActivityFeatures");
Xw = fs.recordingMatrix("WaveformFeatures");
group_defs = struct('Name',    {'Activity', 'Waveform'}, ...
                    'Columns', {string(Xa.Properties.VariableNames), string(Xw.Properties.VariableNames)});
results_by_group = Regressor.regressByGroups(X, Y, group_defs, opts);
results_by_group(1).Result            % RegressionResult array of the "Activity" group
```

## Static methods

### `result = Regressor.regress(X, Y, opts)`

K-fold cross-validated regression on a feature table.

| Argument | Type | Description |
|---|---|---|
| `X` | table | (N × F) feature table from `fs.recordingMatrix()` etc. |
| `Y` | numeric | (N × 1) numeric target vector |
| `opts` | struct | Options (all fields optional; see below) |

Returns a `(1 × K) RegressionResult` array, one object per fold. Rows with a NaN target are left out.

### `[p, null_r2s, obs_r2] = Regressor.permutationTest(X, Y, opts, n_permutations)`

Permutation test for regression significance. Shuffles Y `n_permutations` times (default 100) to build a null R² distribution. With `opts.CVGroups`, the shuffle respects the groups (whole groups when `Y` is constant within them, otherwise within groups).

| Return | Description |
|---|---|
| `p` | One-sided p-value with the +1 correction, `(1 + #{null R² ≥ observed}) / (1 + n_permutations)`; null runs with an undefined R² count as reaching the observed value |
| `null_r2s` | (`n_permutations` × 1) column vector of null R² values |
| `obs_r2` | Observed **pooled out-of-fold R²** (`R2_pooled` of `summarizeFolds`) of the unpermuted run |

The test errors (`Regressor:undefinedR2`) when the observed R² is undefined, i.e. the target is constant.

### `results = Regressor.regressByGroups(X, Y, feature_group_defs, opts)`

Runs a separate regression for each named feature group.

| Argument | Description |
|---|---|
| `feature_group_defs` | Struct array; each element has `.Name` (group label) and `.Columns` (string array of exact column names in `X`) |

Returns a struct array with fields `.Name` and `.Result` (the `RegressionResult` array of that group). A group none of whose columns is found in `X` is skipped with a warning.

## Options (`opts` struct fields)

All fields are optional.

| Field | Default | Description |
|---|---|---|
| `Algorithm` | `"rf"` | `"rf"` (random forest), `"svm"`, `"knn"` |
| `KFold` | `5` | Number of CV folds. `-1` = leave-one-group-out, which **requires `CVGroups`**; without it `regress` errors with `Regressor:kfoldNeedsGroups` |
| `NHyper` | `0` | Bayesian hyperparameter optimisation evaluations (0 = skip). With `CVGroups`, the inner folds of the search also keep groups together |
| `Seed` | `[]` | RNG seed for reproducibility |
| `CVGroups` | `[]` | (N × 1) group IDs for hierarchy-aware CV (e.g. `fs.RecordingTable.ChipID`) |
| `StratificationVar` | `[]` | (N × 1) for stratified CV when `CVGroups` not set |
| `NormalizationGroups` | `[]` | (N × 1) group IDs for per-group z-score |
| `NormalizationPipeline` | `""` | `""` (z-score) or `"combat"` (ComBat batch correction) |
| `CovariateLabels` | `[]` | Biological covariates to protect from ComBat removal |
| `ComputeImportance` | `true` | RF only: compute permutation importance per fold. It dominates RF run time; `permutationTest` switches it off for the null runs |
| `ImportanceMethod` | `"auto"` | RF only. `"auto"`: `"test"` when `CVGroups` is set, otherwise `"oob"`. `"test"`: shuffle one feature at a time in the held-out fold and record the increase in test MSE. `"oob"`: out-of-bag permuted-predictor importance |

With grouped data the out-of-bag rows come from the same groups as the training rows, so `"oob"` favours features that merely identify a group; `"test"` measures what generalises to unseen groups. The test-fold importance uses a single shuffle per feature and fold, so it is noisy for small folds (values can be negative); read the mean and spread over folds from `RegressionResult.summarizeImportance`. The two methods are on different scales and not comparable.

## Notes

- Requires the **Statistics and Machine Learning Toolbox** (`fitrensemble`, `fitrsvm`).
- MATLAB has no `fitrknn`: `"knn"` uses the `KNNRegressionModel` class, an inverse-distance-weighted k-nearest-neighbour regressor (K = 5 by default; with `NHyper > 0`, K is chosen by leave-one-out error, leave-one-group-out when `CVGroups` is set).
- Predictor importance is available for `"rf"` only.
- Report the pooled out-of-fold `R2_pooled` from `RegressionResult.summarizeFolds`. Per-fold R² compares a fold with its own mean and is undefined or unstable when a grouped fold holds one target value or few rows.
- `result(k).train_acc` is the OOB R² for `"rf"` and the resubstitution R² for `"svm"`/`"knn"`; it is optimistic when rows of a group are correlated and is not a test estimate.
- `CVGroups` is passed to `GroupedCV.byGroups()` internally. See [GroupedCV](GroupedCV.md).
- Note: `"cnb"` (naive Bayes) is not available for regression — use [Classifier](Classifier.md) for classification tasks.
