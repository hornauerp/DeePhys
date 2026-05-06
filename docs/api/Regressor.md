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

% Summarise
summary = RegressionResult.summarizeFolds(result);
fprintf('R²: %.3f ± %.3f\n', summary.mean_R2, std([result.computeMetrics().R2]));

% Feature importance
T = RegressionResult.summarizeImportance(result);
disp(T(1:10, :));

% Permutation test
[p, null_r2s, obs_r2] = Regressor.permutationTest(X, Y, opts, 500);

% Per-feature-group regression
group_defs = struct('Activity', ["FiringRate","MeanISI"], 'Network', ["Synchrony","Clustering"]);
results_by_group = Regressor.regressByGroups(X, Y, group_defs, opts);
```

## Static methods

### `result = Regressor.regress(X, Y, opts)`

K-fold cross-validated regression on a feature table.

| Argument | Type | Description |
|---|---|---|
| `X` | table | (N × F) feature table from `fs.recordingMatrix()` etc. |
| `Y` | numeric | (N × 1) numeric target vector |
| `opts` | struct | Options (all fields optional; see below) |

Returns a `(1 × K) RegressionResult` array, one object per fold.

### `[p, null_r2s, obs_r2] = Regressor.permutationTest(X, Y, opts, n_permutations)`

Permutation test for regression significance. Shuffles Y `n_permutations` times to build a null R² distribution.

| Return | Description |
|---|---|
| `p` | p-value: fraction of null R² values ≥ observed |
| `null_r2s` | (1 × n_permutations) null R² distribution |
| `obs_r2` | Observed mean R² |

### `results = Regressor.regressByGroups(X, Y, feature_group_defs, opts)`

Runs separate regression for each named feature group. Returns a struct where each field is a `RegressionResult` array.

## Options (`opts` struct fields)

All fields are optional.

| Field | Default | Description |
|---|---|---|
| `Algorithm` | `"rf"` | `"rf"` (random forest), `"svm"`, `"knn"` |
| `KFold` | `5` | Number of CV folds; `-1` = leave-one-group-out |
| `NHyper` | `0` | Bayesian hyperparameter optimisation evaluations (0 = skip) |
| `Seed` | `[]` | RNG seed for reproducibility |
| `CVGroups` | `[]` | (N × 1) group IDs for hierarchy-aware CV (e.g. `fs.RecordingTable.ChipID`) |
| `StratificationVar` | `[]` | (N × 1) for stratified CV when `CVGroups` not set |
| `NormalizationGroups` | `[]` | (N × 1) group IDs for per-group z-score |
| `NormalizationPipeline` | `""` | `""` (z-score) or `"combat"` (ComBat batch correction) |
| `CovariateLabels` | `[]` | Biological covariates to protect from ComBat removal |

## Notes

- Requires the **Statistics and Machine Learning Toolbox** (`fitrensemble`, `fitrsvm`, `fitrknn`).
- OOB predictor importance is available for `"rf"` only.
- Training accuracy uses OOB R² for `"rf"` and resubstitution R² for `"svm"`/`"knn"`.
- `CVGroups` is passed to `GroupedCV.byGroups()` internally. See [GroupedCV](GroupedCV.md).
- Note: `"cnb"` (complement naive Bayes) is not available for regression — use [Classifier](Classifier.md) for classification tasks.
