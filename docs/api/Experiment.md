# Experiment

Thin orchestration layer connecting `FeatureStore` to the `Classifier`, `Regressor`, and `DimReducer` analysis classes. **Handle class.** Optional — users can call the analysis classes directly with feature tables from `FeatureStore`.

`Experiment` assembles the right feature matrices, label vectors, and CV groupings from the `FeatureStore`, then delegates to the standalone analysis classes and stores results.

## Example

```matlab
% Build from saved processors
procs = RecordingProcessor.loadMany(proc_paths);
exp   = Experiment.fromProcessors(procs);

% Or load processors and build in one step
exp = Experiment.fromPaths(proc_paths);

% Classification (results stored in exp.Results.Classification)
opts = struct('Algorithm', 'rf', 'KFold', 5, 'CVLevel', 'recording');
result = exp.classify('Unit', 'Mutation', opts);

% Regression
exp.regress('Recording', 'DIV', struct('Algorithm', 'rf'));

% Dimensionality reduction
exp.reduce('Unit', struct('Method', 'UMAP'));

% Access results
exp.Results.Classification.Mutation     % (1 x K) ClassificationResult array
exp.Results.Regression.DIV              % (1 x K) RegressionResult array
exp.Results.DimReduction.Unit.UMAP      % DimReductionResult

% Access FeatureStore directly
[X, ids] = exp.FeatureStore.unitMatrix('all');
```

## Static constructors

| Method | Description |
|---|---|
| `Experiment.fromProcessors(proc_array)` | Build from a `RecordingProcessor` array |
| `Experiment.fromProcessors(proc_array, parameters)` | Same, with parameter overrides |
| `Experiment.fromPaths(processor_paths)` | Load processors from file paths and build |
| `Experiment.fromLegacyGroup(rg)` | Migrate from an old `RecordingGroup` object |
| `Experiment.returnDefaultParams()` | Return the default parameter struct |

There is no factory that takes a saved `FeatureStore` alone. To analyse a saved store without the processors, load it with `FeatureStore.load(path)` and call `Classifier` / `Regressor` / `DimReducer` directly on `fs.unitMatrix(...)` etc. (see [Classifier](Classifier.md)).

## Analysis methods

### `result = exp.classify(level, classification_var, opts)`

Cross-validated classification of a metadata variable.

| Argument | Default | Description |
|---|---|---|
| `level` | `"Recording"` | `"Unit"`, `"Recording"`, or `"Culture"` |
| `classification_var` | `"Mutation"` | Name of a metadata column to classify |
| `opts` | `struct()` | Options struct (see below) |

Returns a `(1 x K) ClassificationResult` array (one per CV fold). Also stored in `exp.Results.Classification.(classification_var)`.

### `result = exp.regress(level, regression_var, opts)`

Cross-validated regression of a numeric metadata variable.

| Argument | Default | Description |
|---|---|---|
| `level` | `"Recording"` | `"Unit"`, `"Recording"`, or `"Culture"` |
| `regression_var` | `"DIV"` | Name of a numeric metadata column |
| `opts` | `struct()` | Options struct (see below) |

Returns a `(1 x K) RegressionResult` array. Also stored in `exp.Results.Regression.(regression_var)`.

### `result = exp.reduce(level, opts)`

Dimensionality reduction at unit, recording, or culture level.

| Argument | Default | Description |
|---|---|---|
| `level` | `"Unit"` | `"Unit"`, `"Recording"`, or `"Culture"` |
| `opts` | `struct()` | Options struct (see below) |

Returns a `DimReductionResult`. Also stored in `exp.Results.DimReduction.(level).(method)`.

## Options (`opts`)

All fields are optional. Unspecified fields use the defaults shown.

### Feature selection

| Field | Default | Description |
|---|---|---|
| `FeatureGroups` | `"all"` | Feature groups to include (string array or `"all"`) |
| `FeatureSet` | `"full"` | `"full"` or `"core"` (curated non-redundant subset) |
| `ParentFeatures` | `[]` | Parent feature groups to prefer (e.g. `"ACG"`, `"all"`) |
| `MinFiringRate` | `0` | Exclude units with firing rate below this (Hz) |

### ML algorithm

| Field | Default | Description |
|---|---|---|
| `Algorithm` | `"rf"` | `"rf"`, `"svm"`, `"cnb"`, `"knn"` |
| `KFold` | `5` | Number of CV folds; `-1` = leave-one-group-out (uses the CV groups of the `"Unit"`/`"Recording"` levels, see `CVLevel`). At the `"Culture"` level there are no CV groups: `classify` then does leave-one-out, and `regress` errors (`Regressor:kfoldNeedsGroups`) |
| `NHyper` | `0` | Bayesian hyperparameter optimisation evaluations (0 = skip) |
| `Prior` | `"empirical"` | `"empirical"` or `"uniform"` (balanced class weighting) |
| `Seed` | `[]` | RNG seed for reproducible CV splits |
| `RF` | `struct()` | `classify` only: random forest overrides, e.g. `struct('MinLeafSize', 1, 'MinParentSize', 2)` at the `"Recording"` / `"Culture"` level with few rows per class. See [Classifier](Classifier.md) |
| `ComputeImportance` | `true` | RF only: compute permutation importance per fold (forwarded to `Classifier` / `Regressor`) |
| `ImportanceMethod` | `"auto"` | RF only: `"auto"` (held-out-fold importance when CV groups are set, i.e. at the `"Unit"`/`"Recording"` levels; out-of-bag at the `"Culture"` level), `"test"` or `"oob"`. See [Classifier](Classifier.md) |

### Cross-validation grouping

| Field | Default | Description |
|---|---|---|
| `CVLevel` | `"culture"` | `"recording"` or `"culture"` — keeps all units from the same group in the same fold |

### Normalization

| Field | Default | Description |
|---|---|---|
| `NormalizationVar` | `""` | Metadata field for per-group z-score (e.g. `"PlatingDate"`) |
| `NormalizationPipeline` | `""` | `""` (z-score) or `"combat"` (ComBat batch correction) |
| `CovariateLabels` | `[]` | Covariate for ComBat (defaults to `Y` when pipeline is `"combat"`) |

### Culture-level aggregation

| Field | Default | Description |
|---|---|---|
| `IdentityKeys` | `["ChipID","PlatingDate"]` | Metadata columns that define a unique culture |
| `GroupingVar` | `"DIV"` | Metadata column for the time axis |
| `GroupingValues` | `[7, 14, 21, 28]` | Required time points per culture |
| `Normalization` | `""` | `""`, `"baseline"`, or `"scaled"` for culture feature normalization |

### Dimensionality reduction

| Field | Default | Description |
|---|---|---|
| `Method` | `"UMAP"` | `"UMAP"`, `"PCA"`, or `"tSNE"` |

## Properties

| Name | Type | Description |
|---|---|---|
| `FeatureStore` | FeatureStore | All features and metadata as tables |
| `Processors` | RecordingProcessor | (1 x N) kept for spike time access |
| `Parameters` | struct | Selection and analysis parameters |
| `Results` | struct | Stores `.Classification`, `.Regression`, `.DimReduction` results |

## Notes

- `classify` and `regress` at the `"Unit"` and `"Recording"` levels group CV folds by culture (`CVLevel = "culture"`, default; `IdentityKeys`) or by `RecordingID` (`CVLevel = "recording"`), so rows of one culture/recording never appear in both train and test.
- `classify` at `"Culture"` level requires `IdentityKeys`, `GroupingVar`, and `GroupingValues`. Cultures missing any required time point are dropped.
- Results are stored in `exp.Results` and also returned directly, so you can use either pattern.
- Calling `classify` twice for the same `classification_var` overwrites the previous result (with a warning).
- The `FeatureStore` is accessible via `exp.FeatureStore` for direct feature matrix extraction.
