# NormalizationPipeline

Configurable, serializable normalization pipeline with a fit/transform pattern. Normalizes feature matrices without data leakage: parameters are fitted on the training set and applied to the test set.

Used internally by `Classifier`, `Regressor`, and `DimReducer`. Can also be used standalone for custom pipelines.

## Example

```matlab
% Use a pre-configured pipeline
np = NormalizationPipeline.groupThenGlobal();

% Fit on training data, transform both sets
[X_train_norm, np_fitted] = np.fit_transform(np, X_train, group_labels_train);
X_test_norm = np_fitted.transform(np_fitted, X_test, group_labels_test);

% ComBat batch correction (protects biological signal in Y)
np_combat = NormalizationPipeline.combatThenGlobal();
[X_norm, np_fitted] = np_combat.fit_transform(np_combat, X_train, batch_ids, ...
    CovariateLabels=Y_train);
X_test_norm = np_fitted.transform(np_fitted, X_test, batch_ids_test, ...
    CovariateLabels=Y_test);

% Build a custom pipeline manually
steps = {struct('type','group_zscore'), struct('type','clip'), struct('type','max_abs_scale')};
np_custom = NormalizationPipeline(steps);
```

## Static factory methods

Three pre-configured pipelines are available:

| Method | Steps |
|---|---|
| `NormalizationPipeline.groupThenGlobal()` | per-group z-score → global z-score → clip (±5) → max-abs scale |
| `NormalizationPipeline.globalOnly()` | global z-score → clip (±5) → max-abs scale |
| `NormalizationPipeline.combatThenGlobal()` | ComBat batch correction → global z-score → clip (±5) → max-abs scale |

## Instance methods

### `[X_out, np_fitted] = np.fit_transform(np, X, group_labels, options)`

Fit normalization parameters on `X` and return the transformed matrix and the fitted pipeline. `group_labels` is an (N × 1) string/categorical vector identifying batches or groups.

Optional name-value argument: `CovariateLabels` — (N × 1) biological class labels to protect from ComBat removal.

### `X_out = np.transform(np, X, group_labels, options)`

Apply a previously fitted pipeline to new data. Must be called after `fit_transform`.

## Pipeline steps

Steps are composable. Each is a struct with a `type` field:

| Type | Description | Requires `group_labels` |
|---|---|---|
| `group_zscore` | Per-group z-score (subtract group mean, divide by group std) | Yes |
| `global_zscore` | Global z-score (training-set mean and std) | No |
| `clip` | Clip values to ±5 (after z-scoring) | No |
| `max_abs_scale` | Scale each column to [−1, 1] by dividing by max absolute value | No |
| `combat` | ComBat parametric empirical Bayes batch correction (Johnson et al. 2007) | Yes (batch IDs) |

## ComBat details

ComBat fits batch-specific additive (γ) and multiplicative (δ) effects using parametric empirical Bayes shrinkage. Biological signal is estimated via OLS and protected:

- Pass `CovariateLabels` to prevent removal of biological signal correlated with batch.
- Singleton batches (only one sample) trigger a warning and fall back to global z-score.
- Unseen batches at transform time trigger an error (batch must be seen during fit).
- Fitted parameters are stored in `np.FitParams.combat` for inspection and serialization.

## Properties

| Name | Description |
|---|---|
| `Steps` | Cell array of step structs defining the pipeline |
| `FitParams` | Struct of fitted statistics (populated after `fit_transform`) |
