# RegressionResult

Result container for one cross-validation fold of regression. Inherits shared properties from `MLResult`. You receive a `(1 × K)` array of these from `Regressor.regress()` and `Experiment.regress()`.

## Example

```matlab
result = Regressor.regress(X, Y, opts);

% Aggregate metrics
summary = RegressionResult.summarizeFolds(result);
fprintf('R² (pooled out-of-fold): %.3f  RMSE: %.3f\n', summary.R2_pooled, sqrt(summary.mean_MSE));

% Per-fold flat table
T = RegressionResult.results2table(result);

% Feature importances (RF only)
imp = RegressionResult.summarizeImportance(result);
disp(imp(1:10, :));

% Per-fold metrics
for k = 1:numel(result)
    m = result(k).computeMetrics();
    fprintf('Fold %d — R²: %.3f, RMSE: %.3f\n', k, m.R2, m.RMSE);
end
```

## Instance methods

### `metrics = result.computeMetrics()`

Returns a struct with per-fold regression metrics:

| Field | Description |
|---|---|
| `MSE` | Mean squared error on test set |
| `RMSE` | Root mean squared error |
| `MAE` | Mean absolute error |
| `R2` | Coefficient of determination (1 − SS_res / SS_tot); `NaN` if target is constant |
| `Correlation` | Pearson correlation between Y_pred and Y_test |
| `N` | Number of test objects |
| `mse_train` | Training MSE (OOB for RF, resubstitution otherwise; optimistic for correlated rows) |

## Static methods

### `summary = RegressionResult.summarizeFolds(result_array)`

Aggregate metrics across all folds. Returns a struct:

| Field | Description |
|---|---|
| `per_fold` | Table with columns `Fold`, `MSE`, `R2`, `Correlation` |
| `mean_MSE` | Mean fold MSE |
| `mean_R2` | Mean fold R² (folds with an undefined R² are skipped) |
| `mean_Correlation` | Mean fold Pearson correlation (undefined folds skipped) |
| `R2_pooled` | `1 − SS_res / SS_tot` over the pooled out-of-fold predictions; `NaN` if the pooled target is constant. **The R² to report and to test** (`Regressor.permutationTest` uses it) |

Fold R² compares each fold with its own mean; with grouped CV a fold can hold a single target value (R² undefined) or few rows (unstable), so prefer `R2_pooled` over `mean_R2`.

### `T = RegressionResult.results2table(result_array)`

Convert to a flat table with one row per test object across all folds. Columns: `Fold`, `Y_pred`, `Y_test`, `Residual`.

### `T = RegressionResult.summarizeImportance(result_array)`

Aggregate predictor importance across folds (RF only). With `CVGroups` (default `ImportanceMethod = "auto"`) it is the increase in held-out MSE when a feature is shuffled in the test fold (one shuffle per fold, so noisy for small folds); otherwise out-of-bag permuted-predictor importance. Returns a table sorted by `MeanImportance` descending. See [ClassificationResult.summarizeImportance](ClassificationResult.md) for column descriptions.

## Properties (inherited from `MLResult`)

| Name | Description |
|---|---|
| `Mdl` | Trained compact regression model |
| `Y_pred` | Predicted values (test set) |
| `Y_test` | True values (test set) |
| `train_acc` | Training R² (OOB for RF, resubstitution for others) |
| `predImp` | Permutation importance vector of this fold (RF only; `[]` otherwise or with `ComputeImportance = false`); held-out-fold or out-of-bag, see `ImportanceMethod` |
| `Features` | Feature names used (string array) |
| `Parameters` | Struct with algorithm, normalization_var, K_fold, N_hyper, seed |
| `Preprocessing` | Fitted imputation, normalisation and column mask of this fold |

## Properties (RegressionResult-specific)

| Name | Description |
|---|---|
| `mse_train` | Training mean squared error (OOB for RF, resubstitution otherwise) |
