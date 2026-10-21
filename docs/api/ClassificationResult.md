# ClassificationResult

Result container for one cross-validation fold of classification. Inherits shared properties from `MLResult`. You receive a `(1 × K)` array of these from `Classifier.classify()` and `Experiment.classify()`.

## Example

```matlab
result = Classifier.classify(X, Y, opts);

% Aggregate metrics across folds
summary = ClassificationResult.summarizeFolds(result);
fprintf('Balanced accuracy: %.2f\n', summary.balanced_accuracy);   % pooled out-of-fold
fprintf('Accuracy: %.2f ± %.2f\n', summary.mean_accuracy, summary.std_accuracy);
fprintf('F1:       %.2f\n', summary.mean_F1);

% Per-fold flat table (one row per test object)
T = ClassificationResult.results2table(result);
writetable(T, 'predictions.csv');

% Feature importances (RF only)
imp = ClassificationResult.summarizeImportance(result);
disp(imp(1:10, :));   % top 10 features by mean permutation importance

% Per-fold metrics
for k = 1:numel(result)
    m = result(k).computeMetrics();
    fprintf('Fold %d — Acc: %.3f, F1: %.3f\n', k, m.Accuracy, m.F1_score);
end
```

## Instance methods

### `metrics = result.computeMetrics()`

Returns a struct with per-fold classification metrics:

| Field | Description |
|---|---|
| `Accuracy` | Fraction of correct predictions on test set |
| `N` | Number of test objects |
| `train_acc` | Training accuracy (OOB for RF, resubstitution for others) |
| `Precision` | Macro-averaged precision |
| `Recall` | Macro-averaged recall |
| `F1_score` | Macro-averaged F1 score |
| `ConfusionMatrix` | (C × C) confusion matrix |

Precision/Recall/F1 are always macro-averaged over the classes present in the fold (also for two classes), so fold values of binary and multi-class problems are the same metric. `summarizeFolds` also returns `macro_F1` over the pooled out-of-fold predictions. `train_acc` is optimistic when rows of a group (units of a culture) are correlated; report the held-out `balanced_accuracy`.

## Static methods

### `summary = ClassificationResult.summarizeFolds(result_array)`

Aggregate metrics across all folds. Returns a struct:

| Field | Description |
|---|---|
| `per_fold` | Table with columns `Fold`, `Accuracy`, `F1_score` |
| `mean_accuracy` | Mean fold accuracy |
| `std_accuracy` | Std of fold accuracies |
| `mean_F1` | Mean fold F1 score |
| `std_F1` | Std of fold F1 scores |
| `balanced_accuracy` | Mean per-class recall over the pooled out-of-fold predictions; robust to class imbalance and unequal fold sizes (the number to report and to test) |
| `macro_F1` | Macro-averaged F1 over the same pooled predictions |

### `T = ClassificationResult.results2table(result_array)`

Convert to a flat table with one row per test object across all folds. Columns: `Fold`, `Y_pred`, `Y_test`, `Correct`, `Score_1`, `Score_2`, …

### `T = ClassificationResult.summarizeImportance(result_array)`

Aggregate predictor importance across folds (RF only). Which importance it is depends on `ImportanceMethod` in `Classifier.classify`: with `CVGroups` (default `"auto"`), the increase in held-out classification error when a feature is shuffled in the test fold (one shuffle per fold, so noisy for small folds); otherwise out-of-bag permuted-predictor importance. Returns a table sorted by `MeanImportance` descending:

| Column | Description |
|---|---|
| `Feature` | Feature name |
| `MeanImportance` | Mean importance across folds |
| `StdImportance` | Std across folds |
| `MinImportance` | Min across folds |
| `MaxImportance` | Max across folds |

Returns an empty table when no fold has importance data (e.g., non-RF algorithm).

## Properties (inherited from `MLResult`)

| Name | Description |
|---|---|
| `Mdl` | Trained compact classifier model |
| `Y_pred` | Predicted labels (test set) |
| `Y_test` | True labels (test set) |
| `train_acc` | Training accuracy |
| `predImp` | Permutation importance vector of this fold (RF only; `[]` otherwise or with `ComputeImportance = false`); held-out-fold or out-of-bag, see `ImportanceMethod` |
| `Features` | Feature names used (string array) |
| `Parameters` | Struct with algorithm, normalization_var, K_fold, N_hyper, seed |
| `Preprocessing` | Fitted imputation, normalisation and column mask of this fold; replayed by `Classifier.predict(result(k), X_new)` |

## Properties (ClassificationResult-specific)

| Name | Description |
|---|---|
| `scores` | (N_test × C) posterior probabilities / prediction scores; column j is class j of `GroupLabels` (`NaN` if the fold's training rows held no such class) |
| `GroupLabels` | Class label names; `Y_pred` and `Y_test` hold class indices into this list |
