# ClassificationResult

Result container for one cross-validation fold of classification. Inherits shared properties from `MLResult`. You receive a `(1 × K)` array of these from `Classifier.classify()` and `Experiment.classify()`.

## Example

```matlab
result = Classifier.classify(X, Y, opts);

% Aggregate metrics across folds
summary = ClassificationResult.summarizeFolds(result);
fprintf('Accuracy: %.2f ± %.2f\n', summary.mean_accuracy, summary.std_accuracy);
fprintf('F1:       %.2f\n', summary.mean_F1);

% Per-fold flat table (one row per test object)
T = ClassificationResult.results2table(result);
writetable(T, 'predictions.csv');

% Feature importances (RF only)
imp = ClassificationResult.summarizeImportance(result);
disp(imp(1:10, :));   % top 10 features by mean OOB importance

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

For binary classification, Precision/Recall/F1 are for the first class. For multi-class, they are macro-averaged.

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

### `T = ClassificationResult.results2table(result_array)`

Convert to a flat table with one row per test object across all folds. Columns: `Fold`, `Y_pred`, `Y_test`, `Correct`, `Score_1`, `Score_2`, …

### `T = ClassificationResult.summarizeImportance(result_array)`

Aggregate OOB predictor importance across folds (RF only). Returns a table sorted by `MeanImportance` descending:

| Column | Description |
|---|---|
| `Feature` | Feature name |
| `MeanImportance` | Mean OOB importance across folds |
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
| `predImp` | OOB predictor importance vector (RF only; `[]` otherwise) |
| `Features` | Feature names used (string array) |
| `Parameters` | Struct with algorithm, K_fold, N_hyper, seed |

## Properties (ClassificationResult-specific)

| Name | Description |
|---|---|
| `scores` | (N_test × C) posterior probabilities / prediction scores |
| `GroupLabels` | Class label names used in this fold |
