# Tutorial 3 — Phenotype Analysis

**Source:** [`Tutorials/3_PhenotypeAnalysis/PhenotypeAnalysis.m`](https://github.com/hornauerp/DeePhys/blob/main/Tutorials/3_PhenotypeAnalysis/PhenotypeAnalysis.m)

## What this tutorial covers

Classification, dimensionality reduction, and regression at unit, recording, and culture levels via the `Experiment` orchestration layer and the standalone `Classifier` / `DimReducer` / `Regressor` classes:

1. Build an `Experiment` from saved `RecordingProcessor` files (`Experiment.fromProcessors` / `Experiment.fromPaths`; `Experiment.fromLegacyGroup` migrates an old `RecordingGroup`). There is no factory for a saved `FeatureStore` alone: use `FeatureStore.load` with the direct `Classifier` / `Regressor` API instead
2. Run unit-level classification (e.g., WT vs. mutant)
3. Run culture-level regression (e.g., predict the drug concentration)
4. Run UMAP dimensionality reduction
5. Inspect `ClassificationResult` and `RegressionResult` objects
6. Extract feature importances

## Prerequisites

- Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
- DeePhys on the MATLAB path

## Key objects introduced

| Object | Role |
|---|---|
| `Experiment` | Orchestrates multi-recording ML analyses |
| `Classifier` | Random forest / SVM / KNN classification with grouped CV |
| `Regressor` | Random forest / SVM / KNN regression with grouped CV |
| `DimReducer` | UMAP / PCA / t-SNE dimensionality reduction |
| `ClassificationResult` | Stores accuracy, confusion matrix, feature importances (report `balanced_accuracy` from `summarizeFolds`) |
| `RegressionResult` | Stores R², RMSE, feature importances (report the pooled out-of-fold `R2_pooled` from `summarizeFolds`) |

## Example

```matlab
procs = RecordingProcessor.loadMany(proc_paths);
exp   = Experiment.fromProcessors(procs);

% Unit-level classification
opts = struct('Algorithm', 'rf', 'KFold', 5, 'CVLevel', 'recording');
result = exp.classify('Unit', 'Mutation', opts);
summary = ClassificationResult.summarizeFolds(result);
fprintf('Balanced accuracy: %.2f (fold accuracy %.2f)\n', summary.balanced_accuracy, summary.mean_accuracy);

% Dimensionality reduction
opts_dr = struct('Method', 'UMAP', 'FeatureGroups', 'all');
exp.reduce('Unit', opts_dr);
reduction = exp.Results.DimReduction.Unit.UMAP;
```

## Next step

[Tutorial 4 — Cell Type Classification](04_cell_type_classification.md)
