# MLPipeline

Static ML helpers: classifier/regressor creation, cross-validation splits, and label pooling. **No instantiation needed** — all methods are static.

## Example

```matlab
% Train a random forest classifier
[clf, train_acc] = MLPipeline.createClassifier(X_train, Y_train, 'rf');

% Predict on held-out data
Y_pred = predict(clf, X_test);

% Grouped k-fold CV (keeps recordings from the same chip together)
% Uses the GroupedCV class, not MLPipeline
gcv = GroupedCV.byGroups(group_ids, 5, Y_train);

% Default hyperparameters
params = MLPipeline.returnDefaultParams();
params.RF.NumCycles = 1000;
[clf, ~] = MLPipeline.createClassifier(X_train, Y_train, 'rf', 0, params);

% Hyperparameter search whose inner folds keep groups (e.g. cultures) together
[clf, ~] = MLPipeline.createClassifier(X_train, Y_train, 'rf', 30, params, train_groups);
[mdl, ~] = MLPipeline.createRegressor(X_train, y_train, 'rf', 30, train_groups);
```

## Static methods

### `MLPipeline.returnDefaultParams()`

Returns the default hyperparameter struct. Fields:

| Field | Default | Description |
|---|---|---|
| `RF.NumCycles` | 500 | Number of trees |
| `RF.MinLeafSize` | 5 | Minimum leaf size |
| `RF.NumVariablesToSample` | `'auto'` | sqrt(F) for classification, F/3 for regression |
| `RF.Surrogate` | `'on'` | Surrogate splits for missing values |
| `RF.Reproducible` | `true` | Reproducible tree templates |
| `RF.Prior` | `'empirical'` | Class prior (`'uniform'` for balanced weighting) |
| `UMAP.NNeighbors` | 100 | k-NN for UMAP graph |
| `UMAP.NNeighborsCulture` | 10 | k-NN for culture-level UMAP |
| `UMAP.MinDist` | 1 | Minimum distance in embedding |
| `UMAP.Spread` | 5 | Effective scale of embedded points |
| `UMAP.SGDTasks` | 20 | Stochastic gradient descent tasks |
| `UMAP.ClusterDetail` | `'adaptive'` | Cluster resolution strategy |

### `MLPipeline.createClassifier(X_train, Y_train, alg, N_hyper, params, train_groups)`

Train a classifier with optional hyperparameter optimisation.

| Argument | Description |
|---|---|
| `X_train` | (N × F) feature matrix |
| `Y_train` | (N × 1) class labels |
| `alg` | `"rf"` (default), `"svm"`, `"cnb"` (naive Bayes), `"knn"` |
| `N_hyper` | Hyperparameter optimisation evaluations (0 = skip) |
| `params` | Parameter struct from `returnDefaultParams()` |
| `train_groups` | (N × 1) CV group IDs of the training rows (default `[]`). Only used when `N_hyper > 0`: the inner 5-fold cross-validation of the hyperparameter search then keeps each group together (`MLUtils.groupedPartition`), so the search does not reward hyperparameters that memorise groups. `Classifier.classify` passes the groups of the training rows when `CVGroups` is set |

Returns `[clf, train_acc]`.

### `MLPipeline.createRegressor(X_train, Y_train, alg, N_hyper, train_groups)`

Same interface as `createClassifier` for regression tasks (no `params` argument; `train_groups` as above, used when `N_hyper > 0`). Supported algorithms: `"rf"`, `"svm"`, `"knn"`. Returns `[mdl, train_r2]`. `"knn"` returns a `KNNRegressionModel` (MATLAB has no `fitrknn`): K = 5, or with `N_hyper > 0` chosen by leave-one-out error (leave-one-group-out when `train_groups` is given).

## Notes

- `'rf'` trains a bagged ensemble via `fitcensemble`/`fitrensemble` (MATLAB Statistics and Machine Learning Toolbox).
- Set `N_hyper > 0` to enable Bayesian hyperparameter optimisation (slow; use for final models).
- `RF.Prior = 'uniform'` gives equal weight to each class regardless of class frequency — useful for imbalanced datasets.
- For classification, OOB error is used for `'rf'` training accuracy; resubstitution loss for other algorithms. With `N_hyper > 0`, `train_acc` is instead one minus the cross-validated objective of the search. All of these are optimistic for correlated (grouped) rows.
- For regression, OOB R² is reported for `'rf'`; resubstitution R² for `'svm'`/`'knn'`.
