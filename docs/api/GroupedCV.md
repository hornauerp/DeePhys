# GroupedCV

Hierarchy-aware cross-validation that respects data structure. **Prevents data leakage** when rows are not independent — e.g., multiple units from the same recording, or multiple recordings from the same culture.

Standard `cvpartition` splits at the row level, so units from the same recording can appear in both train and test folds, inflating accuracy estimates. `GroupedCV` splits at the group level first, then expands masks to the object level.

## Example

```matlab
% Units from the same recording are always in the same fold
gcv = GroupedCV.byGroups(fs.UnitTable.RecordingID, 5, Y);
fprintf('%d folds, %d total units\n', gcv.NumFolds, numel(Y));

% Iterate folds
for k = 1:gcv.NumFolds
    [train_idx, test_idx] = gcv.fold(k);
    X_train = X(train_idx, :);
    X_test  = X(test_idx, :);
    Y_train = Y(train_idx);
    Y_test  = Y(test_idx);
    % ... fit and evaluate ...
end

% Leave-one-group-out (one recording per fold)
gcv_lofo = GroupedCV.byGroups(rec_ids, -1);

% Inspect configuration
s = gcv.describe();
fprintf('Groups: %d, mean size: %.1f\n', s.NumGroups, s.MeanGroupSize);
```

## Constructor

`gcv = GroupedCV(group_ids, K)` — split by numeric group index vector.

`gcv = GroupedCV(group_ids, K, class_labels)` — same, with stratification (dominant class per group used via `mode()`).

For raw string/categorical IDs, use the factory method below.

## Static factory

### `gcv = GroupedCV.byGroups(group_ids, K, Y)`

The preferred entry point. Accepts any unique-able vector (strings, integers, categoricals).

| Argument | Default | Description |
|---|---|---|
| `group_ids` | — | (N × 1) group assignment per object (e.g. `fs.UnitTable.RecordingID`) |
| `K` | `5` | Number of folds; `-1` = leave-one-group-out |
| `Y` | `[]` | (N × 1) class labels for stratified group-level splitting (optional) |

## Instance methods

### `[train_idx, test_idx] = gcv.fold(k)`

Returns `(N × 1)` logical masks for fold `k`.

### `summary = gcv.describe()`

Returns a struct with: `NumFolds`, `NumGroups`, `NumObjects`, `Stratified`, `GroupSizes`, `MeanGroupSize`, `MinGroupSize`, `MaxGroupSize`.

## Properties

| Name | Description |
|---|---|
| `NumFolds` | Number of CV folds |
| `GroupIDs` | (N × 1) group assignment per object |
| `UniqueGroups` | Unique group values |
| `GroupPartition` | Underlying `cvpartition` at the group level |
| `ClassLabels` | (N × 1) class labels if stratification was used; `[]` otherwise |

## Notes

- **Stratification caveat:** Class labels are reduced to one per group (the dominant class via `mode()`). When groups contain mixed classes, per-fold class balance is not guaranteed.
- When `K > number of unique groups`, K is automatically reduced to the number of groups (with a warning).
- `Classifier.classify` and `Regressor.regress` use `GroupedCV.byGroups` internally when `opts.CVGroups` is set.
