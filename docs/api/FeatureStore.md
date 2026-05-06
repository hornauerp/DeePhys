# FeatureStore

Table-centric container for features across a collection of recordings. **Handle class.** The central architectural piece of DeePhys v2 — all downstream analysis (ML, statistics, plotting) works with `FeatureStore` tables.

## Example

```matlab
% Build from processed recordings
fs = FeatureStore.fromProcessors([proc1, proc2, proc3]);
fs.save('experiment.mat');

% Load later
fs = FeatureStore.load('experiment.mat');

% Subset to specific conditions
fs_wt = fs.subsetByMetadata('Mutation', {'WT'});

% Feature matrix for classification
[X, unitIDs] = fs.unitMatrix('all');
Y = fs.UnitTable.Mutation;
result = Classifier.classify(X, Y, opts);
```

## Static constructors

| Method | Description |
|---|---|
| `FeatureStore.fromProcessors(proc_array)` | Assemble from a `RecordingProcessor` array |
| `FeatureStore.fromLegacyRecordings(mearec_array)` | Migrate from old `MEArecording` array (optional second and third arguments `unit_features` and `network_features` are string arrays selecting which feature groups to include; default `"all"`) |
| `FeatureStore.load(path)` | Load a saved `FeatureStore.mat` |

## Instance methods

| Method | Description |
|---|---|
| `fs.save(path)` | Save to `.mat` |
| `fs.subset(recording_ids)` | Return a new `FeatureStore` containing only the specified recordings |
| `fs.subsetByMetadata(field, values)` | Return a new `FeatureStore` filtered by a metadata column |
| `[X, ids] = fs.unitMatrix(feature_groups, parent_features, ...)` | Extract (N_units × F) feature matrix |
| `[X, ids] = fs.recordingMatrix(feature_groups, parent_features, ...)` | Extract (N_recordings × F) feature matrix |
| `[X, ids] = fs.cultureMatrix(identity_keys, grouping_var, grouping_values, normalization, ...)` | Extract (N_cultures × F) wide feature matrix |
| `report = fs.summarizeFeatureGroups()` | Print and return a diagnostic report of feature group coverage |

## Properties

| Name | Type | Description |
|---|---|---|
| `UnitTable` | table | One row per unit: UnitID, RecordingID, metadata columns, feature columns |
| `RecordingTable` | table | One row per recording: RecordingID, metadata columns, feature columns |
| `MetadataTable` | table | One row per recording: RecordingID, all metadata fields |

## `unitMatrix` options

```matlab
[X, ids] = fs.unitMatrix('all')                       % all child features
[X, ids] = fs.unitMatrix('all', 'ACG')                % prefer Parent_ACG* over child ACG*
[X, ids] = fs.unitMatrix('all', 'all')                % all parent features where available
[X, ids] = fs.unitMatrix('all', [], FeatureSet="core") % curated non-redundant feature subset
[X, ids] = fs.unitMatrix('all', [], MinFiringRate=0.1) % exclude quiet units
```

`feature_groups` can be `"all"` or a string array like `["ActivityFeatures","WaveformFeatures"]`.

## `cultureMatrix`

Aggregates recording-level features into a wide table with one row per culture (identified by `identity_keys`, e.g. `["ChipID","PlatingDate"]`). Columns are suffixed by the grouping variable value (e.g. `FiringRate_7`, `FiringRate_14`).

```matlab
[X, culture_ids] = fs.cultureMatrix(["ChipID","PlatingDate"], "DIV", [7 14 21 28], "baseline")
```

## Notes

- `unitMatrix('all')` excludes metadata and `Parent_*` columns automatically (only child feature columns are returned).
- `subsetByMetadata` accepts a cell array of values for OR logic.
- The three tables are always kept in sync — never modify them directly.
- `summarizeFeatureGroups()` is useful after `fromProcessors()` to confirm which analyses were run.
