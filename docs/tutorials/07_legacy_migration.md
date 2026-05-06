# Tutorial 7 — Legacy Migration

**Source:** [`Tutorials/7_LegacyMigration/LegacyMigration.m`](https://github.com/hornauerp/DeePhys/blob/dev/Tutorials/7_LegacyMigration/LegacyMigration.m)

## What this tutorial covers

Converts old `MEArecording` / `RecordingGroup` objects saved from DeePhys v1 to the new v2 API. After migration, all downstream analyses (Tutorials 2–6) work without change.

## Migration paths

| Path | Use when |
|---|---|
| A — `RecordingProcessor.fromLegacyMat(path)` | You have a single `.mat` with a `MEArecording` object |
| B — `Experiment.fromLegacyGroup(rg)` | You have a `RecordingGroup` in the workspace |
| C — `FeatureStore.fromLegacyRecordings(mearec_array)` | You have a `MEArecording` array and only need features |
| D — `UnitData.fromLegacyUnit(unit_obj, recording_id)` | You need to migrate a single `Unit` object |
| E — `UnitData.fromLegacyUnitArray(unit_array, recording_id)` | You need to batch-migrate an array of `Unit` objects |

## Prerequisites

- Saved `.mat` files containing old `MEArecording` or `RecordingGroup` objects
- DeePhys on the MATLAB path

## Example

```matlab
% Path A: single recording
proc = RecordingProcessor.fromLegacyMat('/path/to/MEArecording.mat');
proc.save('/path/to/save');

% Path C: batch feature-only migration
fs = FeatureStore.fromLegacyRecordings(mearec_array);
fs.save('migrated_experiment.mat');
```

## What is preserved

- All electrophysiological features present in the old objects
- Metadata fields
- Unit spike times and waveforms (where stored)

## What changes

- Handle-class objects → value classes (no more `saveobj`/`loadobj` hacks)
- `RecordingGroup` analyses → `Experiment` methods
- Per-unit features accessed via `FeatureStore.UnitTable` instead of `Unit.Features`
