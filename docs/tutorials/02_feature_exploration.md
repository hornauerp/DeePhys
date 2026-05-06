# Tutorial 2 — Feature Exploration

**Source:** [`Tutorials/2_FeatureExploration/FeatureExploration.m`](https://github.com/hornauerp/DeePhys/blob/dev/Tutorials/2_FeatureExploration/FeatureExploration.m)

## What this tutorial covers

How to inspect and extract feature matrices from a `FeatureStore` for plotting, quality checking, and preparing ML inputs:

1. Load a saved `FeatureStore`
2. Inspect `UnitTable` and `RecordingTable` columns
3. Subset by metadata fields (e.g., mutation, DIV)
4. Extract feature matrices for ML (`unitMatrix`, `recordingMatrix`)
5. Check for missing values and feature distributions

## Prerequisites

- Tutorial 1 completed — a saved `FeatureStore.mat`
- DeePhys on the MATLAB path

## Key methods introduced

| Method | Description |
|---|---|
| `fs.subsetByMetadata(field, values)` | Return a new FeatureStore filtered by a metadata column |
| `fs.unitMatrix('all')` | Extract (N_units × F) feature matrix, excluding metadata columns |
| `fs.recordingMatrix('all')` | Extract (N_recordings × F) feature matrix |

## Example

```matlab
fs = FeatureStore.load('experiment.mat');

% Subset to WT and HET only
fs_sub = fs.subsetByMetadata('Mutation', {'WT', 'HET'});

% Feature matrix for classification
[X, unitIDs] = fs_sub.unitMatrix('all');
Y = fs_sub.UnitTable.Mutation;
```

## Next step

[Tutorial 3 — Phenotype Analysis](03_phenotype_analysis.md)
