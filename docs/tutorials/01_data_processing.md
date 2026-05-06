# Tutorial 1 — Data Processing

**Source:** [`Tutorials/1_DataProcessing/DataProcessing.m`](https://github.com/hornauerp/DeePhys/blob/dev/Tutorials/1_DataProcessing/DataProcessing.m)

## What this tutorial covers

The full pipeline from raw Kilosort output to a `FeatureStore` ready for analysis:

```
SpikeData → RecordingProcessor → FeatureStore
```

Specifically:

1. Define paths and a metadata struct
2. Load one Kilosort directory into a `SpikeData` object
3. Create a `RecordingProcessor` and run QC (quality-control unit filtering)
4. Compute unit-level and network-level features
5. Compute connectivity (CCG / STTC)
6. Save the `RecordingProcessor` to disk
7. Assemble multiple processors into a `FeatureStore`

## Prerequisites

- DeePhys on the MATLAB path (`startup` from repo root)
- At least one Kilosort output directory

## Key objects introduced

| Object | Role |
|---|---|
| `SpikeData` | Immutable container for raw Kilosort output |
| `RecordingProcessor` | Lazy computation engine for one recording |
| `FeatureStore` | Table-centric container for features across recordings |

## Running the tutorial

Open `Tutorials/1_DataProcessing/DataProcessing.m` in MATLAB, fill in the path placeholders in Section 1, then run section by section with Ctrl+Enter.

```matlab
% Minimal example
sd   = SpikeData.fromKilosort('/path/to/ks', metadata);
proc = RecordingProcessor(sd);
proc.runQC();
proc.computeUnitFeatures();
proc.computeNetworkFeatures();
proc.save('/path/to/save');

fs = FeatureStore.fromProcessors([proc1, proc2]);
fs.save('experiment.mat');
```

## Next step

[Tutorial 2 — Feature Exploration](02_feature_exploration.md)
