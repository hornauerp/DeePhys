# Tutorial 6 — Transfer Learning

**Source:** [`Tutorials/6_TransferLearning/TransferLearning.m`](https://github.com/hornauerp/DeePhys/blob/dev/Tutorials/6_TransferLearning/TransferLearning.m)

## What this tutorial covers

Two scenarios for classifying units when labelled data comes from a different platform or experiment:

**Scenario C2 — External labels → DeePhys units**
Patch-clamp or optogenetic labels from an external recording train a classifier for DeePhys units. No drug-response bootstrap needed.

**Scenario C3 — DeePhys labels → External units**
DeePhys-derived training labels classify units from an external recording (different device, sampling rate, etc.).

Both methods embed DeePhys and external units jointly in a single unsupervised UMAP, then apply graph label propagation — the same mechanism as Tutorial 4.

## Prerequisites

- Tutorial 4 completed (CellTypeClassifier trained)
- External recording data (same or different device)
- DeePhys on the MATLAB path

## Key concepts

| Concept | Description |
|---|---|
| Joint UMAP embedding | DeePhys and external units embedded together without supervision |
| Graph label propagation | Labels spread from seed nodes (known E/I) to unlabelled nodes via UMAP k-NN graph |
| Feature alignment | Common feature subset computed across platforms |

## Example (C2)

```matlab
% External patch-clamp labels classify DeePhys units
% First build a CellTypeClassifier with your DeePhys data
ctc = CellTypeClassifier(fs, ud, params);

% Then classify using external labels as training ground truth
ctc.classifyUnitsWithExternalTrain(external_features, external_labels);
```

## Next step

[Tutorial 7 — Legacy Migration](07_legacy_migration.md)
