# Tutorial 4 — Cell Type Classification

**Source:** [`Tutorials/4_CellTypeClassification/CellTypeClassification.m`](https://github.com/hornauerp/DeePhys/blob/dev/Tutorials/4_CellTypeClassification/CellTypeClassification.m)

## What this tutorial covers

Unsupervised excitatory/inhibitory (E/I) labelling of sorted units using a transductive graph-based pipeline:

1. Embed all units jointly in a single unsupervised UMAP
2. Detect inhibitory communities via Louvain clustering (enriched for drug-responsive units)
3. Propagate E/I labels to all units via graph label propagation
4. Save a trained `CellTypeClassifier`

No separate test-set projection is required — all units are embedded together.

## Prerequisites

- Tutorial 1 completed (RecordingProcessors saved)
- Drug-response recordings (stimulated + baseline) for inhibitory ground-truth detection
- DeePhys on the MATLAB path

## Key objects introduced

| Object | Role |
|---|---|
| `CellTypeClassifier` | Holds the UMAP embedding, Louvain partition, and E/I labels |
| `EIAnalyzer` | Downstream analysis using CellTypeClassifier output |

## Algorithm overview

```
All units → UMAP embedding
         → Louvain community detection
         → Score communities by fraction of drug-responsive units
         → Top community = inhibitory ground truth
         → Graph label propagation → E/I label for every unit
```

Labels: `1` = excitatory, `2` = inhibitory, `NaN` = unclassified.

## Example

```matlab
% Build UnitData array from processors
ud = [];
for i = 1:numel(procs)
    ud = [ud, procs(i).Units];
end

ctc = CellTypeClassifier(fs, ud, params);
ctc.identifyResponsiveUnits();
ctc.generateTrainLabels();
ctc.classifyUnits();

% Save (uses MATLAB's saveobj/loadobj)
save('CellTypeClassifier.mat', 'ctc');

% Inspect
sum(ctc.UnitLabels == 1)  % excitatory count
sum(ctc.UnitLabels == 2)  % inhibitory count
```

## Next step

[Tutorial 5 — E/I Analysis](05_ei_analysis.md)
