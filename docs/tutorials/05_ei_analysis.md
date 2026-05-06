# Tutorial 5 — E/I Analysis

**Source:** [`Tutorials/5_EIAnalysis/EIAnalysis.m`](https://github.com/hornauerp/DeePhys/blob/dev/Tutorials/5_EIAnalysis/EIAnalysis.m)

## What this tutorial covers

Network-level E/I analysis using labels from a trained `CellTypeClassifier`:

1. Load a trained `CellTypeClassifier`
2. Compute excitatory and inhibitory population firing-rate traces
3. Detect burst episodes and extract aligned burst cutouts
4. Quantify unit–population correlations
5. Compare E/I dynamics across conditions

## Prerequisites

- Tutorial 4 completed (CellTypeClassifier trained and saved)
- DeePhys on the MATLAB path

## Key objects introduced

| Object | Role |
|---|---|
| `EIAnalyzer` | Computes E/I activity traces and burst statistics from CellTypeClassifier output |

## What "E/I analysis" means here

Rather than a single ratio, DeePhys quantifies the temporal dynamics of excitatory and inhibitory population activity separately:

- **Population firing rate traces** — smoothed spike density for E and I subpopulations
- **Burst detection** — ISIN algorithm applied to population traces
- **Burst-aligned cutouts** — E and I activity in ±window around each burst onset
- **Unit–population correlations** — how strongly each unit's firing tracks its population

## Example

```matlab
data = load('CellTypeClassifier.mat');
ctc  = data.ctc;

ei = EIAnalyzer(ctc);
ei.computeActivity();
ei.detectBursts();
ei.extractBurstCutouts();
ei.PlotBurstCutouts();
```

## Next step

[Tutorial 6 — Transfer Learning](06_transfer_learning.md)
