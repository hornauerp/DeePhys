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
- **Burst detection** — burst state where the I/E ratio is low and E activity is above its noise floor; burst peaks from the total population trace
- **Burst-aligned cutouts** — E and I activity in a ±window around each burst peak (cross-correlation aligned, outliers rejected)
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


Which units form the E and I populations is set by `Labels.Source`: `"graph"` (default,
`UnitLabels`), `"composition"` (`estimateComposition` labels) or `"confident"` (only units
passing `confidentLabels`; the others are left out of both populations):

```matlab
ctc.estimateComposition();
ei = EIAnalyzer(ctc, struct('Labels', struct('Source', "confident", 'MinProbability', 0.9)));
```

Labels that come from the classifier carry its features: the drug-free ACG encodes how a unit
fires around network bursts. Burst-locked E / I dynamics computed with those labels are partly
built into the labels; responder-defined populations (as in the published analysis) or
waveform-only labels avoid this.
## Next step

[Tutorial 6 — Transfer Learning](06_transfer_learning.md)
