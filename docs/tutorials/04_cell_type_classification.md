# Tutorial 4 — Cell Type Classification

**Source:** [`Tutorials/4_CellTypeClassification/CellTypeClassification.m`](https://github.com/hornauerp/DeePhys/blob/dev/Tutorials/4_CellTypeClassification/CellTypeClassification.m)

## What this tutorial covers

Unsupervised excitatory/inhibitory (E/I) labelling of sorted units using a transductive graph-based pipeline:

1. Build one deterministic k-nearest-neighbour graph of all units (drug-free ACG + waveform)
2. Find communities (Leiden) significantly enriched in drug-responsive units: their responders become the inhibitory seeds, units of depleted communities the excitatory seeds (consensus over restarts)
3. Propagate E/I labels to all units over the graph
4. Check the result with `qualityReport` instead of tuning parameters

No separate test-set projection is required — all units are classified together.

## Prerequisites

- Tutorial 1 completed (RecordingProcessors saved)
- Drug-response recordings (stimulated + baseline) for inhibitory ground-truth detection
- DeePhys on the MATLAB path

## Key objects introduced

| Object | Role |
|---|---|
| `CellTypeClassifier` | Holds the unit graph, the community partition, the training seeds and the E/I labels |
| `EIAnalyzer` | Downstream analysis using CellTypeClassifier output |

## Algorithm overview

```
All units → exact kNN graph (drug-free ACG + waveform)
         → Leiden communities (30 restarts)
         → Hypergeometric test per community (BH-FDR): enriched / depleted in responders
         → Consensus seeds: responders in enriched communities = inhibitory,
                            units of depleted communities = excitatory
         → Graph label propagation → E/I label for every unit
         → qualityReport (seed counts, condition balance, independent validation)
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


## Culture composition and confident labels

The graph labels are counted at the balanced seed prior. For the inhibitory fraction of each
culture, and for labels that use it, run the quantification step (calibrated probabilities
with a per-culture EM prior adjustment):

```matlab
T = ctc.estimateComposition();          % per culture: Fraction, CILow, CIHigh
R = ctc.spikeCountReliability();        % label stability vs spike count (spike-train thinning)
labels = ctc.confidentLabels(MinProbability=0.8);   % NaN where not confident / too few spikes
```

`confidentLabels` keeps units with P(inhibitory) >= 0.8 or <= 0.2 and at least
`R.MinSpikes` spikes (the smallest tested count with >= 90 % label agreement after thinning).
With the composition probabilities, the expected share of wrong labels among the kept units is
printed. Downstream analyses can use only these units, e.g.
`EIAnalyzer(ctc, struct('Labels', struct('Source', "confident")))`.

## Several perturbations of the same cell class

When interneurons were both activated and silenced in different cultures (e.g. hM3Dq and
hM4Di), both responder sets can seed the classifier. Each rule names the cultures and the
response direction of the target cells:

```matlab
params.Bootstrap.ResponseRules = {'Construct', 1, "increase"; 'Construct', 2, "decrease"};
ctc = CellTypeClassifier(fs, ud, params);
ctc.identifyResponsiveUnits();          % responders of all rules, pooled
ctc.discardSegmentUnits();
ctc.generateTrainLabels();
ctc.classifyUnits();
V = ctc.validateResponseRules();        % train on each rule alone, validate on the others
```

Pool the rules only if `validateResponseRules` shows that each one predicts the responders of
the other; otherwise keep one perturbation out of training as the independent check.
## Next step

[Tutorial 5 — E/I Analysis](05_ei_analysis.md)
