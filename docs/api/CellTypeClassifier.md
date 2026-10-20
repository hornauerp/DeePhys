# CellTypeClassifier

Transductive excitatory / inhibitory (E/I) labelling of sorted units. **Handle class.**

All units of an experiment are placed on one deterministic k-nearest-neighbour graph (ACG and waveform features). Ground truth, either drug responders or culture metadata, seeds label propagation over this graph, so every unit is classified together and no separate test set is projected. The default pipeline has no parameters to tune; check a run with `qualityReport` instead.

Labels: `1` = excitatory, `2` = inhibitory, `NaN` = unclassified.

## Example

```matlab
ctc = CellTypeClassifier(fs, ud, params);   % FeatureStore + UnitData array from the same assembly

ctc.identifyResponsiveUnits();     % ground truth: firing-rate response, or culture metadata
ctc.discardSegmentUnits();         % dose-response designs: keep baseline units only
ctc.generateTrainLabels();         % graph, communities, training seeds
ctc.classifyUnits();               % label propagation over the graph

labels = ctc.UnitLabels;           % 1 = excitatory, 2 = inhibitory, NaN = unclassified
report = ctc.qualityReport();      % acceptance checks, no ground truth needed
```

## Ground truth (`Parameters.Bootstrap.GroundTruthMethod`)

| Method | Inhibitory seeds | Excitatory seeds |
|---|---|---|
| `"two_window"` | Units whose rate increases between a pre and a post window of one recording (bootstrap permutation test) | Inferred from the non-responsive pool |
| `"full_curve"` (default) | Units with a monotonic dose-response, a minimum fold change and a confirmed difference between baseline and maximum dose. Needs units tracked across recordings. | Inferred from the non-responsive pool |
| `"control_calibrated"` | As `"full_curve"`, with the response thresholds calibrated on control cultures (recommended when controls exist) | Inferred from the non-responsive pool |
| `"metadata"` | Cell-type column of the `UnitTable` (pure cultures) | Cell-type column of the `UnitTable` |

Several perturbations of the same cell class can be combined with `Bootstrap.ResponseRules`; `validateResponseRules` checks that they identify the same class.

## Pipeline methods

| Method | Role |
|---|---|
| `identifyResponsiveUnits(filter)` | Assign ground truth (`ResponsiveUnitIdx`, `ResponsiveUnitDirection`, `ResponsiveStrength`) |
| `discardSegmentUnits()` | Drop non-baseline segments of dose-response recordings (call before `generateTrainLabels`) |
| `generateTrainLabels()` | Normalisation, unit graph, communities and training seeds. Default seeds: Leiden communities enriched (inhibitory) or depleted (excitatory) in responders, as a consensus over restarts. |
| `classifyUnits()` | Label propagation over the graph (`Classification.Method = "graph"`) or distance-weighted kNN (`"knn"`). Sets `UnitLabels` and `UnitConfidence`. |
| `classifyUnitsEnsemble()` | Majority vote over several seeds (`Parameters.Ensemble`, off by default); `UnitConfidence` is the share of agreeing runs |
| `classifyExternalUnits(wf, acg, sr)` | Label units from other data (waveforms, ACGs, sampling rate) on a joint graph with the DeePhys units. The external cut-off is label-free and opt-in. |
| `classifyUnitsWithExternalTrain(wf, acg, labels, sr)` | The reverse: externally labelled units (e.g. patch-clamp or optogenetically identified) seed the labels of the DeePhys units |
| `graphMatrix()` | The symmetric unit graph used for communities and propagation |
| `clearCache()` | Invalidate cached preprocessing after changing `Parameters` |

## Quantification and confident labels

| Method | Result |
|---|---|
| `estimateComposition(GroupBy=..., NComponents=..., NBoot=...)` | Inhibitory fraction per culture with confidence intervals, and a calibrated probability of being inhibitory per unit. Probabilities come from a ridge logistic regression on leading principal components; each culture's fraction is re-estimated by EM under prior shift, so counting labels is not pulled towards the balance of the seeds. Stored in `Composition`. With both classes labelled (`"metadata"`) it trains on every labelled unit. |
| `confidentLabels(MinProbability=..., MinSpikes=...)` | Labels restricted to units with P(inhibitory) at or above `MinProbability` (or at or below 1 − `MinProbability`); `NaN` for the rest. Uses the calibrated composition probabilities when available and reports the expected share of wrong labels. |
| `spikeCountReliability(Counts=..., NUnits=...)` | How many spikes a unit needs for a stable label, by thinning the spike trains of well-sampled units. `R.MinSpikes` can be passed to `confidentLabels`. |

## Checks

| Method | What it checks |
|---|---|
| `qualityReport(ConditionVar=..., ValidationFilter=...)` | PASS / WARN for training seeds, responder calls (estimated FDR from controls), condition balance, an optional independent check and the inhibitory fraction with a confidence interval. Run after `classifyUnits`. |
| `validateIndependentResponders(filter)` | AUC of P(inhibitory) separating expected interneurons from expected excitatory neurons among the responders of an independent condition (cultures not used for the seeds) |
| `validateResponseRules()` | With several response rules: train on each rule alone and validate on the others |
| `validateTrainingLabels()` | Leave-one-culture-out agreement of the training labels |
| `assessStability(NRuns=...)` | Agreement (adjusted Rand index) of the labels across random seeds |

`optimizeUnsupervisedUMAP()` remains for the optional graph and community settings. `optimizeHyperparameters()` was removed: the default pipeline needs no tuning.

## Key properties

| Property | Content |
|---|---|
| `Parameters` | Defaults merged with the user overrides (`returnDefaultParams`) |
| `ResponsiveUnitIdx`, `ResponsiveUnitDirection` | Ground truth from `identifyResponsiveUnits` |
| `TrainLabels` | Training seeds and diagnostics |
| `UnitLabels`, `UnitConfidence` | Labels and their confidence (1 for training units) |
| `Graph` | Unit graph: `.W` (sparse), `.k`, `.Method` (`"exact"` or `"umap"`) |
| `Composition` | `estimateComposition` result: `.Table` per culture, `.Probability` and `.Labels` per unit |

## Saving

`save('ctc.mat', 'ctc')` keeps everything needed for inference. The UMAP model and the preprocessing caches are not serialisable; they are regenerated on the next pipeline call.

## Migration

`CellTypeClassifier.fromLegacyGroup(rg, params)` builds a classifier from an old `RecordingGroup`.

## See also

- [Tutorial 4 — Cell Type Classification](../tutorials/04_cell_type_classification.md)
- [EIAnalyzer](EIAnalyzer.md): E/I activity of the classified cultures
