# RecordingProcessor

Lazy computation engine for one recording. **Handle class.** Each `compute*` method checks `Status` before running and short-circuits if the step is already done — safe to call repeatedly.

## Example

```matlab
% Full pipeline for one recording
sd   = SpikeData.fromKilosort('/path/to/ks', metadata);
proc = RecordingProcessor(sd);
proc.runQC();
proc.computeUnitFeatures();
proc.computeNetworkFeatures();
proc.computeConnectivity();
proc.save('/path/to/ks/RecordingProcessor.mat');

% Load later
proc = RecordingProcessor.load('/path/to/ks/RecordingProcessor.mat');

% One-step construction from Kilosort output
proc = RecordingProcessor.fromKilosort('/path/to/ks', metadata);

% Batch load
procs = RecordingProcessor.loadMany({path1, path2, path3});
```

## Constructor

`proc = RecordingProcessor(spike_data)` — construct without running any computation.

`proc = RecordingProcessor(spike_data, parameters)` — override default parameters.

## Static constructors

| Method | Description |
|---|---|
| `RecordingProcessor.fromKilosort(path, metadata)` | Construct from Kilosort directory (wraps `SpikeData.fromKilosort`) |
| `RecordingProcessor.fromLegacyMat(path)` | Migrate from old `MEArecording` `.mat` file |
| `RecordingProcessor.load(path)` | Load a saved `RecordingProcessor.mat` |
| `RecordingProcessor.loadMany(paths)` | Load an array of saved processors |
| `RecordingProcessor.returnDefaultParams()` | Return the default parameter struct |

## Computation methods

| Method | What it does |
|---|---|
| `proc.runQC()` | Filter templates by quality criteria; populates `proc.Units` |
| `proc.computeUnitFeatures()` | Extract per-unit electrophysiological features into `UnitFeatureTable` |
| `proc.computeParentFeatures()` | Compute ACGs from the parent (concatenated) recording; stores as `Parent_ACG*` columns in `UnitFeatureTable` |
| `proc.computeNetworkFeatures()` | Extract network-level features into `NetworkFeatureTable` |
| `proc.computeConnectivity()` | Compute CCG / FullCCG / STTC / DDC connectivity; stores in `proc.Connectivity` |
| `proc.computeCellTypeFeatures()` | E/I-stratified graph, balance, burst, and correlation features (requires `CellTypeLabels`) |
| `proc.computeSpatialAnalysis()` | Spatial distribution features for units and cell types |
| `proc.runAll()` | Convenience: run all enabled analyses in sequence |

| Static method | What it does |
|---|---|
| `RecordingProcessor.applyLabelsFromClassifier(proc_array, ctc)` | Copy `UnitLabels` from a `CellTypeClassifier` to processors by UnitID matching |

## Properties

| Name | Type | Description |
|---|---|---|
| `SpikeData` | SpikeData | Raw recording data (value class) |
| `Units` | UnitData | (1 × N) QC-passing units |
| `Parameters` | struct | Analysis parameters (merged with defaults) |
| `UnitFeatureTable` | table | (N_units × F) per-unit feature table |
| `NetworkFeatureTable` | table | (1 × F) network-level feature table |
| `Connectivity` | struct | CCG / STTC / DDC results |
| `Bursts` | struct | Burst detection results |
| `CellTypeLabels` | double | (1 × N_units) 1=exc, 2=inh, `NaN`=unclassified |
| `Status` | struct | Tracks which analyses have been run |

## Default Parameters

Access all defaults with `RecordingProcessor.returnDefaultParams()`. Override by passing a partial struct to the constructor:

```matlab
proc = RecordingProcessor(sd, struct('QC', struct('FiringRate', [0.1 10])));
```

### QC parameters (`params.QC`)

| Field | Default | Description |
|---|---|---|
| `Amplitude` | `[20 1000]` | Accepted peak-to-trough amplitude range (µV) |
| `FiringRate` | `[0.01 10]` | Accepted firing rate range (Hz) |
| `RPV` | `0.02` | Max fraction of refractory period violations (0–1) |
| `RefractoryPeriod` | `0.002` | Refractory period (s) used for RPV computation |
| `Axon` | `0.8` | Waveform axon score threshold (units with score > threshold are excluded) |
| `Noise` | `14` | Noise amplitude threshold |
| `NoiseCutout` | `[-1 1]` | Time window (ms) around trough for noise check |
| `PowerCutoff` | `1.3` | Max high-frequency power ratio |
| `N_Units` | `10` | Minimum units per recording (recordings with fewer are dropped) |
| `MinSpikeCount` | `10` | Minimum spike count per unit |
| `KSLabel.Enable` | `false` | Filter by Kilosort label (`cluster_KSLabel.tsv`) |
| `KSLabel.AcceptedLabels` | `"good"` | Labels to keep (`"good"`, `"mua"`, or both) |
| `Bombcell.Enable` | `false` | Filter using Bombcell quality metrics (`cluster_bc_unitType.tsv`; missing table: warning, standard QC only) |
| `Bombcell.AcceptedTypes` | `[1, 2]` | Bombcell unit types to keep (1=good, 2=MUA) |

### Burst parameters (`params.Bursts`)

| Field | Default | Description |
|---|---|---|
| `Method` | `"ISIN"` | `"ISIN"` (legacy: N by Dunn index, valley threshold, tonic heuristic), `"logISI"`, or `"ISINJoint"` (`detectNetworkBursts`: N and threshold chosen jointly against circular-shift surrogates; distinct maxima of the significance map are levels of grouping (bursts, sub-events, superbursts / epochs); within the burst level the finest scale that keeps bursts whole, so all synchronous events recruiting `MinBurstUnits` units (significantly more than expected by chance) count; superbursts - trains of bursts without silence between them - are split into their bursts (`proc.Bursts.Superburst`, `proc.Bursts.Info.Superbursts`); returns no bursts for cultures without network bursts; diagnostics in `proc.Bursts.Info`, levels in `proc.Bursts.Info.Levels`) |
| `TonicMethod` | `"gini"` | `"ISINJoint"` only: units removed before detection (`"gini"`: not burstier than Poisson, adaptive bins of ~3 spikes; `"coupling"` population coupling; `"regular"` more regular than Poisson; `"coupling+regular"`; `"either"`, `"both"`, `"legacy"`, `"none"`) |
| `TonicMinRate` | `"median"` | `"ISINJoint"` only: units firing below this rate (Hz, or `"median"` of the recording's units) are never removed as tonic; low-rate units carry burst spikes rather than background |
| `ShapeLevel` | `0.25` | Burst onset / offset at this fraction of each burst's peak above baseline (`burstShape`); gives the `BurstWidth`, `BurstRise`, `BurstFall`, `BurstFWHM` features, comparable between detectors and bursts (the legacy `BurstDuration`, `RiseTime`, `FallTime` depend on the detector edges) |
| `MergeFactor` | `0.5` | Legacy methods: merge bursts closer than MergeFactor x threshold |
| `Binning` | `0.1` | Bin size (s) for burst statistics |

### Connectivity parameters

| Field | Default | Description |
|---|---|---|
| `Analyses.Connectivity` | `["CCG","STTC"]` | Enabled methods (`"CCG"`, `"STTC"`, `"DDC"`, `"FullCCG"`) |
| `CCG.BinSize` | `0.001` | Bin width (s) |
| `CCG.Duration` | `0.1` | One-sided lag (s) |
| `CCG.Alpha` | `0.001` | Significance threshold for edge detection |
| `CCG.InhibitoryCriterion` | `"strict"` | Inhibitory links: `"strict"` (causal trough below the pre-spike mean and the slow 10 ms-to-edge flank; no reverse artefact of strong excitation) or `"cellexplorer"` (CellExplorer `ce_MonoSynConvClick`, the previous DeePhys default; many false positives) |
| `STTC.MaxLag` | `0.05` | Synchrony window (s) |
| `STTC.JitterSD` | `0.1` | SD (s) of the spike jitter of the surrogates; must exceed `MaxLag` (2 × `MaxLag`) |
| `STTC.N_Surrogates` | `100` | Surrogate count for significance |
| `STTC.Percentile` | `95` | Significance percentile |
| `DDC.BinSize` | `0.001` | Bin width (s) |
| `DDC.N_Surrogates` | `100` | Surrogate count |
| `DDC.Percentile` | `95` | Significance percentile |
| `GraphFeatures.Types` | `"binary"` | Graph features from the unweighted graphs (`"binary"`), the weighted graphs (`"weighted"`, columns `<Feature>_<Method>_w`), or both |
| `GraphFeatures.SmallWorldnessIterations` | `100` | Randomised graphs for small-worldness |

Connectivity results per method (`proc.Connectivity.<Method>`):

| Method | Unweighted graph | Weighted graph |
|---|---|---|
| CCG | `.bd` (directed; -1 = inhibitory link) | — |
| STTC | `.bu` (undirected, surrogate-significant) | `.wu` (tiling coefficients) |
| DDC | `.bd` (directed, surrogate-significant) | `.wd` (signed DDC coefficients, all pairs), `.wd_sig` (significant links only) |

Weighted graph features use the absolute weights of the significant links, normalised to (0, 1]. For DDC, positive coefficients indicate excitatory and negative ones inhibitory influence; the sign is kept in `.wd` / `.wd_sig`.

## Notes

- `proc.save(path)` serialises to plain MATLAB types — no handle-object hacks needed.
- All `compute*` methods are idempotent: calling them a second time is a no-op. Exceptions: a step whose sub-analysis failed is marked `"failed"` and re-runs on the next call (its own columns are replaced, those of other steps kept); `computeConnectivity(methods)` adds methods not computed before; `computeUnitFeatures(groups)` / `computeParentFeatures(params)` on a finished step warn that the arguments are ignored (set the `Status` field to `"pending"` to recompute).
- Override parameters by passing a nested struct, e.g.: `RecordingProcessor(sd, struct('QC', struct('FiringRate', [0.1 10])))`.
- Default connectivity methods are `["CCG","STTC"]`. DDC and FullCCG are available but not enabled by default.
- `Status` tracks each analysis step with tri-state values: `"pending"`, `"done"`, or `"failed"`.
