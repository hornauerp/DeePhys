# EIAnalyzer

Excitatory / inhibitory population activity of the cultures classified by a `CellTypeClassifier`. **Handle class.**

There are two analysis paths:

- **Legacy path** (`computeActivity` → `detectBursts` → `extractBurstCutouts` → `normalizeBurstCutouts`). These are the published methods and are kept unchanged for reproducibility. `detectBursts` labels burst bins by I/E < `MaxIERatio`, so any statement about the I/E ratio during those bursts is circular.
- **Burst-event path.**
  - Bursts come from `detectNetworkBursts`, using all units and ignoring the labels.
  - Every metric is computed per recording and tested against random "I" sets matched in firing rate.
  - The computations are standalone functions (`burstEIMetrics`, `burstUnitTiming`, `burstRhythm`) that the methods call for each recording.

## Example

```matlab
ctc.identifyResponsiveUnits({'Construct', 1});   % cultures with an activating DREADD in interneurons
ctc.discardSegmentUnits();
ctc.generateTrainLabels();
ctc.classifyUnits();

eia = EIAnalyzer(ctc, struct('Labels', struct('Source', "confident", 'MinProbability', 0.8)));
eia.computeNetworkBursts();
eia.computeBurstEI();
eia.computeUnitTiming();
eia.computeBurstRhythm();
[T, G] = eia.summarizeBurstEI(GroupBy=["DIV", "Construct"]);
eia.PlotBurstEI(1);

% ground truth of the perturbation instead of predicted labels
eia2 = EIAnalyzer(ctc, struct('Labels', struct('Source', "responders", 'ResponderFilter', {{'Construct', 1}})));
```

## Labels (`Parameters.Labels.Source`)

| Source | E / I populations |
|---|---|
| `"graph"` (default) | `UnitLabels` of `classifyUnits` |
| `"composition"` | `Composition.Labels` (`estimateComposition`) |
| `"confident"` | `confidentLabels`: units with P(inhibitory) ≥ `MinProbability` or ≤ 1 − `MinProbability`. All other units are left out of both populations. |
| `"responders"` | I = units whose `ResponsiveUnitDirection` equals `ResponderDirection` (default `"increase"`). E = the other units of the same cultures. The cultures are those matching `ResponderFilter` (e.g. `{'Construct', 1}`), or, if it is empty, those with at least one responder. |

Unlabelled units still take part in burst detection.

## Burst-event methods

| Method | Result |
|---|---|
| `computeNetworkBursts()` | `NetworkBursts(r)`: `RecordingID`, `Culture`, `UnitRows` (FeatureStore rows), `Duration`, `Bursts`, `Info`. Detector options are in `Parameters.NetworkBursts.Options`. |
| `computeBurstEI()` | `BurstEI(r).Result` (see `burstEIMetrics` below) |
| `computeUnitTiming()` | `UnitTiming`: one table row per unit (see `burstUnitTiming` below) |
| `computeBurstRhythm()` | `BurstRhythm(r).Result` (see `burstRhythm` below) |
| `summarizeBurstEI(GroupBy=...)` | `T`: one row per recording with metadata, every metric with its null (`_Null`, `_PLower`, `_PUpper`), unit-timing and rhythm summaries. `G`: group medians, plus the share of recordings with `IEContrastIn_PLower` < 0.05. |
| `PlotBurstEI(r)` | E and I aligned at onset, peak and offset, and I/E with the 95 % band of the null |

## `burstEIMetrics(spike_times, spike_units, labels, bursts, duration, ...)`

### Epoch rates

Each unit's rate is computed in four epochs, taken from each burst's own onset and offset:

| Epoch | Window |
|---|---|
| In | [Onset, Offset] |
| Pre | `PreWindow` before the onset |
| Post | `PostWindow` after the offset |
| Between | outside every burst ± `Margin` |

The function also builds each unit's mean rate aligned at onset, peak and offset.

Every population, and every null draw, is an average over its units. I/E is taken after averaging; bins with E = 0 are not dropped.

### Metrics

| Metric | Meaning |
|---|---|
| `IFoldIn`, `EFoldIn`, `IEContrastIn` | In / Between for I and for E, and their ratio. I/E lower during bursts gives `IEContrastIn` < 1. |
| `...Pre`, `...Post` | The same before onset and after offset |
| `ILeadTime`, `IEMaxTime` | The I maximum and the I/E maximum of the onset-aligned traces, within `LeadWindow` |
| `IEMinTime` | The I/E minimum of the peak-aligned traces |
| `ITerminal`, `ITerminalVsIn` | The offset-aligned I maximum relative to the Between rate and to the In rate |
| `IAfterPeak`, `IAfterPeakTime`, `EAfterPeak` | After-peak: the height of the largest maximum within `AfterWindow` after the burst peak, above the lowest point between the burst peak and it. It is relative to the Between rate, and is 0 if the rate only decays. If both I and E show one, it is network reverberation rather than I alone. |
| `IAfterRel`, `EAfterRel` | The after-peak height as a fraction of the range (maximum − minimum) of the population's peak-aligned trace, from 0 to 1. Unlike the Between-normalised values, this is comparable between I and E. |
| `AntiPhase` | Bursts inside a superburst, after its first: (I in the gaps between them / I inside them) / (the same for E) |

Folds are NaN when a population never fires between bursts, or when the Between epoch is shorter than `MinBetween` (10 s); recordings with frequent bursts leave no interburst reference.

### Null and output

- **Null:** `NNull` random "I" sets. Each has the same number of units per firing-rate stratum (`NullStrata` quantiles), drawn from the labelled units.
- **Output:** `Metrics` table with `Value`, `NullMedian`, `NullLow`/`NullHigh` (2.5 / 97.5 %), `PLower` (share of null draws ≤ the value) and `PUpper`. Also `Rates`, `Traces` (with the null band of I/E) and `Units`.

## `burstUnitTiming(spike_times, spike_units, bursts, duration, ...)`

One row per unit:

| Column | Meaning |
|---|---|
| `Latency` | Median first-spike time relative to onset, within [Onset − `PreWindow`, Offset] |
| `Participation` | Share of bursts with a spike in [Onset, Offset] |
| `PreRatio` | Rate before onset ÷ the unit's mean rate |
| `InRatio` | Rate in the burst ÷ mean rate |
| `PostRatio` | Rate after offset ÷ mean rate |

These measures are defined per unit, so labels (soft, confident or ground truth) can be applied afterwards.

## `burstRhythm(spike_times, spike_units, groups, bursts, ...)`

- **`Locked`** — the burst-locked sub-peaks.
  - The fast part is the burst-averaged rate (2 ms bins, aligned at the peak) minus a 25 ms Gaussian smooth.
  - Reported: its split-half reliability (odd vs even bursts), its P against circular shifts, and its frequency (the highest autocorrelation peak at 25–200 ms lag).
  - The window must hold at least two cycles, so a short burst core only admits higher frequencies.
- **`Shared`** — co-fluctuations between disjoint unit sets that are not locked to the burst.
  - Measured as cross-spectral coherency minus a shift predictor that pairs bursts k and k + j.
  - Pairs are either random halves of one group or two groups (e.g. E vs I). For two groups, `Lag` is the lag of the second group.
  - Reported: the peak against the shifted pairings, and band means.
