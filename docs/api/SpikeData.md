# SpikeData

Pure data container for Kilosort spike-sorted output. **Value class** — no computation, no side effects. Trivially serializable because it holds only plain arrays and structs.

## Example

```matlab
metadata = struct('ChipID', 'Chip001', 'DIV', 21, 'Mutation', 'WT');
sd = SpikeData.fromKilosort('/path/to/ks_output', metadata);

fprintf('Duration: %.1f s, N units: %d\n', sd.Duration, ...
        numel(unique(sd.SpikeUnits)));
```

## Constructor (static)

`sd = SpikeData.fromKilosort(input_path)` — load from Kilosort output directory.

`sd = SpikeData.fromKilosort(input_path, metadata)` — same, with metadata struct.

`sd = SpikeData.fromStruct(s)` — reconstruct from a previously saved plain struct.

`sd = SpikeData.fromLegacyMEArecording(mearec)` — extract SpikeData from an old `MEArecording` object (used internally by the legacy migration pipeline).

## Properties

| Name | Type | Description |
|---|---|---|
| `RecordingID` | string | Deterministic hash of the spike data (the same on every machine and mount point); legacy objects: hash of `InputPath` |
| `SortingID` | string | Sorting the recording belongs to: segments of one sorting (`segment_N` folders, `split_sortings` windows, or an explicit `ParentInputPath`) share it |
| `InputPath` | string | Absolute path to the Kilosort output directory |
| `ParentPath` | string | Path to the parent (concatenated) Kilosort directory; `""` if none |
| `Metadata` | struct | Recording metadata (ChipID, PlatingDate, DIV, Mutation, …) |
| `SamplingRate` | double | Samples per second |
| `Duration` | double | Recording duration in seconds |
| `ElectrodeCoordinates` | double | (N_channels × 2) XY electrode positions |
| `SpikeTimes` | double | (N_spikes × 1) spike times in seconds |
| `SpikeUnits` | double | (N_spikes × 1) 1-indexed template assignment |
| `TemplateWaveforms` | double | (N_templates × N_samples × N_channels) |
| `TemplateAmplitudes` | double | (N_templates × 1) trough amplitude of each template in recording units (ADC): unwhitened template × median spike scaling / `scaleproc`. `[]` when unavailable (see below) |
| `GainToUV` | double | µV per recording unit from SpikeInterface's `spikeinterface_recording.json` (`gain_to_uV`); `NaN` if unknown |

## Notes

- Template amplitudes need `whitening_mat_inv.npy` and `amplitudes.npy` of a Kilosort 1–3 sorting. Kilosort fits the whitened data scaled by `ops.scaleproc` (default 200, or `scaleproc` from `spikeinterface_params.json`), so a spike in recording units is `template × whitening_mat_inv × scaling / scaleproc`. Other sorters (by `sorter_name` in `spikeinterface_params.json`) and phy exports with an identity `whitening_mat_inv` and no sorter metadata get `[]`, since their scaling convention differs.
- `RecordingID` is derived from the spike data, so it is stable across sessions and paths; units of segments with one `SortingID` get the same `UnitID`.
- `DIV` is computed automatically from `PlatingDate` + `RecordingDate` if not provided in metadata.
- If a sibling `qc_output/` folder exists alongside `InputPath`, it is auto-detected as `ParentPath` (split-recording support).
