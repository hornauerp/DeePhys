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
| `RecordingID` | string | Deterministic hash of `InputPath` |
| `InputPath` | string | Absolute path to the Kilosort output directory |
| `ParentPath` | string | Path to the parent (concatenated) Kilosort directory; `""` if none |
| `Metadata` | struct | Recording metadata (ChipID, PlatingDate, DIV, Mutation, …) |
| `SamplingRate` | double | Samples per second |
| `Duration` | double | Recording duration in seconds |
| `ElectrodeCoordinates` | double | (N_channels × 2) XY electrode positions |
| `SpikeTimes` | double | (N_spikes × 1) spike times in seconds |
| `SpikeUnits` | double | (N_spikes × 1) 1-indexed template assignment |
| `TemplateWaveforms` | double | (N_templates × N_samples × N_channels) |

## Notes

- `RecordingID` is derived from `InputPath` so it is stable across sessions.
- `DIV` is computed automatically from `PlatingDate` + `RecordingDate` if not provided in metadata.
- If a sibling `qc_output/` folder exists alongside `InputPath`, it is auto-detected as `ParentPath` (split-recording support).
