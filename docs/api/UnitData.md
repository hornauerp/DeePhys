# UnitData

Per-unit raw data container. **Value class**, no back-reference to the parent recording. References its parent by `RecordingID` (a string key), not an object handle — eliminating the circular-reference problem that required `saveobj`/`loadobj` hacks in the old `Unit` ↔ `MEArecording` design.

## Example

```matlab
% UnitData objects are created internally by RecordingProcessor.runQC()
% and stored in proc.Units. You rarely construct them manually.

proc.runQC();
ud = proc.Units(1);

fprintf('UnitID: %s, N spikes: %d, KS label: %s\n', ...
        ud.UnitID, numel(ud.SpikeTimes), ud.KSLabel);

% Legacy migration
ud = UnitData.fromLegacyUnit(old_unit_obj, recording_id);
```

## Constructor

`ud = UnitData(template_id, ref_electrode, spike_times, ref_waveform, recording_id, sampling_rate, duration)`

All arguments are optional — call with no arguments for array pre-allocation.

`ud = UnitData.fromLegacyUnit(unit_obj, recording_id)` — migrate from old `Unit` handle object.

`ud = UnitData.fromStruct(s)` — reconstruct from saved plain struct.

## Properties

| Name | Type | Description |
|---|---|---|
| `UnitID` | string | Stable deterministic ID (hash of ChipID + PlatingDate + TemplateID + RefElectrode) |
| `RecordingID` | string | Links to parent `SpikeData.RecordingID` |
| `TemplateID` | double | 1-indexed template number from Kilosort |
| `ReferenceElectrode` | double | Channel index of the peak-amplitude electrode |
| `SpikeTimes` | double | (N_spikes × 1) spike times in seconds |
| `ReferenceWaveform` | double | (N_samples × 1) peak-channel waveform |
| `SamplingRate` | double | Copied from parent recording |
| `RecordingDuration` | double | Copied from parent recording (seconds) |
| `KSLabel` | string | Kilosort label (`"good"` / `"mua"`); `""` if not available |
| `BombcellType` | double | BC classification: 0=noise, 1=good, 2=MUA, `NaN`=not run |
| `FullACG` | double | Full-recording autocorrelogram |
| `FullACGBinSize` | double | Bin size (s) of `FullACG` |
| `FullACGLag` | double | One-sided lag (s) of `FullACG` |
