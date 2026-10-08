# Getting Started

## Requirements

- MATLAB R2019b or later
- A spike-sorted dataset in the [phy format](https://github.com/cortex-lab/phy) (produced by [Kilosort](https://github.com/MouseLand/Kilosort) or compatible sorters via [SpikeInterface](https://spikeinterface.readthedocs.io))
- **Statistics and Machine Learning Toolbox** — required for classification (`fitcensemble`, `fitcsvm`, `fitcknn`, `fitcnb`), regression (`fitrensemble`, `fitrsvm`), and statistical tests (`lillietest`)

All third-party dependencies (readNPY, CCG, othercolor, catch22, BCT, UMAP) are bundled in `Functions/` and `Toolboxes/` and added to the path by `startup.m`. No separate installation is required.

## Installation

Clone the repository and run the startup script:

```bash
git clone https://github.com/hornauerp/DeePhys.git
```

```matlab
% In MATLAB, from the repo root:
startup
```

`startup.m` adds all necessary paths (Classes, Functions, Toolboxes) to the MATLAB path for the current session. Add this call to your `startup.m` in the MATLAB preferences folder to make it permanent.

## Input data format

DeePhys reads Kilosort output directories, or the output of any other spike sorter exported to phy format (e.g. SpikeInterface `export_to_phy`). The following files are used:

| File | Content |
|---|---|
| `spike_times.npy` | Spike times in samples (required) |
| `spike_templates.npy` | Template assignment per spike (Kilosort). Without it, `spike_clusters.npy` is used, whose cluster indices must then address `templates.npy` (one template per unit) |
| `templates.npy` | Template waveforms (required) |
| `channel_positions.npy` | Electrode XY coordinates (optional; needed for the spatial features) |
| `params.py` | Must contain `sample_rate`; optionally `n_samples` for an accurate duration (e.g. `sample_rate = 30000.0`; any spacing) |

Label tables are optional. With `QC.KSLabel.Enable = true`, `cluster_KSLabel.tsv` provides `good`/`mua` labels; with `QC.Bombcell.Enable = true`, Bombcell's `cluster_bc_unitType.tsv` does. If an enabled table is missing (e.g. Bombcell was not run), a warning is issued and only the standard QC applies. Units are defined by the templates: merges or splits from manual phy curation (`spike_clusters.npy` alongside `spike_templates.npy`) are not used.

The `RecordingID` of a sorting is a hash of its spike data, so the same sorting gets the same ID on every machine and mount point.

## Metadata

Each recording is described by a MATLAB struct with scalar fields. Any fields included here are stored in the `FeatureStore` and become available for subsetting and as ML labels.

```matlab
metadata = struct();
metadata.ChipID        = 'Chip001';
metadata.PlatingDate   = '2024-01-15';
metadata.RecordingDate = '2024-02-05';  % DIV computed automatically if absent
metadata.DIV           = 21;
metadata.Mutation      = 'WT';
metadata.Concentration = 0;
```

`DIV` (days in vitro) is computed automatically from `PlatingDate` and `RecordingDate` if not provided.

## Next steps

Work through the tutorials in order for a complete walkthrough:

1. [Data Processing](tutorials/01_data_processing.md) — `SpikeData` → `RecordingProcessor` → `FeatureStore`
2. [Feature Exploration](tutorials/02_feature_exploration.md) — inspect feature matrices
3. [Phenotype Analysis](tutorials/03_phenotype_analysis.md) — ML classification and dimensionality reduction
4. [Cell Type Classification](tutorials/04_cell_type_classification.md) — E/I labelling
5. [E/I Analysis](tutorials/05_ei_analysis.md) — network burst analysis
6. [Transfer Learning](tutorials/06_transfer_learning.md) — external data interoperability
7. [Legacy Migration](tutorials/07_legacy_migration.md) — upgrade from old MEArecording objects
