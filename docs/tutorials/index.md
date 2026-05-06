# Tutorials

The tutorials cover the complete DeePhys analysis workflow from raw Kilosort output to publishable results. Each tutorial is a self-contained MATLAB script in the `Tutorials/` folder.

**Work through them in order** — each tutorial builds on the outputs of the previous one.

| # | Tutorial | What you'll learn |
|---|---|---|
| 1 | [Data Processing](01_data_processing.md) | Load spike data, run QC, extract features, build a FeatureStore |
| 2 | [Feature Exploration](02_feature_exploration.md) | Inspect feature matrices, check distributions, prepare ML inputs |
| 3 | [Phenotype Analysis](03_phenotype_analysis.md) | Classification, dimensionality reduction, regression |
| 4 | [Cell Type Classification](04_cell_type_classification.md) | Unsupervised E/I labelling via UMAP + graph label propagation |
| 5 | [E/I Analysis](05_ei_analysis.md) | Burst detection, E/I activity traces, unit–population correlations |
| 6 | [Transfer Learning](06_transfer_learning.md) | Classify units from external recordings or patch-clamp labelled data |
| 7 | [Legacy Migration](07_legacy_migration.md) | Convert old MEArecording / RecordingGroup objects to the new API |

## Prerequisites

All tutorials assume:
- DeePhys is on the MATLAB path (`startup` run from the repo root)
- Spike-sorted data in the phy format (see [Getting Started](../getting_started.md))
- Each tutorial notes additional prerequisites (e.g., Tutorial 5 requires Tutorial 4 output)
