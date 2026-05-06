# DeePhys

**Deep electrophysiological phenotype characterization** — a MATLAB toolbox for analysing extracellular recordings from neuronal cultures on high-density microelectrode arrays (HD-MEAs).

![Analysis schematic](https://github.com/user-attachments/assets/dfa63a20-211d-43e7-85ff-fc12b16bb68d)

## What DeePhys does

Starting from spike-sorted data in the [phy format](https://github.com/cortex-lab/phy), DeePhys lets you:

- **Extract** electrophysiological features at the single-unit and network level
- **Visualize** developmental trajectories across recordings and conditions (plotting utilities available in `EIAnalyzer` and `Functions/`; not yet fully documented)
- **Classify** experimental conditions using machine learning
- **Identify** condition-predictive biomarkers
- **Evaluate** treatment effects
- **Dissect** heterogeneous cell populations at single-cell resolution
- **Label** excitatory and inhibitory neurons without patch-clamp

## Quick start

```matlab
% Run once from the repo root to add all paths
startup

% Load one Kilosort output
sd = SpikeData.fromKilosort('/path/to/ks_output', metadata);

% Process and extract features
proc = RecordingProcessor(sd);
proc.runQC();
proc.computeUnitFeatures();
proc.computeNetworkFeatures();
proc.save('/path/to/save');

% Assemble multi-recording feature store
fs = FeatureStore.fromProcessors([proc1, proc2, proc3]);
```

See [Getting Started](getting_started.md) for installation and prerequisites, then work through the [Tutorials](tutorials/index.md) for the full workflow.

## Citation

If DeePhys contributed to your work, please cite:

> Hornauer et al. (2023). *DeePhys* — deep electrophysiological phenotyping of neuronal cultures.
> [Stem Cell Reports](https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(23)00501-5)

## Links

- [GitHub repository](https://github.com/hornauerp/DeePhys)
- [Zenodo dataset (original paper)](https://doi.org/10.5281/zenodo.7876370)
- [Issue tracker](https://github.com/hornauerp/DeePhys/issues)
