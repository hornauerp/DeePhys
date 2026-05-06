# Welcome to *DeePhys*
The package for **Deep electrophysiological phenotype characterization**:

<img src="https://github.com/user-attachments/assets/dfa63a20-211d-43e7-85ff-fc12b16bb68d" alt="Analysis schematic" style="width:500px;"/>

Created with [BioRender](https://www.biorender.com)

**[Full documentation](https://deephys.readthedocs.io)**

## Overview
*DeePhys* was created to facilitate the analysis of extracellular recordings of neuronal cultures using high-density microelectrode arrays (HD-MEAs). *DeePhys* allows users to easily:
- Extract electrophysiological features from spikesorted HD-MEA recordings
- Visualize differential developmental trajectories 
- Apply machine learning algorithms to classify different conditions
- Obtain biomarkers predictive of the respective condition
- Evaluate the effect of treatments
- Dissect heterogeneous cell populations/cultures on the single-cell level

## Requirements
Currently *DeePhys* is only available on MATLAB, so a recent MATLAB installation (>2019b) is required.

## Installation
The package is ready-to-use right after cloning. Run `startup.m` from the repo root to add all paths.

## Usage
Code requires spikesorted data in the [phy format](https://github.com/cortex-lab/phy). For help with spikesorting check out the [SpikeInterface package](https://spikeinterface.readthedocs.io/en/latest/).

We provide [seven tutorials](Tutorials/) covering the full analysis workflow:

1. [Data Processing](Tutorials/1_DataProcessing/DataProcessing.m) — load Kilosort output, run QC, extract features
2. [Feature Exploration](Tutorials/2_FeatureExploration/FeatureExploration.m) — inspect and prepare feature matrices
3. [Phenotype Analysis](Tutorials/3_PhenotypeAnalysis/PhenotypeAnalysis.m) — classification, dimensionality reduction, regression
4. [Cell Type Classification](Tutorials/4_CellTypeClassification/CellTypeClassification.m) — excitatory/inhibitory labelling
5. [E/I Analysis](Tutorials/5_EIAnalysis/EIAnalysis.m) — network burst detection and E/I quantification
6. [Transfer Learning](Tutorials/6_TransferLearning/TransferLearning.m) — external data interoperability
7. [Legacy Migration](Tutorials/7_LegacyMigration/LegacyMigration.m) — migrate old MEArecording objects to the new API

A dataset from our most recent paper will accompany the tutorials (link to be added upon publication). The [dataset of the original paper](https://doi.org/10.5281/zenodo.7876370) is still available.

## Citation
If you find this package helpful or used in your analyses, please cite the [DeePhys paper](https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(23)00501-5) and link to this GitHub repository.

## Dependencies
This package uses several packages/toolboxes, all of which are **bundled in the repository** (under `Functions/` and `Toolboxes/`) and added to the MATLAB path automatically by `startup.m`. No separate installation is required.

- the `readNPY` function provided by the [npy-matlab package](https://github.com/kwikteam/npy-matlab)
- the `CCG` function provided by the [FMAToolbox](https://github.com/michael-zugaro/FMAToolbox)
- the [`othercolor` function](https://ch.mathworks.com/matlabcentral/fileexchange/30564-othercolor)
- the [`catch22` toolbox](https://github.com/DynamicsAndNeuralSystems/catch22) as published [here](https://doi.org/10.1007/s10618-019-00647-x)
- the ISIN burst detection algorithm as published [here](https://www.frontiersin.org/articles/10.3389/fncom.2013.00193/full)
- the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/home)
- the [UMAP for MATLAB](https://github.com/MalcolmSlaney/StanfordPosner) toolbox

## Issues?
If you face any problems or bugs, or have ideas for additions to this package please open an [issue](https://github.com/hornauerp/DeePhys/issues).
