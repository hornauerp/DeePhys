# Welcome to *DeePhys*
The package for **Deep electrophysiological phenotype characterization**:

<img src="https://github.com/hornauerp/EphysDopa/blob/2e43777e3fd1fe0e7467c4b3bf0aa25afb88b602/Figures/EphysDopaSchematic_v2202222.png" alt="Analysis schematic" style="width:500px;"/>

Created with [BioRender](BioRender.com)

## Overview
*DeePhys* was created to facilitate the analysis of extracellular recordings of neuronal cultures using high-density microelectrode arrays (HD-MEAs). *DeePhys* allows users to easily:
- Extract electrophysiological features from spikesorted HD-MEA recordings
- Visualize differential developmental trajectories 
- Apply machine learning algorithms to classify different conditions
- Obtain biomarkers predictive of the respective condition
- Evaluate the effect of treatments
- Dissect heterogeneous cell populations/cultures on the single-cell level

## Requirements
Currently *DeePhys* is only available on MATLAB, so a recent MATLAB installation (>2019b) is required. We plan on expanding *DeePhys* to Python in the near future.

## Installation
The package is ready-to-use right after cloning. As the download via git-lfs is heavily limited, please download the spike sorted data for the tutorial [here](https://zenodo.org/records/10635138) and replace the folder SortingExamples with the unzipped version from there. Preprocessed data can be downloaded [here](https://zenodo.org/records/7876371).

## Usage
Code requires spikesorted data in the [phy format](https://github.com/cortex-lab/phy). For help with spikesorting check out the [Spikeinterface package](https://spikeinterface.readthedocs.io/en/latest/). 

The analysis pipeline is subdivided into the following modules (links to the tutorials):
- [Feature Extraction](/Tutorials/1_FeatureExtraction)
- [Phenotype Generation](/Tutorials/2_PhenotypeGeneration)
- [Treatment Evaluation](/Tutorials/3_TreatmentEvaluation)
- [Single Cell Analysis](/Tutorials/4_SingleCellEvaluation)


## Citation
The *DeePhys* package was first published on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.03.31.486582v1), but was since heavily updated and is no longer compatible to the [prior version](https://github.com/hornauerp/EphysDopa).

## Disclaimer
This package uses several packages/toolboxes:
- the `readNPY` function provided by the [npy-matlab package](https://github.com/kwikteam/npy-matlab)
- the `CCG` function provided by the [FMAToolbox](https://github.com/michael-zugaro/FMAToolbox)
- the [`othercolor` function](https://ch.mathworks.com/matlabcentral/fileexchange/30564-othercolor).
- the [`catch22` toolbox](https://github.com/DynamicsAndNeuralSystems/catch22) as published [here](https://doi.org/10.1007/s10618-019-00647-x)
- the ISIN burst detection algorithm as published [here](https://www.frontiersin.org/articles/10.3389/fncom.2013.00193/full)
- the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/home)
