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
The package is ready-to-use right after cloning. 

## Usage
Code requires spikesorted data in the [phy format](https://github.com/cortex-lab/phy). For help with spikesorting check out the [Spikeinterface package](https://spikeinterface.readthedocs.io/en/latest/). 

Different parts of the analysis are subdevided into different analysis scripts:


## Citation
This package was published in "*DeePhys*, a machine learning-driven platform for electrophysiological phenotype screening" and additionally contains code to replicate the figures used in the publication.

## Disclaimer
This package uses several packages/toolboxes:
- the `readNPY` function provided by the [npy-matlab package](https://github.com/kwikteam/npy-matlab)
- the `CCG` function provided by the [FMAToolbox](https://github.com/michael-zugaro/FMAToolbox)
- the [`othercolor` function](https://ch.mathworks.com/matlabcentral/fileexchange/30564-othercolor).
- the [`catch22` toolbox](https://github.com/DynamicsAndNeuralSystems/catch22) as published [here](https://doi.org/10.1007/s10618-019-00647-x)
- the ISIN burst detection algorithm as published [here](https://www.frontiersin.org/articles/10.3389/fncom.2013.00193/full)
- the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/home)