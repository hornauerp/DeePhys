# Welcome to *DeePhys*
The package for **Deep electrophysiological phenotype characterization**:

<img src="https://github.com/user-attachments/assets/dfa63a20-211d-43e7-85ff-fc12b16bb68d" alt="Analysis schematic" style="width:500px;"/>

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
Currently *DeePhys* is only available on MATLAB, so a recent MATLAB installation (>2019b) is required.

## Installation
The package is ready-to-use right after cloning. 

## Usage
Code requires spikesorted data in the [phy format](https://github.com/cortex-lab/phy). For help with spikesorting check out the [Spikeinterface package](https://spikeinterface.readthedocs.io/en/latest/). 

We provide [several tutorials](/Tutorials) and a [dataset (to be added)](https://zenodo.org/) from our most recent paper (to be published) to facilitate picking up the analysis workflow:
The analysis pipeline can be subdivided into the [Feature Extraction](/Tutorials/1_FeatureExtraction), which performs a basic quality control and extracts features from the spike-sorted data, and the [Phenotype Generation](/Tutorials/2_PhenotypeGeneration) that runs the core analyses. The [BasicFeatureExtraction tutorial](/Tutorials/1_FeatureExtraction/BasicFeatureExtraction.m) shows how to set up the feature extraction, while the [FeatureInspection tutorial](/Tutorials/1_FeatureExtraction/FeatureInspection.m) shows how to inspect the analysis results and adjust them if necessary. After that, you can follow the [NetworkPhenotype tutorial](/Tutorials/2_PhenotypeGeneration/NetworkPhenotype.mat) or the [SingleCellPhenotype tutorial](/Tutorials/2_PhenotypeGeneration/SingleCellPhenotype.mat), which showcase how to perform the analyses shown in the original [DeePhys paper](https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(23)00501-5) and our latest publication on it (currently in print).

The [dataset of the original paper](https://doi.org/10.5281/zenodo.7876370) and the [old tutorials](/Tutorials/Old) are still available, however, there will probably be compatibility problems due to frequent updates. They might, however, still be useful to understand specific analysis workflows.

## Citation
If you find this package helpful or used in your analyses, please cite the [DeePhys paper](https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(23)00501-5) and link to this github repository.

## Dependencies
This package uses several packages/toolboxes:
- the `readNPY` function provided by the [npy-matlab package](https://github.com/kwikteam/npy-matlab)
- the `CCG` function provided by the [FMAToolbox](https://github.com/michael-zugaro/FMAToolbox)
- the [`othercolor` function](https://ch.mathworks.com/matlabcentral/fileexchange/30564-othercolor).
- the [`catch22` toolbox](https://github.com/DynamicsAndNeuralSystems/catch22) as published [here](https://doi.org/10.1007/s10618-019-00647-x)
- the ISIN burst detection algorithm as published [here](https://www.frontiersin.org/articles/10.3389/fncom.2013.00193/full)
- the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/home)

## Issues?
If you face any problems or bugs, or have ideas for additions to this package please open an [issue](https://github.com/hornauerp/DeePhys/issues). 
