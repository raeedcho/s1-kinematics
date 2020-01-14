# S1 kinematic encoding analysis

This repository contains analysis code relevant to the manuscript "Area 2 of primary somatosensory cortex encodes kinematics of the whole arm", which can be found [here](https://www.biorxiv.org/content/10.1101/643205v2).

## How to get the code

The best way to get the code is to clone this repository. When cloning the repository, make sure to use the -r flag to also clone the submodules in `lib/`.

```
git clone -r https://github.com/raeedcho/s1-kinematics.git
```

I've tested this code most recently on MATLAB 2019b with the following toolboxes installed:

- (Coming soon)

## Execution instructions

The main scripts used to generate the figures from the behavioral and neural data (found [here](https://doi.org/10.5061/dryad.nk98sf7q7)) are `plotRFMaps.m`, `twoworkspace_analysis.m`, and `actpas_analysis.m`.

Before running the main analysis scripts, add `lib/` and all its subdirectories to the MATLAB path. Download the data to disk and modify the `dataroot` at the top of the scripts to point to the unzipped data directory. Then simply run the scripts. Note that because these scripts use extensive cross-validation, the analysis will take a long time to run.
