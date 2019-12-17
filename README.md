# S1 kinematic encoding analysis

This repository contains analysis code relevant to the manuscript "Area 2 of primary somatosensory cortex encodes kinematics of the whole arm", which can be found [here](https://www.biorxiv.org/content/10.1101/643205v2).

## Execution instructions

The main scripts used to generate the figures from the behavioral and neural data (found [here]) are `twoworkspace_analysis.m` and `actpas_analysis.m`. I've tested the code most recently on MATLAB 2019b with most of the official toolboxes installed, but it should work on earlier versions as well.

When cloning the repository, make sure to use the `-r` flag to also clone the submodules in `lib/`. Before running the main analysis scripts, add `lib/` and all its subdirectories to the MATLAB path. Download the data to disk and alter the `datadir` at the top of the scripts accordingly. Then simply run the scripts. Note that because these scripts use extensive cross-validation, the analysis will take a long time to run.
