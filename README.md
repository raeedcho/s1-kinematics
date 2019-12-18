# S1 kinematic encoding analysis

This repository contains analysis code relevant to the manuscript "Area 2 of primary somatosensory cortex encodes kinematics of the whole arm", which can be found [here](https://www.biorxiv.org/content/10.1101/643205v2).

## Execution instructions

The main scripts used to generate the figures from the behavioral and neural data (found [here]) are `twoworkspace_analysis.m` and `actpas_analysis.m`. I've tested the code most recently on MATLAB 2019b with most of the official toolboxes installed, but it should work on earlier versions as well.

When cloning the repository, make sure to use the `-r` flag to also clone the submodules in `lib/`. Before running the main analysis scripts, add `lib/` and all its subdirectories to the MATLAB path. Download the data to disk and alter the `datadir` at the top of the scripts accordingly. Then simply run the scripts. Note that because these scripts use extensive cross-validation, the analysis will take a long time to run.

# MIT License

Copyright 2019 Raeed Chowdhury

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
