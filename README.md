# S1 kinematic encoding analysis

This repository contains analysis code relevant to the manuscript "Area 2 of primary somatosensory cortex encodes kinematics of the whole arm", which can be found [here](https://doi.org/10.7554/eLife.48198).

## How to get the code

### Using `git`

The best way to get the code is to clone this repository with `git`. When cloning the repository, make sure to use the -r flag to also clone the submodules in `lib/`.

```
git clone -r https://github.com/raeedcho/s1-kinematics.git
```

### Downloading the release

You can also download the code from [here](https://github.com/raeedcho/s1-kinematics/releases/download/v1.0/s1-kinematics-v1.0.zip).

## Execution instructions

The main scripts used to generate the figures from the behavioral and neural data (found [here](https://doi.org/10.5061/dryad.nk98sf7q7)) are `plotRFMaps.m`, `twoworkspace_analysis.m`, and `actpas_analysis.m`.

Before running the main analysis scripts, add `lib/` and all its subdirectories to the MATLAB path. Download the data to disk and modify the `dataroot` at the top of the scripts to point to the unzipped data directory. Then simply run the scripts. Note that because these scripts use extensive cross-validation, the analysis will take a long time to run.

## Environmental requirements

I've tested this code most recently on MATLAB 2019a with the following toolboxes installed (note that many of these toolboxes are not actually required for the code to run--this was just my most recent MATLAB environment):

- Simulink                                              Version 9.3         (R2019a)
- Bioinformatics Toolbox                                Version 4.12        (R2019a)
- Communications Toolbox                                Version 7.1         (R2019a)
- Computer Vision Toolbox                               Version 9.0         (R2019a)
- Control System Toolbox                                Version 10.6        (R2019a)
- Curve Fitting Toolbox                                 Version 3.5.9       (R2019a)
- DSP System Toolbox                                    Version 9.8         (R2019a)
- Database Toolbox                                      Version 9.1         (R2019a)
- Datafeed Toolbox                                      Version 5.8.1       (R2019a)
- Deep Learning Toolbox                                 Version 12.1        (R2019a)
- Embedded Coder                                        Version 7.2         (R2019a)
- Filter Design HDL Coder                               Version 3.1.5       (R2019a)
- Fixed-Point Designer                                  Version 6.3         (R2019a)
- Fuzzy Logic Toolbox                                   Version 2.5         (R2019a)
- GPU Coder                                             Version 1.3         (R2019a)
- Global Optimization Toolbox                           Version 4.1         (R2019a)
- HDL Coder                                             Version 3.14        (R2019a)
- HDL Verifier                                          Version 5.6         (R2019a)
- Image Acquisition Toolbox                             Version 6.0         (R2019a)
- Image Processing Toolbox                              Version 10.4        (R2019a)
- Instrument Control Toolbox                            Version 4.0         (R2019a)
- MATLAB Coder                                          Version 4.2         (R2019a)
- MATLAB Compiler                                       Version 7.0.1       (R2019a)
- MATLAB Compiler SDK                                   Version 6.6.1       (R2019a)
- MATLAB Report Generator                               Version 5.6         (R2019a)
- Mapping Toolbox                                       Version 4.8         (R2019a)
- Model Predictive Control Toolbox                      Version 6.3         (R2019a)
- Optimization Toolbox                                  Version 8.3         (R2019a)
- Parallel Computing Toolbox                            Version 7.0         (R2019a)
- Partial Differential Equation Toolbox                 Version 3.2         (R2019a)
- Phased Array System Toolbox                           Version 4.1         (R2019a)
- Reinforcement Learning Toolbox                        Version 1.0         (R2019a)
- Robotics System Toolbox                               Version 2.2         (R2019a)
- Robust Control Toolbox                                Version 6.6         (R2019a)
- Sensor Fusion and Tracking Toolbox                    Version 1.1         (R2019a)
- Signal Processing Toolbox                             Version 8.2         (R2019a)
- SimBiology                                            Version 5.8.2       (R2019a)
- Simulink 3D Animation                                 Version 8.2         (R2019a)
- Simulink Check                                        Version 4.3         (R2019a)
- Simulink Code Inspector                               Version 3.4         (R2019a)
- Simulink Coder                                        Version 9.1         (R2019a)
- Simulink Control Design                               Version 5.3         (R2019a)
- Simulink Coverage                                     Version 4.3         (R2019a)
- Simulink Design Optimization                          Version 3.6         (R2019a)
- Simulink Design Verifier                              Version 4.1         (R2019a)
- Simulink Report Generator                             Version 5.6         (R2019a)
- Simulink Requirements                                 Version 1.3         (R2019a)
- Simulink Test                                         Version 3.0         (R2019a)
- Stateflow                                             Version 10.0        (R2019a)
- Statistics and Machine Learning Toolbox               Version 11.5        (R2019a)
- Symbolic Math Toolbox                                 Version 8.3         (R2019a)
- System Composer                                       Version 1.0         (R2019a)
- System Identification Toolbox                         Version 9.10        (R2019a)
- Text Analytics Toolbox                                Version 1.3         (R2019a)
- Wavelet Toolbox                                       Version 5.2         (R2019a)
