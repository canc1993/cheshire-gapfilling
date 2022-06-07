# CHESHIRE
## Overview

```CHESHIRE``` (Chebyshev Spectral Hyperlink Predictor) is a hyperlink prediction algorithm which can be used to gap-fill genome-scale metabolic models. 
This package contains the source code of our paper:

Can Chen, Chen Liao, and Yang-Yu Liu. "Filling Gaps in Genome-scale Metabolic Models through Deep Learning." BioRxiv (2022) [PDF].

## System Requirements

### Hardware Requirements
```CHESHIRE``` requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB

CPU: 4+ cores, 2+ GHz/core

### OS Requirements
```CHESHIRE``` is supported for macOS, and has been tested on the following systems:

macOS Big Sur (version 11.6.2)

macOS Monterey (version 12.4)


### Python Dependencies
```CHESHIRE``` depends on the Python scientific stack:

```
cobra==0.22.1
joblib==1.1.0
numpy==1.21.2
optlang==1.5.2
pandas==1.3.2
torch==1.9.0
torch_geometric==1.7.2
torch_scatter==2.0.8
tqdm==4.62.1
```



## Demo

The file "main.py" is the main program for running the genome-scale metabolic model gap-filling experiments. 
All the parameters are saved in the files "config.py" and "input_parameters.txt."
The folder "data/gems" contains genome-scale metabolic models for gap-filling, and
the folder "data/pools" contains reaction pools.
The folder "results/scores" contains reaction scores produced by CHESHIRE, and
the folder "results/gaps" contains a list of suggested reactions based on the reaction scores. 
