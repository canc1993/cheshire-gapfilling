# CHESHIRE
## Overview

GEnome-scale Metabolic models (GEMs) are powerful tools to predict cellular metabolism and physiological states in living organisms. However, even highly curated GEMs have gaps (i.e., missing reactions) due to our imperfect knowledge of metabolic processes. Here we present a deep learning-based method -- **CHE**byshev **S**pctral **H**yperl**I**nk p**RE**dictor (```CHESHIRE```) -- to predict missing reactions of GEMs purely from the metabolic network topology. ```CHESHIRE``` takes a metabolic network and a pool of candidate reactions as the input and outputs confidence scores for candidate reactions. This package contains the source code of our paper:

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

## Installation Guide
```CHESHIRE``` can be simply installed from GitHub:

```
git clone https://github.com/canc1993/cheshire-gapfilling.git
```

## Demo

The file "main.py" is the main program for running the genome-scale metabolic model gap-filling experiments. 
All the parameters are saved in the files "config.py" and "input_parameters.txt."
The folder "data/gems" contains genome-scale metabolic models for gap-filling, and
the folder "data/pools" contains reaction pools.
The folder "results/scores" contains reaction scores produced by CHESHIRE, and
the folder "results/gaps" contains a list of suggested reactions based on the reaction scores. 
