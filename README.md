# CHESHIRE
## Overview

GEnome-scale Metabolic models (GEMs) are powerful tools to predict cellular metabolism and physiological states in living organisms. However, even highly curated GEMs have gaps (i.e., missing reactions) due to our imperfect knowledge of metabolic processes. Here we present a deep learning-based method -- **CHE**byshev **S**pctral **H**yperl**I**nk p**RE**dictor (```CHESHIRE```) -- to predict missing reactions of GEMs purely from the metabolic network topology. ```CHESHIRE``` takes a metabolic network and a pool of candidate reactions as the input and outputs confidence scores for candidate reactions. This package contains the source code of our paper:

Can Chen, Chen Liao, and Yang-Yu Liu. "Teasing out Missing Reactions in Genome-scale Metabolic Networks through Deep Learning." BioRxiv (2022) [[PDF](https://www.biorxiv.org/content/10.1101/2022.06.27.497720v1.full.pdf)].

## System Requirements

### Hardware Requirements
The package requires only a standard computer with enough RAM to support the operations defined by users. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB

CPU: 4+ cores, 2+ GHz/core

### OS Requirements
The package is supported for macOS, and has been tested on the following systems:

macOS Big Sur (version 11.6.2)

macOS Monterey (version 12.4)


### Dependencies
The package depends on the Python scientific stack:

```
cobra==0.22.1
joblib==1.1.0
numpy==1.21.2
optlang==1.5.2
pandas==1.3.2
torch==1.9.0
torch_geometric==1.7.2
torch_scatter==2.0.8
torch_sparse==0.6.11 
tqdm==4.62.1
```

Users are required to additionally install the ```cplex``` solver (https://www.ibm.com/analytics/cplex-optimizer) from IBM for running the package. Note that cplex only works with certain python versions (e.g., cplex 12.10 has APIs for python3.6 an python3.7)/

## Usage
 
Step 1. Download the package from GitHub and ente the folder:

```
git clone https://github.com/canc1993/cheshire-gapfilling.git
cd cheshire-gapfilling
```

Step 2. Prepare the input data. We have provided a demonstration that can be run through ```python3 main.py```.

There are three folders within the directory  "cheshire-gapfilling/data".
The folder ```data/gems``` contains GEMs you want to gap-fill. An an example, we included two GEMs, ```GCF_000005845.2``` and ```GCF_000160535.1```, for gap-filling.
The folder ```data/pools``` contains a reaction pool ('bigg_universe.xml') constructed from the BiGG database. The candidate reactions from the pool wil be used for gap-filling. To use your own pool, rename it to "". If this database is constructed from not Bigg database, edit input_parameters.txt to change the suffix of exchange reactions accordingly.
The folder ```data/others``` contain files for gap prediction. The file media.csv specifies the culture medium used to simulate the GEMs. The substrate_exchange_reactions.csv contains the compounds whose fermentation phenotypes are tested. 


## Demo

The file ```main.py``` is the main program for running the GEM gap-filling experiments. 
All the parameters used in the package are saved in the files ```config.py``` (CHESHIRE parameters) and ```input_parameters.txt``` (phenotypic prediction parameters).

The folder ```results/scores``` contains candidate reaction scores produced by ```CHESHIRE```, and
the folder ```results/gaps``` contains a list of suggested reactions based on the reaction scores and phenotypes. 
