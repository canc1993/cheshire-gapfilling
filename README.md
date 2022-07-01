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

Users are required to additionally install the ```cplex``` solver (https://www.ibm.com/analytics/cplex-optimizer) from IBM for running the package. Note that cplex only works with certain python versions (e.g., CPLEX_Studio12.10 has APIs for python3.6 an python3.7)/

## Usage
 
Step 1. Download the package from GitHub and ente the folder:

```
git clone https://github.com/canc1993/cheshire-gapfilling.git
cd cheshire-gapfilling
```

Step 2. Prepare the input data in the directory ```cheshire-gapfilling/data```. We have provided a demonstration that 

There are three folders within the directory ```cheshire-gapfilling/data```.

The folder ```data/gems``` contains GEMs for gapfilling. An an example, we included two GEMs, ```GCF_000005845.2``` and ```GCF_000160535.1```. To use your own GEMs, remember to edit ```GEM_DIRECTORY``` in the input_parameters.txt to specify the directory that contains your GEMs.

The folder ```data/pools``` contains a reaction pool that provides candidates of missing reactions. As an example, we included a pool constructed from the BiGG database (```bigg_universe.xml```). To use your own pool, rename it to "". If your pool is constructed from a biochemical reaction database other than BiGG (e.g., ModelSeed), edit ```REACTION_POOL``` and ```EX_SUFFIX``` in the input_parameters.txt to specify the path to the pool file and the suffix of exchange reactions used in your pool respectively.

The folder ```data/others``` contain files for gap prediction. For all fermentation compounds listed in ```substrate_exchange_reactions.csv```, we will test whether this phenotype is potentially missing from the input GEMs and can be recovered by adding reactions predicted by CHESHIRE. To use your own list of fermentation compounds, edit ```SUBSTRATE_EX_RXNS``` in the input_parameters.txt to specify the path to your file. Additionally, The file ```media.csv``` specifies the culture medium (compound ID and maximum uptake flux) used to simulate the GEMs. To use your own culture medium file, edit ```CULTURE_MEDIUM``` in the input_parameters.txt to specify the file location.

Step 3. Setup parameters. Firs of all, CHESHIRE has two main programs: (1) ranking the candidate reactions in the reaction pool based on the predicted scores (i.e., likelihoo of being ) and (2) among the predicted reactions, find out the key reactions which lead to positive phenotypic predictiosn that are potentially missing in the input GEMs. The second program is time-consuming if too many reactions are added. To run th first step only, comment out ```xxx``` in the main.py.

Step 4. Run gapfilling by ```python3 main.py```.

Step 5. Interpret the results. 

