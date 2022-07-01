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

The folder ```data/pools``` contains a reaction pool that provides candidates of missing reactions. As an example, we included a pool constructed from the BiGG database (```bigg_universe.xml```). To use your own pool, rename it to "". If your pool is constructed from a biochemical reaction database other than BiGG (e.g., ModelSeed), edit ```REACTION_POOL``` and ```EX_SUFFIX``` in the input_parameters.txt to specify the path to the pool file and the suffix of exchange reactions used in your pool respectively. Also edit ```NAMESPACE``` in the same file to indicate which namespace of biochemical reaction database is used. We currently only support ```bigg``` and ```modelseed```.

The folder ```data/others``` contain files for gap prediction. For all fermentation compounds listed in ```substrate_exchange_reactions.csv```, we will test whether this phenotype is potentially missing from the input GEMs and can be recovered by adding reactions predicted by CHESHIRE. To use your own list of fermentation compounds, edit ```SUBSTRATE_EX_RXNS``` in the input_parameters.txt to specify the path to your file. Additionally, The file ```media.csv``` specifies the culture medium (compound ID and maximum uptake flux) used to simulate the GEMs. To use your own culture medium file, edit ```CULTURE_MEDIUM``` in the input_parameters.txt to specify the file location. Make sure that the medium file have columns ```compound``` and the namespace (contains maximum uptake flux) you specified in ```NAMESPACE``` in the input_parameters.txt. 

Step 3. Setup parameters. Firs of all, CHESHIRE has two main programs: (1) ranking the candidate reactions in the reaction pool based on the predicted scores (i.e., likelihoo of being ) and (2) among the predicted reactions, find out the key reactions which lead to positive phenotypic predictiosn that are potentially missing in the input GEMs. The second program is time-consuming if too many reactions are added. To run th first step only, comment out ```validate()``` in the main.py.

Other parameters that are not mentioned in Step 2 are discussed below:
```GAPFILLED_RXNS_DIRECTORY```: directory which contains scores of candidate reactions in the pool. DO NOT CHANGE the default folder.
```NUM_GAPFILLED_RXNS_TO_ADD```: number of top reactions predicted by CHESHIRE will be added to the GEMs for fermentation test.
```ADD_RANDOM_RXNS```: a boolean value (0/1) indicating whether random reactions from the pool are added to the GEMs for fermentation test.
```NUM_CPUS```: number of CPUs used for adding reactions. Note that the first program that predicts scores is not parallelized.
```RESOLVE_EGC```: a boolean value (0/1) indicating whether it takes extra time to resolve energy-generating cycles.
```OUTPUT_DIRECTORY```: the directory that suggested gaps will be outputed to.
```OUTPUT_FILENAME```: file name of suggested gaps.
```FLUX_CUTOFF```: a fermentation phenotype is considered as positive if the maximum secretion flux is larger than the cutoff value.
```ANAEROBIC```: a boolean value indicating whether the fermentation test is performed in anaerobic condition. If yes, then reactions involving oxygen molecules will be skipped during reaction addition.
```BATCH_SIZE```: An integer to indicate how many reactions are added in a batch (together). If RESOLVE_EGC is true, reaction needs to be added one by one because the occurance of EGC may depends on multiple reactions that are added. By setting the batch size to a small number (e.g., 10), we wil add 10 reactions in a batch and add the next batch if this batch does not result in EGC. If an EGC is found, then we go back to add the reactions 1by1 . 

Step 4. Run gapfilling by ```python3 main.py```.

Step 5. Interpret the results. A folder ```results``` will be generated automatically. It contains three subfolders:

1. scores: Predicted scores for each GEM. Rows are reaction IDs from the pool and columns are each individual prediction. By default, we run it twice. To change this number, edit ```config.py``` and change the default of parameter ```xxx``` to the number you want.
2. gaps: suggested gaps include fermentation simualtions for all GEMs. The columns are 
minimum__no_gapfill: minimum flux of the input GEM (lower bound of flux variability analysis)
maximum__no_gapfill: maximum flux of the input GEM (upper bound of flux variability analysis)
biomass__no_gapfill: biomas production rate of the input GEM
normalized_maximum__no_gapfill: maximum flux divide by biomass production rate for the input GEM
phenotype__no_gapfill: a binary value indicating whether normalzied_maximum__no_gapfill exceeds the flux cutoff (```FLUX_CUTOFF```) specified in the input_parameters.txt
minimum__w_gapfill: minimum flux of the gapfilled GEM (lower bound of flux variability analysis)
maximum__w_gapfill: maximum flux of the gapfilled GEM (upper bound of flux variability analysis)
biomass__w_gapfill: biomas production rate of the gapfilled GEM
normalized_maximum__w_gapfill: maximum flux divide by biomass production rate for the gapfilled GEM
phenotype__w_gapfill: a binary value indicating whether normalzied_maximum__w_gapfill exceeds the flux cutoff (```FLUX_CUTOFF```) specified in the input_parameters.txt
gem_file: the input GEM ID
random_rxns: whether this test is added by random reactions (basically the value of ```ADD_RANDOM_REACTIONS```)
num_rxns_to_add: number of reactions added during gapfilling
rxn_ids_added: which reactions have been added
key reactions: if we found a fermentation phenotypic change from 0 (input GEM) to 1 (gapfilled GEM), we used linear-mixed effects model to determine the minimum number of reactions that are necessary to achieve this phenotypic transition.

num_rxns_to_add: number of reactions added during gapfilling
