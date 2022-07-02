# CHESHIRE
## Overview

GEnome-scale Metabolic models (GEMs) are powerful tools to predict cellular metabolism and physiological states in living organisms. However, even highly curated GEMs have gaps (i.e., missing reactions) due to our imperfect knowledge of metabolic processes. Here we present a deep learning-based method -- **CHE**byshev **S**pctral **H**yperl**I**nk p**RE**dictor (```CHESHIRE```) -- to predict missing reactions of GEMs purely from the metabolic network topology. ```CHESHIRE``` takes a metabolic network and a pool of candidate reactions as the input and outputs confidence scores for candidate reactions. Among the top candidates, we further identify key reactions that lead to secretion of fermentation compounds in the gapfilled GEMs. This package contains the source code of our paper:

Can Chen, Chen Liao, and Yang-Yu Liu. "Teasing out Missing Reactions in Genome-scale Metabolic Networks through Deep Learning." BioRxiv (2022) [[PDF](https://www.biorxiv.org/content/10.1101/2022.06.27.497720v1.full.pdf)].

## System Requirements

### Hardware Requirements
The package requires only a standard computer with enough RAM to support the operations defined by users. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB

CPU: 4+ cores, 2+ GHz/core

### OS Requirements
The package has been tested on MacOS Big Sur (version 11.6.2) and Monterey (version 12.3, 12.4).

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

Users are required to additionally install the ```cplex``` solver (https://www.ibm.com/analytics/cplex-optimizer) from IBM to run the package. Note that cplex only works with certain python versions (e.g., CPLEX_Studio12.10 has APIs for python3.6 an python3.7).

## Usage

To run the demonstration, type "python3 main.py" in your terminal. Follow the following steps to identify missing reactions in your own GEMs:
 
**Step 1. Download the package**

```
git clone https://github.com/canc1993/cheshire-gapfilling.git
cd cheshire-gapfilling
```

**Step 2. Prepare input files** 

All input files should be deposited to the directory ```cheshire-gapfilling/data```. There are three folders:

1. The folder ```data/gems``` contains GEMs to be tested. Each GEM is a xml file.

2. The folder ```data/pools``` contains a reaction pool under name ```universe.xml```. Each pool is a GEM that has the extension ```.xml```. To use your own pool, remember to rename it to ```universe.xml```. Also remember to edit ```EX_SUFFIX``` and ```NAMESPACE``` in the input_parameters.txt to specify the suffix of exchange reactions and which namespace of biochemical reaction database is used. For ```NAMESPACE```, we currently only support ```bigg``` and ```modelseed```.

3. The folder ```data/fermentation``` contains two files for GEM simulations. For each pair of input and gapfilled GEMs, our algorithm simulates their fermentation phenotypes and, if the phenotype is positiv in the gapfilled GEM but negative in the input GEM, we identity and output the minimum number of reactions among the to candidates that allow the phenotypic change (column ```key_rxns``` in ```results/gaps/suggested_gaps/csv```). The file ```substrate_exchange_reactions.csv``` contains a list of fermentation compounds that we will search for missing phenotypes in the input GEMs. The file requires at least two columns, including a column ```compound``` to specify the conventional compound  names (e.g., Sucrose) and a column named by ```NAMESPACE``` in the input_parameters.txt to specify the compound IDs (e.g., sucr) in the GEMs. To use your own list of fermentation compounds, remembe to rename it to ```substrate_exchange_reactions.csv```. Additionally, the file ```media.csv``` specifies the culture medium used to simulate the GEMs. This file also requires at leas two columns, including a column named by ```NAMESPACE``` in the input_parameters.txt to specify the compound IDs in the GEMs and another column ```flux``` to specify the maximum uptake flux for each culture medium component.

**Step 3. Setup simulation parameters**

Firs of all, CHESHIRE has two main programs: (1) score the candidate reactions in the pool for their likelihood of being missing in the input GEMs; and (2) among the top candidate reactions with the highest likelihood, find out the minimum set that leads to new metabolic secretions that are potentially missing in the input GEMs. The second program is time-consuming if the number of top candidates added to the input GEMs for simulations is too large (this parameter is controlled by ```NUM_GAPFILLED_RXNS_TO_ADD``` in the input_parameters.txt). ***If you only want the scores and rankings of candidate reactions, comment out ```validate()``` in the main.py***.

Other parameters that are not mentioned in Step 2 are discussed below:
```GEM_DIRECTORY```: directory where GEMs are stored. DO NOT CHANGE.
```REACTION_POOL```: directory of reaction pool. DO NOT CHANGE.
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

**Step 4. Run CHESHIRE by ```python3 main.py```**

**Step 5. Interpret the results**

A folder ```results``` will be generated automatically. It contains three subfolders:

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
