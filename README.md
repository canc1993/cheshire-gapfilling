# CHESHIRE
## Overview

GEnome-scale Metabolic models (GEMs) are powerful tools to predict cellular metabolism and physiological states in living organisms. However, even highly curated GEMs have gaps (i.e., missing reactions) due to our imperfect knowledge of metabolic processes. Here we present a deep learning-based method -- **CHE**byshev **S**pctral **H**yperl**I**nk p**RE**dictor (```CHESHIRE```) -- to predict missing reactions of GEMs purely from the metabolic network topology. ```CHESHIRE``` takes a metabolic network and a pool of candidate reactions as the input and outputs confidence scores for candidate reactions. Among the top candidates, we further identify key reactions that lead to secretion of fermentation compounds in the gap-filled GEMs. This package contains the source code of our paper:

Can Chen*, Chen Liao*, and Yang-Yu Liu. "Teasing out Missing Reactions in Genome-scale Metabolic Networks through Hypergraph Learning." Nature Communications, vol. 14, p.2375, 2023 [[PDF](https://www.nature.com/articles/s41467-023-38110-7)].

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

1. The folder ```data/zimmermann``` contains 2 GEMs from Zimmermann et al. as examples. Each GEM is a xml file.

Zimmermann, J., Kaleta, C. & Waschina, S. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. Genome Biol 22, 81 (2021). https://doi.org/10.1186/s13059-021-02295-1

2. The folder ```data/pools``` contains a reaction pool under name ```bigg_universe.xml```. Each pool is a GEM that has the extension ```.xml```. To use your own pool, remember to rename it to ```universe.xml```. Also remember to edit ```EX_SUFFIX``` and ```NAMESPACE``` in the input_parameters.txt to specify the suffix of exchange reactions and which namespace of biochemical reaction database is used. For ```NAMESPACE```, we currently only support ```bigg``` and ```modelseed```.

3. The folder ```data/fermentation``` contains two files for GEM simulations. For each pair of input and gap-filled GEMs, our algorithm simulates their fermentation phenotypes and, if the phenotype is positive in the gap-filled GEM but negative in the input GEM, we identify and output the minimum number of reactions among the top candidates that allow the phenotypic change (column ```key_rxns``` in ```results/gaps/suggested_gaps.csv```). The file ```substrate_exchange_reactions.csv``` contains a list of fermentation compounds that we will search for missing phenotypes in the input GEMs. The file requires at least two columns, including a column ```compound``` to specify the conventional compound  names (e.g., sucrose) and a column named by ```NAMESPACE``` in the input_parameters.txt to specify the compound IDs (e.g., sucr) in the GEMs. To use your own list of fermentation compounds, remember to rename it to ```substrate_exchange_reactions.csv```. Additionally, the file ```media.csv``` specifies the culture medium used to simulate the GEMs. This file also requires at least two columns, including a column named by ```NAMESPACE``` in the input_parameters.txt to specify the compound IDs in the GEMs and another column ```flux``` to specify the maximum uptake flux for each culture medium component.

**Step 3. Setup simulation parameters**

First of all, CHESHIRE has three main programs: (1) score the candidate reactions in the pool for their likelihood of being missing in the input GEMs (function ```get_predicted_score()``` in main.py); (2) score the mean similarity of the candidate reactions to the existing reactions in the input GEMs (function ```get_similarity_score()``` in main.py); and (3) among the top candidate reactions with the highest likelihood, find out the minimum set that leads to new metabolic secretions that are potentially missing in the input GEMs (function ```validate()``` in main.py). The last program is time-consuming if the number of top candidates added to the input GEMs for simulations is too large (this parameter is controlled by ```NUM_GAPFILLED_RXNS_TO_ADD``` in the input_parameters.txt). ***If you only want the scores and rankings of candidate reactions, comment out ```validate()``` in main.py***.

All simulation parameters are defined in the input_parameters.txt:

1. ```CULTURE_MEDIUM``` (mandatory): filepath of culture medium. For the moment, use ```./data/fermentation/media.csv``` always.

2. ```REACTION_POOL``` (mandatory): filepath of reaction pool. Note that the reaction pool should use the same namespace (e.g., BiGG, ModelSeed) as the GEM files. For the moment, use ```./data/pools/universe.xml``` always.

3. ```GEM_DIRECTORY``` (mandatory): directory of input GEMs. For the moment, use ```./data/gems/``` always.

4. ```GAPFILLED_RXNS_DIRECTORY``` (mandatory): filepath of candidate reaction scores. For the moment, use ```./results/scores``` always.

5. ```NUM_GAPFILLED_RXNS_TO_ADD``` (mandatory): number of top candidate reactions predicted by CHESHIRE to be added to the input GEMs for the fermentation test.

6. ```ADD_RANDOM_RXNS``` (mandatory): a boolean value (0/1). If true, we will randomly select ```NUM_GAPFILLED_RXNS_TO_ADD``` reactions from the reaction pool, instead of using reactions with highest CHESHIRE scores.

7. ```SUBSTRATE_EX_RXNS``` (mandatory): filepath of fermentation compounds to be tested. For the moment, use ```./data/fermentation/substrate_exchange_reactions.csv``` always.

8. ```NUM_CPUS``` (optional, default = 1): number of CPUs used for simulations in ```validate()```. Note that the first program ```predict()``` that outputs reaction scores is not parallelized.

9. ```EX_SUFFIX``` (optional, default = "_e"): suffix of exchange reactions.

10. ```RESOLVE_EGC``` (optional, default = 1): a boolean value (0/1). If true, it takes extra time to resolve energy-generating cycles (see our manuscript for how we resolve the EGCs).

11. ```OUTPUT_DIRECTORY``` (optional, default = "./results/gaps"): output directory of simulation results.

12. ```OUTPUT_FILENAME``` (optional, default = "suggested_gaps.csv"): filename of simulation results to output.

13. ```FLUX_CUTOFF``` (optional, default = 1e-5): a fermentation phenotype is considered as positive if the maximum secretion flux is larger than the cutoff value.

14. ```ANAEROBIC``` (optional, default = 1): a boolean value (0/1). If true, candidate reactions involving oxygen molecules will be skipped during gap-filling and simulations.

15. ```BATCH_SIZE``` (optional, default = 10): An integer to indicate how many reactions are added in a batch (together) during gap-filling. We test for EGC after each batch and if EGC is found, reactions in the batch will be added one by one.

16. ```NAMESPACE``` (optional, default = "bigg"): Namespace of GEMs and reaction pool. Currently, we only support BiGG (```bigg```) and ModelSeed (```modelseed```).

17. ```MIN_PREDICTED_SCORES``` (optional, default = 0.9995): Cnadidate reactions with predicted scores below this cutoff will be discarded and not used for gap-filling.

**Step 4. Run CHESHIRE by ```python3 main.py```**

**Step 5. Interpret the results**

The output files will be saved to the folder ```cheshire-gapfilling/results```. The directory contains three subfolders:

1. ```universe```: A merged pool that combines the reactions in the user-provided pool (```./data/pools/universe.xml```) and all reactions in the input GEMs (```./data/gems```).

2. ```scores```: Predicted reaction scores for each GEM. Rows are reaction IDs from the pool and columns are each individual Monte-Carlo simulation run. To rank the reactions, we use the mean scores across all runs. By default, we run it once. To change this number, edit ```config.py``` and change the default of parameter ```num_iter``` to increase the prediction robustness.

3. ```gaps```: Simulations of metabolic fermentation for all input GEMs and their corresponding gap-filled models (i.e., after adding top candidate reactions). Each row is an exchange reaction (i.e., a compound that can be secreted) and columns are explained as follows:

* ```minimum__no_gapfill```: minimum secretion flux of the input GEM (lower bound of flux variability analysis)

* ```maximum__no_gapfill```: maximum secretion flux of the input GEM (upper bound of flux variability analysis)

* ```biomass__no_gapfill```: biomass production rate of the input GEM

* ```normalized_maximum__no_gapfill```: ```maximum__no_gapfill``` divided by ```biomass__no_gapfill``` for the input GEM

* ```phenotype__no_gapfill```: a binary value (0/1) indicating whether ```normalzied_maximum__no_gapfill``` >= ```FLUX_CUTOFF``` (specified in the input_parameters.txt)

* ```minimum__w_gapfill```: minimum secretion flux of the gap-filled GEM (lower bound of flux variability analysis)

* ```maximum__w_gapfill```: maximum flux of the gap-filled GEM (upper bound of flux variability analysis)

* ```biomass__w_gapfill```: biomass production rate of the gap-filled GEM

* ```normalized_maximum__w_gapfill```: ```maximum__w_gapfill``` divided by ```biomass__w_gapfill``` for the gap-filled GEM

* ```phenotype__w_gapfill```: a binary value (0/1) indicating whether ```normalzied_maximum__w_gapfill``` >= ```FLUX_CUTOFF``` (specified in the 
input_parameters.txt)

* ```gem_file```: the GEM file name

* ```random_rxns```: a binary value (0/1) indicating whether ```gem_file``` is gap-filled by adding random reactions (equal to the value of ```ADD_RANDOM_RXNS``` in the input_parameters.txt)

* ```num_rxns_to_add```: number of reactions added during gap-filling.

* ```rxn_ids_added```: IDs of candidate reactions that have been added

* ```essential reactions```: If we found a fermentation phenotypic change from 0 (input GEM) to 1 (gap-filled GEM), we used mixed-integer linear programming to determine the minimum number of reactions that are necessary to achieve this phenotypic transition. Otherwise this field is left empty.
