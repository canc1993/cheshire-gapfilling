# CHESHIRE
## About

CHESHIRE (Chebyshev Spectral Hyperlink Predictor) is a hyperlink prediction algorithm which can be used to gap-fill genome-scale metabolic models. 
This toolbox contains the source code of our paper:

Can Chen, Chen Liao, and Yang-Yu Liu. "Filling Gaps in Genome-scale Metabolic Models through Deep Learning." BioRxiv (2022) [PDF].

## How To Use

The file "main.py" is the main program for running the genome-scale metabolic model gap-filling experiments. 
All the parameters are saved in the files "config.py" and "input_parameters.txt."
The folder "data/gems" contains genome-scale metabolic models for gap-filling, and
the folder "data/pools" contains reaction pools.
The folder "results/scores" contains reaction scores produced by CHESHIRE, and
the folder "results/gaps" contains a list of suggested reactions based on the reaction scores. 
