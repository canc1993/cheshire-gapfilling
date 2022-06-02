import os
import torch
import pandas as pd
import cobra
import math
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.solver import linear_reaction_coefficients
import numpy as np


def get_filenames(path):
    return sorted(os.listdir(path))


def get_data(path, sample):
    model = cobra.io.read_sbml_model(path + '/' + sample)
    biomass = linear_reaction_coefficients(model)
    model.remove_reactions(biomass, remove_orphans=True)
    stoichiometric_matrix = create_stoichiometric_matrix(model)
    incidence_matrix = np.abs(stoichiometric_matrix) > 0
    remove_rxn_index = np.sum(incidence_matrix, axis=0) <= 1
    model.remove_reactions(model.reactions[remove_rxn_index], remove_orphans=True)
    return model


def create_pool():
    print('-------------------------------------------------------')
    print('merging GEMs with reaction pool...')
    path = 'data/gems'
    namelist = get_filenames(path)
    model_pool = cobra.io.read_sbml_model('data/pools/bigg_universe.xml')
    for sample in namelist:
        if sample.endswith('xml'):
            model = get_data(path, sample)
            model_pool.merge(model)
    cobra.io.write_sbml_model(model_pool, 'data/pools/comb_universe.xml')


def get_data_from_pool(path, sample, model_pool_df):
    model = get_data(path, sample)
    rxns_df = pd.read_csv(path + '/reactions_w_gene_reaction_rule.csv')
    rxns = rxns_df.reaction[rxns_df.id == sample[:-4]].to_numpy()
    model_df = model_pool_df[rxns]
    cols2use = model_pool_df.columns.difference(model_df.columns)
    return model_df, model_pool_df[cols2use]


def create_neg_incidence_matrix(incidence_matrix):
    incidence_matrix_neg = torch.zeros(incidence_matrix.shape)
    for i, edge in enumerate(incidence_matrix.T):
        nodes = torch.where(edge)[0]
        nodes_comp = torch.tensor(list(set(range(len(incidence_matrix))) - set(nodes.tolist())))
        edge_neg_l = torch.tensor(np.random.choice(nodes, math.floor(len(nodes) * 0.5), replace=False))
        edge_neg_r = torch.tensor(np.random.choice(nodes_comp, len(nodes) - math.floor(len(nodes) * 0.5), replace=False))
        edge_neg = torch.cat((edge_neg_l, edge_neg_r))
        incidence_matrix_neg[edge_neg, i] = 1
    return incidence_matrix_neg


def hyperlink_score_loss(y_pred, y):
    negative_score = torch.mean(y_pred[y == 0])
    logistic_loss = torch.log(1 + torch.exp(negative_score - y_pred[y == 1]))
    loss = torch.mean(logistic_loss)
    return loss


def create_label(incidence_matrix_pos, incidence_matrix_neg):
    y_pos = torch.ones(len(incidence_matrix_pos.T))
    y_neg = torch.zeros(len(incidence_matrix_neg.T))
    return torch.cat((y_pos, y_neg))