import os
import torch
import numpy as np
import pandas as pd
import networkx as nx
import cobra
import math
from node2vec import Node2Vec
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.solver import linear_reaction_coefficients
import warnings
import re

warnings.filterwarnings("ignore")


def get_filenames(path):
    return sorted(os.listdir(path))


def get_data(path, sample):
    model = cobra.io.read_sbml_model(path + '/' + sample)
    if path[-4:] == 'bigg':
        biomass_rxns = np.array(pd.read_csv('./data/pools/bigg_biomass_reactions.csv').bigg_id.to_list())
        rxns = np.array([rxn.id for rxn in model.reactions])
        biomass = model.reactions[np.isin(rxns, biomass_rxns)]
    else:
        biomass = linear_reaction_coefficients(model)
    model.remove_reactions(biomass, remove_orphans=True)
    ## modelseed
    if path[-9:] == 'modelseed':
        rxns = [rxn.id for rxn in model.reactions]
        r = re.compile("EX_.*")
        ex_rxns = list(filter(r.match, rxns))
        remove_ex_rxn_index = [rxns.index(ex_rxn) for ex_rxn in ex_rxns]
        ex_rxns = [model.reactions[index] for index in remove_ex_rxn_index]
        model.remove_reactions(ex_rxns, remove_orphans=True)
    stoichiometric_matrix = create_stoichiometric_matrix(model)
    incidence_matrix = np.abs(stoichiometric_matrix) > 0
    remove_rxn_index = np.sum(incidence_matrix, axis=0) <= 1
    model.remove_reactions(model.reactions[remove_rxn_index], remove_orphans=True)
    return model, np.abs(create_stoichiometric_matrix(model)) > 0


def create_pool(name):
    path = './data/'+ name
    namelist = get_filenames(path)
    model_pool = cobra.io.read_sbml_model('./data/pools/bigg_universe.xml')
    for sample in namelist:
        if sample.endswith('xml'):
            model = get_data(path, sample)[0]
            model_pool.merge(model)
    cobra.io.write_sbml_model(model_pool, './data/pools/' + name + '_universal_model_pool.xml')


def get_data_from_pool(path, sample, model_pool_df, rxns_w_genes):
    added_rxns = None
    if rxns_w_genes:
        model = get_data(path, sample)[0]
        all_rxns = [rxn.id for rxn in model.reactions]
        rxns_df = pd.read_csv(path + '/reactions_w_gene_reaction_rule.csv')
        rxns = rxns_df.reaction[rxns_df.id == sample[:-4]].to_numpy()
        added_rxns = np.array(list(set(all_rxns) - set(list(rxns))))
    else:
        model = get_data(path, sample)[0]
        rxns = np.array([rxn.id for rxn in model.reactions])
    model_df = model_pool_df[rxns]
    cols2use = model_pool_df.columns.difference(model_df.columns)
    return model_df, model_pool_df[cols2use], added_rxns


def get_data_from_pool2(path, sample, model_pool):
    model = get_data(path, sample)[0]
    rxns = np.array([rxn.id for rxn in model.reactions])
    model_pool.merge(model)
    model_pool_df = create_stoichiometric_matrix(model_pool, array_type='DataFrame')
    model_df = model_pool_df[rxns]
    cols2use = model_pool_df.columns.difference(model_df.columns)
    return model_df, model_pool_df[cols2use]


def create_hyperedge_index(incidence_matrix):
    row, col = torch.where(incidence_matrix.T)
    hyperedge_index = torch.cat((col.view(1, -1), row.view(1, -1)), dim=0)
    return hyperedge_index


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


# def train_test_split(incidence_matrix, y_label, train_size):
#     size = len(y_label)
#     shuffle_index = torch.randperm(size)
#     incidence_matrix = incidence_matrix[:, shuffle_index]
#     y_label = y_label[shuffle_index]
#     split_index = round(train_size * size)
#     return incidence_matrix[:, :split_index], y_label[:split_index], incidence_matrix[:, split_index:], y_label[split_index:]


#def node2vec(incidence_matrix, emb_dim):
#    size = len(incidence_matrix)
#    hyperedge_index = create_hyperedge_index(incidence_matrix)
#    node_set = hyperedge_index[0]
#    b = hyperedge_index[1].int()
#    edgelist = []
#    for i in range(b[-1] + 1):
#        nodes = node_set[b == i]
#        num_nodes = len(nodes)
#        adj_matrix = torch.triu(torch.ones(num_nodes, num_nodes) - torch.eye(num_nodes))
#        row, col = torch.where(adj_matrix)
#        row, col = nodes[row], nodes[col]
#        for j in range(len(row)):
#            edgelist.append((int(row[j]), int(col[j])))
#    G = nx.Graph()
#    G.add_nodes_from(range(size))
#    G.add_edges_from(edgelist)
#    node2vec_emb = Node2Vec(G, dimensions=emb_dim, quiet=True)
#    emb = node2vec_emb.fit()
#    return torch.tensor(emb.wv.vectors)
