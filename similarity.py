import pandas as pd
from utils import get_data
import numpy as np
import cobra
from cobra.util.array import create_stoichiometric_matrix
import scipy.spatial.distance as distance
import os


def get_similarity_score(name, top_N):
    path = './results/predicted_scores'
    all_files = sorted(os.listdir(path))
    model_pool = cobra.io.read_sbml_model('./data/pools/bigg_universe.xml')

    for sample in all_files:
        if sample.endswith('csv'):
            model_pool_copy = model_pool.copy()
            scores = pd.read_csv(path +'/' + sample, index_col=0).mean(axis=1).sort_values(ascending=False)
            model = get_data('./data/' + name, sample[:-4] + '.xml')[0]
            rxns = [rxn.id for rxn in model.reactions]
            model_pool_copy.merge(model)
            model_pool_df = create_stoichiometric_matrix(model_pool_copy, array_type='DataFrame')

            candidate_rxns = scores.index.tolist()[:top_N]
            candidate_matrix = model_pool_df[candidate_rxns].values
            model_matrix = model_pool_df[rxns].values
            feature = np.concatenate((candidate_matrix, model_matrix), axis=1)
            feature = feature[np.sum(np.abs(feature), axis=1) > 0, :]
            corr_matrix = np.abs(1 - distance.cdist(feature.T, feature.T, 'correlation'))
            similarity_max = corr_matrix[:candidate_matrix.shape[1], candidate_matrix.shape[1]:].max(axis=1)

            predicted_scores = scores.loc[candidate_rxns].values
            all_scores = np.concatenate((predicted_scores.reshape(-1, 1), similarity_max.reshape(-1, 1)), axis=1)
            all_scores_df = pd.DataFrame(data=all_scores, index=candidate_rxns, columns=['predicted_scores', 'similarity_scores'])
            all_scores_df.to_csv('./results/similarity_scores/' + sample)


#if __name__ == "__main__":
#    similarity(name='zimmermann', top_N=2000)
