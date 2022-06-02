from utils import *
import torch
from CHESHIRE import CHESHIRE
from tqdm import tqdm
import config
import pandas as pd
import cobra

args = config.parse()


def train(feature, y, incidence_matrix, model, optimizer):
    model.train()
    optimizer.zero_grad()
    y_pred = model(feature, incidence_matrix)
    loss = hyperlink_score_loss(y_pred, y)
    loss.backward()
    optimizer.step()


def test(feature, incidence_matrix, model):
    model.eval()
    with torch.no_grad():
        y_pred = model(feature, incidence_matrix)
    return torch.squeeze(y_pred)


def predict():
    print('-------------------------------------------------------')
    path = 'data/gems'
    namelist = get_filenames(path)
    # read reaction pool
    bigg_pool = cobra.io.read_sbml_model('data/pools/bigg_universe.xml')
    model_pool = cobra.io.read_sbml_model('data/pools/comb_universe.xml')
    model_pool_df = create_stoichiometric_matrix(model_pool, array_type='DataFrame')
    for sample in namelist:
        if sample.endswith('.xml'):
            print('training CHESHIRE and predicting reaction scores: ' + sample[:-4] + '...')
            # read the model and reaction pool
            rxn_df, rxn_pool_df = get_data_from_pool(path, sample, model_pool_df)
            incidence_matrix_pos = np.abs(rxn_df.to_numpy()) > 0
            incidence_matrix_pos = torch.tensor(incidence_matrix_pos, dtype=torch.float)
            incidence_matrix_pos = torch.unique(incidence_matrix_pos, dim=1)
            incidence_matrix_cand = np.abs(rxn_pool_df.to_numpy()) > 0
            incidence_matrix_cand = torch.tensor(incidence_matrix_cand, dtype=torch.int64)
            score = torch.empty((incidence_matrix_cand.shape[1], args.num_iter))
            for i in range(args.num_iter):
                # create negative reactions
                incidence_matrix_neg = create_neg_incidence_matrix(incidence_matrix_pos)
                incidence_matrix_neg = torch.unique(incidence_matrix_neg, dim=1)
                incidence_matrix = torch.cat((incidence_matrix_pos, incidence_matrix_neg), dim=1)
                y = create_label(incidence_matrix_pos, incidence_matrix_neg)
                model = CHESHIRE(input_dim=incidence_matrix_pos.shape, emb_dim=args.emb_dim, conv_dim=args.conv_dim, k=args.k, p=args.p)
                optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
                for _ in tqdm(range(args.max_epoch)):
                    # training
                    train(incidence_matrix_pos, y, incidence_matrix, model, optimizer)
                # predicting reaction scores
                score[:, i] = test(incidence_matrix_pos, incidence_matrix_cand, model)
            score_df = pd.DataFrame(data=score.detach().numpy(), index=rxn_pool_df.columns)
            bigg_rxns = set([item.id for item in bigg_pool.reactions])
            common_rxns = list(bigg_rxns & set(score_df.index))
            common_score_df = score_df.T[common_rxns].T
            common_score_df.to_csv('./results/scores/' + sample[:-4] + '.csv')