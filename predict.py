from utils import *
import torch
from CHESHIRE import CHESHIRE
from tqdm import tqdm
import config
import pandas as pd
import cobra
from os.path import exists

args = config.parse()


def train(feature, y, incidence_matrix, model, optimizer):
    model.train()
    optimizer.zero_grad()
    y_pred = model(feature, incidence_matrix)
    loss = hyperlink_score_loss(y_pred, y)
    loss.backward()
    optimizer.step()


def predict(feature, incidence_matrix, model):
    model.eval()
    with torch.no_grad():
        y_pred = model(feature, incidence_matrix)
    return torch.squeeze(y_pred)


def get_prediction_score(name):
    path = './data/' + name
    namelist = get_filenames(path)
    universe_pool = cobra.io.read_sbml_model('./data/pools/bigg_universe.xml')
    for sample in namelist:
        if sample.endswith('.xml'):
            universe_pool_copy = universe_pool.copy()
            rxn_df, rxn_pool_df = get_data_from_pool2(path, sample, universe_pool_copy)
            incidence_matrix_pos = np.abs(rxn_df.to_numpy()) > 0
            incidence_matrix_pos = torch.tensor(incidence_matrix_pos, dtype=torch.float)
            incidence_matrix_pos = torch.unique(incidence_matrix_pos, dim=1)
            incidence_matrix_cand = np.abs(rxn_pool_df.to_numpy()) > 0
            incidence_matrix_cand = torch.tensor(incidence_matrix_cand, dtype=torch.int64)

            incidence_matrix_neg = create_neg_incidence_matrix(incidence_matrix_pos)
            incidence_matrix_neg = torch.unique(incidence_matrix_neg, dim=1)
            incidence_matrix = torch.cat((incidence_matrix_pos, incidence_matrix_neg), dim=1)
            y = create_label(incidence_matrix_pos, incidence_matrix_neg)
            model = CHESHIRE(input_dim=incidence_matrix_pos.shape, emb_dim=args.emb_dim, conv_dim=args.conv_dim, k=args.k, p=args.p)
            optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
            for _ in tqdm(range(args.max_epoch)):
                train(incidence_matrix_pos, y, incidence_matrix, model, optimizer)
            score = predict(incidence_matrix_pos, incidence_matrix_cand, model)
            score_df = pd.DataFrame(data=score.detach().numpy(), index=rxn_pool_df.columns)
            if exists('./results/predicted_scores/' + sample[:-4] + '.csv'):
                exist_score_df = pd.read_csv('./results/predicted_scores/' + sample[:-4] + '.csv', index_col=0)
                score_df = pd.concat([exist_score_df, score_df], axis=1)
                score_df.to_csv('./results/predicted_scores/' + sample[:-4] + '.csv')
            else:
                score_df.to_csv('./results/predicted_scores/' + sample[:-4] + '.csv')


#if __name__ == "__main__":
#    repeat = 5
#    for i in range(repeat):
#        main(name='zimmermann')
