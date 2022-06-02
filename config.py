import argparse


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--emb_dim', type=int, default=256)
    parser.add_argument('--conv_dim', type=int, default=128)
    parser.add_argument('--k', type=int, default=3)
    parser.add_argument('--p', type=float, default=0.1)
    parser.add_argument('--num_iter', type=int, default=1)
    parser.add_argument('--max_epoch', type=int, default=30)
    parser.add_argument('--lr', type=float, default=0.01)
    parser.add_argument('--weight_decay', type=float, default=5e-4)
    return parser.parse_args()