import torch
import torch.nn as nn
import torch_geometric.nn as gnn
from torch_scatter import scatter


class CHESHIRE(nn.Module):
    def __init__(self, input_dim, emb_dim, conv_dim, k, p):
        super(CHESHIRE, self).__init__()
        self.linear_encoder = nn.Linear(input_dim[1], emb_dim)
        self.tanh = nn.Hardtanh()
        self.norm_emb = gnn.GraphNorm(emb_dim)
        self.dropout = nn.AlphaDropout(p)
        self.graph_conv = gnn.ChebConv(emb_dim, conv_dim, k)
        self.max_pool = gnn.global_max_pool
        self.linear = nn.Linear(2 * conv_dim, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, feature, incidence_matrix):
        x = self.tanh(self.linear_encoder(feature))
        x, hyperedge_index = self.partition(x, incidence_matrix)
        edge_index, batch = self.expansion(hyperedge_index)
        x = self.dropout(self.norm_emb(x, batch))
        x = self.tanh(self.graph_conv(x, edge_index))
        y_maxmin = self.max_pool(x, batch) - self.min_pool(x, batch)
        y_norm = self.norm_pool(x, batch)
        y = torch.cat((y_maxmin, y_norm), dim=1)
        return self.sigmoid(self.linear(y))

    @staticmethod
    def norm_pool(x, batch):
        size = int(batch.max().item() + 1)
        counts = torch.unique(batch, sorted=True, return_counts=True)[1]
        return torch.sqrt(scatter(x ** 2, batch, dim=0, dim_size=size, reduce='sum') / counts.view(-1, 1))

    @staticmethod
    def min_pool(x, batch):
        size = int(batch.max().item() + 1)
        return scatter(x, batch, dim=0, dim_size=size, reduce='min')

    @staticmethod
    def expansion(hyperedge_index):
        node_set = hyperedge_index[0]
        index = hyperedge_index[1].int()
        edge_index = torch.empty((2, 0), dtype=torch.int64)
        batch = torch.empty(len(node_set), dtype=torch.int64)
        for i in range(index[-1] + 1):
            nodes = node_set[index == i]
            batch[nodes.long()] = i
            num_nodes = len(nodes)
            adj_matrix = torch.ones(num_nodes, num_nodes) - torch.eye(num_nodes)
            row, col = torch.where(adj_matrix)
            row, col = nodes[row], nodes[col]
            edge = torch.cat((row.view(1, -1), col.view(1, -1)), dim=0)
            edge_index = torch.cat((edge_index, edge), dim=1)
        return edge_index, batch

    @staticmethod
    def partition(x, incidence_matrix):
        row, col = torch.where(incidence_matrix.T)
        hyperedge_index = torch.cat((col.view(1, -1), row.view(1, -1)), dim=0)
        node_set, sort_index = torch.sort(hyperedge_index[0])
        hyperedge_index[1] = hyperedge_index[1][sort_index]
        x = x[node_set.long(), :]
        hyperedge_index[0] = torch.arange(0, len(hyperedge_index[0]))
        index_set, sort_index = torch.sort(hyperedge_index[1])
        hyperedge_index[1] = index_set
        hyperedge_index[0] = hyperedge_index[0][sort_index]
        return x, hyperedge_index
