from torch_geometric.nn import GCNConv, global_mean_pool, Linear
import torch.nn.functional as F
import torch


class GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        torch.manual_seed(1234)
        self.conv1 = GCNConv(6, 32)
        self.conv2 = GCNConv(32, 32)
        self.conv3 = GCNConv(32, 32)
        self.classifier = Linear(32, 6)

    def forward(self, x, edge_index, batch):
        h = self.conv1(x, edge_index)
        h = h.tanh()
        h = self.conv2(h, edge_index)
        h = h.tanh()
        h = self.conv3(h, edge_index)
        h = h.tanh()  # Final GNN embedding space.
        # Apply a final (linear) classifier.
        h = global_mean_pool(h, batch)
        
        h = F.dropout(h, p=0.5, training = self.training)
        out = self.classifier(h)

        return out

