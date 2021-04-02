import params
import torch
from torch import nn
import torch.nn.functional as F


class FFNN(torch.nn.Module):
    def __init__(self, embeddingSize, nSe, nD, nLayer=params.N_LAYER, device=torch.device('cpu')):
        super(FFNN, self).__init__()

        self.model = nn.Sequential()
        if nLayer == 0:
            mm = nn.Linear(nD, nSe).to(device)
            torch.nn.init.xavier_uniform_(mm.weight.data)
            self.model.add_module('io', mm)
            self.model.add_module('reluio', nn.ReLU())
            return

        mm = nn.Linear(nD, embeddingSize).to(device)
        # mm.weight.data.uniform_(0.001, 1)
        torch.nn.init.xavier_uniform_(mm.weight.data)
        self.model.add_module('inp', mm)
        self.model.add_module('relu1', nn.ReLU())
        for i in range(nLayer):
            mm = nn.Linear(embeddingSize, embeddingSize).to(device)
            # mm.weight.data.uniform_(0.001, 1)
            torch.nn.init.xavier_uniform_(mm.weight.data)

            self.model.add_module('layer_%s' % i, mm)
            self.model.add_module('relu_%s' % i, nn.ReLU())


        mm = nn.Linear(embeddingSize, nSe).to(device)
        # mm.weight.data.uniform_(0.001, 1)
        torch.nn.init.xavier_uniform_(mm.weight.data)


        self.model.add_module('out', mm)
        self.model.add_module('reluo', nn.ReLU())

    def forward(self, inp):
        return self.model(inp)
