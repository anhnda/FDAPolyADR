import params
import torch
from torch import nn
import torch.nn.functional as F


class FFNN(torch.nn.Module):
    def __init__(self, embeddingSize, nSe, nD, nLayer=params.N_LAYER, device=torch.device('cpu')):
        super(FFNN, self).__init__()

        self.model = nn.Sequential()
        self.model.add_module('inp', nn.Linear(nD, embeddingSize).to(device))
        self.model.add_module('relu1', nn.ReLU())
        for i in range(nLayer):
            self.model.add_module('layer_%s' % i, nn.Linear(embeddingSize, embeddingSize).to(device))
            self.model.add_module('relu_%s' % i, nn.ReLU())
        self.model.add_module('out', nn.Linear(embeddingSize, nSe).to(device))
        self.model.add_module('reluo',  nn.ReLU())

    def forward(self, inp):
        return self.model(inp)
