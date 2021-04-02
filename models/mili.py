import params
import torch
from torch import nn
import torch.nn.functional as F


class MILI(torch.nn.Module):
    def __init__(self, embeddingSize, inputDim, C, DK, nSe, nD, nLayer=params.N_LAYER, device=torch.device('cpu')):
        super(MILI, self).__init__()
        self.li = nn.Linear(inputDim, embeddingSize).to(device)
        self.WW = []
        self.VV = []
        self.device = device
        for i in range(nLayer):
            vi = torch.rand((C, embeddingSize, DK), requires_grad=True).to(device)
            wi = torch.rand((C, DK, 1), requires_grad=True).to(device)
            self.register_parameter('w_%s' % i, wi)
            self.register_parameter('v_%s' % i, vi)

            self.WW.append(wi)
            self.VV.append(vi)

        outSize = C * nSe

        self.lo = nn.Linear(outSize, nSe).to(device)
        self.nLayer = nLayer
    def forward(self, inp, mask):
        x = self.li(inp)
        for i in range(self.nLayer):
            re = torch.matmul(x, self.VV[i])
            if i == 0:
                re[~mask] = float('-inf')
            re = torch.softmax(re, dim=-1)

        return
