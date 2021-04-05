import params
import torch
from torch import nn
import torch.nn.functional as F


class MIL(torch.nn.Module):
    def __init__(self, embeddingSize, inputDim, C, DK, nSe, nD, nLayer=1, device=torch.device('cpu')):
        super(MIL, self).__init__()
        self.li = nn.Linear(inputDim, embeddingSize).to(device)
        self.WW = []
        self.VV = []
        self.device = device
        self.nSe = nSe
        self.nD = nD
        self.C = C
        self.DK = DK
        self.embeddingSize = embeddingSize
        nLayer = 1
        self.nLayer = nLayer

        for i in range(nLayer):
            vi = torch.nn.Parameter(torch.rand((embeddingSize, DK), requires_grad=True).to(device))
            wi = torch.nn.Parameter(torch.rand((DK, 1), requires_grad=True).to(device))
            self.register_parameter('w_%s' % i, wi)
            self.register_parameter('v_%s' % i, vi)

            self.WW.append(wi)
            self.VV.append(vi)

        outSize = embeddingSize

        self.lo = nn.Linear(outSize, nSe).to(device)
    def forward(self, inp, mask):
        re = self.li(inp)
        # print(torch.isnan(inp).any())
        # print(torch.isnan(inp).any())
        # print(inp.shape)
        # print(re.shape)
        # print(mask.shape)
        # re = re.unsqueeze(1)
        # print(re.shape)
        mask = mask[:, 0, :]
        for i in range(self.nLayer):

            re1 = torch.matmul(re, self.VV[i])
            # print("MATMUL VV", re1.shape)
            re1 = torch.tanh(re1)
            # print(re1.shape, self.WW[i].shape)
            re1 = torch.matmul(re1, self.WW[i])
            re1 = re1.squeeze()
            # print(re1.shape, mask.shape)

            if i == 0:
                re1[~mask] = float('-inf')
            # print(torch.max(re1), torch.min(re1), re1.shape, mask.shape)
            re1 = torch.softmax(re1, dim=-1)
            # print(torch.isnan(re1).any())
            re1 = torch.unsqueeze(re1, -1)

            # print(re.shape, re1.shape)
            re = torch.mul(re, re1)

            # print(torch.isnan(re).any())

            # print("DOne...1")
            # print("Before sum: ",re.shape)
            re = torch.sum(re, dim=1)


            re = torch.relu(re)
            # print(torch.isnan(re).any())

            # print("Done...2")
            # re = torch.unsqueeze(re, 1)
            # print("Done...3")
            # print(re.shape)
            # exit(-1)

        # print("Last layer", re.shape)
        # re = re.squeeze()
        # print(re.shape)
        # print("II")
        out = F.relu(self.lo(re))

        return out
