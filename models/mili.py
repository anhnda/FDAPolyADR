import torch
import torch.nn.functional as F


class WHGNN(torch.nn.Module):
    def __init__(self, featureSize, embeddingSize, nSe, nD, nLayer=params.N_LAYER, device=torch.device('cpu')):
        super(WHGNN, self).__init__()
        self.nSe = nSe
        self.nD = nD
        self.nV = nSe + nD
        self.embeddingSize = embeddingSize
        self.device = device
        self.feature2EmbedLayer1 = torch.nn.Linear(featureSize, embeddingSize)
        self.feature2EmbedLayer2 = torch.nn.Linear(embeddingSize, embeddingSize)
        self.embeddingSe = torch.nn.Embedding(nSe, embeddingSize)
        self.embeddingSe.weight.data.uniform_(0.001, 0.3)
        self.layerWeightList = []
        self.dimWeightList = []
        self.nLayer = nLayer
