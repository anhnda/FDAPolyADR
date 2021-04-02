import torch

ar1 = torch.ones(3,5)
indices = torch.tensor([[0,1, 2], [0, 2, 3], [0, 3, 4]]).long()
print(ar1)
print(indices)


def convertToIndices(indices):
    nD, nC = indices.shape
    ar = torch.arange(0, nD)
    ar = ar.repeat(nC, 1).t()
    ar = torch.vstack([ar.reshape(-1), indices.reshape(-1)])
    return ar

def setValue(tensor, indices, value):
    tensor[tuple(indices)] = value

indices  = convertToIndices(indices)
print(indices)
ar1[tuple(indices)] = 0
print(ar1)