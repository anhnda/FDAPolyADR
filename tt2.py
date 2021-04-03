import torch
import numpy as np
ar1 = torch.ones(3,5)
indices = torch.tensor([[0,1, 2], [0, 2, 3], [0, 3, 4]]).long()



def convertToIndices(indices):
    nD, nC = indices.shape
    ar = torch.arange(0, nD)
    ar = ar.repeat(nC, 1).t()
    ar = torch.vstack([ar.reshape(-1), indices.reshape(-1)])
    return ar

def setValue(tensor, indices, value):
    indices = convertToIndices(indices)
    tensor[tuple(indices)] = value


v  = []
for i in range(10000):
    vv = []
    for ii in range(4):
        vi = np.zeros(2000)
        vv.append(vi)
    vv = np.asarray(vv)
    v.append(vv)
v = np.asarray(v)
print(v.shape)