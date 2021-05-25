import torch
import numpy as np
ar1 = torch.ones(3,5)
indices = torch.tensor([[0,1, 2], [0, 2, 3], [0, 3, 4]]).long()

d1 = "1"
d2 = "2"
v = [d1, d2]
vv = sorted(v)
print(vv)
