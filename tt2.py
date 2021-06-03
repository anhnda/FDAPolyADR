import torch
import numpy as np
ar = [i for i in range(10)]
ar = np.tile(ar,100)
v  = np.random.choice(ar, (100, 3), replace=False)
print(v.shape)
print(v)