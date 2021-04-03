import torch

B = 5  # Batch size
C = 6  # Num chanels
D = 4  # W dim
N = 4  # Input size
K = 3  # Feature size
PD = 2  # PADING SIZE
NP = N + PD  # INPUT WITH PADDING

inp = torch.randn(B, 1, N, K)
pad = torch.nn.ZeroPad2d((0, 0, 0, PD))
xLen = torch.LongTensor([2, 3, 4, 4, 4])

ipad = pad(inp)
print(ipad.shape)

w = torch.rand((C, D, 1), requires_grad=True)
V = torch.rand((C, K, D), requires_grad=True)

re = torch.matmul(ipad, V)
# B * C * NP * D

print(re.shape)
re = torch.tanh(re)

re = torch.matmul(re, w)

idx = torch.arange(NP).unsqueeze(0).unsqueeze(0).repeat((B, C, 1))
print(idx.shape)
# print(idx)
len_expanded = xLen.unsqueeze(0).unsqueeze(0)
print(len_expanded.shape)

re = re.squeeze()
print(re.shape)

len_expanded = xLen.unsqueeze(-1).unsqueeze(-1).expand(re.size())

# print(len_expanded.shape)
# print(len_expanded)

mask = len_expanded > idx
# print(mask)

re[~mask] = float('-inf')
# print(re)




re2 = torch.softmax(re, dim=-1)

print(re2.shape)
print(re2)

re2 = torch.unsqueeze(re2, -1)
print(re2.shape)
v = torch.mul(ipad, re2)
print(v)
print(v.shape)
v = torch.sum(v, dim=2)
print(v)
print(v.shape)
v = torch.unsqueeze(v, 1)
print(v.shape)

v = torch.squeeze(v)
v = v.reshape(v.size(0), -1)
print(v.shape)