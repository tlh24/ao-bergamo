from __future__ import division, print_function
from time import time
from os.path import join
from itertools import product
import unittest
import sys
import os
import tempfile
import shutil
import gc   

import numpy as np
import h5py
import matplotlib.pyplot as plt
import math
from math import sqrt
import pdb
import pickle

import torch

from torch import nn, optim
from torch.nn import init
from torch.nn import functional as F
from torch.autograd import Function

# from pympler.tracker import SummaryTracker

class EqualLR:
	def __init__(self, name):
		self.name = name

	def compute_weight(self, module):
		weight = getattr(module, self.name + '_orig')
		fan_in = weight.data.size(1) * weight.data[0][0].numel()
		# print("compute_weight ", weight.size(), 'fan_in', fan_in)

		return weight * sqrt(2 / fan_in)
	
	# Idea: rather than changing the intialization of the weights, 
	# here they just use N(0, 1), and dynamically scale the weights by 
	# sqrt(2 / fan_in), per He 2015. 
	# gradient updates are applied to the unscaled weights (checked), 
	# which means that the gradient updates are also scaled by sqrt(2/fan_in)

	@staticmethod
	def apply(module, name):
		fn = EqualLR(name)

		weight = getattr(module, name)
		del module._parameters[name]
		module.register_parameter(name + '_orig', nn.Parameter(weight.data))
		module.register_forward_pre_hook(fn)

		return fn

	def __call__(self, module, input):
		weight = self.compute_weight(module)
		setattr(module, self.name, weight)

def equal_lr(module, name='weight'):
    EqualLR.apply(module, name)
    return module

class EqualLinear(nn.Module):
	def __init__(self, in_dim, out_dim):
		super().__init__()

		linear = nn.Linear(in_dim, out_dim)
		linear.weight.data.normal_()
		linear.bias.data.zero_()

		self.linear = equal_lr(linear)

	def forward(self, input):
		return self.linear(input)
	  
	def compute_weight(self):
		w = self.linear.state_dict()['weight_orig']
		fan_in = w.data.size(1) * w.data[0][0].numel() 
		return w * (2 / sqrt(fan_in))

device = torch.device('cuda:0')
batch_size = 16 # 64 converges more slowly, obvi
# batch_size = 32, after 50k steps validation loss is 0.000206
# 16, validation loss is 0.000191
# layers: in, 1024, 512, 256, 97: validation loss 7e-5, 1.2e-4, 8.7e-5
# layers: in, 1024, 1024 512, 256, 97: validation loss 1.3e-4
# layers: in, 512, 256, 97: 1e-4
# straight nn.Linear, no EqualLinear, 500k train steps before: 1.7e-3
# EqualLinear, LR = 1e-4, decay = 2e-4, 2M steps, betas (0.0 0.9), validation loss = 3.4e-5.
# LR = 5e-5, weight_decay = 1e-4, 2.5M steps, betas (0.0 0.9), validation loss 3.5e-5, 
# -- EqualLinear really works! 
ncentroids = 1075
analyzer = nn.Sequential(
	EqualLinear(ncentroids, 1024), 
	nn.LeakyReLU(0.2), 
	EqualLinear(1024, 512), 
	nn.LeakyReLU(0.2), 
	EqualLinear(512, 256), 
	nn.LeakyReLU(0.2), 
	EqualLinear(256, 97)).cuda(device)

optimizer = optim.AdamW(analyzer.parameters(), lr=1e-4, betas=(0.9, 0.999), weight_decay=2e-4)
lossfunc = torch.nn.SmoothL1Loss() # mean reduction

niters = 2500000
slowloss = 0.0
losses = np.zeros((2,niters)) 
nvalidate = 100000

f = h5py.File('../rundata/centroids_full_cleaned.mat', 'r')
A = f.get('A')[:]
v = f.get('v')[:]
VS = f.get('VS')[:]
flatwf = f.get('flatwf')[:]
A = torch.from_numpy(A.astype(np.single)) 
v = torch.from_numpy(v.astype(np.single))
VS = torch.from_numpy(VS.astype(np.single))
flatwf = torch.from_numpy(flatwf.astype(np.single))
(nsamp, ncentroids) = A.shape;
print('VS shape' , VS.shape)
print('flatwf shape' , flatwf.shape)
f.close()
print(f'ncentroids:{ncentroids}'); 

# save the input matrix to convert from SVD modes to wavefronts. 
# (this still will rely on the system flat !!! )
VS = torch.transpose(VS, 0, 1) # pytorch seems to be left-multiply
# (input is the second dimension)
flatwf = torch.squeeze(flatwf) # matlab seems to add that second dimension.
np.save('dmcontrolnet/VS.weight.npy', VS.numpy(), allow_pickle=False)
np.save('dmcontrolnet/VS.bias.npy', flatwf.numpy(), allow_pickle=False)


for k in range(niters):
	r = torch.randint(nvalidate, nsamp, (batch_size,)) # on the cpu
	x = A[r, :].cuda(device) # copy to gpu
	y = v[r, :].cuda(device)
	x.grad = None; 
	y.grad = None; 
	y = y / 0.15 # scale to +- 1.0; don't forget to un-scale later!
	analyzer.zero_grad()
	predict = analyzer(x)
	loss = lossfunc(y, predict)
	loss.backward()
	optimizer.step()
	slowloss = 0.99*slowloss + 0.01 * loss.detach()
	if k % 200 == 0 : 
		print(f'{k} loss: {loss}; slowloss {slowloss}')
	losses[0,k] = loss
	losses[1,k] = slowloss

sumloss = 0.0
with torch.no_grad():
	for k in range(0, nvalidate, batch_size):
		r = torch.arange(k, k+batch_size, 1, dtype=torch.int64)
		x = A[r, :].cuda(device) # copy to gpu
		y = v[r, :].cuda(device)
		y = y / 0.15 # scale to +- 1.0
		predict = analyzer(x)
		sumloss = sumloss + lossfunc(y, predict)

sumloss = sumloss / (nvalidate / batch_size)
print(f'validation loss: {sumloss}')

# this muse be called after a gradient-free forward pass
# which call the __call__ function of EqualLR class, 
# thereby properly scaling the weights. 
d = {}
d['0.weight'] = analyzer[0].linear.weight
print('0.weight shape', analyzer[0].linear.weight.shape)
d['0.bias'] = analyzer[0].linear.bias
d['2.weight'] = analyzer[2].linear.weight
d['2.bias'] = analyzer[2].linear.bias
d['4.weight'] = analyzer[4].linear.weight
d['4.bias'] = analyzer[4].linear.bias
d['6.weight'] = analyzer[6].linear.weight
d['6.bias'] = analyzer[6].linear.bias

for key in d:
	val = d[key].cpu().detach()
	np.save('dmcontrolnet/' + key + '.npy', val.numpy(), allow_pickle=False)
	#torch.save(d[key], 'dmcontrolnet/' + key + '.pt')

dmcontrolnet = nn.Sequential(
	nn.Linear(ncentroids, 1024), 
	nn.LeakyReLU(0.2), 
	nn.Linear(1024, 512), 
	nn.LeakyReLU(0.2), 
	nn.Linear(512, 256), 
	nn.LeakyReLU(0.2), 
	nn.Linear(256, 97)).cuda(device)

dmcontrolnet.load_state_dict(d)
# dmcontrolnet.load_state_dict(analyzer.state_dict())

# check this simplified network is working.
sumloss = 0.0
with torch.no_grad():
	for k in range(0, nvalidate, batch_size):
		r = torch.arange(k, k+batch_size, 1, dtype=torch.int64)
		x = A[r, :].cuda(device) # copy to gpu
		y = v[r, :].cuda(device)
		y = y / 0.15 # scale to +- 1.0
		dmcontrolnet.zero_grad()
		predict = dmcontrolnet(x)
		l = lossfunc(y, predict)
		#print(predict)
		sumloss = sumloss + l

sumloss = sumloss / (nvalidate / batch_size)
print(f'validation loss: {sumloss}')
