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
        # print("compute_weight ", weight.size())

        return weight * sqrt(2 / fan_in)
    
    # Idea: rather than changing the intialization of the weights, 
    # here they just use N(0, 1), and dynamically scale the weights by 
    # sqrt(2 / fan_in), per He 2015. 
    # updates are applied to the unscaled weights (checked), 
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

device = torch.device('cuda:0')
batch_size = 16 # 64 converges more slowly, obvi
niters = 500000

f = h5py.File('../rundata/centroids_cleaned.mat')
A = f.get('A')[:]
v = f.get('v')[:]
A = torch.from_numpy(A) 
v = torch.from_numpy(v)
(nsamp, ncentroids) = A.shape;
f.close()

analyzer = nn.Sequential(
    EqualLinear(ncentroids, 1024), 
    nn.LeakyReLU(0.2), 
    EqualLinear(1024, 512), 
    nn.LeakyReLU(0.2), 
    EqualLinear(512, 256), 
    nn.LeakyReLU(0.2), 
    EqualLinear(256, 97)).cuda(device)

optimizer = optim.Adam(analyzer.parameters(), lr=0.001, betas=(0.0, 0.99))
lossfunc = torch.nn.SmoothL1Loss() # mean reduction

slowloss = 0.0
losses = np.zeros((2,niters)) 


for k in range(niters):
    r = torch.randint(0, nsamp, (batch_size,)) # on the cpu
    x = A[r, :].cuda(device) # copy to gpu
    y = v[r, :].cuda(device)
    x.grad = None; 
    y.grad = None; 
    y = y / 0.15 # scale to +- 1.0
    analyzer.zero_grad()
    predict = analyzer(x)
    loss = lossfunc(y, predict)
    loss.backward()
    optimizer.step()
    slowloss = 0.99*slowloss + 0.01 * loss.detach()
    if k % 20 == 0 : 
        print(f'{k} loss: {loss}; slowloss {slowloss}')
        gc.collect()
    losses[0,k] = loss
    losses[1,k] = slowloss


# save the results in a file! 
#scipy.io.savemat('forward_nn.mat', mdict=
                 #{'l1': analyzer.)

#for name, param in analyzer.named_parameters():
    #print('name: ', name)
    #print(type(param))
    #print('param.shape: ', param.shape)
    #print('param.requires_grad: ', param.requires_grad)
    #print('=====')

# wfile = h5py.File('../data/Cforward_nn.mat', analyzer.state_dict)
