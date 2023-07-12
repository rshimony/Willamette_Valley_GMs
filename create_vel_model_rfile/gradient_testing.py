#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 16:20:38 2023

@author: rshimony
"""

import matplotlib.pyplot as plt
import sys
import os

flush = sys.stdout.flush()

from pySW4.utils import geo
from pySW4.prep import rfileIO
import numpy as np
import time
import math
from warnings import warn
import xarray as xr
from scipy.interpolate import griddata
import pandas as pd
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point,LineString

#%%

vel_model = np.array([\
[0,   -999,  -999,   -999, -999,-999],
[1,  1.987e3, 0.6e3, 1.899e3, 78.34,  39.17],
[2,  1e20, 1e20, 1e20, 1e20,  1e20]])

vel_model_grad = np.array([\
[0,   0,  0,   0, 0,0],
[1,  0.614, 0.487, 0.119, 0,  0],
[2,  0, 0, 0, 0, 0]])
    
hvi_ = 100

layr_block_index = np.array([0,0,1,1,1,1,2,2,2,2,2,2,2,2])



unq_counts = np.unique(layr_block_index,return_counts=True)[1]

idx_ar = np.array([])
for i in unq_counts:
    id_ar = np.arange(i)
    idx_ar = np.concatenate((idx_ar,id_ar*hvi_))

grad_ar = np.zeros([len(layr_block_index),5])
j=0
for i in [3,1,2,4,5]:
    grad_ar[:,j] = vel_model_grad[layr_block_index,i]
    j=j+1
    
k_array_grad = grad_ar*idx_ar[:, np.newaxis]   
    
    
    
# np.unique(vs_depth,return_counts=True)
# Out[64]: 
# (array([  400.        ,   720.22766667,  1040.45533333,  1360.683     ,
#          1680.91066667,  2001.13833333,  2321.366     ,  2641.59366667,
#          2961.82133333,  3282.049     ]),
#  array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))

# np.unique(vs_depth,return_counts=True)[1]
# Out[65]: array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# np.unique(vs_depth,return_counts=True)[1][0]
# Out[66]: 1
