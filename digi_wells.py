#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:36:46 2023

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

#%%

ira_baker_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/Ira_Baker/IraB_digitized.csv')
porter_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/Porter/Porter1_digitized.csv')
mpf_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/MP_Farms/MPF_digitized.csv')
hjm_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/HJ_Miller/HJM_digitized.csv')
wolv_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/Wolverton/Wolv_digitized.csv')
ind_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/Independence/Ind_digitized.csv')
bruer_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/Bruer/Bruer_digitized.csv')
klohs_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Wells_Data/Klohs/Klohs_digitized.csv')

ira_baker_depth = ira_baker_f['m']
porter_depth = porter_f['m']
mpf_depth = mpf_f['m']
hjm_depth = hjm_f['m']
wolv_depth = wolv_f['m']
ind_depth = ind_f['m']
bruer_depth = bruer_f['m']
klohs_depth = klohs_f['m']

ira_baker_val = ira_baker_f['Vs']
porter_val = porter_f['Vs']
mpf_val = mpf_f['Vs']
hjm_val = hjm_f['Vs']
wolv_val = wolv_f['Vs']
ind_val = ind_f['Vs']
bruer_val = bruer_f['Vs']
klohs_val = klohs_f['Vs']


#%%

print(len(ira_baker_depth))
print(len(porter_depth))
print(len(mpf_depth))
print(len(hjm_depth))
print(len(wolv_depth))
print(len(ind_depth))
print(len(bruer_depth))
print(len(klohs_depth))
print(len(wolv_depth))

#%%

ira_baker_depth_min = ira_baker_f['m'][:len(wolv_depth)]
porter_depth_min = porter_f['m'][:len(wolv_depth)]
mpf_depth_min = mpf_f['m'][:len(wolv_depth)]
hjm_depth_min = hjm_f['m'][:len(wolv_depth)]
wolv_depth_min = wolv_f['m'][:len(wolv_depth)]
ind_depth_min = ind_f['m'][:len(wolv_depth)]
bruer_depth_min = bruer_f['m'][:len(wolv_depth)]
klohs_depth_min = klohs_f['m'][:len(wolv_depth)]

ira_baker_val_min = ira_baker_f['Vs'][:len(wolv_depth)]
porter_val_min = porter_f['Vs'][:len(wolv_depth)]
mpf_val_min = mpf_f['Vs'][:len(wolv_depth)]
hjm_val_min = hjm_f['Vs'][:len(wolv_depth)]
wolv_val_min = wolv_f['Vs'][:len(wolv_depth)]
ind_val_min = ind_f['Vs'][:len(wolv_depth)]
bruer_val_min = bruer_f['Vs'][:len(wolv_depth)]
klohs_val_min = klohs_f['Vs'][:len(wolv_depth)]

wells_val_avg = (ira_baker_val_min + porter_val_min + mpf_val_min + hjm_val_min + wolv_val_min + 
                 ind_val_min + bruer_val_min + klohs_val_min)/8

plt.figure(figsize = (4,8))
plt.plot(wells_val_avg,ira_baker_depth_min,linewidth=0.4,color='k')
plt.gca().invert_yaxis()
plt.show()

#%%

sort_vals = np.sort(wells_val_avg)

coef = np.polyfit(sort_vals,ira_baker_depth_min,1)
poly1d_fn = np.poly1d(coef)

y = coef[0]* wells_val_avg + coef[1]

plt.figure(figsize = (4,8))
plt.plot(wells_val_avg,ira_baker_depth_min,linewidth=0.4,color='k')
# plt.plot(poly1d_fn(sort_vals), '--r')
plt.gca().invert_yaxis()
plt.show()

#'m*x+b'

#%%

res = stats.linregress(wells_val_avg, ira_baker_depth_min)
y = res.slope*wells_val_avg + res.intercept

plt.figure(figsize = (4,8))
plt.plot(wells_val_avg,y,linewidth=0.4,color='k')
# plt.plot(poly1d_fn(sort_vals), '--r')
plt.gca().invert_yaxis()
plt.show()

vals_arr = np.arange(0.9,3.9,0.01)
y_long = res.slope*vals_arr + res.intercept

plt.figure(figsize = (4,8))
plt.plot(vals_arr,y_long,linewidth=0.4,color='k')
# plt.plot(poly1d_fn(sort_vals), '--r')
plt.gca().invert_yaxis()
plt.show()

#%%

plt.figure(figsize = (6,9))
plt.plot(ira_baker_val,ira_baker_depth,linewidth=0.4,alpha=0.6,label='Stacked Wells Data')
plt.plot(porter_val,porter_depth,linewidth=0.4,alpha=0.6)
plt.plot(mpf_val,mpf_depth,linewidth=0.8,alpha=0.6)
plt.plot(hjm_val,hjm_depth,linewidth=0.8,alpha=0.6)
plt.plot(wolv_val,wolv_depth,linewidth=0.8,alpha=0.6)
plt.plot(ind_val,ind_depth,linewidth=0.8,alpha=0.6)
plt.plot(bruer_val,bruer_depth,linewidth=0.8,alpha=0.6)
plt.plot(klohs_val,klohs_depth,linewidth=0.8,alpha=0.6)
plt.plot(vals_arr,y_long,'r--',linewidth=2.0,label='Data Fit')
plt.gca().invert_yaxis()

plt.xlabel('Vs [m/s]',size=32)
plt.ylabel('Depth [m]',size=32)
plt.xticks([1,2.5,4])
plt.tick_params(labelsize=32)
# plt.legend()

plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/digi_wells.jpg'
                    ,dpi=300,bbox_inches='tight')




















