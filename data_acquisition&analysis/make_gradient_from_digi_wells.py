#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:39:58 2024

@author: rshimony
"""
'''
Reading in digitized wells .csv files, 
calculating a linear velocity gradient of the averaged velocity profile,
extrapolating the gradient to the deepest well log depth,
ploting:
    1) A stacked wells data figure with all indivdual profiles, averaged profile and linear gradient.
    2) A single figure for each ditized velocity profile.
    3) A subplot figure with all digitized profiles together.
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

#%%
## Reading in digitized wells files
ira_baker_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/IraB_digitized.csv')
porter_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/Porter1_digitized.csv')
mpf_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/MPF_digitized.csv')
hjm_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/HJM_digitized.csv')
wolv_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/Wolv_digitized.csv')
ind_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/Ind_digitized.csv')
bruer_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/Bruer_digitized.csv')
klohs_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/Klohs_digitized.csv')

#%%
## Inserting a "wells name" column to each well file
ira_baker_f.insert(0, "well_name", ['Ira Baker']*len(ira_baker_f), True)
porter_f.insert(0, "well_name", ['Porter']*len(porter_f), True)
mpf_f.insert(0, "well_name", ['M & P Farms']*len(mpf_f), True)
hjm_f.insert(0, "well_name", ['HJ Miller']*len(hjm_f), True)
wolv_f.insert(0, "well_name", ['Wolverton']*len(wolv_f), True)
ind_f.insert(0, "well_name", ['Independence']*len(ind_f), True)
bruer_f.insert(0, "well_name", ['Bruer']*len(bruer_f), True)
klohs_f.insert(0, "well_name", ['Klohs']*len(klohs_f), True)

## Creating a concatinated flatfile with all digitized well log data.
concat_wells_ff = pd.concat([ira_baker_f, porter_f, mpf_f, hjm_f, wolv_f, ind_f, bruer_f, klohs_f])
concat_wells_ff.to_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/concat_wells_ff.csv')

## Making lists of depths, Vs values and well names for further analysis
depth_ls = [ira_baker_f['m'] , porter_f['m'] , mpf_f['m'] , hjm_f['m'] , wolv_f['m'] , ind_f['m'] , bruer_f['m'] , klohs_f['m']]
val_ls = [ira_baker_f['Vs'] , porter_f['Vs'] , mpf_f['Vs'] , hjm_f['Vs'] , wolv_f['Vs'] , ind_f['Vs'] , bruer_f['Vs'] , klohs_f['Vs']]
name_ls = ['Ira Baker' , 'Porter' , 'M & P Farms' , 'HJ Miller' , 'Wolverton' , 'Independence' , 'Bruer' , 'Klohs']

#%%
## Getting maximum and minimum depths and values each well
max_depth_ls = []
max_val_ls = []
min_val_ls = []
depth_len_ls = []
for i in range(len(depth_ls)):
    max_depth = max(depth_ls[i])
    max_depth_ls.append(max_depth)
    max_val = max(val_ls[i])
    max_val_ls.append(max_val)  
    min_val = min(val_ls[i])
    min_val_ls.append(min_val) 
    depth_len = len(depth_ls[i])
    depth_len_ls.append(depth_len)
min_len_depth_well = depth_ls[np.argmin(depth_len_ls)]
max_depth_well = depth_ls[np.argmax(max_depth_ls)]

## Getting the part, from each well, up to the shallowest well depth to create an averaged profile.
# Depths
ira_baker_depth_min = ira_baker_f['m'][:len(min_len_depth_well)]
porter_depth_min = porter_f['m'][:len(min_len_depth_well)]
mpf_depth_min = mpf_f['m'][:len(min_len_depth_well)]
hjm_depth_min = hjm_f['m'][:len(min_len_depth_well)]
wolv_depth_min = wolv_f['m'][:len(min_len_depth_well)]
ind_depth_min = ind_f['m'][:len(min_len_depth_well)]
bruer_depth_min = bruer_f['m'][:len(min_len_depth_well)]
klohs_depth_min = klohs_f['m'][:len(min_len_depth_well)]
# Vs values
ira_baker_val_min = ira_baker_f['Vs'][:len(min_len_depth_well)]
porter_val_min = porter_f['Vs'][:len(min_len_depth_well)]
mpf_val_min = mpf_f['Vs'][:len(min_len_depth_well)]
hjm_val_min = hjm_f['Vs'][:len(min_len_depth_well)]
wolv_val_min = wolv_f['Vs'][:len(min_len_depth_well)]
ind_val_min = ind_f['Vs'][:len(min_len_depth_well)]
bruer_val_min = bruer_f['Vs'][:len(min_len_depth_well)]
klohs_val_min = klohs_f['Vs'][:len(min_len_depth_well)]
# Shallowest well part from the deepest well 
max_depth_well_norm = max_depth_well[:len(min_len_depth_well)]

## Calculating the averaged profile
wells_val_avg = (ira_baker_val_min + porter_val_min + mpf_val_min + hjm_val_min + wolv_val_min + 
                 ind_val_min + bruer_val_min + klohs_val_min)/len(depth_ls)


#%%
## Calculating the gradient
res = stats.linregress(wells_val_avg*1000, max_depth_well_norm)
vals_arr = np.arange(900,3900,10)
y = res.slope*(wells_val_avg*1000) + res.intercept
y_long = res.slope*vals_arr + res.intercept

#%%
## Plotting stacked wells data figure with all indivdual profiles, averaged profile and linear gradient.
plt.figure(figsize = (6,9))
plt.plot(val_ls[0],depth_ls[0],linewidth=0.4,alpha=0.6,color='gray',label='Stacked Wells Data')
plt.plot(val_ls[1],depth_ls[1],linewidth=0.4,alpha=0.6,color='gray')
plt.plot(val_ls[2],depth_ls[2],linewidth=0.8,alpha=0.6,color='gray')
plt.plot(val_ls[3],depth_ls[3],linewidth=0.8,alpha=0.6,color='gray')
plt.plot(val_ls[4],depth_ls[4],linewidth=0.8,alpha=0.6,color='gray')
plt.plot(val_ls[5],depth_ls[5],linewidth=0.8,alpha=0.6,color='gray')
plt.plot(val_ls[6],depth_ls[6],linewidth=0.8,alpha=0.6,color='gray')
plt.plot(val_ls[7],depth_ls[7],linewidth=0.8,alpha=0.6,color='gray')
plt.plot(wells_val_avg,max_depth_well_norm,linewidth=1.4,color='k',label='Mean Data Profile')
plt.plot(wells_val_avg,y,'r',linewidth=2.0,label='Data Fit')
plt.plot((vals_arr/1000),y_long,'r:',linewidth=2.0,label='Data Fit Extrapolated')
plt.gca().invert_yaxis()
plt.xlabel('Vs [km/s]',size=32)
plt.ylabel('Depth [m]',size=32)
plt.xticks([1,2.5,4])
plt.tick_params(labelsize=32)
plt.legend(fontsize="24",bbox_to_anchor=(1.0, -0.14))

plt.savefig('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/digi_wells/digi_wells_stacked.jpg'
                    ,dpi=300,bbox_inches='tight')

#%%
## Plotting single figure for each ditized velocity profile.
highest_val = max(max_val_ls)
lowest_val = min(min_val_ls)
highest_depth = max(max_depth_well)
lowest_depth = min(max_depth_well)

for i in range(len(depth_ls)):
    
    plt.figure(figsize = (6,9))
    plt.plot(val_ls[i],depth_ls[i],linewidth=0.8,color='k')
    plt.title(name_ls[i] , fontsize = 28)
    plt.xlabel('Vs [m/s]',size=24)
    plt.ylabel('Depth [m]',size=24)
    plt.xticks([1,2,3,4])
    plt.tick_params(labelsize=24)
    plt.xlim(lowest_val-0.25 , highest_val+0.05)
    plt.ylim(lowest_depth-200 , highest_depth+100)
    plt.gca().invert_yaxis()

    plt.savefig('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/digi_wells/'+name_ls[i]+ '.jpg'
                        ,dpi=300,bbox_inches='tight')
    
#%%
## Plotting subplot figure with all digitized profiles together.
fig,axs = plt.subplots(figsize=(16, 12),nrows=2, ncols=4,sharex=True,sharey=True)
axs = axs.flatten()
    
for i in range(len(depth_ls)):
    axs[i].plot(val_ls[i],depth_ls[i],linewidth=0.8,color='k')
    axs[i].set_title(name_ls[i] , fontsize = 28)
    axs[i].set_xticks([1,2,3,4])
    axs[i].tick_params(labelsize=24)
    axs[i].set_xlim(lowest_val-0.25 , highest_val+0.05)
    axs[i].set_ylim(lowest_depth-200 , highest_depth+100)
    axs[i].invert_yaxis()
    
fig.supxlabel('Vs [m/s]',size=26)
fig.supylabel('Depth [m]',size=26,x=0.0)
fig.tight_layout()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/digi_wells/all_digi_wells.jpg'
                        ,dpi=300,bbox_inches='tight')


































