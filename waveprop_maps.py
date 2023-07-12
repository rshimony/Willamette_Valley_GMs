#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 13:53:38 2023

@author: rshimony
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from glob import glob

#%%
plot_region = [-124.19,-121.51,43.4,46.1]

outputfilename= '/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/map4waveprop.png'

image_data_ls = sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/gp1_deep_1d/mag*.csv'))

image_data_ls_short = [image_data_ls[10],image_data_ls[36],image_data_ls[86],image_data_ls[114]]

# cdict = {'red':  ((0., 1, 1), (0.03, 1, 1), (0.20, 0, 0), (0.66, 1, 1), (0.89, 1, 1), (1, 0.5, 0.5)),
#          'green':((0., 1, 1), (0.03, 1, 1), (0.20, 0, 0), (0.375, 1, 1), (0.64, 1, 1), (0.91, 0, 0), (1, 0, 0)),
#          'blue': ((0., 1, 1), (0.08, 1, 1), (0.20, 1, 1), (0.34, 1, 1), (0.65, 0, 0), (1, 0, 0))}
# whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)

#%%   
 
for i in range(len(image_data_ls_short)):
    time = image_data_ls_short[i].split('/')[-1].split('_')[1][:4]
    
    print(time)
    
    fig, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True
    im = plt.imread(outputfilename)

    im = ax.imshow(im, extent=plot_region,alpha=1)

    # get colormap and create a colormap object
    # ncolors = 256
    # color_array = plt.get_cmap(whitejet)(range(ncolors))
    # color_array[:,-1] = np.linspace(0.0,1.0,ncolors) # change alpha values
    # map_object = LinearSegmentedColormap.from_list(name='whitejet_alpha',colors=color_array)
    # plt.register_cmap(cmap=map_object) # register this new colormap with matplotlib

    image_data = pd.read_csv(image_data_ls_short[i])
    im=ax.imshow(image_data,extent=plot_region, cmap='whitejet_alpha',vmin=0 ,vmax=0.00003)#0.1,vmax=0.00015, cmap=cmap) #2
    
    cbar=fig.colorbar(im,fraction=0.03, pad=0.02) #0.02
    cbar.set_label('vel [m/s]') #cbar.ax.tick_params(labelsize=30)
    
    fontsize = 12
    ax.set_aspect(1.4) 
    ax.set_xticks([-124,-123,-122])
    ax.set_yticks([44,45,46])
    ax.set_ylabel('latitude',fontsize=fontsize)
    ax.set_xlabel('longitude',fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax.set_title('time = '+time+' s',fontsize=fontsize+5)
    ax.axis('off')
    #ax.text(134, 44.3, 'time = '+t+' s', color='r', fontsize=fontsize, weight='bold')
    
    plt.close(fig)
    fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/pngs/mag_gp1_'+time+'s.png',dpi=300,\
        bbox_inches='tight',facecolor='white', edgecolor='none')
        
 #%%

time = image_data_ls[50].split('/')[-1].split('_')[1][:3]

print(time)

fig, ax = plt.subplots()
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
im = plt.imread(outputfilename)

im = ax.imshow(im, extent=plot_region,alpha=1)

# get colormap and create a colormap object
# ncolors = 256
# color_array = plt.get_cmap(whitejet)(range(ncolors))
# color_array[:,-1] = np.linspace(0.0,1.0,ncolors) # change alpha values
# map_object = LinearSegmentedColormap.from_list(name='whitejet_alpha',colors=color_array)
# plt.register_cmap(cmap=map_object) # register this new colormap with matplotlib

image_data = pd.read_csv(image_data_ls[50])
im=ax.imshow(image_data,extent=plot_region, cmap='whitejet_alpha', vmin=0, vmax=0.00003)#0.1,vmax=0.00015, cmap=cmap) #2

cbar=fig.colorbar(im,fraction=0.03, pad=0.02) #0.02
cbar.set_label('mag displ') #cbar.ax.tick_params(labelsize=30)

fontsize = 12
ax.set_aspect(1.5) 
ax.set_xticks([-124,-123,-122])
ax.set_yticks([44,45,46])
ax.set_ylabel('latitude',fontsize=fontsize)
ax.set_xlabel('longitude',fontsize=fontsize)
ax.tick_params(labelsize=fontsize)
ax.set_title('time = '+time+' s',fontsize=fontsize+5)
ax.axis('off')
#ax.text(134, 44.3, 'time = '+t+' s', color='r', fontsize=fontsize, weight='bold')

# plt.close(fig)
plt.show()
# fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/gp1_1d/png/mag_'+time+'s.png',dpi=500,\
    # bbox_inches='tight',facecolor='white', edgecolor='none')       
        
#%%
fig, ax = plt.subplots()
image_data = pd.read_csv(image_data_ls[50])
im=ax.imshow(image_data,extent=plot_region, cmap='jet', vmin=0, vmax=0.08,alpha=0.8)      
plt.show()       
        
        
        
        
        
        
        
        
        
        
        
        
        