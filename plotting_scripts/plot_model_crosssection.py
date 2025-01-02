#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:06:35 2024

@author: rshimony
"""
'''
This script takes SW4 cross section images of the x axis (E-W) from the 4 new WV models,
Creates longitude ticks to replace distance ticks,
Plots E-W cross sections.
This script takes ONLY ONE MODEL at a time.
'''

import matplotlib.pyplot as plt
from glob import glob
import pySW4 as sw4
import numpy as np

#%%
###INPUTS###
# Cross section images
cs_images=sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sw4_ims_for_cs/ims4cs_salem_srv_sm/*=0000.x*.s.*'))
# Output directory
cs_fig_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/models_cs/srv_sm/'

#%%
## Get longitude ticklabels
# Reading a cross section image
filename_ticks=cs_images[1]
image_ticks=sw4.read_image(filename_ticks)
layer_ticks = image_ticks.patches[0]

# Setting longitude array
model_max_lon = -121.3924
model_min_lon = -124.19
model_npts_x = layer_ticks.ni
lon_arr = np.linspace(model_min_lon,model_max_lon,model_npts_x)

# Setting cross section distances array
x_arr = layer_ticks.x

# Longitude ticks for cross section plot
lon_ticks = [-124.0,-123.5,-123.0,-122.5,-122.0,-121.5]

# Replacing corresponding distamces with lons
tick_idx = []
x_ticks = []
for i in lon_ticks:
    tick_id = np.argwhere(lon_arr>i)[0][0]
    tick_idx.append(tick_id)
    x_tick = x_arr[tick_id]
    x_ticks.append(x_tick)
    
# Setting longitude ticklabels
lon_ticklabels = ['124' + u'\N{DEGREE SIGN}' + 'W' , '123.5' + u'\N{DEGREE SIGN}' + 'W' , '123' + u'\N{DEGREE SIGN}' + 'W',
                  '122.5' + u'\N{DEGREE SIGN}' + 'W' , '122' + u'\N{DEGREE SIGN}' + 'W' , '121.5' + u'\N{DEGREE SIGN}' + 'W']

#%%
## Plotting E-W cross sections

for i in range(len(cs_images)):

    cs_ax_name = cs_images[i].split('/')[-1].split('.')[2]
    cs_ax = cs_images[i].split('/')[-1].split('.')[2].split('=')[0]
    
    filename=cs_images[i]
    image=sw4.read_image(filename)
    layer = image.patches[0]
    model=image.patches[0].data
    extent=([layer.extent[0],layer.extent[1],layer.extent[3],layer.extent[2]])
    
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(13,4))

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(lon_ticklabels)
    
    ax.plot([150000,200000],[4500,4500],'k')
    ax.plot([150000,150000],[4500,5000],'k')
    ax.plot([200000,200000],[4500,5000],'k')
    ax.text(164000,5200,'50 km',fontsize=25)
    
    ax.set_xlabel('W-E (km)',fontsize=30)

    ax.set_yticks([0,1000,2000,3000,4000,5000,6000])
    ax.set_yticklabels(['0','1','2','3','4','5','6'])
    ax.tick_params(axis='both', which='major', labelsize=28)
    ax.set_ylabel('Depth (km)',fontsize=30)
    im=ax.imshow(model,extent=extent,cmap='gist_earth',vmin=0,vmax=3000)
    ax.set_ylim(6000,0)
    ax.set_aspect('auto')
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins1 = inset_axes(ax, width="35.%",height='7%',loc='lower left',bbox_to_anchor=(0.01, 0.2, 1, 1),
                       bbox_transform=ax.transAxes, borderpad=2)
    cb=fig.colorbar(im, cax=axins1, orientation="horizontal")
    cb.set_label('Vs [m/sec]',size=22)
    cb.ax.tick_params(labelsize=22) 
    # plt.show()
    fig.savefig(cs_fig_dir + cs_ax_name + '.jpg',dpi=300,bbox_inches='tight')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    