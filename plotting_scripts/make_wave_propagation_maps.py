#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:42:41 2024

@author: rshimony
"""
'''
Making ground velocity maps taken every 0.5s of SW4 simulation.
These images show the wave propagation pattern in the simulation
Script steps:
    1) Creating WV ploygon
    2) Creating a base map (empty shaded relief region map with WV polygon outlines)
    * Might need environment change at this point
    3) Saving SW4 images data to .csv files
    4) Creating new colormap for plotting ground velocity data
    5) Plotting ground velocity figures
'''
import pygmt
import pandas as pd
import numpy as np
from shapely.geometry.polygon import Polygon

## Creating WV basin polygon
inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/basin_depth_poly_outer.csv')
lon_outer = outer_f['wv_poly_lons']
lat_outer = outer_f['wv_poly_lats']   
  
inner1_poly = np.zeros((len(lon_inner1),2))
for i in range(len(lon_inner1)):
    inner1_poly[i,0] = lon_inner1[i]
    inner1_poly[i,1] = lat_inner1[i]
    
inner2_poly = np.zeros((len(lon_inner2),2))
for i in range(len(lon_inner2)):
    inner2_poly[i,0] = lon_inner2[i]
    inner2_poly[i,1] = lat_inner2[i]
    
outer_poly = np.zeros((len(lon_outer),2))
for i in range(len(lon_outer)):
    outer_poly[i,0] = lon_outer[i]
    outer_poly[i,1] = lat_outer[i]

full_poly = Polygon(outer_poly, [inner1_poly , inner2_poly])

xpn_m, ypn_m = full_poly.exterior.xy 
xph1_m, yph1_m = full_poly.interiors[0].xy
xph2_m, yph2_m = full_poly.interiors[1].xy

#%%
## Creating a base map the SW4 images to be plotted on top
# This is an empty shaded relief regional map, with only the outline of the WV plotted.
plot_region = [-124.19,-121.38,43.4,46.8]
fig = pygmt.Figure()
topo_data = '@earth_relief_03s'

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):       
    fig.grdimage(topo_data,region=plot_region, projection="M15c", cmap="grey",shading=True,frame=None)
    
    fig.plot(x=xpn_m, y=ypn_m,pen="0.8p,black")
    fig.plot(x=xph1_m, y=yph1_m,pen="0.8p,black")
    fig.plot(x=xph2_m, y=yph2_m,pen="0.8p,black")

    fig.coast(water="whitesmoke")

    fig.show()
    fig.savefig('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/map4waveprop_final.png',dpi=600)

#%%
'''
MIGHT NEED AN EVIRONMENT CHANGE BEFORE RUNNING THIS SECTION
Since PyGMT often can't live in the same environment as pySW4
(pySW4 is the code package used to read and analyze SW4 outputs)
'''
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from glob import glob
import pySW4 as sw4

#%%
## Saving SW4 images data to .csv files for easier manipulation of data for plotting
# Reading SW4 images (only '.z=0' images. These are 'map view' images recording ground velocity every 0.5 second)
vel_im_list = np.array(sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/springfield_steph_freq8_mags/*.z=0.magdudt.sw4img')))
## Output csv files directory
csv_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/springfield_steph_freq8_csv_data/'
## Converting SW4 images to csv files and saving
for i in range(len(vel_im_list)):
    imagefile = sw4.read_image(vel_im_list[i])
    
    clcle_no = vel_im_list[i].split('.')[-4][-4:]
    t=str(round((float(clcle_no)*0.0184366),1))
    
    print(clcle_no,t)
    
    patch=imagefile.patches[0]
    matrix=patch.data
  
    pd.DataFrame(matrix[::-1]).to_csv(csv_dir+t+'s_velmag.csv', header=None, index=None)

#%%
## Creating new colormap for plotting SW4 ground velocity data.
## This colormap is transparent when data value is 0.
cdict = {'red':  ((0., 1, 1), (0.03, 1, 1), (0.20, 0, 0), (0.66, 1, 1), (0.89, 1, 1), (1, 0.5, 0.5)),
          'green':((0., 1, 1), (0.03, 1, 1), (0.20, 0, 0), (0.375, 1, 1), (0.64, 1, 1), (0.91, 0, 0), (1, 0, 0)),
          'blue': ((0., 1, 1), (0.08, 1, 1), (0.20, 1, 1), (0.34, 1, 1), (0.65, 0, 0), (1, 0, 0))}
whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)

# get colormap and create a colormap object
ncolors = 256
color_array = plt.get_cmap(whitejet)(range(ncolors))
color_array[:,-1] = np.linspace(0.0,1.0,ncolors) # change alpha values
map_object = LinearSegmentedColormap.from_list(name='whitejet_alpha',colors=color_array)
plt.register_cmap(cmap=map_object) # register this new colormap with matplotlib


#%%  
## Plotting SW4 ground velocity images (image taken every 0.5s throughout the simulation)
# Setting plotting region
plot_region = [-124.19,-121.38,43.4,46.8]
# Setting base map image
base_image= '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/map4waveprop_final.png'
# Reading images .csv files (holds SW4 ground velocity data)
image_data_ls = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/springfield_steph_freq8_csv_data/*'))
# Setting ground velocity figures directory
mag_figs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/springfield_steph_freq8_mags_figs/'
# Plotting images and saving
for i in range(len(image_data_ls)):
    time = image_data_ls[i].split('/')[-1].split('_')[0][:-1]
    if len(time) == 3:
        time = '0' + time
    
    print(time)
    
    fig, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True
    im = plt.imread(base_image)

    im = ax.imshow(im, extent=plot_region,alpha=1)

    image_data = pd.read_csv(image_data_ls[i])
    im=ax.imshow(image_data[:-50],extent=plot_region, cmap='whitejet_alpha',vmin=0 ,vmax=0.0001)
    
    cbar=fig.colorbar(im,fraction=0.03, pad=0.02)
    cbar.set_label('vel [m/s]') 
    
    fontsize = 12
    ax.set_aspect(1.4) 
    ax.set_xticks([-124,-123,-122])
    ax.set_yticks([44,45,46])
    ax.set_ylabel('latitude',fontsize=fontsize)
    ax.set_xlabel('longitude',fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax.set_title('time = '+time+' s',fontsize=fontsize+5)
    ax.axis('off')
    
    plt.close(fig)
    fig.savefig(mag_figs_dir+time+'s.png',dpi=300, bbox_inches='tight',facecolor='white', edgecolor='none')
        

















        
        
        
        
        
        
        
        
        
        
        















        
        
        
        
        
        
        
        
        
        
        