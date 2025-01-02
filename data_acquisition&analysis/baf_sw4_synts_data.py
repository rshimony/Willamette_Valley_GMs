#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 12:11:33 2024

@author: rshimony
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import pySW4 as sw4
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

'''
SW4 synthetics data Basin Amplification Factor (BAF) calculation and plotting.
This script reads in PGV images from SW4 simulations (hmax files),
calculates PGV ratios between the selected new WV model (SRV 1D) and the USGS CVM,
calculate BAF value from these and plots map of PGV (shaking) ratio with BAF value.
'''

#%%

## Creating WV polygon - cartesian, to fit SW4 output data
inner1_cart_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/poly_cart_inner1.csv')
x_inner1 = inner1_cart_f['inner1_x']
y_inner1 = inner1_cart_f['inner1_y']

inner2_cart_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/poly_cart_inner2.csv')
x_inner2 = inner2_cart_f['inner2_x']
y_inner2 = inner2_cart_f['inner2_y']

outer_cart_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/poly_cart_outer.csv')
x_outer = outer_cart_f['outer_x']
y_outer = outer_cart_f['outer_y']   
  
inner1_cart_poly = np.zeros((len(x_inner1),2))
for i in range(len(x_inner1)):
    inner1_cart_poly[i,0] = x_inner1[i]
    inner1_cart_poly[i,1] = y_inner1[i]
    
inner2_cart_poly = np.zeros((len(x_inner2),2))
for i in range(len(x_inner2)):
    inner2_cart_poly[i,0] = x_inner2[i]
    inner2_cart_poly[i,1] = y_inner2[i]
    
outer_cart_poly = np.zeros((len(x_outer),2))
for i in range(len(x_outer)):
    outer_cart_poly[i,0] = x_outer[i]
    outer_cart_poly[i,1] = y_outer[i]

full_cart_poly = Polygon(outer_cart_poly, [inner1_cart_poly , inner2_cart_poly])


## Reading in events catalog and extracting event location
events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/event_catalog_valevents.csv')
scottsmills_lon = events['longitude'][0]
scottsmills_lat = events['latitude'][0]
salem_lon = events['longitude'][1]
salem_lat = events['latitude'][1]
springfield_lon = events['longitude'][2]
springfield_lat = events['latitude'][2]

## Setting plot region
plot_region = [-124.19,-121.38,43.4,46.8]

## Base map image for the SW4 images to be plotted on top
# This is an empty shaded relief regional map, with only the outline of the WV plotted.
base_image= '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/map4waveprop_final.png'

## Reading SW4 images. These images (.cycle=4014.z=0.hmaxdudt.sw4img) are taken by SW4 at the last cycle of the simulation
#   and showing the maximun velocity recorded at each gridpoint thoughout the simulation (basically acting as a PGV image)
# Reading data from the SRV 1D model (the selected WV model) and USGS CVM, from all events:
salem_srv_1d_hmax = sw4.read_image('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/salem_srv_1d_maxdudt/image.cycle=4014.z=0.hmaxdudt.sw4img')
salem_steph_hmax = sw4.read_image('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/salem_steph_maxdudt/image.cycle=4014.z=0.hmaxdudt.sw4img')

scottmills_srv_1d_hmax = sw4.read_image('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/scottsmills_srv_1d_maxdudt/image.cycle=4014.z=0.hmaxdudt.sw4img')
scottmills_steph_hmax = sw4.read_image('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/scottsmills_steph_maxdudt/image.cycle=4014.z=0.hmaxdudt.sw4img')

springfield_srv_1d_hmax = sw4.read_image('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/springfield_srv_1d_maxdudt/image.cycle=4014.z=0.hmaxdudt.sw4img')
springfield_steph_hmax = sw4.read_image('/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/images/springfield_steph_maxdudt/image.cycle=4014.z=0.hmaxdudt.sw4img')

## Setting output directory
baf_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/baf_shakemaps/sw4_synts/'

#%%
## Extracting images data
salem_srv_1d_hmax_data = salem_srv_1d_hmax.patches[0].data
salem_steph_hmax_data = salem_steph_hmax.patches[0].data

scottmills_srv_1d_hmax_data = scottmills_srv_1d_hmax.patches[0].data
scottmills_steph_hmax_data = scottmills_steph_hmax.patches[0].data

springfield_srv_1d_hmax_data = springfield_srv_1d_hmax.patches[0].data
springfield_steph_hmax_data = springfield_steph_hmax.patches[0].data

## Calculating shaking ratio (for BAF)
baf_salem = salem_srv_1d_hmax_data/salem_steph_hmax_data
baf_scottsmills = scottmills_srv_1d_hmax_data/scottmills_steph_hmax_data
baf_springfield = springfield_srv_1d_hmax_data/springfield_steph_hmax_data

## Saving shaking ratios to .csv files
pd.DataFrame(baf_salem[::-1]).to_csv(baf_dir+'salem_baf.csv', header=None, index=None)
pd.DataFrame(baf_scottsmills[::-1]).to_csv(baf_dir+'scottsmills_baf.csv', header=None, index=None)
pd.DataFrame(baf_springfield[::-1]).to_csv(baf_dir+'springfield_baf.csv', header=None, index=None)

#%%
## Calculating WV BAF value from SW4 images
def sw4_baf_calculate(baf_array):
    '''
    This function takes in SW4 shaking ratio (baf) array,
     extracts ONLY WV data points and calculates the mean shaking ratio value of the WV,
     this value is BAF of the WV.
    '''
    y_baf_data , x_baf_data = baf_array.shape
    
    counter = 1
    
    baf_data_poly_vel = []
    baf_data_poly_x = []
    baf_data_poly_y = []
    
    for y in range(y_baf_data):
        if counter % 100 == 0:
            print(str(counter) + ' from ' + str(y_baf_data))
        for x in range(x_baf_data):
            
            grd_point = Point(x*100,y*100)
        
            if full_cart_poly.contains(grd_point) == True:
                baf_data_poly_vel.append(baf_array[y,x])
                baf_data_poly_x.append(x)
                baf_data_poly_y.append(y)
                
        counter = counter + 1 
    
    mean_baf = np.array(baf_data_poly_vel).mean()
    
    return mean_baf

## Initianting the function for all events:
salem_baf = sw4_baf_calculate(baf_salem)
scottsmills_baf = sw4_baf_calculate(baf_scottsmills)
springfield_baf = sw4_baf_calculate(baf_springfield)

#%%
## Plotting SW4 BAF maps:
def plot_sw4_baf_map(baf_data_file,baf_value,eq_name,evlon,evlat):
    '''
    SW4 BAF maps plotting function. Plots BAF for a SINGLE EVENT at a time.
    Inputs:
        - baf_data_file: .csv file that holds the shaking ratio values of all grid points
        - baf_value: BAF (mean shaking ratio value) of the WV (calculated in previous section).
        - eq_name: Event name (string)
        - evlon, evlat: Event location
    '''
    fontsize = 8
    fig, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True
    im = plt.imread(base_image)
    
    im = ax.imshow(im, extent=plot_region,alpha=1)
    
    baf_data = pd.read_csv(baf_data_file)
    im=ax.imshow(baf_data[:-50],extent=plot_region, cmap='Reds',vmin=1, vmax=6,alpha=0.6)
    
    cbar=fig.colorbar(im,fraction=0.03, pad=0.02)
    cbar.set_label('Max Vel Ratio [SRV 1D/USGS CVM]',fontsize=fontsize+4) 
    cbar.ax.tick_params(labelsize=fontsize+5)
    
    ax.set_aspect(1.4) 
    ax.set_xticks([-124,-123,-122])
    ax.set_yticks([44,45,46])
    ax.tick_params(labelsize=fontsize+8)
    ax.text(-123.0, 43.6, 'BAF = ' + str(baf_value)[:4] , color='k', fontsize=fontsize+3, weight='bold', style='italic')
    ax.text(-124.0, 46.6, eq_name+' EQ' , color='k', fontsize=fontsize+4, weight='bold')
    ax.scatter(evlon,evlat,s=60,c='gold',edgecolor='k',marker='*', linewidths=0.5)
    
    plt.close(fig)
    fig.savefig(baf_dir+eq_name+'_baf.png',dpi=300, bbox_inches='tight',facecolor='white', edgecolor='none')

## Making plotts:
# Salem
salem_baf_data_file = baf_dir + 'salem_baf.csv'
plot_sw4_baf_map(salem_baf_data_file,salem_baf,'Salem',salem_lon,salem_lat)
# Scotts Mills
scottsmills_baf_data_file = baf_dir + 'scottsmills_baf.csv'
plot_sw4_baf_map(scottsmills_baf_data_file,scottsmills_baf,'Scotts Mills',scottsmills_lon,scottsmills_lat)
# Springfield
springfield_baf_data_file = baf_dir + 'springfield_baf.csv'
plot_sw4_baf_map(springfield_baf_data_file,springfield_baf,'Springfield',springfield_lon,springfield_lat)










































