#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:30:36 2024

@author: rshimony
"""
'''
Plot maps of wells location. 
This script creates maps of wells used in constructing the new WV model.
Creates two maps: 
    1) Wells that penetrate the Eugene-Spencer formation
    2) Wells that penetrate the SRV formation
Input is "depth to units" flatfile, that holds the wells data.
'''

import pandas as pd
import pygmt
import numpy as np
from shapely.geometry.polygon import Polygon

#%%
## Input flatfile
depth2units_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/depth2units_wells.csv')

## Extracting wells data from flatfile
srv_lon_units = depth2units_f['Longitude'][depth2units_f['srv_depth'].notna()]
srv_lat_units = depth2units_f['Latitude'][depth2units_f['srv_depth'].notna()]
srv_depth = depth2units_f['srv_depth'][depth2units_f['srv_depth'].notna()]
srv_well_name = depth2units_f['SiteName'][depth2units_f['srv_depth'].notna()]

eugspe_lon_units = depth2units_f['Longitude'][depth2units_f['eugspe_depth'].notna()]
eugspe_lat_units = depth2units_f['Latitude'][depth2units_f['eugspe_depth'].notna()]
eugspe_depth = depth2units_f['eugspe_depth'][depth2units_f['eugspe_depth'].notna()]
eugspe_well_name = depth2units_f['SiteName'][depth2units_f['eugspe_depth'].notna()]

## Output directory
output_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wells_maps/'

#%%
## Creating WV polygon
inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/wv_poly.csv')
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
buff_poly = full_poly.buffer(-0.05)

xpn_m, ypn_m = full_poly.exterior.xy 
xph1_m, yph1_m = full_poly.interiors[0].xy
xph2_m, yph2_m = full_poly.interiors[1].xy

#%%

## Plotting EugSpe wells map
topo_data = '@earth_relief_03s'
fig = pygmt.Figure()
region=[-124.28,-121.4,43.4,46.3]
with pygmt.config(FONT_ANNOT_PRIMARY=24):
    pygmt.makecpt(
        cmap='gray',series=[-5000,6000])
    fig.basemap(region=region, projection='M0/0/10c')
    fig.grdimage(topo_data, cmap=True,shading=True,frame=True)
fig.coast(water="dimgray")
#basin outline
fig.plot(x=xpn_m, y=ypn_m,pen="1p,black", label="WV Basin Model Outline")
fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
#wells
pygmt.makecpt(
    cmap='abyss',reverse=True,
    series=[0,4000],no_bg=True)
fig.plot(x=eugspe_lon_units, y=eugspe_lat_units, 
          style="c0.25c", 
          fill=eugspe_depth,
          pen="black",
          transparency=10,
          cmap=True,
          no_clip='r')
with pygmt.config(FONT_ANNOT_PRIMARY=22,FONT_LABEL=22):
    fig.colorbar(position="JMB+w8c/0.7c",frame=["af+lDepth to Unit (m)"]) 
fig.savefig(output_dir + 'eugspe_wells_map.png',dpi=300)



## Plotting SRV wells map
topo_data = '@earth_relief_03s'
fig = pygmt.Figure()
region=[-124.28,-121.4,43.4,46.3]
with pygmt.config(FONT_ANNOT_PRIMARY=24):
    pygmt.makecpt(
        cmap='gray',series=[-5000,6000])
    fig.basemap(region=region, projection='M0/0/10c')
    fig.grdimage(topo_data, cmap=True,shading=True,frame=True)
fig.coast(water="dimgray")
#basin outline
fig.plot(x=xpn_m, y=ypn_m,pen="1p,black", label="WV Basin Model Outline")
fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
#wells
pygmt.makecpt(
    cmap='abyss',reverse=True,
    series=[0,4000],no_bg=True)
fig.plot(x=srv_lon_units, y=srv_lat_units, 
          style="c0.25c", 
          fill=srv_depth,
          pen="black",
          transparency=10,
          cmap=True,
          no_clip='r')
with pygmt.config(FONT_ANNOT_PRIMARY=22,FONT_LABEL=22):
    fig.colorbar(position="JMB+w8c/0.7c",frame=["af+lDepth to Unit (m)"]) 
fig.savefig(output_dir + 'srv_wells_map.png',dpi=300)





















































