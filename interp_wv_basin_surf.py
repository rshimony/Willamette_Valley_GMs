#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 11:58:56 2023

@author: rshimony
"""

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

#%%

def extract_ifile_basin_data(ifile_path):
    ifile_dists_file = np.genfromtxt(ifile_path)
    
    ifile_dists_lons = []
    ifile_dists_lats = []
    ifile_dists = []
    
    for i in range(1,len(ifile_dists_file)):
    
        ifile_dist = ifile_dists_file[i][2]
        ifile_dist_lons = ifile_dists_file[i][0]
        ifile_dist_lats = ifile_dists_file[i][1]
        
        ifile_dists.append(ifile_dist)
        ifile_dists_lons.append(ifile_dist_lons)
        ifile_dists_lats.append(ifile_dist_lats)
    
    ifile_dists_arr = np.array(ifile_dists)
    ifile_dists_basin = ifile_dists_arr[ifile_dists_arr>0]
    
    ifile_dists_lons_arr = np.array(ifile_dists_lons)
    ifile_dists_lons_basin = ifile_dists_lons_arr[ifile_dists_arr>0]
    
    ifile_dists_lats_arr = np.array(ifile_dists_lats)
    ifile_dists_lats_basin = ifile_dists_lats_arr[ifile_dists_arr>0]
    
    return ifile_dists_lons_basin , ifile_dists_lats_basin , ifile_dists_basin

# ifile_grav_top = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist_grav_top.txt'
ifile_grav_top = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_gp1_deep_singmat.txt'
ifile_grav_bot = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist_grav_bot.txt'
ifile_gp2 = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_gp2_singmat.txt'
ifile_dist2edge = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist.txt'

grav_top_lon , grav_top_lat , grav_top_depth = extract_ifile_basin_data(ifile_grav_top)
grav_bot_lon , grav_bot_lat , grav_bot_depth = extract_ifile_basin_data(ifile_grav_bot)
gp2_lon , gp2_lat , gp2_depth = extract_ifile_basin_data(ifile_gp2)
dist2edge_lon , dist2edge_lat , dist2edge_depth = extract_ifile_basin_data(ifile_dist2edge)

#%%

inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_outer.csv')
lon_outer = outer_f['outer_lon']
lat_outer = outer_f['outer_lat']   
  
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

xpn, ypn = full_poly.exterior.xy
xph1, yph1 = full_poly.interiors[0].xy
xph2, yph2 = full_poly.interiors[1].xy

east_patch_max_lat = 45.0
east_patch_min_lat = 44.6
east_patch_max_lon = -122.1
east_patch_min_lon = -122.65
east_patch_dom = Polygon([[east_patch_min_lon,east_patch_min_lat] , [east_patch_max_lon,east_patch_min_lat] , 
                          [east_patch_max_lon,east_patch_max_lat] , [east_patch_min_lon,east_patch_max_lat]])

xpn_c = []
ypn_c = []

for i in range(len(xpn)):
    grd_pnt = Point(xpn[i],ypn[i])
    if east_patch_dom.contains(grd_pnt) == False:
        xpn_c.append(xpn[i])
        ypn_c.append(ypn[i])

#%%
    
ifile_dists_lons = []
ifile_dists_lats = []
ifile_dists = []

ifile_dists_file = np.genfromtxt(ifile_dist2edge)

for i in range(1,len(ifile_dists_file)):

    ifile_dist = ifile_dists_file[i][2]
    ifile_dist_lons = ifile_dists_file[i][0]
    ifile_dist_lats = ifile_dists_file[i][1]
    
    ifile_dists.append(ifile_dist)
    ifile_dists_lons.append(ifile_dist_lons)
    ifile_dists_lats.append(ifile_dist_lats)
    
#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_dists_lons, ifile_dists_lats, c=ifile_dists)
plt.colorbar()
plt.scatter(xpn_c, ypn_c, c='red' , s=20)
plt.scatter(xph1, yph1, c='red' , s=20)
plt.scatter(xph2, yph2, c='red' , s=20)
plt.show() 

#%%

depth2units_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/depth2units_wells/depth2units_wells.csv')

lon_units = depth2units_f['Longitude']
lat_units = depth2units_f['Latitude']
srv_depth = depth2units_f['srv_depth']
yamhill_depth = depth2units_f['yamhill_depth']
spencer_depth = depth2units_f['Spencer_depth']
eugene_depth = depth2units_f['eugene_depth']

nonan_depth = spencer_depth[spencer_depth.notna()]
nonan_lons = lon_units[spencer_depth.notna()]
nonan_lats = lat_units[spencer_depth.notna()]

#%%
# utme_low = X_t_flat[(vp_vs < min_vp_vs) & (vp_flat < 1e4)]

pn_c_zeros = np.zeros_like(xpn_c)
ph1_zeros = np.zeros_like(xph1)
ph2_zeros = np.zeros_like(xph2)

pn_c_zeros[(np.array(xpn_c)>-122.9) & (np.array(xpn_c)<=-122.8)] = 1500
pn_c_zeros[(np.array(xpn_c)>-122.8) & (np.array(xpn_c)<=-122.7)] = 2500
pn_c_zeros[(np.array(xpn_c)>-122.7) & (np.array(xpn_c)<=-122.6)] = 3200
pn_c_zeros[np.array(xpn_c) > -122.6] = 4000

poly_conc_lons = np.concatenate((xpn_c, xph1 , xph2 , lon_units), axis=None)
poly_conc_lats = np.concatenate((ypn_c, yph1 , yph2 , lat_units), axis=None)
poly_conc_depths = np.concatenate((pn_c_zeros, ph1_zeros , ph2_zeros , yamhill_depth), axis=None)

plt.figure(figsize=(8,10))
plt.scatter(poly_conc_lons, poly_conc_lats, c=poly_conc_depths)
plt.colorbar()
plt.show() 

wv_poly_interp = griddata((poly_conc_lons, poly_conc_lats), poly_conc_depths, (dist2edge_lon[0::200], dist2edge_lat[0::200]), method='linear')

plt.figure(figsize=(8,10))
plt.scatter(dist2edge_lon[0::200], dist2edge_lat[0::200], c=wv_poly_interp)
plt.colorbar()
plt.show() 

#%%

ifile_lons_arr = np.array(ifile_dists_lons)
ifile_lats_arr = np.array(ifile_dists_lats)
ifile_dists_arr = np.array(ifile_dists)

ifile_lons_sliced = ifile_lons_arr[0::500]
ifile_lats_sliced = ifile_lats_arr[0::500]
ifile_dists_sliced = ifile_dists_arr[0::500]

ifile_dists_sliced_cal = (ifile_dists_sliced**0.5)*17

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons_sliced, ifile_lats_sliced, c=ifile_dists_sliced_cal)
plt.colorbar()
plt.show() 

#%%

# conc_lons = np.concatenate((ifile_lons_sliced, nonan_lons), axis=None)
# conc_lats = np.concatenate((ifile_lats_sliced, nonan_lats), axis=None)
# conc_depths = np.concatenate((ifile_dists_sliced_cal, nonan_depth), axis=None)

# plt.figure(figsize=(8,10))
# plt.scatter(conc_lons, conc_lats, c=conc_depths)
# plt.colorbar()
# plt.show() 

conc_lons = np.concatenate((ifile_lons_sliced, dist2edge_lon[0::200]), axis=None)
conc_lats = np.concatenate((ifile_lats_sliced, dist2edge_lat[0::200]), axis=None)
conc_depths = np.concatenate((ifile_dists_sliced_cal, wv_poly_interp), axis=None)

plt.figure(figsize=(8,10))
plt.scatter(conc_lons, conc_lats, c=conc_depths)
plt.colorbar()
plt.show() 

#%%

wv_depth_interp = griddata((conc_lons, conc_lats), conc_depths, (ifile_dists_lons, ifile_dists_lats), method='linear')

nonan_depth_interp = wv_depth_interp[pd.notna(wv_depth_interp)]
nonan_depth_idx = np.argwhere(pd.notna(wv_depth_interp))
nonan_depth_idx_flat = nonan_depth_idx.flatten()

wv_interp_surf = np.array(ifile_dists)
wv_interp_surf[nonan_depth_idx_flat] = nonan_depth_interp


# plt.figure(figsize=(8,10))
# plt.scatter(ifile_dists_lons, ifile_dists_lats, c=ifile_dists)
# plt.colorbar()
# plt.scatter(nonan_lons, nonan_lats, c='red' , s=25)
# plt.show() 

plt.figure(figsize=(8,10))
plt.scatter(ifile_dists_lons, ifile_dists_lats, c=wv_depth_interp)
plt.colorbar()
plt.show() 

plt.figure(figsize=(8,10))
plt.scatter(ifile_dists_lons, ifile_dists_lats, c=wv_interp_surf)
plt.colorbar()
plt.show() 



































