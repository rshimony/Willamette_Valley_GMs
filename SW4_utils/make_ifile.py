#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:09:07 2023

@author: rshimony
"""

import numpy as np
from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from geopy import distance
from scipy.ndimage import gaussian_filter
import pandas as pd

#%%

lons = np.arange(-124.19,-121.51,0.01)
lats = np.arange(43.4,46.7,0.01)

f_path = '/Users/rshimony/Desktop/zzthick4.axyz'

f = np.genfromtxt(f_path)

f_lat = []
f_lon = []
f_z = []

for i in range(len(f)):
    f_lon.append(f[i][0])
    f_lat.append(f[i][1])
    f_z.append(f[i][2])
    
uniq_f_lat = np.unique(f_lat)
uniq_f_lon = np.unique(f_lon)

#%%

f_z_filt = gaussian_filter(f_z, sigma=30)

plt.figure(figsize=(8,10))
plt.scatter(f_lon, f_lat, c=f_z_filt)
plt.colorbar()
plt.show()

#%%

ifile = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_gp1_deep_singmat.txt')

ifile_lons = []
ifile_lats = []
ifile_dists = []

for i in range(1,len(ifile)):
    ifile_lon = ifile[i][0]
    ifile_lat = ifile[i][1]
    ifile_dist = ifile[i][2]
    
    ifile_lons.append(ifile_lon)
    ifile_lats.append(ifile_lat)
    ifile_dists.append(ifile_dist)
    
#%%

ifile_bot = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist_grav_bot.txt')

ifile_dists_bot = []

for i in range(1,len(ifile_bot)):
    ifile_dist_bot = ifile_bot[i][2]
    
    ifile_dists_bot.append(ifile_dist_bot)

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons, ifile_lats, c=ifile_dists)
plt.colorbar()
plt.show()

#%%

max_f_lat = max(uniq_f_lat)
min_f_lat = min(uniq_f_lat)
max_f_lon = max(uniq_f_lon)
min_f_lon = min(uniq_f_lon)

grav_dom = Polygon([[min_f_lon,min_f_lat] , [max_f_lon,min_f_lat] , [max_f_lon,max_f_lat] , [min_f_lon,max_f_lat]])
grav_dom_north = Polygon([[min_f_lon,45.27] , [max_f_lon,45.27] , [max_f_lon,max_f_lat] , [min_f_lon,max_f_lat]])

#%%

ifile_lons_gdom = []
ifile_lats_gdom = []
index = []

for i in range(len(ifile_lons)):
    grd_pnt = Point(ifile_lons[i],ifile_lats[i])
    if grav_dom.contains(grd_pnt) == True:
        ifile_lons_gdom.append(ifile_lons[i])
        ifile_lats_gdom.append(ifile_lats[i])
        index.append(i)
        
ifile_lons_gdom_n = []
ifile_lats_gdom_n = []
index_n = []

for i in range(len(ifile_lons)):
    grd_pnt = Point(ifile_lons[i],ifile_lats[i])
    if grav_dom_north.contains(grd_pnt) == True:
        ifile_lons_gdom_n.append(ifile_lons[i])
        ifile_lats_gdom_n.append(ifile_lats[i])
        index_n.append(i)   
   
#%%

interp = NearestNDInterpolator(list(zip(f_lon, f_lat)), f_z)
interp_depth = interp(ifile_lons_gdom, ifile_lats_gdom)
interp_depth_m = interp_depth*1000

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons_gdom, ifile_lats_gdom, c=interp_depth_m)
plt.colorbar()
plt.show()

f_z_int_filt = gaussian_filter(interp_depth_m, sigma=50)

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons_gdom, ifile_lats_gdom, c=f_z_int_filt)
plt.colorbar()
plt.show() 


#%%

east_patch_max_lat = 45.0
east_patch_min_lat = 44.6
east_patch_max_lon = -122.1
east_patch_min_lon = -122.65

east_patch_dom = Polygon([[east_patch_min_lon,east_patch_min_lat] , [east_patch_max_lon,east_patch_min_lat] , 
                          [east_patch_max_lon,east_patch_max_lat] , [east_patch_min_lon,east_patch_max_lat]])



ifile_lons_east_patch_dom = []
ifile_lats_east_patch_dom = []
index_east_patch_dom = []

for i in range(len(ifile_lons)):
    grd_pnt = Point(ifile_lons[i],ifile_lats[i])
    if east_patch_dom.contains(grd_pnt) == True:
        ifile_lons_east_patch_dom.append(ifile_lons[i])
        ifile_lats_east_patch_dom.append(ifile_lats[i])
        index_east_patch_dom.append(i)
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

south_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_south_poly.csv')
lon_south = south_f['lon']
lat_south = south_f['lat']  
   
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
    
south_poly = np.zeros((len(lon_south),2))
for i in range(len(lon_south)):
    south_poly[i,0] = lon_south[i]
    south_poly[i,1] = lat_south[i]


full_poly = Polygon(outer_poly, [inner1_poly , inner2_poly])

#%%

inpoly_lons = []
inpoly_lats = []
inpoly_index = []

for i in range(len(ifile_lons)):
    grd_point = Point(ifile_lons[i],ifile_lats[i])
    if full_poly.contains(grd_point) == True:
        inpoly_lons.append(ifile_lons[i])
        inpoly_lats.append(ifile_lats[i])
        inpoly_index.append(i)

#%%

inpoly_gdom_lons = []
inpoly_gdom_lats = []
index_gdom_inpoly = []

counter_igp = 1

for i in index:
    if counter_igp % 1000 == 0:
        print(str(counter_igp) + ' from ' + str(len(index)))
    if i in inpoly_index:
        index_gdom_inpoly.append(i)
        inpoly_gdom_lons.append(ifile_lons[i])
        inpoly_gdom_lats.append(ifile_lats[i])
    
    counter_igp = counter_igp + 1
    

inpoly_gdom_n_lons = []
inpoly_gdom_n_lats = []
index_gdom_n_inpoly = []

counter_igp_n = 1

for i in index_n:
    if counter_igp_n % 1000 == 0:
        print('n ' + str(counter_igp_n) + ' from ' + str(len(index_n)))
    if i in inpoly_index:
        index_gdom_n_inpoly.append(i)
        inpoly_gdom_n_lons.append(ifile_lons[i])
        inpoly_gdom_n_lats.append(ifile_lats[i])
    
    counter_igp_n = counter_igp_n + 1
   
#%%

idx_dict = {'inpoly_lons':inpoly_lons, 'inpoly_lats':inpoly_lats , 'inpoly_index':inpoly_index}
idx_df = pd.DataFrame(idx_dict)
idx_df.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_idx.csv',index=False)
   
gdom_dict = {'index_gdom_inpoly':index_gdom_inpoly, 'inpoly_gdom_lats':inpoly_gdom_lats , 'inpoly_gdom_lons':inpoly_gdom_lons}
gdom_df = pd.DataFrame(gdom_dict)
gdom_df.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_gdom_idx.csv',index=False)   

gdom_n_dict = {'index_gdom_n_inpoly':index_gdom_n_inpoly, 'inpoly_gdom_n_lats':inpoly_gdom_n_lats , 'inpoly_gdom_n_lons':inpoly_gdom_n_lons}
gdom_n_df = pd.DataFrame(gdom_dict)
gdom_n_df.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_gdom_idx.csv',index=False)  

#%%

idx_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_idx.csv')
gdom_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_gdom_idx.csv')

inpoly_lons = idx_f['inpoly_lons']
inpoly_lats = idx_f['inpoly_lats']
inpoly_index = idx_f['inpoly_index']

inpoly_gdom_index = gdom_f['index_gdom_inpoly']
# inpoly_gdom_n_index = gdom_f['index_gdom_inpoly']

#%%

ifile_gbase = np.zeros_like(ifile_dists)
ifile_gbase[index] = f_z_int_filt
gbase_depths = ifile_gbase[inpoly_gdom_index]

# gbase_depths_n = ifile_gbase[index_gdom_n_inpoly]

# ifile_dists_grav = (np.array(ifile_dists)**0.5)*18.5 ### this calibration is for the eug_spe version of the basin
ifile_dists_grav = (np.array(ifile_dists))  ### this is for the SRV. (with the deep ifile as input for ifile dists).
ifile_dists_grav[inpoly_gdom_index] = gbase_depths

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons, ifile_lats, c=ifile_dists_grav)
plt.colorbar()
plt.show() 

#%%

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('cal_depth_ifile_gp1_deep_1d.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_dists_grav)):
        
        layer_1 = ifile_dists_grav[i]
        layer_2 = ifile_dists_bot[i]
        
        f.write("%s %s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1,layer_2))

#%%

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('cal_depth_ifile_gp1_deep_singmat.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,1)
    f.write(header_line)
    for i in range(len(ifile_dists_grav)):
        
        layer_1 = ifile_dists_grav[i]
        
        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1))
        
#%%

ifile_original = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_gp1_deep_singmat.txt')

ifile_original_lons = []
ifile_original_lats = []
ifile_original_dists = []

for i in range(1,len(ifile_original)):
    ifile_original_lon = ifile_original[i][0]
    ifile_original_lat = ifile_original[i][1]
    ifile_original_dist = ifile_original[i][2]
    
    ifile_original_lons.append(ifile_original_lon)
    ifile_original_lats.append(ifile_original_lat)
    ifile_original_dists.append(ifile_original_dist)
    
#%%
topo_layer = np.zeros_like(ifile_dists_grav)

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('ifile_gp1_srv_singmat_4rfile.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_dists_grav)):
        
        layer_0 = topo_layer[i]
        layer_1 = ifile_dists_grav[i]
        
        f.write("%s %s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_0,layer_1))
        
#%%

#dists_cal_eug_spe = (np.array(ifile_original_dists)**0.5)*18.5 #### this needs to be done previously

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('ifile_gp1_eug_spe.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_dists_grav)):
        
        layer_1 = ifile_dists_grav[i]
        
        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1))
        
#%%
topo_layer = np.zeros_like(ifile_dists_grav)
ifile_dists_grav_100add = ifile_dists_grav + 150

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('ifile_gp1_eug_spe_4rfile_100add.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_dists_grav)):
        
        layer_0 = topo_layer[i]
        layer_1 = ifile_dists_grav_100add[i]
        
        f.write("%s %s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_0,layer_1))        
    
#%%
nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('ifile_d2e_eug_spe.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,1)
    f.write(header_line)
    for i in range(len(ifile_dists_grav)):
        
        layer_1 = ifile_dists_grav[i]
        
        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1)) 
#%%
     
plt.figure(figsize=(8,10))
plt.scatter(ifile_lons, ifile_lats, c=ifile_dists_grav)
plt.colorbar()      
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/ifile_basindepth.jpg',dpi=300,bbox_inches='tight')
  
plt.figure(figsize=(8,10))
plt.scatter(ifile_lons, ifile_lats, c=ifile_dists_grav_100add)
plt.colorbar()      
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/ifile_100_basindepth.jpg',dpi=300,bbox_inches='tight')   
        
#%%

topo_layer = np.zeros_like(ifile_dists)
basin_layer = topo_layer + 4000
ifile_notaper = np.zeros_like(ifile_dists)

notaper_depths = basin_layer[inpoly_index]
ifile_notaper[inpoly_index] = notaper_depths
        
plt.figure(figsize=(8,10))
plt.scatter(ifile_lons, ifile_lats, c=ifile_notaper)
plt.colorbar()      
plt.show()
        
 #%%

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('ifile_notaper_4rfile.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_notaper)):
        
        layer_0 = topo_layer[i]
        layer_1 = ifile_notaper[i]
        
        f.write("%s %s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_0,layer_1))        
    
nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('ifile_notaper.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,1)
    f.write(header_line)
    for i in range(len(ifile_notaper)):
        
        layer_1 = ifile_notaper[i]
        
        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1))        
        
#%%

idx_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_idx.csv')
gdom_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_gdom_idx.csv')

inpoly_lons = idx_f['inpoly_lons']
inpoly_lats = idx_f['inpoly_lats']
inpoly_index = idx_f['inpoly_index']

inpoly_gdom_index = gdom_f['index_gdom_inpoly']  

cart_ifile_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/basin_depth_poly_cartesian.csv')
depth_poly = cart_ifile_df.depth    
        
       
        
        
        
        
        
        
        
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
    
    
    
    
    
    
    
    
    
    
    