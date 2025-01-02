#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:12:40 2024

@author: rshimony
"""

'''
Makes a Z2.5 and Z1.0 surface data file from data contained in an rfile (SW4 simulations data file).
Script outline:
    1) Creating WV polygon
    2) Reading rfiles (USGS CVM and WV model (SRV 1D))
    3) Setting extents and grid sizes for sampling rfiles and output z_surface files
    4) Writing z_surface files from both rfiles and saving to text file.
    5) Reading in the newly made z_surface files
    6) Making a "Valley Only" depth z_surface files. 
    * These files have z depth values only for grid point within the WV polygon, while the rest of the domain (outside the WV) is 0.
    7) Plotting maps of the new z_surfaces (8 maps: z1.0 and z2.5 for cvm and wv model, and for full and valley z_surface files.)
'''

from pySW4.prep import rfileIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

#%%
## Crating WV polygon
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

## rfile input and read
rfile_path_cvm = '/Users/rshimony/Desktop/WillametteValley/wv_project/rfiles/usgs_cvm/stephenson_rfile2.rfile'
rfile_cvm = rfileIO.read(rfile_path_cvm, 'all', verbose=True)

rfile_path_wv_model = '/Users/rshimony/Desktop/WillametteValley/wv_project/rfiles/wv_basin_model_1d/srv/rfile_srv_1d.rfile'
rfile_wv_model = rfileIO.read(rfile_path_wv_model, 'all', verbose=True)

#%%

## Set extent and grid sampling size
# Cartesian grid - for sampling the rfile
e = 216 # E-W length [km]
n = 378 # N-S length [km]
h_diff = 0.2 # [km]

# Geographic grid - for the output z_surface files
lon_min, lon_max = [-124.19,-121.38]
lat_min, lat_max = [43.4,46.8]
lon_num, lat_num = [1080,1890]

e_arr = np.arange(0,e,h_diff)
n_arr = np.arange(0,n,h_diff)

lon = np.linspace(lon_min,lon_max,lon_num)
lat = np.linspace(lat_min,lat_max,lat_num)

#%%

# Initiate the z_surface file
def write_z_sur_from_rfile(outfile_path, rfile):
    z_surface_path = outfile_path
    z_surface_file=open(z_surface_path,'w+')
    header_line="lon lat x y z1.0 z2.5\n"
    z_surface_file.write(header_line)
    
    # Writing z_surface file
    counter = 1
    for i in range(len(e_arr)):
        if counter % 10 == 0:
            print(str(counter) + ' from ' + str(len(e_arr)))
        for j in range(len(n_arr)):
            
            vel_z_prof = rfile.get_z_profile(n_arr[j],e_arr[i], property='vs')[1]
            depth_z_prof = rfile.get_z_profile(n_arr[j],e_arr[i], property='vs')[0]
    
            z_1_idx = np.argwhere(vel_z_prof>=1000)[0][0]
            z_25_idx = np.argwhere(vel_z_prof>=2500)[0][0]
    
            z_1 = depth_z_prof[z_1_idx]
            z_25 = depth_z_prof[z_25_idx]
    
            line="%s %s %s %s %s %s\n"%(lon[i],lat[j],e_arr[i],n_arr[j],z_1,z_25)
            z_surface_file.write(line)
        counter = counter + 1
    
    z_surface_file.close()

write_z_sur_from_rfile('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_cvm', rfile_cvm)
write_z_sur_from_rfile('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_wv_model', rfile_wv_model)

#%%

## Reading z_surface file
def reading_z_surface(z_surface_file_path):
    z_sur = np.genfromtxt(z_surface_file_path)
    
    z_sur_lon_ls = []
    z_sur_lat_ls = []
    z_sur_x_ls = []
    z_sur_y_ls = []
    z_sur_z1_ls = []
    z_sur_z25_ls = []
    
    for i in range(1,len(z_sur)):
        z_sur_lon = z_sur[i][0]
        z_sur_lat = z_sur[i][1]
        z_sur_x = z_sur[i][2]
        z_sur_y = z_sur[i][3]
        z_sur_z1 = z_sur[i][4]
        z_sur_z25 = z_sur[i][5]
        
        z_sur_lon_ls.append(z_sur_lon)
        z_sur_lat_ls.append(z_sur_lat)
        z_sur_x_ls.append(z_sur_x)
        z_sur_y_ls.append(z_sur_y)
        z_sur_z1_ls.append(z_sur_z1)
        z_sur_z25_ls.append(z_sur_z25)
    
    return z_sur_lon_ls,z_sur_lat_ls,z_sur_x_ls,z_sur_y_ls,z_sur_z1_ls,z_sur_z25_ls

z_sur_cvm = reading_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_cvm')
z_sur_wv_model = reading_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_wv_model')

#%%

z_sur_lons = z_sur_wv_model[0]
z_sur_lats = z_sur_wv_model[1]
z_sur_x = z_sur_wv_model[2]
z_sur_y = z_sur_wv_model[3]

def make_valley_z_sur(full_z_sur_array):
    '''
    Making a new z_surface list where the depth data is only for grid point within the WV polygon.
    The rest of the domain (outside the WV) is 0.
    '''
    counter = 1

    valley_z_sur = []
    
    for i in range(len(full_z_sur_array)):
        if counter % 100 == 0:
            print(str(counter) + ' from ' + str(len(full_z_sur_array)))
            
        grd_point = Point(z_sur_lons[i],z_sur_lats[i])
    
        if full_poly.contains(grd_point) == False:
            full_z_sur_array[i] = 0
            
        valley_z_sur.append(full_z_sur_array[i])
                
        counter = counter + 1
    
    return valley_z_sur

z1_cvm = z_sur_cvm[4]
z2pt5_cvm = z_sur_cvm[5]
z1_wv_model = z_sur_wv_model[4]
z2pt5_wv_model = z_sur_wv_model[5]

valley_z1_cvm = make_valley_z_sur(z1_cvm)
valley_z2pt5_cvm = make_valley_z_sur(z2pt5_cvm)
valley_z1_wv_model = make_valley_z_sur(z1_wv_model)
valley_z2pt5_wv_model = make_valley_z_sur(z2pt5_wv_model)

np.savetxt('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_valley_cvm.xyz', 
           np.array([np.array(z_sur_lons), np.array(z_sur_lats), np.array(z_sur_x), np.array(z_sur_y), np.array(valley_z1_cvm),np.array(valley_z2pt5_cvm)]).T,
           delimiter='\t', header="lon,lat,x,y,z1.0,z2.5", fmt="%s")

np.savetxt('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_valley_wv_model.xyz', 
           np.array([np.array(z_sur_lons), np.array(z_sur_lats), np.array(z_sur_x), np.array(z_sur_y), np.array(valley_z1_wv_model),np.array(valley_z2pt5_wv_model)]).T,
           delimiter='\t', header="lon,lat,x,y,z1.0,z2.5", fmt="%s")
#%%

## Reading and plotting z_surface

def plot_z_surface(surface_file, surface_type, outdir, fig_name):
    z_sur = np.genfromtxt(surface_file)
    
    z_sur_lon_ls = []
    z_sur_lat_ls = []
    z_sur_x_ls = []
    z_sur_y_ls = []
    z_sur_z1_ls = []
    z_sur_z25_ls = []
    
    for i in range(1,len(z_sur)):
        z_sur_lon = z_sur[i][0]
        z_sur_lat = z_sur[i][1]
        z_sur_x = z_sur[i][2]
        z_sur_y = z_sur[i][3]
        z_sur_z1 = z_sur[i][4]
        z_sur_z25 = z_sur[i][5]
        
        z_sur_lon_ls.append(z_sur_lon)
        z_sur_lat_ls.append(z_sur_lat)
        z_sur_x_ls.append(z_sur_x)
        z_sur_y_ls.append(z_sur_y)
        z_sur_z1_ls.append(z_sur_z1)
        z_sur_z25_ls.append(z_sur_z25)
    
    if surface_type == 'z1pt0':
        plt.figure(figsize=(8,10))
        plt.scatter(z_sur_lon_ls,z_sur_lat_ls,c=z_sur_z1_ls)
        plt.colorbar()
        plt.savefig(outdir + fig_name ,dpi=300,bbox_inches='tight')
        
    if surface_type == 'z2pt5':
        plt.figure(figsize=(8,10))
        plt.scatter(z_sur_lon_ls,z_sur_lat_ls,c=z_sur_z25_ls)
        plt.colorbar()
        plt.savefig(outdir + fig_name ,dpi=300,bbox_inches='tight')
    
out_fig_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/'

plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_cvm' , 'z1pt0' , out_fig_dir , 'z1.0_cvm_full.png')
plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_cvm' , 'z2pt5' , out_fig_dir , 'z2.5_cvm_full.png')
plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_wv_model' , 'z1pt0' , out_fig_dir , 'z1.0_wv_model_full.png')
plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_wv_model' , 'z2pt5' , out_fig_dir , 'z2.5_wv_model_full.png')

plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_valley_cvm.xyz' , 'z1pt0' , out_fig_dir , 'z1.0_cvm_valley.png')
plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_valley_cvm.xyz' , 'z2pt5' , out_fig_dir , 'z2.5_cvm_valley.png')
plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_valley_wv_model.xyz' , 'z1pt0' , out_fig_dir , 'z1.0_wv_model_valley.png')
plot_z_surface('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_valley_wv_model.xyz' , 'z2pt5' , out_fig_dir , 'z2.5_wv_model_valley.png')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



































