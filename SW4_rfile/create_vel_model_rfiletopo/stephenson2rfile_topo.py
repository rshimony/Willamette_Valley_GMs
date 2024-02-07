#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 16:30:51 2023

@author: rshimony
"""

import matplotlib.pyplot as plt
import sys
import os

flush = sys.stdout.flush()

from pySW4.utils import geo
from pySW4.prep import rfileIO
import numpy as np
import time
import math
from warnings import warn
import xarray as xr
from scipy.interpolate import griddata

wk_dir = '/Users/rshimony/Documents/My_Notebooks/misc/'    
os.chdir(wk_dir)
 
#%%

def llz2utm(lon,lat,projection_zone='None'):
    '''
    Convert lat,lon to UTM
    '''
    from numpy import zeros,where,chararray
    import utm
    from pyproj import Proj
    from scipy.stats import mode
    
    x=zeros(lon.shape)
    y=zeros(lon.shape)
    zone=zeros(lon.shape)
    b=chararray(lon.shape)
    if projection_zone==None:
        #Determine most suitable UTM zone
        for k in range(len(lon)):
            #x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k]-360)
            x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k])
        zone_mode=mode(zone)
        i=where(zone==zone_mode)[0]
        letter=b[i[0]]
        z=str(int(zone[0]))+letter
    else:
        z=projection_zone
    p = Proj(proj='utm',zone=z,ellps='WGS84')
    x,y=p(lon,lat)
    return x,y

def roundup(x):
      return int(math.ceil(x / 100.0)) * 100

#%% Inputs and options
start = time.time()

hhi = np.array([200, 300, 900])
hvi = np.array([100,  300,  900])
block_layers = np.array([0, 1200, 9900, 59400]) 

outpath = '/Users/rshimony/Desktop/WillametteValley/Models/steph_rfile/'
os.chdir(outpath)

rfileoutput = 'stephenson_rfile_topo.rfile'

# rfile input
CVM_L1_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vp_16l1_my_model.nc'    
CVM_L1_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vs_16l1_my_model.nc'
CVM_L2_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vp_16l2_my_model.nc'
CVM_L2_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vs_16l2_my_model.nc'
CVM_L3_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vp_16l3_my_model.nc'
CVM_L3_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vs_16l3_my_model.nc'

# CVM_L1_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vp_16l1.nc'    
# CVM_L1_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vs_16l1.nc'
# CVM_L2_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vp_16l2.nc'
# CVM_L2_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vs_16l2.nc'
# CVM_L3_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vp_16l3.nc'
# CVM_L3_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vs_16l3.nc'

'''
coordinated of SW corner - full Stephenson model

lat = 40.2
lon = -129
'''

filenames = [CVM_L1_P, CVM_L1_S, CVM_L1_P, CVM_L1_S, CVM_L2_P, CVM_L2_S, CVM_L3_P, CVM_L3_S]

block1_vp = xr.open_dataset(CVM_L1_P)
block1_vs = xr.open_dataset(CVM_L1_S)

block2_vp = xr.open_dataset(CVM_L2_P)
block2_vs = xr.open_dataset(CVM_L2_S)

block3_vp = xr.open_dataset(CVM_L3_P)
block3_vs = xr.open_dataset(CVM_L3_S)
    
#%% helper functions
def calcQ(vs3d):

    # print('Setting minimum Qs=60\n')

    vs = vs3d.flatten()
    qs = np.zeros(len(vs))

    for i in range(len(vs)):
        if (vs[i] < 1000):
            qs[i] = (0.23*vs[i]) - 80
        elif (vs[i] >= 1000 and vs[i] <= 2500):
            qs[i] = 150 + 0.1*(vs[i] - 1000)
        else:
            qs[i] = 300.

        if qs[i] < 60:
                qs[i] = 60

    qp = 2*qs

    return qp.reshape(vs3d.shape), qs.reshape(vs3d.shape)


def calcRho(vp3d):

    # From Stpehnson et al., 2017 - Using the empirical relations of Brocher (2005)
    vp = vp3d.flatten()/1000. # Put in km/s
    rho = np.zeros(len(vp))

    for i in range(len(vp)):
        rho[i] = 1.6612*vp[i] - 0.4721*vp[i]**2 + 0.0671*vp[i]**3 - 0.0043*vp[i]**4 + 0.000106*vp[i]**5
        if rho[i] < 2.0:
            rho[i] = 2.0
        elif rho[i] > 3.5:
            rho[i] = 3.5

    rho = rho * 1000 # Back in kg/m3

    return rho.reshape(vp3d.shape)
    
#%% Writing file header
rfilename = os.path.join(outpath, rfileoutput)

proj4 = ('+proj=tmerc +datum=NAD83 +units=m +lon_0=-124.19 +lat_0=43.4 +no_defs') ##

magic        = np.int32(1)
precision    = np.int32(4)
attenuation  = np.int32(1)
az           = np.float64(0.0)
lon0         = np.float64(-124.19) 
lat0         = np.float64(43.4) 
nb           = np.int32(len(block_layers)) 

header = (magic, precision, attenuation, az, lon0, lat0, proj4.encode(), nb)

print('\nrfile initialized at {}.\n'
      'Writing file header...\n'
      '======================'.format(rfilename))
flush
with open(rfilename, 'wb') as f:
    rfileIO.write_hdr(f, *header)
    
print('Done!')


#%% Writing block headers
print('Writing block headers')

# Writing block headers
for i in range(nb):
    
    if i == 0:

        # Block No.0 (topo)
        hh0 = np.float64(hhi[i])    # grid size in hori direction
        hv0 = np.float64(hvi[i])    # grid size in vertical direction
        z00 = np.float64(0)         # base z-level for block b.
        nc0 = np.int32(1)           # no of components in block b (in this case 1, just elevation)
        ni0 = np.int32(len(block1_vp.utmn))   # no of grid points along y axis
        nj0 = np.int32(len(block1_vp.utme))   # no of grid points along x axis
        nk0 = np.int32(1)           # no of grid points along the z axis (in this case 1)

        x_extent = (ni0 - 1) * hh0
        y_extent = (nj0 - 1) * hh0

        print ('\nrfile horizontal dimentions are '
               '{:.2f} km x {:.2f} km\n'
               u'with origin at {:.4f}°N, {:.4f}°E\n'
               'ni0: {}, nj0: {}'.format(x_extent * 1e-3, y_extent * 1e-3,
                                         lat0, lon0,
                                         ni0, nj0))

        block0_hdr = [hh0, hv0, z00, nc0, ni0, nj0, nk0]
        print(block0_hdr)

    else:
        # the data block
        if i == 1:
            zmin_blk = block_layers[i-1]
        else:
            zmin_blk = globals()['z'+str(i-1)].max() 
            
        zmax_blk = block_layers[i]

        globals()['z'+str(i)] = np.arange(zmin_blk, zmax_blk + hvi[i-1], hvi[i-1])

        hh_blk = np.float64(hhi[i-1])
        hv_blk = np.float64(hvi[i-1])
        z0_blk = np.float64(zmin_blk)
        nc_blk = np.int32(5)
        ni_blk = np.int32(len(globals()['block'+str(i)+'_vp'].utmn))
        nj_blk = np.int32(len(globals()['block'+str(i)+'_vp'].utme))
        nk_blk = np.int32(globals()['z'+str(i)].size)

        x_extent = (ni_blk - 1) * hh_blk
        y_extent = (nj_blk - 1) * hh_blk

        print ('\nBlock No.{} horizontal dimentions are '
                '{:.2f} km x {:.2f} km\n'
                'with a vertical extent of {} m high to {} m low.\n'
                'ni1: {}, nj1: {}, nk1: {}'.format(i, x_extent * 1e-3, y_extent * 1e-3,
                                          z0_blk, globals()['z'+str(i)].max(),
                                          ni_blk, nj_blk, nk_blk))

        globals()['block'+str(i)+'_hdr'] = [hh_blk, hv_blk, z0_blk, \
                      nc_blk, ni_blk, nj_blk, nk_blk]
        print(globals()['block'+str(i)+'_hdr'])

# Write block headers
print('\nWriting block headers...\n'
        '========================')
flush

with open(rfilename, 'ab') as f:
    for i in range(nb):
        rfileIO.write_block_hdr(f, *globals()['block'+str(i)+'_hdr'])

print('Done!')


#%% Write topo block data section
print('\nWriting topography data...\n'
        '==========================')
flush


# download topo from NOAA
# convert grd to dat
# interpolate dat at required located to create a 2D topo 
# multiply it by -1 since topo is positive downwards



# import subprocess
# f = open("stephenson_topo2.csh", "w")
# f.write("#!/bin/bash\n")
# f.write("\n")
# f.write("# To download Topography & bathymetry data (from OpenSWPC)\n")
# f.write("#if [ ! -e ETOPO1_Bed_g_gmt4.grd.gz ]; then\n")
# f.write(" #   curl -o ETOPO1_Bed_g_gmt4.grd.gz https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd.gz\n")
# f.write("#fi\n")
# f.write("#gunzip ETOPO1_Bed_g_gmt4.grd.gz\n")
# f.write("region=-124.19/-121.51/43.4/46.1\n")
# f.write("\n")
# f.write("# Topography grd file\n")
# f.write("topo=./ETOPO1_Bed_g_gmt4.grd\n")
# f.write("gmt grdcut -R$region $topo -Gtopo.cut_csz.grd\n")
# f.write("\n")
# f.write("gmt grd2xyz topo.cut_csz.grd > topo.cut_csz.dat\n")
# f.write("\n")
# f.close()
# subprocess.run('conda activate base', shell=True)
# subprocess.call(['sh', './stephenson_topo.csh'])
# subprocess.run('conda activate branch-env', shell=True)


f_path = 'topo.cut_csz.dat'

f = np.genfromtxt(f_path)

f_lat = f[:,1]
f_lon = f[:,0]
f_z = f[:,2]

topo_utme , topo_utmn = llz2utm(f_lon,f_lat,projection_zone=10)

xi_blk1, yi_blk1 = np.meshgrid(block1_vp.utme.values, block1_vp.utmn.values)
topo_grid = griddata((topo_utme, topo_utmn), f_z, (xi_blk1, yi_blk1), method='linear')

min_topo = roundup(-topo_grid.min())
max_topo = roundup(-topo_grid.max())

print('maximum elevation = ',topo_grid.max(), 'and minimum elev = ',topo_grid.min())

assert min_topo < 1200, 'The minimum elevation must be less than the base of the first block. (Atleast for now)'

if max_topo < block_layers[0]:
    print('Changing the highest topography to', str(max_topo),'m')
    block_layers[0] = max_topo
    
    
# plt.figure()
# plt.scatter(topo_utme,topo_utmn)
# plt.scatter(xi_blk1,yi_blk1)
# plt.show()

# plt.figure()
# plt.scatter(f_lon,f_lat,c=f_z)
# plt.show()
# plt.figure(figsize=(8,10))
# plt.scatter(xi_blk1,yi_blk1,c=vp_prev_base)
# plt.show()

with open(rfilename, 'ab') as f:
    (-topo_grid.astype(np.float32)).tofile(f)
print('Done!')


# #%% write block No.1 data section
# for blk in range(nb-1):
    
#     #blk = 1
#     vp_blk = globals()['block'+str(blk+1)+'_vp']
#     vs_blk = globals()['block'+str(blk+1)+'_vs']
    
#     z_blk = globals()['z'+str(blk+1)]
#     ni_blk = np.int32(len(vp_blk.utmn))
#     nj_blk = np.int32(len(vp_blk.utme))
#     hhi_blk = hhi[blk]
    
#     # init the k_array
#     k_array = np.ones((z_blk.size, 5))
    
#     print('\nWriting data section of Block No.'+str(blk+1)+'\n'
#             '==================================')
#     print('z_blk =')
#     print(z_blk)
    
    
#     # Add repeat layers for the second and third block
#     if blk == 1 or blk == 2:
#         # Open earlier files
#         vp_blk_prev = globals()['block'+str(blk)+'_vp']
#         vs_blk_prev = globals()['block'+str(blk)+'_vs']

#         vp_prev = vp_blk_prev.vp.values
#         vs_prev = vs_blk_prev.vs.values

#         # Interpolate last layer of previous block onto current spacing
#         xi_new, yi_new = np.meshgrid(vp_blk.utme.values, vp_blk.utmn.values)
#         xi_prev, yi_prev = np.meshgrid(vp_blk_prev.utme.values, vp_blk_prev.utmn.values)
        
#         vp_prev_base = griddata((xi_prev.T.flatten(), yi_prev.T.flatten()), vp_prev[:,:,-1].flatten(), (xi_new, yi_new), method='linear')
#         vs_prev_base = griddata((xi_prev.T.flatten(), yi_prev.T.flatten()), vs_prev[:,:,-1].flatten(), (xi_new, yi_new), method='linear')

#         if np.isnan(vs_prev_base[0,0]):
#             vp_prev_base[0,:] = vp_prev_base[1,:]
#         if np.isnan(vs_prev_base[0,0]):
#             vs_prev_base[0,:] = vs_prev_base[1,:]
        
#         # Add additional layer at the top
#         vp_corr = np.dstack((vp_prev_base.T, vp_blk.vp.values))
#         vs_corr = np.dstack((vs_prev_base.T, vs_blk.vs.values))

#         vp_blk_prev.close()
#         vs_blk_prev.close()
    
    
#     with open(rfilename, 'ab') as f:
#         for i in range(ni_blk):
#             for j in range(nj_blk):
                
#                 # Determing k_array at location index [i,j]  #rho,vp,vs,qp,qs
#                 if blk == 0:
#                     vp = vp_blk.vp.values[j,i,:]
#                     vs = vs_blk.vs.values[j,i,:]
#                 else:
#                     vp = vp_corr[j,i,:]
#                     vs = vs_corr[j,i,:]
                    
#                 rho = calcRho(vp)
#                 qp, qs = calcQ(vs)
                
#                 # remove layers with 1e20 values
#                 index = np.where(vp>1e7)
#                 vp[index] = -999
#                 vs[index] = -999
#                 rho[index] = -999
#                 qp[index] = -999
#                 qs[index] = -999
                
                
#                 k_array[:,0]= rho
#                 k_array[:,1]= vp
#                 k_array[:,2]= vs
#                 k_array[:,3]= qp
#                 k_array[:,4]= qs
                
#                 if k_array.min() < 50 and k_array.min() > -999:
#                     print(i,j,k_array.min())
            
                
                
            
            
#                 # crate anothe karray with -999 of the same length (10 km to 1.2 km)
#                 # replicate the first element of the karray at position i,j to fill the karray from say -10km up to 1.2 km
#                 # What is the topo at position i,j
#                 # how many elements between the topo and the base of the block 1 (1.2 km)
#                 # extract these elemtns from the extended karray (i,j) from the end or below
#                 # insert it to the -999 array from behind
            
            
            
#                 if blk == 0:
#                     # crate anothe karray with -999 of the same length (10 km to 1.2 km)
#                     len_karray = int(1+(1200-max_topo)/hvi[0])
#                     tmp_karray = -999*np.ones((len_karray, 5))
                
#                     # replicate the first element of the karray at position i,j to fill the karray from say -10km up to 1.2 km
#                     repetitions = len_karray- len(k_array)
#                     k_array_n = np.concatenate((np.tile(k_array[0,:], (repetitions, 1)), k_array), axis=0)
                
            
            
            
            
            
            
#                 # write k_array to file
#                 k_array.astype(np.float32).tofile(f)
#                 if k_array.size != (len(z_blk)*5):
#                     print(k_array.size)
                    
                
                    
#                 if j == 0 and i%100 == 0:
#                     print('\rWriting row {} of {}, '
#                           'column {} of {} of block 1...'.format(i + 1, ni_blk,
#                                                                   j + 1, nj_blk)),
#                 flush
    
#     print('\nDone!\n')

# # ####################################################################
# end = time.time()
# time_elaps = end - start
# if time_elaps < 60:
#     print(f'Duration: {round(time_elaps)} seconds')
# else:
#     print(f'Duration: {round(time_elaps/60)} minutes')




# #%% Reading and plotting the rfile
# model = rfileIO.read(rfileoutput, 'all', verbose=True)

# #%% ---------------------------------------------------------------------------
# model.plot_topography(cmap='terrain')
# figpath = os.getcwd() +'/fig.rfiletopo.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 
# plt.show()

# #%%

# cs = model.get_cross_section(0, 370, 100, 100) # lat1,lat2,lon1,lon2 in km
# fig, ax, cb = cs.plot(property='vp',aspect=10, cmap='jet') #vmin=0, 
# figpath = os.getcwd() +'/fig.crossection2.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 
# plt.show()


# #%%
# # ---------------------------------------------------------------------------
# # Vs profile at different locations (x,y) within the domain
# lat_bin = np.linspace(100.0, 300, 37)
# lon_bin = np.linspace(100.0, 200, 21)

# fig, ax = plt.subplots()

# for i in range(len(lat_bin)):
#     for j in range(len(lon_bin)):
#         z, properties = model.get_z_profile(lon_bin[j], lat_bin[i])

#         p = np.ma.masked_equal(properties.T[2], -999)
        
    
#         ax.plot(p, z,linewidth=0.5) 

# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.set_xlabel('Vs (m/s)', fontsize=20)
# ax.xaxis.set_label_position('top')
# ax.set_ylabel('Z, km', fontsize=20)
# ax.tick_params(axis='x',labelsize=15)
# ax.tick_params(axis='y',labelsize=15)

# figpath = os.getcwd() +'/fig.profile_vs.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 

# #%%
# cs = model.get_cross_section(100, 100, 0, 215) # lat1,lat2,lon1,lon2 in km
# fig, ax, cb = cs.plot(property='vs',aspect=1, cmap='jet') #vmin=0, 
# figpath = os.getcwd() +'/fig.crossection2.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 
# plt.show()

# #%%
# cs = model.get_cross_section(0, 375, 100, 100) # lat1,lat2,lon1,lon2 in km
# fig, ax, cb = cs.plot(property='vs',aspect=2, cmap='jet') #vmin=0, 
# figpath = os.getcwd() +'/fig.crossection2.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 
# plt.show()
















