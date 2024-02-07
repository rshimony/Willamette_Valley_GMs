#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:41:16 2023

@author: rshimony
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pygmt
import pandas as pd
#%%

vs30_file = '/Users/rshimony/Documents/Vs30_global/earthquake-global_vs30-master/PNW/pnw.grd'
water_file = '/Users/rshimony/Documents/Vs30_global/earthquake-global_vs30-master/PNW/landmask_water.grd'
land_file = '/Users/rshimony/Documents/Vs30_global/earthquake-global_vs30-master/PNW/landmask_land.grd'

vs30_grd = Dataset(vs30_file, 'r', format='NETCDF4')
water_grd = Dataset(water_file, 'r', format='NETCDF4')
land_grd = Dataset(land_file, 'r', format='NETCDF4')


lon_grd=vs30_grd.variables['lon'][:]
lat_grd=vs30_grd.variables['lat'][:]
vs30=vs30_grd.variables['z'][:,:]
water=water_grd.variables['z'][:]
land=land_grd.variables['z'][:]

LON,LAT = np.meshgrid(lon_grd,lat_grd)

lon_val_idx = np.argwhere((lon_grd>-124.19) & (lon_grd< -121.51))[:,0]
lat_val_idx = np.argwhere((lat_grd>43.4) & (lat_grd< 46.7))[:,0]

lon_grd_val = lon_grd[lon_val_idx]
lat_grd_val = lat_grd[lat_val_idx]

# LON_val,LAT_val = np.meshgrid(lon_grd_val,lat_grd_val)
# vs30_val = vs30[lon_val_idx,lat_val_idx]

# flat_lon = np.array(LON.flatten())
# flat_lat = np.array(LAT.flatten())
# flat_vs30 = np.array(vs30.flatten())


# plt.figure(figsize=(8,10))
# plt.scatter(flat_lon,flat_lat,c=flat_vs30)
# plt.colorbar()
# plt.show()

#lons = np.arange(-124.19,-121.51,0.01)
# lats = np.arange(43.4,46.7,0.01)

#%%
st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/domain_stations.csv')
st_names = np.array(st_file['st_name_dom'])
st_lons = np.array(st_file['st_lon_dom'])
st_lats = np.array(st_file['st_lat_dom']) 
st_names_unq = pd.unique(st_names)
st_lons_unq = pd.unique(st_lons)
st_lats_unq = pd.unique(st_lats)

# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

# make color pallets
# pygmt.makecpt(
#     cmap='geo',
#     series='-2000/4000/100',
#     continuous=True)

## topo:
topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=18):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=['WSne', "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(vs30_file, cmap='jet',shading=False,frame=True)
    fig.colorbar(frame="af+lVs30 (m/s)")
    fig.coast(water="skyblue3")
    
    # fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    # fig.plot(x=xph1, y=yph1,pen="1p,black")
    # fig.plot(x=xph2, y=yph2,pen="1p,black")
    
    # fig.plot(x=valst_lons, y=valst_lats, 
    #           style="t0.42c", 
    #           fill="navy",
    #           pen="black",
    #           transparency=30,
    #           label = 'Stations') 

    fig.plot(x=st_lons_unq, y=st_lats_unq, 
              style="t0.4c", 
              fill="white",
              pen="black",
              transparency=50,
              )

    # fig.text(text=st_names[20], x=st_lons[20]-0.02, y=st_lats[20]+0.08 , font='10p,Helvetica-Bold,magenta4')
    # fig.text(text=st_names[64], x=st_lons[64]-0.06, y=st_lats[64]-0.06, font='10p,Helvetica-Bold,magenta4')
    # fig.text(text=st_names[89], x=st_lons[89]+0.05, y=st_lats[89]+0.06, font='10p,Helvetica-Bold,magenta4')
    # fig.legend(position="JTR+jTR+o0.2c", box='+gWhite+p1p',)

fig.show()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/vs30_stations.png',dpi=600)