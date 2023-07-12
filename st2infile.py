#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 14:23:17 2022

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry.polygon import Polygon
import cartopy.crs as ccrs
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

xpnn = np.array([xpn])
xph1n = np.array([xph1])
xph2n = np.array([xph2])

ypnn = np.array([ypn])
yph1n = np.array([yph1])
yph2n = np.array([yph2])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

#%%

stations = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/salem_eq/metadata/snr_stations.csv')
st_lat = stations['st_lat']
st_lon = stations['st_lon']
st_name = stations['st_nm']

st_lat_dom = []
st_lon_dom = []
st_name_dom = []

for i in range(len(st_lat)):
    if 43.5 < st_lat[i] < 45.9 and -123.9 < st_lon[i] < -122.1:
        st_lat_dom.append(st_lat[i])
        st_lon_dom.append(st_lon[i])
        st_name_dom.append(st_name[i])
        
st_dom_dict = {'st_lon_dom':st_lon_dom,'st_lat_dom':st_lat_dom, 'st_name_dom':st_name_dom,}
st_dom_df = pd.DataFrame(st_dom_dict)
st_dom_df.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/domain_stations.csv',index=False)

#%%

from cartopy.io.img_tiles import GoogleTiles
class ShadedReliefESRI(GoogleTiles):
    # shaded relief
    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
               'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
               z=z, y=y, x=x)
        return url

#%%
plt.figure(figsize=(8, 12))

ax = plt.axes(projection=ShadedReliefESRI().crs)
ax.set_extent([-124.19,-121.51,43.4,46.1])
ax.add_image(ShadedReliefESRI(), 10)

gl_major = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
gl_major.xlabels_top = False
gl_major.ylabels_right = False
gl_major.xlabel_style = {'size': 13}
gl_major.ylabel_style = {'size': 13}

# scat = ax.scatter(st_lon,st_lat, s=200 , c='orange' , transform=ccrs.PlateCarree(),edgecolors='k', marker='^')
# for i in range(len(st_name)):
#     ax.text(st_lon[i]+0.02, st_lat[i],str(st_name[i]), fontsize=12, color='k', transform=ccrs.PlateCarree())
    
scat1 = ax.scatter(st_lon_dom,st_lat_dom, s=100 , c='orange' , transform=ccrs.PlateCarree(),edgecolors='k', marker='^')
for i in range(len(st_name_dom)):
    ax.text(st_lon_dom[i]+0.02, st_lat_dom[i],str(st_name_dom[i]), fontsize=8, color='k', transform=ccrs.PlateCarree())

ax.plot(xpn,ypn, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
ax.plot(xph1,yph1, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
ax.plot(xph2,yph2, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)

plt.show()

#%%

# Attach REC to exsisting SW4 infile

infile = '/Users/rshimony/Desktop/WillametteValley/Models/dist2depth_s_grav_n/salem/infile_dist2depth_s_grav_n_salem.txt'

with open(infile, 'a') as fout:
    
    fout.write('\n')
    fout.write('**REC stations**\n')

    for i in range(len(st_name_dom)):
        name = st_name_dom[i]
        lat = st_lat_dom[i]
        lon = st_lon_dom[i]
    
        fout.write('rec lat=%e lon=%e depth=0 file=%s variables=velocity\n' % (lat,lon,name))
    fout.close()
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    