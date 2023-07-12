#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:14:32 2023

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

# inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
# lon_inner1 = inner1_f['inner1_lon']
# lat_inner1 = inner1_f['inner1_lat']

# inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
# lon_inner2 = inner2_f['inner2_lon']
# lat_inner2 = inner2_f['inner2_lat']

# outer_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_outer.csv')
# lon_outer = outer_f['outer_lon']
# lat_outer = outer_f['outer_lat']   

# south_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_south_poly.csv')
# lon_south = south_f['lon']
# lat_south = south_f['lat']  
   
# inner1_poly = np.zeros((len(lon_inner1),2))
# for i in range(len(lon_inner1)):
#     inner1_poly[i,0] = lon_inner1[i]
#     inner1_poly[i,1] = lat_inner1[i]
    
# inner2_poly = np.zeros((len(lon_inner2),2))
# for i in range(len(lon_inner2)):
#     inner2_poly[i,0] = lon_inner2[i]
#     inner2_poly[i,1] = lat_inner2[i]
    
# outer_poly = np.zeros((len(lon_outer),2))
# for i in range(len(lon_outer)):
#     outer_poly[i,0] = lon_outer[i]
#     outer_poly[i,1] = lat_outer[i]
    
# south_poly = np.zeros((len(lon_south),2))
# for i in range(len(lon_south)):
#     south_poly[i,0] = lon_south[i]
#     south_poly[i,1] = lat_south[i]


# full_poly = Polygon(outer_poly, [inner1_poly , inner2_poly])

inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/wv_poly.csv')
lon_outer = outer_f['wv_poly_lons']
lat_outer = outer_f['wv_poly_lats'] 

# utme_inner1,utmn_inner1 = llz2utm(np.array(lon_inner1),np.array(lat_inner1),projection_zone=10)
# utme_inner2,utmn_inner2 = llz2utm(np.array(lon_inner2),np.array(lat_inner2),projection_zone=10)
# utme_outer,utmn_outer = llz2utm(np.array(lon_outer),np.array(lat_outer),projection_zone=10)

# inner1_poly_utm = np.vstack((utme_inner1, utmn_inner1)).T
# inner2_poly_utm = np.vstack((utme_inner2, utmn_inner2)).T
# outer_poly_utm  = np.vstack((utme_outer, utmn_outer)).T

# full_poly_utm = Polygon(outer_poly_utm, [inner1_poly_utm , inner2_poly_utm])



#%%

lat0 = 43.4
lon0 = -124.19

diff_200_lon = 0.0026
diff_200_lat = 0.0018


    
delta_lon = lon_outer - lon0
delta_lat = lat_outer - lat0

dec_idx_x = delta_lon/diff_200_lon
dec_idx_y = delta_lat/diff_200_lat
    

idx_poly_points = np.zeros((len(dec_idx_x),2))
for i in range(len(dec_idx_x)):
    idx_poly_points[i,0] = dec_idx_x[i]
    idx_poly_points[i,1] = dec_idx_y[i]

idx_poly = Polygon(idx_poly_points)

ind_x = np.arange(0,1080,1)
ind_y = np.arange(0,1890,1)

ind_x_grd , ind_y_grd = np.meshgrid(ind_x,ind_y)

xpn, ypn = idx_poly.exterior.xy

plt.figure(figsize=(8,10))
plt.scatter(ind_x_grd,ind_y_grd)
plt.plot(xpn,ypn,c='r')
plt.show()

#%%

def poly_latlon2idx(lons,lats,lon0,lat0,diff_lon,diff_lat,diff_m):
    
    lat0 = lat0
    lon0 = lon0

    diff_lon = diff_lon
    diff_lat = diff_lat


        
    delta_lon = lons - lon0
    delta_lat = lats - lat0

    dec_idx_x = delta_lon/diff_lon
    dec_idx_y = delta_lat/diff_lat
    
    x_m = dec_idx_x * diff_m
    y_m = dec_idx_y * diff_m

    idx_poly_points = np.zeros((len(dec_idx_x),2))
    for i in range(len(dec_idx_x)):
        idx_poly_points[i,0] = dec_idx_x[i]
        idx_poly_points[i,1] = dec_idx_y[i]
        
    m_poly_points = np.zeros((len(dec_idx_x),2))
    for i in range(len(dec_idx_x)):
        m_poly_points[i,0] = x_m[i]
        m_poly_points[i,1] = y_m[i]
   
    return(idx_poly_points,m_poly_points)



# lat0 = 43.4
# lon0 = -124.19

# diff_200_lon = 0.0026
# diff_200_lat = 0.0018

# # idx_poly_points = poly_latlon2idx(lon_outer,lat_outer,lon0,lat0,diff_200_lon,diff_200_lat)

# outer_poly_points = poly_latlon2idx(lon_outer,lat_outer,lon0,lat0,diff_200_lon,diff_200_lat)
# inner1_poly_points = poly_latlon2idx(lon_inner1,lat_inner1,lon0,lat0,diff_200_lon,diff_200_lat)
# inner2_poly_points = poly_latlon2idx(lon_inner2,lat_inner2,lon0,lat0,diff_200_lon,diff_200_lat)

# full_poly_idx = Polygon(outer_poly_points, [inner1_poly_points , inner2_poly_points])

# xpn, ypn = full_poly_idx.exterior.xy 
# xph1, yph1 = full_poly_idx.interiors[0].xy
# xph2, yph2 = full_poly_idx.interiors[1].xy



# idx_poly = Polygon(idx_poly_points)

# ind_x = np.arange(0,1080,1)
# ind_y = np.arange(0,1890,1)

# ind_x_grd , ind_y_grd = np.meshgrid(ind_x,ind_y)

# xpn, ypn = idx_poly.exterior.xy

# plt.figure(figsize=(8,10))
# plt.scatter(ind_x_grd,ind_y_grd)
# plt.plot(xpn,ypn,c='r')
# plt.show()


lat0 = 43.4
lon0 = -124.19

diff_200_lon = 0.0013
diff_200_lat = 0.0009

diff_m = 100

outer_poly_points_idx , outer_poly_points_m = poly_latlon2idx(lon_outer,lat_outer,lon0,lat0,diff_200_lon,diff_200_lat,diff_m)
inner1_poly_points_idx , inner1_poly_points_m = poly_latlon2idx(lon_inner1,lat_inner1,lon0,lat0,diff_200_lon,diff_200_lat,diff_m)
inner2_poly_points_idx , inner2_poly_points_m = poly_latlon2idx(lon_inner2,lat_inner2,lon0,lat0,diff_200_lon,diff_200_lat,diff_m)

full_poly_idx = Polygon(outer_poly_points_idx, [inner1_poly_points_idx , inner2_poly_points_idx])
full_poly_m = Polygon(outer_poly_points_m, [inner1_poly_points_m , inner2_poly_points_m])

xpn_idx, ypn_idx = full_poly_idx.exterior.xy 
xph1_idx, yph1_idx = full_poly_idx.interiors[0].xy
xph2_idx, yph2_idx = full_poly_idx.interiors[1].xy

xpn_m, ypn_m = full_poly_m.exterior.xy 
xph1_m, yph1_m = full_poly_m.interiors[0].xy
xph2_m, yph2_m = full_poly_m.interiors[1].xy

plt.figure(figsize=(8,10))
plt.plot(xpn_idx,ypn_idx,c='r')
plt.plot(xph1_idx,yph1_idx,c='g')
plt.plot(xph2_idx,yph2_idx,c='b')
plt.show()

plt.figure(figsize=(8,10))
plt.plot(xpn_m,ypn_m,c='r')
plt.plot(xph1_m,yph1_m,c='g')
plt.plot(xph2_m,yph2_m,c='b')
plt.show()

#%%




