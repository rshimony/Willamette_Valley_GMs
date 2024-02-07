#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:35:21 2023

@author: rshimony
"""

import pandas as pd
import numpy as np
import glob
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point,LineString
import matplotlib.pyplot as plt

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

#%%

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

#%%

inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/wv_poly.csv')
lon_outer = outer_f['wv_poly_lons']
lat_outer = outer_f['wv_poly_lats']

utme_inner1,utmn_inner1 = llz2utm(np.array(lon_inner1),np.array(lat_inner1),projection_zone=10)
utme_inner2,utmn_inner2 = llz2utm(np.array(lon_inner2),np.array(lat_inner2),projection_zone=10)
utme_outer,utmn_outer = llz2utm(np.array(lon_outer),np.array(lat_outer),projection_zone=10)

inner1_poly_utm = np.vstack((utme_inner1, utmn_inner1)).T
inner2_poly_utm = np.vstack((utme_inner2, utmn_inner2)).T
outer_poly_utm  = np.vstack((utme_outer, utmn_outer)).T

full_poly_utm = Polygon(outer_poly_utm, [inner1_poly_utm , inner2_poly_utm])


grd_pnt = Point(404100.0,4806600.0)
if full_poly_utm.contains(grd_pnt) == False:
    print('not inside')
























