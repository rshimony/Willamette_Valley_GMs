#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:15:44 2022

@author: rshimony
"""
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import pandas as pd
import matplotlib.pyplot as plt
from geopy import distance
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

#%%

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
# plt.figure(figsize=(8,10))
# plt.plot(xpn, ypn, c="red")
# plt.plot(xph1, yph1, c="red")
# plt.plot(xph2, yph2, c="red")

xpnn = np.array([xpn])
xph1n = np.array([xph1])
xph2n = np.array([xph2])

ypnn = np.array([ypn])
yph1n = np.array([yph1])
yph2n = np.array([yph2])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)


#%%
nmat=1

lons = np.arange(-124.19,-121.51,0.01)
lats = np.arange(43.4,46.7,0.01)

inpoly_lons = []
inpoly_lats = []

for i in lons:
    for j in lats:
        grd_point = Point(i,j)
        if full_poly.contains(grd_point) == True:
            inpoly_lons.append(i)
            inpoly_lats.append(j)
            
#%%

f_path = '/Users/rshimony/Desktop/zzthick4.axyz'

f = np.genfromtxt(f_path)

f_lat = []
f_lon = []
f_z = []

for i in range(len(f)):
    f_lon.append(f[i][0])
    f_lat.append(f[i][1])
    f_z.append(f[i][2])

# plt.figure()
# plt.scatter(f_lon , f_lat , c=f_z)

#%%

uniq_f_lat = np.unique(f_lat)
uniq_f_lon = np.unique(f_lon)

tual_points = np.zeros((len(f_lat),2))

for i in range(len(f_lat)):
    tual_points[i,0] = f_lon[i]
    tual_points[i,1] = f_lat[i]
    
tp = Point(tual_points[2])

#%%

mt_points = [tual_points[0] , tual_points[1] , tual_points[2] , tual_points[3]]

#%%

f_lats_round_str = [ '%.4f' % elem for elem in f_lat ]
f_lats_round = np.float64(f_lats_round_str)

f_lons_round_str = [ '%.4f' % elem for elem in f_lon ]
f_lons_round = np.float64(f_lons_round_str)

lats_round_str = [ '%.4f' % elem for elem in lats ]
lats_round = np.float64(lats_round_str)

lons_round_str = [ '%.4f' % elem for elem in lons ]
lons_round = np.float64(lons_round_str)
#%%

# matched_lats_str = []
# matched_lons_str = []
matched_points = []
matched_depths = []

# for i in f_lons_round_str:
#     for j in lats_round_str:
#         if i == j:
#             matched_lats_str.append(i)
        
# for i in f_lons_round_str:
#     for j in lons_round_str:
#         if i == j:
#             matched_lons_str.append(i)
counter_2 = 1

for i in lons:
    if counter_2 % 10 == 0:
        print(str(counter_2) + ' from ' + str(len(lons)))
    for j in lats:
        gpnt = np.array([i,j])
        for p in range(len(tual_points)):
            if (gpnt.all()) == (tual_points[p].all()):
                matched_points.append(tual_points[p])
                matched_depths.append(f_z[p])
    counter_2 = counter_2 + 1                
            
#%%

# f_lats_round = [ round(elem, 4) for elem in f_lat ]
# f_lats_round_str = [ '%.4f' % elem for elem in f_lat ]
# f_lats_round = np.float64(f_lats_round_str)

# diffs = []
# for i in range(1,len(f_lat)):
#     diff = f_lat[i] - f_lat[i-1]
#     diffs.append(diff)
    
# diffs_avg = np.average(diffs)
    

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


lons = np.arange(-124.19,-121.51,0.01)
lats = np.arange(43.4,46.7,0.01)

nx=len(lons)
ny=len(lats)

ifile=open('ifile_poly_depth_tual','w+')
header_line="%d %d %d\n"%(nx,ny,nmat)
ifile.write(header_line)

counter = 1

for y in lats:
    if counter % 1000 == 0:
        print(str(counter) + ' from ' + str(len(lats)))
    for x in lons:
        
        grd_point = Point(x,y)
    
        if full_poly.contains(grd_point) == True:
            ds = []
            for c in range(len(xpnt[0])):
                coord_inpoly = (y , x)
                coord_bpoly = (ypnt[0][c] , xpnt[0][c])
                d = distance.distance(coord_inpoly, coord_bpoly).km
                ds.append(d)
            d_min = min(ds)
            z = d_min
        
        # gp = np.array([x,y])
        # for p in range(len(tual_points)):
        #     if (gp.all()) == (tual_points[p].all()):
        #         z = f_z[p]
                
        else: 
            z = 0
                                   
        line="%s %s %s\n"%(x,y,z)
        ifile.write(line)
    counter = counter + 1

ifile.close()























