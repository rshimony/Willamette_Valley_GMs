#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 16:52:02 2023

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

inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/wv_poly.csv')
lon_outer = outer_f['wv_poly_lons']
lat_outer = outer_f['wv_poly_lats'] 

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

xpnn = np.array([xpn_m])
xph1n = np.array([xph1_m])
xph2n = np.array([xph2_m])

ypnn = np.array([ypn_m])
yph1n = np.array([yph1_m])
yph2n = np.array([yph2_m])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

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

plt.figure(figsize=(8,10))
plt.scatter(xpnt,ypnt,c='r',s=2)
plt.show()

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
    
uniq_f_lat = np.unique(f_lat)
uniq_f_lon = np.unique(f_lon)

f_lat_arr = np.array(f_lat)
f_lon_arr = np.array(f_lon)

diff_f = 0.0001

f_points_idx , f_points_m = poly_latlon2idx(f_lon_arr,f_lat_arr,lon0,lat0,diff_200_lon,diff_200_lat,diff_m)

grav_x = (f_points_m[:,0])
grav_y = (f_points_m[:,1])


#%%

northing = np.arange(0,378000,100)
easting = np.arange(0,216000,100)

inpoly_n = []
inpoly_e = []

counter = 1

for j in northing:
    if counter % 100 == 0:
        print(str(counter) + ' from ' + str(len(northing)))
    for i in easting:
        grd_point = Point(i,j)
        if full_poly_m.contains(grd_point) == True:
            inpoly_e.append(i)
            inpoly_n.append(j)

    counter = counter + 1



# coord_inpoly = (inpoly_n[100] , inpoly_e[100])
# coord_bpoly = (ypnt[0][50] , xpnt[0][50])

# dist = np.sqrt( (coord_inpoly[0] - coord_bpoly[0])**2 + (coord_inpoly[1] - coord_bpoly[1])**2 )
# d = distance.distance(coord_inpoly, coord_bpoly).m


dists = []
counter_d = 1

for i in range(len(inpoly_n)):
    
    if counter_d % 100 == 0:
        print(str(counter_d) + ' from ' + str(len(inpoly_n)))
    ds = []
    
    for j in range(len(xpnt[0])):
    
        coord_inpoly = (inpoly_n[i] , inpoly_e[i])
        coord_bpoly = (ypnt[0][j] , xpnt[0][j])

        d = np.sqrt( (coord_inpoly[0] - coord_bpoly[0])**2 + (coord_inpoly[1] - coord_bpoly[1])**2 )
        # d = distance.distance(coord_inpoly, coord_bpoly).m
        ds.append(d)
        
    d_min = min(ds)
    dists.append(d_min)
    counter_d = counter_d + 1
    
    
  
evdict = {'northing_poly':inpoly_n, 'easting_poly':inpoly_e , 'depth':dists}
evdf = pd.DataFrame(evdict)
evdf.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/basin_depth_poly_cartesian.csv',index=False)


#%%

plt.figure(figsize=(8,10))
plt.scatter(inpoly_e,inpoly_n,c=dists)
plt.colorbar()
plt.show()


#%%

cart_ifile_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/basin_depth_poly_cartesian.csv')

n_poly = cart_ifile_df.northing_poly
e_poly = cart_ifile_df.easting_poly
depth_poly = cart_ifile_df.depth

# ifile_depths = []
# counter_i = 1

# for i in easting:
#     if counter_i % 100 == 0:
#         print(str(counter_i) + ' from ' + str(len(easting)))
#     for j in northing:
#         ifile_point = (i,j)
#         for k in range(len(n_poly)):
#             poly_point = (e_poly[k] , n_poly[k])
#             if ifile_point == poly_point:
#                 z = depth_poly[k]
#             else:
#                 z = 0
            
#             ifile_depths.append(z)


#     counter_i = counter_i + 1

#%%

nx=len(easting)
ny=len(northing)

nmat = 1

ifile_cart=open('ifile_cartesian_4rfile','w+')
header_line="%d %d %d\n"%(nx,ny,nmat)
ifile_cart.write(header_line)

counter = 1

for j in northing:
    if counter % 1000 == 0:
        print(str(counter) + ' from ' + str(len(northing)))
    for i in easting:

        line="%s %s %s\n"%(i,j,0)
        ifile_cart.write(line)
    counter = counter + 1

ifile_cart.close()

#%%

ifile_cart = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_cartesian_4rfile')

ifile_cart_e_ls = []
ifile_cart_n_ls = []
ifile_cart_dists = []

for i in range(1,len(ifile_cart)):
    ifile_cart_e = ifile_cart[i][0]
    ifile_cart_n = ifile_cart[i][1]
    ifile_cart_dist = ifile_cart[i][2]
    
    ifile_cart_e_ls.append(ifile_cart_e)
    ifile_cart_n_ls.append(ifile_cart_n)
    ifile_cart_dists.append(ifile_cart_dist)

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_cart_e_ls,ifile_cart_n_ls,c=ifile_cart_dists)
plt.colorbar()
plt.show()

#%%

counter = 1
inpoly_cart_idx = []    

for i in range(len(ifile_cart_e_ls)):

    if counter % 100 == 0:
        print(str(counter) + ' from ' + str(len(ifile_cart_e_ls)))
        
    grd_point = Point(ifile_cart_e_ls[i],ifile_cart_n_ls[i])
    if full_poly_m.contains(grd_point) == True:
        inpoly_cart_idx.append(i)


    counter = counter + 1

#%%

idxdict = {'inpoly_cart_idx':inpoly_cart_idx}
idxdf = pd.DataFrame(idxdict)
idxdf.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/inpoly_cart_idx.csv',index=False)

#%%

inpoly_cart_idx_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/inpoly_cart_idx.csv')

inpoly_cart_idx = inpoly_cart_idx_df.inpoly_cart_idx


ifile_cart_dists_arr = np.array(ifile_cart_dists)

ifile_cart_dists_arr[inpoly_cart_idx] = depth_poly

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_cart_e_ls,ifile_cart_n_ls,c=ifile_cart_dists_arr)
plt.colorbar()
plt.scatter(grav_x,grav_y,c='r')
plt.show()

#%%

uniq_grav_x = np.unique(grav_x)
uniq_grav_y = np.unique(grav_y)

max_grav_x = max(uniq_grav_x)
min_grav_x = min(uniq_grav_x)
max_grav_y = max(uniq_grav_y)
min_grav_y = min(uniq_grav_y)

grav_dom = Polygon([[min_grav_x,min_grav_y] , [max_grav_x,min_grav_y] , [max_grav_x,max_grav_y] , [min_grav_x,max_grav_y]])

#%%

ifile_x_gdom = []
ifile_y_gdom = []
idx_gdom = []

for i in range(len(ifile_cart_e_ls)):
    grd_pnt = Point(ifile_cart_e_ls[i],ifile_cart_n_ls[i])
    if grav_dom.contains(grd_pnt) == True:
        ifile_x_gdom.append(ifile_cart_e_ls[i])
        ifile_y_gdom.append(ifile_cart_n_ls[i])
        idx_gdom.append(i)

interp = NearestNDInterpolator(list(zip(grav_x, grav_y)), f_z)
interp_depth = interp(ifile_x_gdom, ifile_y_gdom)
interp_depth_m = interp_depth*1000

plt.figure(figsize=(8,10))
plt.scatter(ifile_x_gdom, ifile_y_gdom, c=interp_depth_m)
plt.colorbar()
plt.show()

f_z_int_filt = gaussian_filter(interp_depth_m, sigma=50)

plt.figure(figsize=(8,10))
plt.scatter(ifile_x_gdom, ifile_y_gdom, c=f_z_int_filt)
plt.colorbar()
plt.show() 

#%%

ifile_cart_e_ls_arr = np.array(ifile_cart_e_ls)
ifile_cart_n_ls_arr = np.array(ifile_cart_n_ls)

plt.figure(figsize=(8,10))
plt.scatter(ifile_cart_e_ls,ifile_cart_n_ls,c=ifile_cart_dists_arr)
plt.colorbar()
plt.scatter(ifile_cart_e_ls_arr[inpoly_cart_idx],ifile_cart_n_ls_arr[inpoly_cart_idx],c='r')
plt.scatter(ifile_cart_e_ls_arr[idx_gdom],ifile_cart_n_ls_arr[idx_gdom],c='g')
plt.show()

#%%

inpoly_gdom_xs = []
inpoly_gdom_ys = []
index_gdom_inpoly = []

counter_igp = 1

for i in idx_gdom:
    if counter_igp % 1000 == 0:
        print(str(counter_igp) + ' from ' + str(len(idx_gdom)))
    if i in np.array(inpoly_cart_idx):
        index_gdom_inpoly.append(i)
        inpoly_gdom_xs.append(ifile_cart_e_ls[i])
        inpoly_gdom_ys.append(ifile_cart_n_ls[i])
    
    counter_igp = counter_igp + 1
    
# counter_igp = 1

# for i in idx_gdom:
#     if counter_igp % 1000 == 0:
#         print(str(counter_igp) + ' from ' + str(len(idx_gdom)))
#     for j in inpoly_cart_idx:
        
#         if i == j:
#             index_gdom_inpoly.append(i)
#             inpoly_gdom_xs.append(ifile_cart_e_ls[i])
#             inpoly_gdom_ys.append(ifile_cart_n_ls[i])
    
#     counter_igp = counter_igp + 1

#%%

idx_g_dict = {'inpoly_gdom_xs':inpoly_gdom_xs, 'inpoly_gdom_ys':inpoly_gdom_ys , 'index_gdom_inpoly':index_gdom_inpoly}
idx_g_df = pd.DataFrame(idx_g_dict)
idx_g_df.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/inpoly_gdom_cart_idx.csv',index=False)

#%%

ifile_gbase = np.zeros_like(ifile_cart_dists_arr)
ifile_gbase[idx_gdom] = f_z_int_filt
gbase_depths = ifile_gbase[index_gdom_inpoly]

# gbase_depths_n = ifile_gbase[index_gdom_n_inpoly]

# ifile_dists_grav = (np.array(ifile_cart_dists_arr)**0.5)*18.5 ### this calibration is for the eug_spe version of the basin
ifile_dists_grav = (np.array(ifile_cart_dists_arr)**0.5)*33 ### this calibration is for the srv version of the basin
ifile_dists_grav[index_gdom_inpoly] = gbase_depths

plt.figure(figsize=(8,10))
plt.scatter(ifile_cart_e_ls,ifile_cart_n_ls, c=ifile_dists_grav)
plt.colorbar()
plt.show() 

#%%

topo_layer = np.zeros_like(ifile_dists_grav)

nx=len(easting)
ny=len(northing)

nmat = 1

with open('ifile_cartesian_4rfile_cal_srv.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,nmat)
    f.write(header_line)
    for i in range(len(ifile_dists_grav)):
        
        layer_0 = topo_layer[i] 
        layer_1 = ifile_dists_grav[i]
        
        f.write("%s %s %s %s\n" % (ifile_cart_e_ls[i],ifile_cart_n_ls[i],layer_0,layer_1))

#%%

ifile_cart = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_cartesian_4rfile_cal.txt',skip_header=1)

ifile_cart_e_ls = ifile_cart[:,0]
ifile_cart_n_ls = ifile_cart[:,1]
ifile_cart_dists = ifile_cart[:,3]


ifile_cart_srv = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_cartesian_4rfile_cal_srv.txt',skip_header=1)

ifile_cart_dists_srv = ifile_cart_srv[:,3]

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_cart_e_ls,ifile_cart_n_ls, c=ifile_cart_dists)
plt.colorbar()
plt.show() 

plt.figure(figsize=(8,10))
plt.scatter(ifile_cart_e_ls,ifile_cart_n_ls, c=ifile_cart_dists_srv)
plt.colorbar()
plt.show() 


#%%

# Making lat/lon ifile with the same dimensions as cartesian ifile

ifile_cart_lons = (ifile_cart_e_ls/100)*0.0013 + (-124.19)
ifile_cart_lats = (ifile_cart_n_ls/100)*0.0009 + (43.4)

plt.figure(figsize=(8,10))
plt.scatter(ifile_cart_lons,ifile_cart_lats, c=ifile_cart_dists_srv)
plt.colorbar()
plt.show() 

#%%

nx=len(easting)
ny=len(northing)

nmat = 1

with open('ifile_latlon_srv.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,nmat)
    f.write(header_line)
    for i in range(len(ifile_cart_dists_srv)):
 
        layer_1 = ifile_cart_dists_srv[i]
        
        f.write("%s %s %s\n" % (ifile_cart_lons[i],ifile_cart_lats[i],layer_1))


with open('ifile_latlon_eugspe.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,nmat)
    f.write(header_line)
    for i in range(len(ifile_cart_dists_srv)):
 
        layer_1 = ifile_cart_dists[i]
        
        f.write("%s %s %s\n" % (ifile_cart_lons[i],ifile_cart_lats[i],layer_1))


#%%

ifile_latlon_eugspe = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_latlon_eugspe.txt',skip_header=1)

ifile_latlon_eugspe_lons = ifile_latlon_eugspe[:,0]
ifile_latlon_eugspe_lats = ifile_latlon_eugspe[:,1]
ifile_latlon_eugspe_depths = ifile_latlon_eugspe[:,2]


ifile_latlon_srv = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_latlon_srv.txt',skip_header=1)

ifile_latlon_srv_lons = ifile_latlon_srv[:,0]
ifile_latlon_srv_lats = ifile_latlon_srv[:,1]
ifile_latlon_srv_depths = ifile_latlon_srv[:,2]


plt.figure(figsize=(8,10))
plt.scatter(ifile_latlon_eugspe_lons,ifile_latlon_eugspe_lats, c=ifile_latlon_eugspe_depths)
plt.colorbar()
plt.show() 

plt.figure(figsize=(8,10))
plt.scatter(ifile_latlon_srv_lons,ifile_latlon_srv_lats, c=ifile_latlon_srv_depths)
plt.colorbar()
plt.show() 






















