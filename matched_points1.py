#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:00:29 2022

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

tual_points = np.zeros((len(f_lat),2))

for i in range(len(f_lat)):
    tual_points[i,0] = f_lon[i]
    tual_points[i,1] = f_lat[i]

#%%

plt.figure(figsize=(8,10))
plt.scatter(f_lon, f_lat, c=f_z)
plt.colorbar()
plt.show()

#%%

f_z_filt = gaussian_filter(f_z, sigma=30)

plt.figure(figsize=(8,10))
plt.scatter(f_lon, f_lat, c=f_z_filt)
plt.colorbar()
plt.show()

#%%
matched_points = []
matched_depths = []

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

ifile = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/dist2depth_model_novs30_nograv/Springfield/cal_depth_ifile_m.txt')

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

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons, ifile_lats, c=ifile_dists , cmap='gist_earth_r')
plt.show() 

#%%

ifile_dist_grav = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/dist2depth_s_grav_n/cal_depth_ifile_dist_grav.txt')

ifile_dist_grav_lons = []
ifile_dist_grav_lats = []
ifile_dist_grav_dists = []

for i in range(1,len(ifile_dist_grav)):
    ifile_dist_grav_lon = ifile[i][0]
    ifile_dist_grav_lat = ifile[i][1]
    ifile_dist_grav_dist = ifile[i][2]
    
    ifile_dist_grav_lons.append(ifile_lon)
    ifile_dist_grav_lats.append(ifile_lat)
    ifile_dist_grav_dists.append(ifile_dist)

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
        
#%%

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

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons_gdom, ifile_lats_gdom, c='r')
plt.show()     
            

#%%

interp = NearestNDInterpolator(list(zip(f_lon, f_lat)), f_z)
interp_depth = interp(ifile_lons_gdom, ifile_lats_gdom)
interp_depth_m = interp_depth*1000

interp_n = NearestNDInterpolator(list(zip(f_lon, f_lat)), f_z)
interp_depth_n = interp_n(ifile_lons_gdom_n, ifile_lats_gdom_n)
interp_depth_m_n = interp_depth_n*1000

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons_gdom, ifile_lats_gdom, c=interp_depth_m)
plt.colorbar()
plt.show()

#%%

f_z_int_filt = gaussian_filter(interp_depth_m, sigma=50)

f_z_int_filt_n = gaussian_filter(interp_depth_m_n, sigma=50)

# plt.figure(figsize=(8,10))
# plt.scatter(ifile_lons_gdom, ifile_lats_gdom, c=f_z_int_filt)
# plt.colorbar()
# plt.show()
#%%
counter_g = 1

for i in range(len(ifile_lons)):
    if counter_g % 10 == 0:
        print(str(counter_g) + ' from ' + str(len(ifile_lons)))
    grd_pnt = Point(ifile_lons[i],ifile_lats[i])
    gpn_t = np.array([ifile_lons[i],ifile_lats[i]])
    if grav_dom.contains(grd_pnt) == True:
        for p in range(len(ifile_lons_gdom)):
            gpn_grv = np.array([ifile_lons_gdom[p],ifile_lats_gdom[p]])
            if (gpn_t.all()) == (gpn_grv.all()):
                ifile_dists[i] = interp_depth[p]
    counter_g = counter_g + 1                

#%%

# ifile_point = np.array([ifile_lons_gdom[0],ifile_lats_gdom[0]])

# ifile_lon_arr = np.array(ifile_lons)
# ifile_lat_arr = np.array(ifile_lats)

# ifile_lon_gdom_arr = np.array(ifile_lons_gdom)
# ifile_lat_gdom_arr = np.array(ifile_lats_gdom)

# rep_arr_lon = np.tile(ifile_lon_arr[0], (len(ifile_lons_gdom),1))

# subt_arr = np.abs(rep_arr_lon - ifile_lon_gdom_arr)



#%%

samp_arr = np.array([1,2,3,4,5])

rep_ind = [1,2]

samp_arr[rep_ind] = [0,10]


#%%

ifile_dists_arr = np.array(ifile_dists)

# ifile_dists_arr[index] = interp_depth_m

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons, ifile_lats, c=ifile_dists_arr)
plt.colorbar()
plt.show()


#%%

with open('cal_depth_ifile_grav.txt', 'w') as f:
    for i in range(len(ifile_dists)):
        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],ifile_dists_arr[i])) 

#%%
topo_layer = np.zeros_like(ifile_dists_arr)
domain_layer = (np.zeros_like(ifile_dists_arr))+11000

with open('cal_depth_ifile_grav_4rfile.txt', 'w') as f:
    for i in range(len(ifile_dists)):
        f.write("%s %s %s %s %s\n" % (ifile_lons[i],ifile_lats[i],topo_layer[i],ifile_dists_arr[i],domain_layer[i])) 


#%%

uniq_ifile_lat = np.unique(ifile_lats)
uniq_ifile_lon = np.unique(ifile_lons)

max_db_lat = min(uniq_f_lat)
min_db_lat = min(uniq_ifile_lat)
max_db_lon = -122.0
min_db_lon = -123.2

db_dom = Polygon([[min_db_lon,min_db_lat] , [max_db_lon,min_db_lat] , [max_db_lon,max_db_lat] , [min_db_lon,max_db_lat]])

#%%

db_lons = []
db_lats = []
db_depths=[]
index_db = []

for i in range(len(ifile_lons)):
    grd_pnt = Point(ifile_lons[i],ifile_lats[i])
    if db_dom.contains(grd_pnt) == True:
        db_lons.append(ifile_lons[i])
        db_lats.append(ifile_lats[i])
        db_depths.append(ifile_dists[i])
        index_db.append(i)



#%%

w_edge_lats = []
w_edge_lons = []
w_edge_depths = []

for i in range(len(db_lons)):
    grd_pnt = Point(db_lons[i],db_lats[i])
    if db_lons[i] == min(db_lons):
        w_edge_lons.append(db_lons[i])
        w_edge_lats.append(db_lats[i])
        w_edge_depths.append(db_depths[i])

#%%

e_edge_lats = []
e_edge_lons = []

for i in range(len(db_lons)):
    grd_pnt = Point(db_lons[i],db_lats[i])
    if db_lons[i] == max(db_lons):
        e_edge_lons.append(db_lons[i])
        e_edge_lats.append(db_lats[i])

e_edge_lats_arr = np.array(e_edge_lats)        
e_edge_depth = np.zeros_like(e_edge_lats_arr)
e_edge_depth = e_edge_depth + 4000
#%%

uniq_db_lons = np.unique(db_lons)

coord_w_edge_t = (w_edge_lats[1754] , w_edge_lons[1754])
coord_e_edge_t = (e_edge_lats[1754] , e_edge_lons[1754])
d_t = distance.distance(coord_w_edge_t, coord_e_edge_t).m
coord_w_edge_b = (w_edge_lats[0] , w_edge_lons[0])
coord_e_edge_b = (e_edge_lats[0] , e_edge_lons[0])
d_b = distance.distance(coord_w_edge_b, coord_e_edge_b).m

D = (d_b + d_t)/2

lon_dists_b = []
lon_dists_t = []
for i in range(len(uniq_db_lons)):
    d_b_lons = distance.distance(coord_w_edge_b, (w_edge_lats[0] , uniq_db_lons[i])).m
    d_t_lons = distance.distance(coord_w_edge_t, (w_edge_lats[1745] , uniq_db_lons[i])).m
    lon_dists_b.append(d_b_lons)
    lon_dists_t.append(d_t_lons)
    
lon_dists = (np.array(lon_dists_b) + np.array(lon_dists_t))/2
#%%
counter_d = 1

db_depths_dip=[]
for i in range(len(e_edge_lats)):
    print(str(counter_d) + ' of ' + str(len(e_edge_lats)))
    m = (e_edge_depth[i] - w_edge_depths[i])/D
    for j in range(len(uniq_db_lons)):
        depth = m * lon_dists[j] + w_edge_depths[i]
        db_depths_dip.append(depth)
        
    counter_d = counter_d + 1   
    

#%%

ifile_grav = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/rfile_combine/cal_depth_ifile_grav.txt')

ifile_grav_lons = []
ifile_grav_lats = []
ifile_grav_dists = []

for i in range(len(ifile_grav)):
    ifile_lon = ifile_grav[i][0]
    ifile_lat = ifile_grav[i][1]
    ifile_dist = ifile_grav[i][2]
    
    ifile_grav_lons.append(ifile_lon)
    ifile_grav_lats.append(ifile_lat)
    ifile_grav_dists.append(ifile_dist)

#%%

ifile_grav_dists_arr = np.array(ifile_grav_dists)

ifile_grav_dists_arr[index_db] = db_depths_dip

#%%

with open('cal_depth_ifile_grav_dip.txt', 'w') as f:
    for i in range(len(ifile_dists)):
        f.write("%s %s %s\n" % (ifile_grav_lons[i],ifile_grav_lats[i],ifile_grav_dists_arr[i])) 

#%%

ifile_grav_dip = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/rfile_combine/cal_depth_ifile_grav_dip.txt')

ifile_grav_dip_lons = []
ifile_grav_dip_lats = []
ifile_grav_dip_dists = []

for i in range(len(ifile_grav_dip)):
    ifile_lon = ifile_grav_dip[i][0]
    ifile_lat = ifile_grav_dip[i][1]
    ifile_dist = ifile_grav_dip[i][2]
    
    ifile_grav_dip_lons.append(ifile_lon)
    ifile_grav_dip_lats.append(ifile_lat)
    ifile_grav_dip_dists.append(ifile_dist)

#%%

uniq_ifile_gd_lons = np.unique(ifile_grav_dip_lons)
uniq_ifile_gd_lats = np.unique(ifile_grav_dip_lats)

ifile_gd_lon_diffs = []

for i in range(1,len(uniq_ifile_gd_lons)):
    lon_diff = uniq_ifile_gd_lats[i-1] - uniq_ifile_gd_lats[i]
    ifile_gd_lon_diffs.append(lon_diff)

#%%

plt.figure(figsize=(8,10))
# plt.scatter(ifile_lons, ifile_lats, c=ifile_dists_arr)
# plt.scatter(ifile_grav_lons, ifile_grav_lats, c=ifile_grav_dists_arr)
plt.scatter(ifile_grav_dip_lons, ifile_grav_dip_lats, c=ifile_grav_dip_dists)
# plt.scatter(db_lons, db_lats, c='r')
# plt.scatter(w_edge_lons, w_edge_lats, c=w_edge_depths, cmap='jet')
# plt.scatter(e_edge_lons, e_edge_lats, c=e_edge_depth, cmap='jet')
# plt.scatter(db_lons, db_lats, c=db_depths_dip, cmap='jet')
plt.colorbar()
# plt.savefig('/Users/rshimony/Desktop/map_ifile_grav_dip.png',dpi=200,bbox_inches='tight')
plt.show()

#%%

topo_layer = np.zeros_like(ifile_grav_dip_dists)
domain_layer = (np.zeros_like(ifile_grav_dip_dists))+11000

nx=len(np.unique(ifile_grav_dip_lons))
ny=len(np.unique(ifile_grav_dip_lats))



with open('cal_depth_ifile_grav_dip.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,3)
    f.write(header_line)
    for i in range(len(ifile_grav_dip_dists)):
        # f.write("%s %s %s %s %s\n" % (ifile_grav_dip_lons[i],ifile_grav_dip_lats[i],topo_layer[i],ifile_grav_dip_dists[i],domain_layer[i])) 
        f.write("%s %s %s\n" % (ifile_grav_dip_lons[i],ifile_grav_dip_lats[i],ifile_grav_dip_dists[i]))


#%%

from PIL import Image

img = np.array(Image.open('/Users/rshimony/Desktop/WillametteValley/Models/rfile_combine/Grids_2_hh=100m/elev_layer2_hh=100m.tif'))

print(type(img))
print(img.shape)

Image.fromarray(np.flipud(img)).save('/Users/rshimony/Desktop/WillametteValley/Models/rfile_combine/Grids_2_hh=100m/elev_layer2_hh=100m_fliped.tif')

#%%

sing_lon_depth = []
sing_lon_depth_filt = []

for i in range(len(f_z)):
    if -123.0098 <= f_lon[i] <= -122.997:
        sing_lon_depth.append(f_z[i])
        sing_lon_depth_filt.append(f_z_filt[i])
        
#%%

plt.figure(figsize=(8,10))
plt.plot(sing_lon_depth)
plt.plot(sing_lon_depth_filt)
plt.show()
        
#%%
sing_lon_depth_int = []
sing_lon_depth_filt_int = []

for i in range(len(f_z_int_filt)):
    if -123.0096 <= ifile_lons_gdom[i] <= -122.9979:
        sing_lon_depth_int.append(interp_depth_m[i])
        sing_lon_depth_filt_int.append(f_z_int_filt[i])
        
#%%

plt.figure(figsize=(8,10))
plt.plot(sing_lon_depth_int)
plt.plot(sing_lon_depth_filt_int)
plt.show()

#%%
counter_g = 1

bot_dom_grav = []
top_dom_grav = []
e_dom_grav = []
w_dom_grav = []

for i in range(len(f_lon)):
    if counter_g % 100 == 0:
        print(str(counter_g) + ' of ' + str(len(f_lon)))
    if f_lat[i] == min(np.unique(f_lat)):
        bot_dom_grav.append(f_z[i])
    if f_lat[i] == max(np.unique(f_lat)):
        top_dom_grav.append(f_z[i])
    if f_lon[i] == max(np.unique(f_lon)):
        e_dom_grav.append(f_z[i])
    if f_lon[i] == min(np.unique(f_lon)):
        w_dom_grav.append(f_z[i])

    counter_g = counter_g + 1
    
#%%

plt.figure(figsize=(8,10))
plt.subplot(2, 2, 1)
plt.plot(bot_dom_grav)
plt.title('Bottom')
plt.subplot(2, 2, 2)
plt.plot(top_dom_grav)
plt.title('Top')
plt.subplot(2, 2, 3)
plt.plot(e_dom_grav)
plt.title('East')
plt.subplot(2, 2, 4)
plt.plot(w_dom_grav)
plt.title('West')

#%%
#bl = base layer

uniq_ifile_lat = np.unique(ifile_lats)
uniq_ifile_lon = np.unique(ifile_lons)

max_bl_lat = max(uniq_ifile_lat)
min_bl_lat = 44.0
max_bl_lon = -122.2
min_bl_lon = min(uniq_f_lon)

bl_dom = Polygon([[min_bl_lon,min_bl_lat] , [max_bl_lon,min_bl_lat] , [max_bl_lon,max_bl_lat] , [min_bl_lon,max_bl_lat]])

#%%

ifile_lons_bl_dom = []
ifile_lats_bl_dom = []
index_bl_dom = []

for i in range(len(ifile_lons)):
    grd_pnt = Point(ifile_lons[i],ifile_lats[i])
    if bl_dom.contains(grd_pnt) == True:
        ifile_lons_bl_dom.append(ifile_lons[i])
        ifile_lats_bl_dom.append(ifile_lats[i])
        index_bl_dom.append(i)

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c='b')
plt.scatter(ifile_lons_bl_dom,ifile_lats_bl_dom,c='r')
plt.show()

#%%

points = np.array([(44.0, 0),(44.1, 1600),(44.2, 3100),(44.3, 4600), (44.5, 5500), (44.7, 4600), (44.9, 3300) , (45.0, 2300), (45.1, 1450),
                   (45.2, 900),(45.3, 500),(45.5, 0),(45.8, 0),(46.2, 0),(46.5, 0),(46.7, 0)])
# get x and y vectors
x = points[:,0]
y = points[:,1]

# calculate polynomial
z = np.polyfit(x, y, 8)
f = np.poly1d(z)

# calculate new x's and y's
x_new = np.unique(ifile_lats_bl_dom)
y_new = f(x_new)

plt.plot(x,y,'o', x_new, y_new)
# plt.xlim([x[0]-1, x[-1] + 1 ])
plt.show()

#%%
spine_depth = np.abs(y_new)

unique_lons_bl_dom = np.unique(ifile_lons_bl_dom)
unique_lats_bl_dom = np.unique(ifile_lats_bl_dom)
    
def calc_parabola_vertex(x1, y1, x2, y2, x3, y3):
		'''
		Adapted and modifed to get the unknowns for defining a parabola:
		http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
		'''

		denom = (x1-x2) * (x1-x3) * (x2-x3);
		A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom;
		B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom;
		C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom;

		return A,B,C
    

plt.plot(x,y,'o', x_new, spine_depth)
# plt.xlim([x[0]-1, x[-1] + 1 ])
plt.show()
#%%
x1,y1=[unique_lons_bl_dom.max(),0]
x2,y2=[-123.0,spine_depth[0]]
x3,y3=[unique_lons_bl_dom.min(),0]

a,b,c = calc_parabola_vertex(x1, y1, x2, y2, x3, y3)

x = unique_lons_bl_dom
y = (a*(x**2))+(b*x)+c

plt.plot(x,y)
# plt.xlim([x[0]-1, x[-1] + 1 ])
plt.show()

#%%

points_bl = np.array([(unique_lons_bl_dom.max(),0) , (-123.3,spine_depth[500]/2.2), (-123.0,(spine_depth[500]/1.3)) ,
                      (-122.5,(spine_depth[500])) , (unique_lons_bl_dom.min(),0)])

x_bl = points_bl[:,0]
y_bl = points_bl[:,1]

# calculate polynomial
z_bl = np.polyfit(x_bl, y_bl, 4)
f_bl = np.poly1d(z_bl)

# calculate new x's and y's
x_new_bl = unique_lons_bl_dom
y_new_bl = f_bl(x_new_bl)

plt.plot(x_bl,y_bl,'o', x_new_bl, y_new_bl)
# plt.xlim([x[0]-1, x[-1] + 1 ])
plt.show()

#%%

'''
loop over unique lats. for each lat, loop over unique lons, calculate parabula with the vertex from the depth from spine_depth for Y and -123.0 for X, 
and the two edge points will be the max and min lons for Xs and 0s for Ys. 
Then take the resultant ABC and use them to calculate depth (Y) at each point (loop over all lons with current lat)
'''    

# counter_bl = 1
# bl_depths = []    

# for i in range(len(unique_lats_bl_dom)):
#     if counter_bl % 100 == 0:
#         print(str(counter_bl) + ' of ' + str(len(unique_lats_bl_dom)))
#     for j in range(len(unique_lons_bl_dom)):
        
#         x1,y1=[unique_lons_bl_dom.max(),0]
#         x2,y2=[-123.0,spine_depth[i]]
#         x3,y3=[unique_lons_bl_dom.min(),0]
        
#         a,b,c = calc_parabola_vertex(x1, y1, x2, y2, x3, y3)

#         x_i = unique_lons_bl_dom[j]
#         y_i = (a*(x_i**2))+(b*x_i)+c
        
#         bl_depths.append(y_i)
        
#     counter_bl = counter_bl + 1
    
    
counter_bl = 1
bl_depths = []    

for i in range(len(unique_lats_bl_dom)):
    if counter_bl % 100 == 0:
        print(str(counter_bl) + ' of ' + str(len(unique_lats_bl_dom)))
    for j in range(len(unique_lons_bl_dom)):
        
        points_bl = np.array([(unique_lons_bl_dom.max(),0) , (-123.3,spine_depth[i]/2.2), (-123.0,(spine_depth[i]/1.3)) ,
                              (-122.5,(spine_depth[i])) , (unique_lons_bl_dom.min(),0)])
        
        x_bl = points_bl[:,0]
        y_bl = points_bl[:,1]
        
        # calculate polynomial
        z_bl = np.polyfit(x_bl, y_bl, 4)
        f_bl = np.poly1d(z_bl)
        
        # calculate new x's and y's
        x_i_bl = unique_lons_bl_dom[j]
        y_i_bl = f_bl(x_i_bl)
        
        bl_depths.append(y_i_bl)
        
    counter_bl = counter_bl + 1
    

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c='b')
plt.scatter(ifile_lons_bl_dom,ifile_lats_bl_dom,c=bl_depths)
plt.colorbar()
plt.show()   

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

s_poly = Polygon(south_poly)

# full_poly = poly.union(s_poly)
# 
xpn, ypn = full_poly.exterior.xy
xph1, yph1 = full_poly.interiors[0].xy
xph2, yph2 = full_poly.interiors[1].xy

xpns, ypns = s_poly.exterior.xy


plt.figure(figsize=(8,10))
plt.plot(xpn,ypn,c='b')
plt.plot(xph1,yph1,c='b')
plt.plot(xph2,yph2,c='b')
plt.plot(xpns,ypns,c='r')
plt.show()  


#%%

w_pol_edge_lons = []
w_pol_edge_lats = []

for i in range((len(xpn))):
    if xpn[i] <= -123.15 and ypn[i] <= 44.9:
        w_pol_edge_lons.append(xpn[i])
        w_pol_edge_lats.append(ypn[i])
        
w_pol_edge_lons.append(max_bl_lon)
w_pol_edge_lats.append(min_bl_lat)

w_pol_edge_lons.insert(0, max_bl_lon)
w_pol_edge_lats.insert(0, 44.9)

      
w_pol_edge = np.zeros((len(w_pol_edge_lats),2))
for i in range(len(w_pol_edge_lats)):
    w_pol_edge[i,0] = w_pol_edge_lons[i]
    w_pol_edge[i,1] = w_pol_edge_lats[i]
         
w_bl_dom = Polygon(w_pol_edge)

xpnw, ypnw = w_bl_dom.exterior.xy

plt.figure(figsize=(8,10))
plt.plot(xpnw,ypnw,c='b')
plt.show()  

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
        
inpoly_s_lons = []
inpoly_s_lats = []
inpoly_s_index = []

for i in range(len(ifile_lons)):
    grd_point = Point(ifile_lons[i],ifile_lats[i])
    if s_poly.contains(grd_point) == True:
        inpoly_s_lons.append(ifile_lons[i])
        inpoly_s_lats.append(ifile_lats[i])
        inpoly_s_index.append(i)
        
plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c='b')
plt.scatter(inpoly_lons,inpoly_lats,c='r')
plt.scatter(inpoly_s_lons,inpoly_s_lats,c='gold')
plt.show()  

#%%

inpoly_dict = {'inpoly_lons':inpoly_lons,'inpoly_lats':inpoly_lats, 'inpoly_index':inpoly_index,}
inpoly_df = pd.DataFrame(inpoly_dict)
inpoly_df.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_index.csv',index=False)

#%%

dists_s_poly = []
counter = 1

for i in range(len(inpoly_s_lats)):
    
    if counter % 100 == 0:
        print(str(counter) + ' from ' + str(len(inpoly_s_lats)))
    ds = []
    
    for j in range(len(xpns)):
    
        coord_inpoly = (inpoly_s_lats[i] , inpoly_s_lons[i])
        coord_bpoly = (ypns[j] , xpns[j])

        d = (distance.distance(coord_inpoly, coord_bpoly).m)*0.2
        
        ds.append(d)
        
    d_min = min(ds)
    dists_s_poly.append(d_min)
    counter = counter + 1
    
    
    
# nx=len(np.unique(ifile_lons))
# ny=len(np.unique(ifile_lats))



# with open('cal_depth_ifile_dists.txt', 'w') as f:
#     header_line="%d %d %d\n"%(nx,ny,1)
#     f.write(header_line)
#     for i in range(len(ifile_lons)):
#         # f.write("%s %s %s %s %s\n" % (ifile_grav_dip_lons[i],ifile_grav_dip_lats[i],topo_layer[i],ifile_grav_dip_dists[i],domain_layer[i])) 
#         f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],dists[i]))
    
#%%

dists_s_poly_arr = (np.array(dists_s_poly)**0.5)*33

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_dists)
plt.scatter(inpoly_s_lons,inpoly_s_lats,c=dists_s_poly_arr)
plt.show()     

#%%

inpoly_depths = []
counter_ip = 1

for i in range(len(ifile_lons)):
    if counter_ip % 100 == 0:
        print(str(counter_ip) + ' from ' + str(len(ifile_lons)))
    for j in inpoly_index:
        if i == j:
            inpoly_depths.append(ifile_dists[i])
    
    counter_ip = counter_ip + 1
    
    
#%%

inpoly_dists = []
for i in range(len(inpoly_index)):
    inpoly_dist = ifile_dists[inpoly_index[i]]
    inpoly_dists.append(inpoly_dist)
    
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

south_patch_max_lat = 44.2
south_patch_min_lat = 44.1
south_patch_max_lon = -122.95
south_patch_min_lon = -123.2

south_patch_dom = Polygon([[south_patch_min_lon,south_patch_min_lat] , [south_patch_max_lon,south_patch_min_lat] , 
                          [south_patch_max_lon,south_patch_max_lat] , [south_patch_min_lon,south_patch_max_lat]])



ifile_lons_south_patch_dom = []
ifile_lats_south_patch_dom = []
index_south_patch_dom = []

for i in range(len(ifile_lons)):
    grd_pnt = Point(ifile_lons[i],ifile_lats[i])
    if south_patch_dom.contains(grd_pnt) == True:
        ifile_lons_south_patch_dom.append(ifile_lons[i])
        ifile_lats_south_patch_dom.append(ifile_lats[i])
        index_south_patch_dom.append(i)    
#%%

# ifile_grav_dists_arr = np.array(ifile_grav_dists)

# ifile_grav_dists_arr[index_db] = db_depths_dip  

inpoly_dists_arr = np.array(inpoly_dists)
inpoly_dists_norm = inpoly_dists_arr + 1000
    
ifile_bl = np.zeros_like(ifile_dists)

ifile_bl[index_bl_dom] = bl_depths

ifile_bl[inpoly_index] = inpoly_dists_norm

ifile_bl[index] = (f_z_int_filt + 500) 

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_bl)
plt.colorbar()
plt.show()   

#%%

index_gdom_inpoly = []

counter_igp = 1

for i in index:
    if counter_igp % 1000 == 0:
        print(str(counter_igp) + ' from ' + str(len(index)))
    if i in inpoly_index:
        index_gdom_inpoly.append(i)
    
    counter_igp = counter_igp + 1
    
#%%

index_gdom_n_inpoly = []

counter_igp_n = 1

for i in index_n:
    if counter_igp_n % 1000 == 0:
        print(str(counter_igp_n) + ' from ' + str(len(index_n)))
    if i in inpoly_index:
        index_gdom_n_inpoly.append(i)
    
    counter_igp_n = counter_igp_n + 1
#%%

ifile_pbase = np.zeros_like(ifile_dists)

ifile_pbase[index_bl_dom] = bl_depths

#%%
# bl_depths_inpoly = []

# for i in inpoly_index:
#     bl_depths_inpoly.append(ifile_pbase[i])


# ifile_poly_bl = np.zeros_like(ifile_dists)

# ifile_poly_bl[inpoly_index] = bl_depths_inpoly

ifile_gbase = np.zeros_like(ifile_dists)

ifile_gbase[index] = f_z_int_filt
# ifile_gbase[index] = interp_depth_m
  
gbase_depths = ifile_gbase[index_gdom_inpoly]

ifile_dists_gbase = np.array(ifile_dists)[index_gdom_inpoly]

ifile_dists_edge_indx = []
for i in range(len(ifile_dists_gbase)):
    if ifile_dists_gbase[i] < 500:
        ifile_dists_edge_indx.append(index_gdom_inpoly[i])
        
#%%

ifile_gbase_n = np.zeros_like(ifile_dists)

ifile_gbase_n[index_n] = f_z_int_filt_n
  
gbase_depths_n = ifile_gbase_n[index_gdom_n_inpoly]

ifile_dists_gbase_n = np.array(ifile_dists)[index_gdom_n_inpoly]

ifile_dists_edge_indx_n = []
for i in range(len(ifile_dists_gbase_n)):
    if ifile_dists_gbase_n[i] < 500:
        ifile_dists_edge_indx_n.append(index_gdom_n_inpoly[i])

#%%        
ifile_dists_edge = np.array(ifile_dists)[ifile_dists_edge_indx]

ifile_dists_edge_norm = ifile_dists_edge / (max(ifile_dists_edge))

ifile_gbase[ifile_dists_edge_indx] = ifile_gbase[ifile_dists_edge_indx] * ifile_dists_edge_norm

gbase_depths_edge_norm = ifile_gbase[index_gdom_inpoly]
     

# ifile_dists_norm_max = ifile_dists / (max(ifile_dists))

# ifile_dists_norm_max_gbase = ifile_dists_norm_max[index_gdom_inpoly]

# gbase_depths_norm = gbase_depths * ifile_dists_norm_max_gbase


# ifile_poly_bl[index_gdom_inpoly] = gbase_depths

# ifile_poly_bl[index_east_patch_dom] = 0

# plt.figure(figsize=(8,10))
# plt.scatter(ifile_lons,ifile_lats,c=ifile_poly_bl)
# plt.colorbar()
# plt.show()   


#%%        
ifile_dists_edge_n = np.array(ifile_dists)[ifile_dists_edge_indx_n]

ifile_dists_edge_norm_n = ifile_dists_edge_n / (max(ifile_dists_edge_n))

ifile_gbase[ifile_dists_edge_indx_n] = ifile_gbase[ifile_dists_edge_indx_n] * ifile_dists_edge_norm_n

gbase_depths_edge_norm_n = ifile_gbase[index_gdom_n_inpoly]    
#%%

ifile_dists_norm = ((np.array(ifile_dists))**0.5)*20

#%%

# ifile_gbase = np.zeros_like(ifile_dists)

# ifile_gbase[index] = f_z_int_filt
# ifile_gbase[index] = interp_depth_m
  
# gbase_depths = ifile_gbase[index_gdom_inpoly]


# ifile_poly_bl[index_gdom_inpoly] = gbase_depths


# ifile_dist_grav = np.array(ifile_dists)
ifile_dist_grav = np.zeros_like(ifile_dists)

ifile_dist_grav[index_south_patch_dom] = 500
#%%
dists_norm = ifile_dists_norm[inpoly_index]
#%%
ifile_dist_grav[inpoly_index] = dists_norm

ifile_dist_grav[index_gdom_inpoly] = gbase_depths_edge_norm

ifile_dist_grav[inpoly_s_index] = dists_s_poly_arr

ifile_dist_grav[index_east_patch_dom] = 0

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_dist_grav,cmap='gist_earth_r')
plt.colorbar()
plt.show()

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_pbase,cmap='gist_earth_r')
plt.colorbar()
plt.show()

#%%
ifile_dist_grav[inpoly_index] = dists_norm

ifile_dist_grav[index_gdom_n_inpoly] = gbase_depths_edge_norm_n

ifile_dist_grav[inpoly_s_index] = dists_s_poly_arr

ifile_dist_grav[index_east_patch_dom] = 0

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_dist_grav,cmap='gist_earth_r')
plt.colorbar()
plt.show()

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_pbase,cmap='gist_earth_r')
plt.colorbar()
plt.show()

#%%
ifile_dis2depth = ifile_dists_norm*3
ifile_dis2depth[inpoly_s_index] = dists_s_poly_arr
ifile_dis2depth[index_east_patch_dom] = 0

plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_dis2depth,cmap='gist_earth_r')
plt.colorbar()
plt.show()

#%%

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))

with open('cal_depth_ifile_dist.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_dist_grav)):
        
        layer_1 = ifile_dis2depth[i]
       

        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1))

#%%        

with open('cal_depth_ifile_dist_grav_top.txt', 'w') as f:
    for i in range(len(ifile_dist_grav)):
        
        layer_1 = ifile_dist_grav[i]
       

        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1))
 
        
with open('cal_depth_ifile_dist_grav_bot.txt', 'w') as f:
    for i in range(len(ifile_dist_grav)):
        
        layer_1 = ifile_pbase[i]
       

        f.write("%s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1))

#%%

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))



with open('cal_depth_ifile_dist_grav.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_dist_grav)):
        
        layer_1 = ifile_dist_grav[i]
        layer_2 = ifile_pbase[i]
        
        if ifile_dist_grav[i] > ifile_pbase[i]:
            layer_2 = layer_1
        

        f.write("%s %s %s %s\n" % (ifile_lons[i],ifile_lats[i],layer_1,layer_2))
        
#%%

ifile_dist_grav = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/dist2depth_s_grav_n/salem/short/cal_depth_ifile_dist_grav.txt')

ifile_dist_grav_lons = []
ifile_dist_grav_lats = []
ifile_dist_grav_dists = []
# ifile_dist_grav_dips = []

for i in range(len(ifile_dist_grav)):
    ifile_lon = ifile_dist_grav[i][0]
    ifile_lat = ifile_dist_grav[i][1]
    ifile_dist = ifile_dist_grav[i][2]
    # ifile_dip = ifile_dist_grav[i][3]
    
    ifile_dist_grav_lons.append(ifile_lon)
    ifile_dist_grav_lats.append(ifile_lat)
    ifile_dist_grav_dists.append(ifile_dist)
    # ifile_dist_grav_dips.append(ifile_dip)
    
#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_dist_grav_lons,ifile_dist_grav_lats,c=ifile_dist_grav_dists)
plt.colorbar()
plt.show()

plt.figure(figsize=(8,10))
plt.scatter(ifile_dist_grav_lons,ifile_dist_grav_lats,c=ifile_dist_grav_dips)
plt.colorbar()
plt.show()

#%%
single_lon_cs_lats = []
single_lon_cs_depths = []

counter_cs = 1

for i in range(len(ifile_lats)):
    if counter_cs % 1000 == 0:
        print(str(counter_cs) + ' from ' + str(len(ifile_lats)))
    if 123.0 < ifile_lons[i] < 123.05:
        single_lon_cs_lats.append(ifile_lats[i])
        single_lon_cs_depths.append(ifile_dist_grav)
        
    counter_cs = counter_cs + 1

#%%           

index_cs =  123.0 < np.array(ifile_lons) < 123.05



    
