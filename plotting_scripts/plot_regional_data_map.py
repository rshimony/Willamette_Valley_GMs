#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:11:59 2024

@author: rshimony
"""
'''
Plotting a regional map with data used in the WV GMs work.
'''
import pandas as pd
import pygmt
import numpy as np
from shapely.geometry.polygon import Polygon

#%%

## Reading in events catalog and extracting event location
events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/event_catalog_valevents.csv')
## idx   EQ ##
##  0    Scotts Mills ##
##  1    Salem ##
##  2    Springfield ##
evlon = events['longitude']
evlat = events['latitude']
evcode = events['Code']

#%% 

## Reading in all events stations diretories
# Salem
st_file_salem = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/salem_eq/metadata/snr_stations.csv').drop_duplicates()
st_names_salem = np.array(st_file_salem['st_nm'])
st_lons_salem = np.array(st_file_salem['st_lon'])
st_lats_salem = np.array(st_file_salem['st_lat']) 
#Scotts Mills
st_file_scottsmills = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/scottsmills_eq/metadata/snr_stations.csv').drop_duplicates()
st_names_scottsmills = np.array(st_file_scottsmills['st_nm'])
st_lons_scottsmills = np.array(st_file_scottsmills['st_lon'])
st_lats_scottsmills = np.array(st_file_scottsmills['st_lat']) 
# Springfield
st_file_springfield = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/metadata/snr_stations.csv').drop_duplicates()
st_names_springfield = np.array(st_file_springfield['st_nm'])
st_lons_springfield = np.array(st_file_springfield['st_lon'])
st_lats_springfield = np.array(st_file_springfield['st_lat']) 

all_eq_st_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/all_eq_st.csv')
st_names_all_eq = np.array(all_eq_st_f['st_nm'])
st_lons_all_eq = np.array(all_eq_st_f['st_lon'])
st_lats_all_eq = np.array(all_eq_st_f['st_lat']) 

#%%
## Extracting wells data from flatfile
depth2units_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wells_data/depth2units_wells.csv')
site_name = depth2units_f['SiteName']
lon_units = depth2units_f['Longitude']
lat_units = depth2units_f['Latitude']
# Getting locations of DIGITIZED wells
digi_wells_names = ['Ira Baker','Porter','M & P Farms','HJ Miller ','Wolverton','Independence','Bruer ','Klohs ']
digi_mask = np.isin(site_name,digi_wells_names)
digi_lons = lon_units[digi_mask]
digi_lats = lat_units[digi_mask]

#%%

## Creating WV basin polygon
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

def extract_ifile_basin_data(ifile_path):
    '''
    Extract WV basin depth and grid points location 
    from ifile. 
    Outputs depths, lats and lons arrays
    '''
    
    ifile_f = np.genfromtxt(ifile_path,skip_header=1)
    
    ifile_lons = ifile_f[:,0]
    ifile_lats = ifile_f[:,1]
    ifile_depths = ifile_f[:,2]
    
    ifile_depths_basin = ifile_depths[ifile_depths>0]
    
    ifile_lons_basin = ifile_lons[ifile_depths>0]

    ifile_lats_basin = ifile_lats[ifile_depths>0]
    
    return ifile_lons_basin , ifile_lats_basin , ifile_depths_basin

ifile_eugspe = '/Users/rshimony/Desktop/WillametteValley/wv_project/ifiles/ifile_latlon_eugspe.txt'
ifile_srv = '/Users/rshimony/Desktop/WillametteValley/wv_project/ifiles/ifile_latlon_srv.txt'

eugspe_lons , eugspe_lats , eugspe_depths = extract_ifile_basin_data(ifile_eugspe)
srv_lons , srv_lats , srv_depths = extract_ifile_basin_data(ifile_srv)

#%%
## Creating outlines polygon of gravity data from McPhee et al.
# Data file
grav_data_file = '/Users/rshimony/Desktop/WillametteValley/wv_project/grav_interp_depths/zzthick4.axyz'
# Reading data file and extracting unique locations
grav_data = np.loadtxt(grav_data_file)
g_lons, g_lats = grav_data[:, 0].astype('float64'), grav_data[:, 1].astype('float64')
uniq_g_lat = np.unique(g_lats)
uniq_g_lon = np.unique(g_lons)
# Creating the polygon and extracting its outlines
g_lat_min = min(uniq_g_lat)
g_lat_max = max(uniq_g_lat)
g_lon_min = min(uniq_g_lon)
g_lon_max = max(uniq_g_lon)

grav_dom = Polygon([[g_lon_min,g_lat_min] , [g_lon_max,g_lat_min] , [g_lon_max,g_lat_max] , [g_lon_min,g_lat_max]])
xpg, ypg = grav_dom.exterior.xy

#%%
## Ploting regional data map
# topo
topo_data = '@earth_relief_03s'

fig = pygmt.Figure()
region=[-125.4,-121.4,43.4,46.3]
with pygmt.config(FONT_ANNOT_PRIMARY=14):
    fig.basemap(region=region, projection='M0/0/10c',rose="JTR+o-1c+w0.7c+f1+l")
    pygmt.makecpt(
        cmap="gray",
        series = (-14000,5000))
    fig.grdimage(topo_data, cmap=True,shading=True,frame=True)
fig.coast(water="whitesmoke")

#basin outline
fig.plot(x=xpn_m, y=ypn_m,pen="0.6p,black", label="WV Basin Model Outline")
fig.plot(x=xph1_m, y=yph1_m,pen="0.6p,black")
fig.plot(x=xph2_m, y=yph2_m,pen="0.6p,black")

#McPhee outline
fig.plot(x=xpg, y=ypg,pen="1p,255/127/0,-", label="McPhee et al. Data Region")

#cities location and names
fig.plot(x=[-123.37, -123.1],y=[44.0, 44.05],pen="1p,25/25/112")
fig.plot(x=[-123.42, -123.0],y=[44.95, 44.95],pen="1p,25/25/112")
fig.plot(x=[-122.4, -122.76],y=[45.4, 45.45],pen="1p,25/25/112")

fig.text(text='Eugene', x=-123.59, y=44.0 , font='8p,Bookman-Demi,25/25/112')
fig.text(text='Salem', x=-123.6, y=44.95 , font='8p,Bookman-Demi,25/25/112')
fig.text(text='Portland', x=-122.15, y=45.4 , font='8p,Bookman-Demi,25/25/112')

#sub-basins location names
fig.plot(x=[-123.45, -123.1],y=[44.6, 44.5],pen="1p,25/25/112")
fig.plot(x=[-122.4, -122.9],y=[45.1, 45.05],pen="1p,25/25/112")
fig.plot(x=[-123.23, -122.85],y=[45.73, 45.52],pen="1p,25/25/112")

fig.text(text='SWV', x=-123.58, y=44.6 , font='8p,Bookman-Demi,25/25/112')
fig.text(text='NWV', x=-122.25, y=45.08 , font='8p,Bookman-Demi,25/25/112')
fig.text(text='TB', x=-123.27, y=45.78 , font='8p,Bookman-Demi,25/25/112')

#stations
fig.plot(x=st_lons_salem , y=st_lats_salem , style='t0.2c',fill='77/175/74',pen='black',transparency=30, label="Stations")
fig.plot(x=st_lons_scottsmills , y=st_lats_scottsmills , style='t0.2c',fill='77/175/74',pen='black',transparency=30)
fig.plot(x=st_lons_springfield , y=st_lats_springfield , style='t0.2c',fill='77/175/74',pen='black',transparency=30)
# All EQ stations
fig.plot(x=st_lons_all_eq , y=st_lats_all_eq , style='t0.2c',fill='255/255/51',pen='black',transparency=30, label="All EQ Stations")

#wells
fig.plot(x=lon_units , y=lat_units , style='c0.2c',fill='152/78/163',pen='black',transparency=30, label="Wells")
fig.plot(x=digi_lons , y=digi_lats , style='h0.2c',fill='55/126/184',pen='black',transparency=30, label="Digitized Wells")

#events location
fig.plot(x=evlon , y=evlat , style='a0.3c',fill='228/26/28',pen='black', label="EQ Epicenter")
fig.text(text=evcode, x=evlon+0.08, y=evlat-0.01 , font='8p,Helvetica-Bold,brown')


# Inset map
with fig.inset(position="jBL+o0.01c",
    box="+gwhite+p1p",
    region=[-126.5,-119.2,40,48.2],
    projection="M3c",
):

    fig.coast(
        land="gray",
        borders=[1, 2],
        shorelines="1/thin",
        water="white",
        area_thresh=10000,
    )
    
    fig.plot(x=xpn_m, y=ypn_m,pen="0.2p,red")
    fig.plot(x=xph1_m, y=yph1_m,pen="0.2,red")
    fig.plot(x=xph2_m, y=yph2_m,pen="0.2,red")
    
    rectangle = [[region[0], region[2], region[1], region[3]+0.05]]
    fig.plot(data=rectangle, style="r+s", pen="2p,blue")

#CS lines
fig.plot(x=[-124.19, -121.51],y=[st_lats_salem[20]+0.03, st_lats_salem[20]+0.03],pen="1p,red,-.",transparency=30, label="Cross Section")
fig.plot(x=[-124.19, -121.51],y=[st_lats_salem[63]-0.09, st_lats_salem[63]-0.09],pen="1p,red,-.",transparency=30)
fig.plot(x=[-124.19, -121.51],y=[st_lats_salem[88]+0.03, st_lats_salem[88]+0.03],pen="1p,red,-.",transparency=30)

fig.text(text="A'", x=-121.53, y=st_lats_salem[88]+0.08 , font='10p,Helvetica-Bold,royalblue4')
fig.text(text="B'", x=-121.53, y=st_lats_salem[63]-0.035 , font='10p,Helvetica-Bold,royalblue4')
fig.text(text="C'", x=-121.53, y=st_lats_salem[20]+0.08 , font='10p,Helvetica-Bold,royalblue4')

fig.text(text='A', x=-124.17, y=st_lats_salem[88]+0.08 , font='10p,Helvetica-Bold,royalblue4')
fig.text(text='B', x=-124.17, y=st_lats_salem[63]-0.035 , font='10p,Helvetica-Bold,royalblue4')
fig.text(text='C', x=-124.17, y=st_lats_salem[20]+0.08 , font='10p,Helvetica-Bold,royalblue4')

# North arrow and scale
fig.basemap(rose="JTR+o-1.4c+w0.9c+f1+l",map_scale="JBR+o-2.9/-0.6+w50")

# Legend
with pygmt.config(FONT_ANNOT_PRIMARY="8p"):
    fig.legend(position="jTL")
    
fig.savefig('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/region_data_map.png',dpi=800)

















