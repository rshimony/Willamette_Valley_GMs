#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 15:43:35 2024

@author: rshimony
"""

import pandas as pd
# import pygmt
import numpy as np
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

#%%

events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/event_catalog_valevents.csv')

evlon = events['longitude']
evlat = events['latitude']
evcode = events['Code'] #B-ScottstMills, D-Springfield, E-Salem

#%%

st_file_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/springfield_eq/metadata/valley_snr_st.csv')
st_file = st_file_f.drop_duplicates()
st_names = np.array(st_file['dom_valst_names'])
st_lons = np.array(st_file['dom_valst_lons'])
st_lats = np.array(st_file['dom_valst_lats']) 

st_file_f_full = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/salem_eq/metadata/domain_stations.csv')
st_file_full = st_file_f_full.drop_duplicates()
st_names_full = np.array(st_file_full['st_name_dom'])
st_lons_full = np.array(st_file_full['st_lon_dom'])
st_lats_full = np.array(st_file_full['st_lat_dom']) 

mask_notin = ~np.isin(st_names_full,st_names)
not_val_st = st_file_full[mask_notin]
st_names_not_val = np.array(not_val_st['st_nm'])
st_lons_not_val = np.array(not_val_st['st_lon'])
st_lats_not_val = np.array(not_val_st['st_lat']) 

#%%

salem_loc = Point(evlon[4],evlat[4])
salem_buffer = salem_loc.buffer(0.5, join_style=2)

st_in_buffer = []
for i in range(len(st_file_full)):
    st_loc = Point(st_lons_full[i],st_lats_full[i])
    if salem_buffer.contains(st_loc) == True:
        st_in_buffer.append(st_names_full[i]) 

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

xpn_m, ypn_m = full_poly.exterior.xy 
xph1_m, yph1_m = full_poly.interiors[0].xy
xph2_m, yph2_m = full_poly.interiors[1].xy

xpnn = np.array([xpn_m])
xph1n = np.array([xph1_m])
xph2n = np.array([xph2_m])

ypnn = np.array([ypn_m])
yph1n = np.array([yph1_m])
yph2n = np.array([yph2_m])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

#%%

def extract_ifile_basin_data(ifile_path):
    ifile_f = np.genfromtxt(ifile_path,skip_header=1)
    
    ifile_lons = ifile_f[:,0]
    ifile_lats = ifile_f[:,1]
    ifile_depths = ifile_f[:,2]
    
    ifile_depths_basin = ifile_depths[ifile_depths>0]
    
    ifile_lons_basin = ifile_lons[ifile_depths>0]

    ifile_lats_basin = ifile_lats[ifile_depths>0]
    
    return ifile_lons_basin , ifile_lats_basin , ifile_depths_basin

ifile_eugspe = '/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_latlon_eugspe.txt'
ifile_srv = '/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_latlon_srv.txt'

eugspe_lons , eugspe_lats , eugspe_depths = extract_ifile_basin_data(ifile_eugspe)
srv_lons , srv_lats , srv_depths = extract_ifile_basin_data(ifile_srv)

#%%

outdir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/single_st_map/springfield/'

# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)

for i in range(len(st_names_not_val)):
    with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
        fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
        fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
        
        # make color pallets
        pygmt.makecpt(
            cmap='abyss',reverse=True,
            series=[0,5000])
        
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        fig.plot(x=eugspe_lons, y=eugspe_lats,
                      style="c0.01c", 
                      fill=eugspe_depths,
                      cmap=True,
                      transparency=50)
        
        fig.coast(water="whitesmoke")
        fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")
        
        fig.plot(x=evlon[3] , y=evlat[3] , style='a0.8c',fill='yellow',pen='black')
        
        #basin outline
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        
        #stations
        fig.plot(x=st_lons_not_val[i] , y=st_lats_not_val[i] , style='t0.6c',fill='red',pen='black',transparency=30)
    
    # Save the plot
    fig.savefig(outdir+ st_names_not_val[i] +'.png',dpi=200)