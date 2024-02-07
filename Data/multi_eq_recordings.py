#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 13:50:58 2024

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from shapely.geometry.polygon import Polygon
import pygmt

#%%

st_file_f_salem = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/salem_eq/metadata/domain_stations.csv')
st_file_salem = st_file_f_salem.drop_duplicates()
st_names_salem = np.array(st_file_salem['st_name_dom'])
st_lons_salem = np.array(st_file_salem['st_lon_dom'])
st_lats_salem = np.array(st_file_salem['st_lat_dom']) 

st_file_f_scottsmills = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/snr_stations_fixed.csv')
st_file_scottsmills = st_file_f_scottsmills.drop_duplicates()
st_names_scottsmills = np.array(st_file_scottsmills['st_nm'])
st_lons_scottsmills = np.array(st_file_scottsmills['st_lon'])
st_lats_scottsmills = np.array(st_file_scottsmills['st_lat']) 

st_file_f_springfield = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/springfield_eq/metadata/snr_stations_fixed.csv')
st_file_springfield = st_file_f_springfield.drop_duplicates()
st_names_springfield = np.array(st_file_springfield['st_nm'])
st_lons_springfield = np.array(st_file_springfield['st_lon'])
st_lats_springfield = np.array(st_file_springfield['st_lat']) 

#%%

salem_scotml_mask = np.isin(st_names_salem,st_names_scottsmills)
salem_scotml_st = st_names_salem[salem_scotml_mask]
sal_scotm_spr_mask = np.isin(salem_scotml_st,st_names_springfield)
sal_scotm_spr_st = salem_scotml_st[sal_scotm_spr_mask]

all_eq_mask = np.isin(st_names_salem,sal_scotm_spr_st)
all_eq_st = st_file_salem[all_eq_mask]

ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/'

df_all_eq_st = pd.DataFrame(data=all_eq_st)
df_all_eq_st.to_csv(ff_dir+'all_eq_st.csv')

#%%

st_file_f_all_eq = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/all_eq_st.csv')
st_file_all_eq = st_file_f_all_eq.drop_duplicates()
st_names_all_eq = np.array(st_file_all_eq['st_name_dom'])
st_lons_all_eq = np.array(st_file_all_eq['st_lon_dom'])
st_lats_all_eq = np.array(st_file_all_eq['st_lat_dom']) 

#%%

events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/event_catalog_valevents.csv')

evlon = events['longitude']
evlat = events['latitude']
evcode = events['Code'] #B-ScottstMills, D-Springfield, E-Salem

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

outdir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/single_st_map/all_eq/'

# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)

for i in range(len(st_names_all_eq)):
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
        
        fig.plot(x=evlon[[1,3,4]] , y=evlat[[1,3,4]] , style='a0.8c',fill='yellow',pen='black')
        fig.text(text='Scotts Mills', x=evlon[1]+0.38, y=evlat[1]-0.03 , font='12p,Helvetica-Bold,firebrick4')
        fig.text(text='Springfield', x=evlon[3]+0.38, y=evlat[3]-0.03, font='12p,Helvetica-Bold,firebrick4')
        fig.text(text='Salem', x=evlon[4]+0.28, y=evlat[4]-0.03, font='12p,Helvetica-Bold,firebrick4')
        
        #basin outline
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        
        #stations
        fig.plot(x=st_lons_all_eq[i] , y=st_lats_all_eq[i] , style='t0.6c',fill='red',pen='black',transparency=30)
    
    # Save the plot
    fig.savefig(outdir+ st_names_all_eq[i] +'.png',dpi=200)







































