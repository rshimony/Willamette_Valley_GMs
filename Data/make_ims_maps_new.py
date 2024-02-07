#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 12:07:07 2023

@author: rshimony
"""

import obspy as obs
import pandas as pd
import pygmt
import numpy as np
import glob
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point,LineString
from shapely.ops import split
from obspy.core.utcdatetime import UTCDateTime
import matplotlib.pyplot as plt
import geopandas as gpd
import os

#%%
events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/event_catalog_valevents.csv')

evlon = events['longitude']
evlat = events['latitude']
evcode = events['Code']

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

def plot_shade_basin_map(basin_model_lon,basin_model_lat,basin_model_depth,fig_dir,fig_name,res_st_lons,res_st_lats,
                          plot_res=False ,plot_res_ratio=False):
    # # Plot the map using pyGMT
    plot_region = [-124.19,-121.51,43.4,46.1]
    
    # Create the first map
    
    fig = pygmt.Figure()
    
    ## topo:
    topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)
    
    with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
        fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
        fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
        # fig.coast(water='skyblue3')
        
        # make color pallets
        pygmt.makecpt(
            cmap='abyss',reverse=True,
            series=[0,5000])
        
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        fig.plot(x=basin_model_lon, y=basin_model_lat,
                      style="c0.01c", 
                      fill=basin_model_depth,
                      cmap=True,
                      transparency=50)
        
        fig.coast(water="whitesmoke")
        fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")
        
        fig.plot(x=evlon[1] , y=evlat[1] , style='a0.8c',fill='yellow',pen='black')
        
        if plot_res is not False:
            pygmt.makecpt(
                cmap='bam',
                series = (-2,2))
            fig.plot(x=res_st_lons, y=res_st_lats, 
                      style="t0.55c", 
                      fill=plot_res,
                      pen="black",
                      cmap=True,
                      transparency=0)
            fig.colorbar(position="JMR+o0.5c/0c+w12c/0.8c",frame=["af+lResiduals (Observed - Predicted)"])
            
        if plot_res_ratio is not False:
            pygmt.makecpt(
                cmap='bam',
                series = (0,2))
            fig.plot(x=res_st_lons, y=res_st_lats, 
                      style="t0.55c", 
                      fill=plot_res_ratio,
                      pen="black",
                      cmap=True,
                      transparency=0)
            fig.colorbar(position="JMR+o0.5c/0c+w12c/0.8c",frame=["af+lResiduals Ratio (Stephenson/Model)"])    
            
    
    # fig.show()
    fig.savefig(fig_dir+fig_name+'.png',dpi=300)

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

fas_obs_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/fas_ff_obs_scottsmills.csv')
fas_steph_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/fas_ff_steph_scottsmills.csv')
fas_val_obs_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/val_fas_ff_obs_scottsmills.csv')
fas_val_steph_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/val_fas_ff_steph_scottsmills.csv')

fas_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/fas_ff_eugspe_1d_scottsmills.csv')
fas_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/fas_ff_eugspe_sm_scottsmills.csv')
fas_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/fas_ff_srv_1d_scottsmills.csv')
fas_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/fas_ff_srv_sm_scottsmills.csv')

fas_val_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/val_fas_ff_eugspe_1d_scottsmills.csv')
fas_val_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/val_fas_ff_eugspe_sm_scottsmills.csv')
fas_val_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/val_fas_ff_srv_1d_scottsmills.csv')
fas_val_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/scottsmills/val_fas_ff_srv_sm_scottsmills.csv')

ims_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/eugspe_1d/ims_ff_eugspe_1d.csv')
ims_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/eugspe_sm/ims_ff_eugspe_sm.csv')
ims_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/srv_1d/ims_ff_srv_1d.csv')
ims_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/srv_sm/ims_ff_srv_sm.csv')

ims_val_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/eugspe_1d/val_ims_ff_eugspe_1d.csv')
ims_val_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/eugspe_sm/val_ims_ff_eugspe_sm.csv')
ims_val_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/srv_1d/val_ims_ff_srv_1d.csv')
ims_val_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/srv_sm/val_ims_ff_srv_sm.csv')
#%%

def find_vec_norm(Im_value_list):
    '''

    '''
    import numpy as np

    comp_1 = Im_value_list[0]**2;
    comp_2 = Im_value_list[1]**2;
    comp_3 = Im_value_list[2]**2;

    vec_norm = np.sqrt(comp_1 + comp_2 + comp_3)

    return(vec_norm)

st_lons = np.array(fas_val_eugspe_1d_file['st_lon'])
st_lats = np.array(fas_val_eugspe_1d_file['st_lat'])
st_name = np.array(fas_val_eugspe_1d_file['st_name'])
f_2 = np.array(fas_val_eugspe_1d_file['f=0.2'])
f_3 = np.array(fas_val_eugspe_1d_file['f=0.3'])
f_4 = np.array(fas_val_eugspe_1d_file['f=0.4'])
f_5 = np.array(fas_val_eugspe_1d_file['f=0.5'])
f_6 = np.array(fas_val_eugspe_1d_file['f=0.6'])
f_7 = np.array(fas_val_eugspe_1d_file['f=0.7'])

f_2_3 = []
f_3_3 = []
f_4_3 = []
f_5_3 = []
f_6_3 = []
f_7_3 = []
st_names_unq = pd.unique(st_name)
st_lons_unq = pd.unique(st_lons)
st_lats_unq = pd.unique(st_lats)

for i in range(len(st_names_unq)):
    same_st = st_name == st_names_unq[i]
    same_st_idx = np.where(same_st)[0]
    f_2_3.append(f_2[same_st_idx])
    f_3_3.append(f_3[same_st_idx])
    f_4_3.append(f_4[same_st_idx])
    f_5_3.append(f_5[same_st_idx])
    f_6_3.append(f_6[same_st_idx])
    f_7_3.append(f_7[same_st_idx])
    

f_2_t = []
f_3_t = []
f_4_t = []
f_5_t = []
f_6_t = []
f_7_t = []
for j in range(len(f_2_3)):
    st_f_2_t = find_vec_norm(f_2_3[j])
    f_2_t.append(st_f_2_t)
    
    st_f_3_t = find_vec_norm(f_3_3[j])
    f_3_t.append(st_f_3_t)
    
    st_f_4_t = find_vec_norm(f_4_3[j])
    f_4_t.append(st_f_4_t)
    
    st_f_5_t = find_vec_norm(f_5_3[j])
    f_5_t.append(st_f_5_t)
    
    st_f_6_t = find_vec_norm(f_6_3[j])
    f_6_t.append(st_f_6_t)
    
    st_f_7_t = find_vec_norm(f_7_3[j])
    f_7_t.append(st_f_7_t)
    

#%%

def find_vec_norm(Im_value_list):
    '''

    '''
    import numpy as np

    comp_1 = Im_value_list[0]**2;
    comp_2 = Im_value_list[1]**2;
    comp_3 = Im_value_list[2]**2;

    vec_norm = np.sqrt(comp_1 + comp_2 + comp_3)

    return(vec_norm)

def extract_spec_ims_props(spec_ims_file):
    
    st_lons = np.array(spec_ims_file['st_lon'])
    st_lats = np.array(spec_ims_file['st_lat'])
    st_name = np.array(spec_ims_file['st_name'])
    f_2 = np.array(spec_ims_file['f=0.2'])
    f_3 = np.array(spec_ims_file['f=0.3'])
    f_4 = np.array(spec_ims_file['f=0.4'])
    f_5 = np.array(spec_ims_file['f=0.5'])
    f_6 = np.array(spec_ims_file['f=0.6'])
    f_7 = np.array(spec_ims_file['f=0.7'])
    
    f_2_3 = []
    f_3_3 = []
    f_4_3 = []
    f_5_3 = []
    f_6_3 = []
    f_7_3 = []
    st_names_unq = pd.unique(st_name)
    st_lons_unq = pd.unique(st_lons)
    st_lats_unq = pd.unique(st_lats)
    
    for i in range(len(st_names_unq)):
        same_st = st_name == st_names_unq[i]
        same_st_idx = np.where(same_st)[0]
        f_2_3.append(f_2[same_st_idx])
        f_3_3.append(f_3[same_st_idx])
        f_4_3.append(f_4[same_st_idx])
        f_5_3.append(f_5[same_st_idx])
        f_6_3.append(f_6[same_st_idx])
        f_7_3.append(f_7[same_st_idx])
        
    
    f_2_t = []
    f_3_t = []
    f_4_t = []
    f_5_t = []
    f_6_t = []
    f_7_t = []
    for j in range(len(f_2_3)):
        st_f_2_t = find_vec_norm(f_2_3[j])
        f_2_t.append(st_f_2_t)
        st_f_3_t = find_vec_norm(f_3_3[j])
        f_3_t.append(st_f_3_t)
        st_f_4_t = find_vec_norm(f_4_3[j])
        f_4_t.append(st_f_4_t)
        st_f_5_t = find_vec_norm(f_5_3[j])
        f_5_t.append(st_f_5_t)
        st_f_6_t = find_vec_norm(f_6_3[j])
        f_6_t.append(st_f_6_t)
        st_f_7_t = find_vec_norm(f_7_3[j])
        f_7_t.append(st_f_7_t)
        
    return st_names_unq,st_lons_unq,st_lats_unq,f_2_t,f_3_t,f_4_t,f_5_t,f_6_t,f_7_t
    
def extract_ims_pgv(ims_file):
    
    st_lons = np.array(ims_file['st_lons'])
    st_lats = np.array(ims_file['st_lats'])
    st_name = np.array(ims_file['st_names'])
    pgv_t_synt = np.array(ims_file['pgv_t_synt'])
    pgv_t_obs = np.array(ims_file['pgv_t_obs'])
    pgv_t_steph = np.array(ims_file['pgv_t_steph'])
    
    return st_name,st_lons,st_lats,pgv_t_synt,pgv_t_obs,pgv_t_steph    


'''
FAS:
 0     1     2      3      4      5       6     7     8
name  lons  lats  f=0.2  f=0.3  f=0.4  f=0.5  f=0.6  f=0.7

PGV:
 0     1     2       3         4         5    
name  lons  lats  pgv_synt  pgv_obs  pgv_steph
'''

fas_obs = extract_spec_ims_props(fas_obs_file)
fas_val_obs = extract_spec_ims_props(fas_val_obs_file)
fas_steph = extract_spec_ims_props(fas_steph_file)
fas_val_steph = extract_spec_ims_props(fas_val_steph_file)

fas_eugspe_1d = extract_spec_ims_props(fas_eugspe_1d_file) 
fas_eugspe_sm = extract_spec_ims_props(fas_eugspe_sm_file)
fas_srv_1d = extract_spec_ims_props(fas_srv_1d_file)
fas_srv_sm = extract_spec_ims_props(fas_srv_sm_file)

fas_val_eugspe_1d = extract_spec_ims_props(fas_val_eugspe_1d_file) 
fas_val_eugspe_sm = extract_spec_ims_props(fas_val_eugspe_sm_file)
fas_val_srv_1d = extract_spec_ims_props(fas_val_srv_1d_file)
fas_val_srv_sm = extract_spec_ims_props(fas_val_srv_sm_file)

pgv_eugspe_1d = extract_ims_pgv(ims_eugspe_1d_file)
pgv_eugspe_sm = extract_ims_pgv(ims_eugspe_sm_file)
pgv_srv_1d = extract_ims_pgv(ims_srv_1d_file)
pgv_srv_sm = extract_ims_pgv(ims_srv_sm_file)

pgv_val_eugspe_1d = extract_ims_pgv(ims_val_eugspe_1d_file)
pgv_val_eugspe_sm = extract_ims_pgv(ims_val_eugspe_sm_file)
pgv_val_srv_1d = extract_ims_pgv(ims_val_srv_1d_file)
pgv_val_srv_sm = extract_ims_pgv(ims_val_srv_sm_file)

#%%

# diff_syn = np.array(t2_psa_obs) - np.array(t2_psa_eugspe_1d)
# diff_steph = np.array(t2_psa_obs) - np.array(t2_psa_steph)


#%%

f2_fas_steph_res = np.log(fas_obs[3]) - np.log(fas_steph[3])
f3_fas_steph_res = np.log(fas_obs[4]) - np.log(fas_steph[4])
f4_fas_steph_res = np.log(fas_obs[5]) - np.log(fas_steph[5])
f5_fas_steph_res = np.log(fas_obs[6]) - np.log(fas_steph[6])
f6_fas_steph_res = np.log(fas_obs[7]) - np.log(fas_steph[7])
f7_fas_steph_res = np.log(fas_obs[8]) - np.log(fas_steph[8])

f2_fas_val_steph_res = np.log(fas_val_obs[3]) - np.log(fas_val_steph[3])
f3_fas_val_steph_res = np.log(fas_val_obs[4]) - np.log(fas_val_steph[4])
f4_fas_val_steph_res = np.log(fas_val_obs[5]) - np.log(fas_val_steph[5])
f5_fas_val_steph_res = np.log(fas_val_obs[6]) - np.log(fas_val_steph[6])
f6_fas_val_steph_res = np.log(fas_val_obs[7]) - np.log(fas_val_steph[7])
f7_fas_val_steph_res = np.log(fas_val_obs[8]) - np.log(fas_val_steph[8])


f2_fas_eugspe_1d_res = np.log(fas_obs[3]) - np.log(fas_eugspe_1d[3])
f3_fas_eugspe_1d_res = np.log(fas_obs[4]) - np.log(fas_eugspe_1d[4])
f4_fas_eugspe_1d_res = np.log(fas_obs[5]) - np.log(fas_eugspe_1d[5])
f5_fas_eugspe_1d_res = np.log(fas_obs[6]) - np.log(fas_eugspe_1d[6])
f6_fas_eugspe_1d_res = np.log(fas_obs[7]) - np.log(fas_eugspe_1d[7])
f7_fas_eugspe_1d_res = np.log(fas_obs[8]) - np.log(fas_eugspe_1d[8])

f2_fas_val_eugspe_1d_res = np.log(fas_val_obs[3]) - np.log(fas_val_eugspe_1d[3])
f3_fas_val_eugspe_1d_res = np.log(fas_val_obs[4]) - np.log(fas_val_eugspe_1d[4])
f4_fas_val_eugspe_1d_res = np.log(fas_val_obs[5]) - np.log(fas_val_eugspe_1d[5])
f5_fas_val_eugspe_1d_res = np.log(fas_val_obs[6]) - np.log(fas_val_eugspe_1d[6])
f6_fas_val_eugspe_1d_res = np.log(fas_val_obs[7]) - np.log(fas_val_eugspe_1d[7])
f7_fas_val_eugspe_1d_res = np.log(fas_val_obs[8]) - np.log(fas_val_eugspe_1d[8])


f2_fas_eugspe_sm_res = np.log(fas_obs[3]) - np.log(fas_eugspe_sm[3])
f3_fas_eugspe_sm_res = np.log(fas_obs[4]) - np.log(fas_eugspe_sm[4])
f4_fas_eugspe_sm_res = np.log(fas_obs[5]) - np.log(fas_eugspe_sm[5])
f5_fas_eugspe_sm_res = np.log(fas_obs[6]) - np.log(fas_eugspe_sm[6])
f6_fas_eugspe_sm_res = np.log(fas_obs[7]) - np.log(fas_eugspe_sm[7])
f7_fas_eugspe_sm_res = np.log(fas_obs[8]) - np.log(fas_eugspe_sm[8])

f2_fas_val_eugspe_sm_res = np.log(fas_val_obs[3]) - np.log(fas_val_eugspe_sm[3])
f3_fas_val_eugspe_sm_res = np.log(fas_val_obs[4]) - np.log(fas_val_eugspe_sm[4])
f4_fas_val_eugspe_sm_res = np.log(fas_val_obs[5]) - np.log(fas_val_eugspe_sm[5])
f5_fas_val_eugspe_sm_res = np.log(fas_val_obs[6]) - np.log(fas_val_eugspe_sm[6])
f6_fas_val_eugspe_sm_res = np.log(fas_val_obs[7]) - np.log(fas_val_eugspe_sm[7])
f7_fas_val_eugspe_sm_res = np.log(fas_val_obs[8]) - np.log(fas_val_eugspe_sm[8])


f2_fas_srv_1d_res = np.log(fas_obs[3]) - np.log(fas_srv_1d[3])
f3_fas_srv_1d_res = np.log(fas_obs[4]) - np.log(fas_srv_1d[4])
f4_fas_srv_1d_res = np.log(fas_obs[5]) - np.log(fas_srv_1d[5])
f5_fas_srv_1d_res = np.log(fas_obs[6]) - np.log(fas_srv_1d[6])
f6_fas_srv_1d_res = np.log(fas_obs[7]) - np.log(fas_srv_1d[7])
f7_fas_srv_1d_res = np.log(fas_obs[8]) - np.log(fas_srv_1d[8])

f2_fas_val_srv_1d_res = np.log(fas_val_obs[3]) - np.log(fas_val_srv_1d[3])
f3_fas_val_srv_1d_res = np.log(fas_val_obs[4]) - np.log(fas_val_srv_1d[4])
f4_fas_val_srv_1d_res = np.log(fas_val_obs[5]) - np.log(fas_val_srv_1d[5])
f5_fas_val_srv_1d_res = np.log(fas_val_obs[6]) - np.log(fas_val_srv_1d[6])
f6_fas_val_srv_1d_res = np.log(fas_val_obs[7]) - np.log(fas_val_srv_1d[7])
f7_fas_val_srv_1d_res = np.log(fas_val_obs[8]) - np.log(fas_val_srv_1d[8])


f2_fas_srv_sm_res = np.log(fas_obs[3]) - np.log(fas_srv_sm[3])
f3_fas_srv_sm_res = np.log(fas_obs[4]) - np.log(fas_srv_sm[4])
f4_fas_srv_sm_res = np.log(fas_obs[5]) - np.log(fas_srv_sm[5])
f5_fas_srv_sm_res = np.log(fas_obs[6]) - np.log(fas_srv_sm[6])
f6_fas_srv_sm_res = np.log(fas_obs[7]) - np.log(fas_srv_sm[7])
f7_fas_srv_sm_res = np.log(fas_obs[8]) - np.log(fas_srv_sm[8])

f2_fas_val_srv_sm_res = np.log(fas_val_obs[3]) - np.log(fas_val_srv_sm[3])
f3_fas_val_srv_sm_res = np.log(fas_val_obs[4]) - np.log(fas_val_srv_sm[4])
f4_fas_val_srv_sm_res = np.log(fas_val_obs[5]) - np.log(fas_val_srv_sm[5])
f5_fas_val_srv_sm_res = np.log(fas_val_obs[6]) - np.log(fas_val_srv_sm[6])
f6_fas_val_srv_sm_res = np.log(fas_val_obs[7]) - np.log(fas_val_srv_sm[7])
f7_fas_val_srv_sm_res = np.log(fas_val_obs[8]) - np.log(fas_val_srv_sm[8])


pgv_steph_res = np.log(pgv_eugspe_1d[4]) - np.log(pgv_eugspe_1d[5])
pgv_steph_res_val = np.log(pgv_val_eugspe_1d[4]) - np.log(pgv_val_eugspe_1d[5])

pgv_eugspe_1d_res = np.log(pgv_eugspe_1d[4]) - np.log(pgv_eugspe_1d[3])
pgv_eugspe_sm_res = np.log(pgv_eugspe_sm[4]) - np.log(pgv_eugspe_sm[3])
pgv_srv_1d_res = np.log(pgv_srv_1d[4]) - np.log(pgv_srv_1d[3])
pgv_srv_sm_res = np.log(pgv_srv_sm[4]) - np.log(pgv_srv_sm[3])

pgv_eugspe_1d_res_val = np.log(pgv_val_eugspe_1d[4]) - np.log(pgv_val_eugspe_1d[3])
pgv_eugspe_sm_res_val = np.log(pgv_val_eugspe_sm[4]) - np.log(pgv_val_eugspe_sm[3])
pgv_srv_1d_res_val = np.log(pgv_val_srv_1d[4]) - np.log(pgv_val_srv_1d[3])
pgv_srv_sm_res_val = np.log(pgv_val_srv_sm[4]) - np.log(pgv_val_srv_sm[3])


f2_fas_eugspe_1d_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_eugspe_1d_res)
f3_fas_eugspe_1d_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_eugspe_1d_res)
f4_fas_eugspe_1d_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_eugspe_1d_res)
f5_fas_eugspe_1d_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_eugspe_1d_res)
f6_fas_eugspe_1d_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_eugspe_1d_res)
f7_fas_eugspe_1d_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_eugspe_1d_res)

f2_fas_eugspe_1d_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_eugspe_1d_res)
f3_fas_eugspe_1d_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_eugspe_1d_res)
f4_fas_eugspe_1d_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_eugspe_1d_res)
f5_fas_eugspe_1d_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_eugspe_1d_res)
f6_fas_eugspe_1d_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_eugspe_1d_res)
f7_fas_eugspe_1d_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_eugspe_1d_res)

f2_fas_eugspe_sm_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_eugspe_sm_res)
f3_fas_eugspe_sm_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_eugspe_sm_res)
f4_fas_eugspe_sm_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_eugspe_sm_res)
f5_fas_eugspe_sm_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_eugspe_sm_res)
f6_fas_eugspe_sm_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_eugspe_sm_res)
f7_fas_eugspe_sm_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_eugspe_sm_res)

f2_fas_eugspe_sm_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_eugspe_sm_res)
f3_fas_eugspe_sm_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_eugspe_sm_res)
f4_fas_eugspe_sm_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_eugspe_sm_res)
f5_fas_eugspe_sm_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_eugspe_sm_res)
f6_fas_eugspe_sm_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_eugspe_sm_res)
f7_fas_eugspe_sm_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_eugspe_sm_res)

f2_fas_srv_1d_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_srv_1d_res)
f3_fas_srv_1d_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_srv_1d_res)
f4_fas_srv_1d_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_srv_1d_res)
f5_fas_srv_1d_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_srv_1d_res)
f6_fas_srv_1d_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_srv_1d_res)
f7_fas_srv_1d_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_srv_1d_res)

f2_fas_srv_1d_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_srv_1d_res)
f3_fas_srv_1d_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_srv_1d_res)
f4_fas_srv_1d_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_srv_1d_res)
f5_fas_srv_1d_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_srv_1d_res)
f6_fas_srv_1d_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_srv_1d_res)
f7_fas_srv_1d_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_srv_1d_res)

f2_fas_srv_sm_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_srv_sm_res)
f3_fas_srv_sm_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_srv_sm_res)
f4_fas_srv_sm_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_srv_sm_res)
f5_fas_srv_sm_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_srv_sm_res)
f6_fas_srv_sm_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_srv_sm_res)
f7_fas_srv_sm_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_srv_sm_res)

f2_fas_srv_sm_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_srv_sm_res)
f3_fas_srv_sm_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_srv_sm_res)
f4_fas_srv_sm_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_srv_sm_res)
f5_fas_srv_sm_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_srv_sm_res)
f6_fas_srv_sm_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_srv_sm_res)
f7_fas_srv_sm_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_srv_sm_res)


pgv_eugspe_1d_ratio = np.abs(pgv_steph_res)/np.abs(pgv_eugspe_1d_res)
pgv_eugspe_sm_ratio = np.abs(pgv_steph_res)/np.abs(pgv_eugspe_sm_res)
pgv_srv_1d_ratio = np.abs(pgv_steph_res)/np.abs(pgv_srv_1d_res)
pgv_srv_sm_ratio = np.abs(pgv_steph_res)/np.abs(pgv_srv_sm_res)

pgv_eugspe_1d_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_eugspe_1d_res_val)
pgv_eugspe_sm_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_eugspe_sm_res_val)
pgv_srv_1d_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_srv_1d_res_val)
pgv_srv_sm_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_srv_sm_res_val)


grav_lat_idx = np.argwhere(fas_obs[2] > 44.98)[:,0]
lat_region = np.zeros_like(fas_obs[2])
lat_region[grav_lat_idx] = 1

lat_region = np.where(lat_region == 1, 'North', 'South')

srv_sm_nm_ls = ["srv_sm"]*len(lat_region)
srv_1d_nm_ls = ["srv_1d"]*len(lat_region)
eugspe_sm_nm_ls = ["eugspe_sm"]*len(lat_region)
eugspe_1d_nm_ls = ["eugspe_1d"]*len(lat_region)
steph_nm_ls = ["steph"]*len(lat_region)


grav_lat_idx_val = np.argwhere(fas_val_obs[2] > 44.98)[:,0]
lat_region_val = np.zeros_like(fas_val_obs[2])
lat_region_val[grav_lat_idx_val] = 1

lat_region_val = np.where(lat_region_val == 1, 'North', 'South')

srv_sm_nm_ls_val = ["srv_sm"]*len(lat_region_val)
srv_1d_nm_ls_val = ["srv_1d"]*len(lat_region_val)
eugspe_sm_nm_ls_val = ["eugspe_sm"]*len(lat_region_val)
eugspe_1d_nm_ls_val = ["eugspe_1d"]*len(lat_region_val)
steph_nm_ls_val = ["steph"]*len(lat_region_val)

#%%

res_map_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/scottsmills/'
try:
    os.mkdir(res_map_dir)
except FileExistsError:
    print('Directory already exists')
    
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f2_fas_steph_res' , fas_obs[1] , fas_obs[2] , plot_res=f2_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f3_fas_steph_res' , fas_obs[1] , fas_obs[2] , plot_res=f3_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f4_fas_steph_res' , fas_obs[1] , fas_obs[2] , plot_res=f4_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f5_fas_steph_res' , fas_obs[1] , fas_obs[2] , plot_res=f5_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f6_fas_steph_res' , fas_obs[1] , fas_obs[2] , plot_res=f6_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f7_fas_steph_res' , fas_obs[1] , fas_obs[2] , plot_res=f7_fas_steph_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f2_fas_val_steph_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f2_fas_val_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f3_fas_val_steph_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f3_fas_val_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f4_fas_val_steph_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f4_fas_val_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f5_fas_val_steph_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f5_fas_val_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f6_fas_val_steph_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f6_fas_val_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f7_fas_val_steph_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f7_fas_val_steph_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'pgv_steph_res' , fas_obs[1] , fas_obs[2] , plot_res=pgv_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'pgv_steph_res_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res=pgv_steph_res_val)

    
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_eugspe_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f2_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_eugspe_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f3_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_eugspe_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f4_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_eugspe_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f5_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_eugspe_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f6_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_eugspe_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f7_fas_eugspe_1d_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_val_eugspe_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f2_fas_val_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_val_eugspe_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f3_fas_val_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_val_eugspe_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f4_fas_val_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_val_eugspe_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f5_fas_val_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_val_eugspe_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f6_fas_val_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_val_eugspe_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f7_fas_val_eugspe_1d_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=pgv_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_1d_res_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res=pgv_eugspe_1d_res_val)


plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_eugspe_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f2_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_eugspe_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f3_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_eugspe_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f4_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_eugspe_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f5_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_eugspe_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f6_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_eugspe_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f7_fas_eugspe_sm_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_val_eugspe_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f2_fas_val_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_val_eugspe_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f3_fas_val_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_val_eugspe_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f4_fas_val_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_val_eugspe_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f5_fas_val_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_val_eugspe_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f6_fas_val_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_val_eugspe_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f7_fas_val_eugspe_sm_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=pgv_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_sm_res_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res=pgv_eugspe_sm_res_val)


plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_srv_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f2_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_srv_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f3_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_srv_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f4_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_srv_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f5_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_srv_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f6_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_srv_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=f7_fas_srv_1d_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_val_srv_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f2_fas_val_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_val_srv_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f3_fas_val_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_val_srv_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f4_fas_val_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_val_srv_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f5_fas_val_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_val_srv_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f6_fas_val_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_val_srv_1d_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f7_fas_val_srv_1d_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_1d_res' , fas_obs[1] , fas_obs[2] , plot_res=pgv_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_1d_res_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res=pgv_srv_1d_res_val)


plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_srv_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f2_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_srv_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f3_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_srv_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f4_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_srv_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f5_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_srv_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f6_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_srv_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=f7_fas_srv_sm_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_val_srv_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f2_fas_val_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_val_srv_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f3_fas_val_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_val_srv_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f4_fas_val_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_val_srv_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f5_fas_val_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_val_srv_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f6_fas_val_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_val_srv_sm_res' , fas_val_obs[1] , fas_val_obs[2] , plot_res=f7_fas_val_srv_sm_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_sm_res' , fas_obs[1] , fas_obs[2] , plot_res=pgv_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_sm_res_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res=pgv_srv_sm_res_val)

################

ratio_map_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/ratio_maps/scottsmills/'
try:
    os.mkdir(ratio_map_dir)
except FileExistsError:
    print('Directory already exists')
    
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f2_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f3_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f4_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f5_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f6_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f7_fas_eugspe_1d_ratio)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f2_fas_eugspe_1d_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f3_fas_eugspe_1d_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f4_fas_eugspe_1d_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f5_fas_eugspe_1d_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f6_fas_eugspe_1d_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f7_fas_eugspe_1d_ratio_val)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=pgv_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=pgv_eugspe_1d_ratio_val)


plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f2_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f3_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f4_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f5_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f6_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f7_fas_eugspe_sm_ratio)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f2_fas_eugspe_sm_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f3_fas_eugspe_sm_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f4_fas_eugspe_sm_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f5_fas_eugspe_sm_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f6_fas_eugspe_sm_ratio_val)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f7_fas_eugspe_sm_ratio_val)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=pgv_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=pgv_eugspe_sm_ratio_val)


plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f2_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f3_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f4_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f5_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f6_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f7_fas_srv_1d_ratio)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f2_fas_srv_1d_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f3_fas_srv_1d_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f4_fas_srv_1d_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f5_fas_srv_1d_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f6_fas_srv_1d_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f7_fas_srv_1d_ratio_val)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_1d_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=pgv_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_1d_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=pgv_srv_1d_ratio_val)


plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f2_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f3_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f4_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f5_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f6_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=f7_fas_srv_sm_ratio)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f2_fas_srv_sm_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f3_fas_srv_sm_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f4_fas_srv_sm_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f5_fas_srv_sm_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f6_fas_srv_sm_ratio_val)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=f7_fas_srv_sm_ratio_val)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_sm_ratio' , fas_obs[1] , fas_obs[2] , plot_res_ratio=pgv_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_sm_ratio_val' , fas_val_obs[1] , fas_val_obs[2] , plot_res_ratio=pgv_srv_sm_ratio_val)


#%%

data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'lat_region':lat_region,'model':steph_nm_ls,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res , 
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             #'f2_fas_val_steph_res':f2_fas_val_steph_res , 'f3_fas_val_steph_res':f3_fas_val_steph_res, 'f4_fas_val_steph_res':f4_fas_val_steph_res , 
             #'f5_fas_val_steph_res':f5_fas_val_steph_res , 'f6_fas_val_steph_res':f6_fas_val_steph_res , 'f7_fas_val_steph_res':f7_fas_val_steph_res ,
             
             'pgv_steph_res':pgv_steph_res} #,'pgv_steph_res_val':pgv_steph_res_val}

df_res = pd.DataFrame(data=data_res)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res.to_csv(resrat_ff_dir+'steph_res_ff.csv')
              
###################

data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'lat_region':lat_region,'model':eugspe_1d_nm_ls,
             'f2_fas_eugspe_1d_res':f2_fas_eugspe_1d_res , 'f3_fas_eugspe_1d_res':f3_fas_eugspe_1d_res, 'f4_fas_eugspe_1d_res':f4_fas_eugspe_1d_res , 
             'f5_fas_eugspe_1d_res':f5_fas_eugspe_1d_res , 'f6_fas_eugspe_1d_res':f6_fas_eugspe_1d_res , 'f7_fas_eugspe_1d_res':f7_fas_eugspe_1d_res ,
             'pgv_eugspe_1d_res':pgv_eugspe_1d_res,
             # 'f2_fas_val_eugspe_1d_res':f2_fas_val_eugspe_1d_res , 'f3_fas_val_eugspe_1d_res':f3_fas_val_eugspe_1d_res, 'f4_fas_val_eugspe_1d_res':f4_fas_val_eugspe_1d_res , 
             # 'f5_fas_val_eugspe_1d_res':f5_fas_val_eugspe_1d_res , 'f6_fas_val_eugspe_1d_res':f6_fas_val_eugspe_1d_res , 'f7_fas_val_eugspe_1d_res':f7_fas_val_eugspe_1d_res ,
             
             'f2_fas_eugspe_1d_ratio':f2_fas_eugspe_1d_ratio , 'f3_fas_eugspe_1d_ratio':f3_fas_eugspe_1d_ratio, 'f4_fas_eugspe_1d_ratio':f4_fas_eugspe_1d_ratio , 
             'f5_fas_eugspe_1d_ratio':f5_fas_eugspe_1d_ratio , 'f6_fas_eugspe_1d_ratio':f6_fas_eugspe_1d_ratio , 'f7_fas_eugspe_1d_ratio':f7_fas_eugspe_1d_ratio ,
             # 'f2_fas_eugspe_1d_ratio_val':f2_fas_eugspe_1d_ratio_val , 'f3_fas_eugspe_1d_ratio_val':f3_fas_eugspe_1d_ratio_val, 'f4_fas_eugspe_1d_ratio_val':f4_fas_eugspe_1d_ratio_val , 
             # 'f5_fas_eugspe_1d_ratio_val':f5_fas_eugspe_1d_ratio_val , 'f6_fas_eugspe_1d_ratio_val':f6_fas_eugspe_1d_ratio_val , 'f7_fas_eugspe_1d_ratio_val':f7_fas_eugspe_1d_ratio_val ,
             
              'pgv_eugspe_1d_ratio':pgv_eugspe_1d_ratio,
              'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
              'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
              'pgv_steph_res':pgv_steph_res} #, 'pgv_eugspe_1d_ratio_val':pgv_eugspe_1d_ratio_val}

df_res = pd.DataFrame(data=data_res)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res.to_csv(resrat_ff_dir+'eugspe_1d_res_ff.csv')

###################

data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'lat_region':lat_region,'model':eugspe_sm_nm_ls,
             'f2_fas_eugspe_sm_res':f2_fas_eugspe_sm_res , 'f3_fas_eugspe_sm_res':f3_fas_eugspe_sm_res, 'f4_fas_eugspe_sm_res':f4_fas_eugspe_sm_res , 
             'f5_fas_eugspe_sm_res':f5_fas_eugspe_sm_res , 'f6_fas_eugspe_sm_res':f6_fas_eugspe_sm_res , 'f7_fas_eugspe_sm_res':f7_fas_eugspe_sm_res ,
             'pgv_eugspe_sm_res':pgv_eugspe_sm_res,
             # 'f2_fas_val_eugspe_sm_res':f2_fas_val_eugspe_sm_res , 'f3_fas_val_eugspe_sm_res':f3_fas_val_eugspe_sm_res, 'f4_fas_val_eugspe_sm_res':f4_fas_val_eugspe_sm_res , 
             # 'f5_fas_val_eugspe_sm_res':f5_fas_val_eugspe_sm_res , 'f6_fas_val_eugspe_sm_res':f6_fas_val_eugspe_sm_res , 'f7_fas_val_eugspe_sm_res':f7_fas_val_eugspe_sm_res ,
             
             'f2_fas_eugspe_sm_ratio':f2_fas_eugspe_sm_ratio , 'f3_fas_eugspe_sm_ratio':f3_fas_eugspe_sm_ratio, 'f4_fas_eugspe_sm_ratio':f4_fas_eugspe_sm_ratio , 
             'f5_fas_eugspe_sm_ratio':f5_fas_eugspe_sm_ratio , 'f6_fas_eugspe_sm_ratio':f6_fas_eugspe_sm_ratio , 'f7_fas_eugspe_sm_ratio':f7_fas_eugspe_sm_ratio ,
             # 'f2_fas_eugspe_sm_ratio_val':f2_fas_eugspe_sm_ratio_val , 'f3_fas_eugspe_sm_ratio_val':f3_fas_eugspe_sm_ratio_val, 'f4_fas_eugspe_sm_ratio_val':f4_fas_eugspe_sm_ratio_val , 
             # 'f5_fas_eugspe_sm_ratio_val':f5_fas_eugspe_sm_ratio_val , 'f6_fas_eugspe_sm_ratio_val':f6_fas_eugspe_sm_ratio_val , 'f7_fas_eugspe_sm_ratio_val':f7_fas_eugspe_sm_ratio_val ,
             
             'pgv_eugspe_sm_ratio':pgv_eugspe_sm_ratio,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             'pgv_steph_res':pgv_steph_res} #, 'pgv_eugspe_sm_ratio_val':pgv_eugspe_sm_ratio_val}

df_res = pd.DataFrame(data=data_res)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res.to_csv(resrat_ff_dir+'eugspe_sm_res_ff.csv')

###################

data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'lat_region':lat_region,'model':srv_1d_nm_ls,
             'f2_fas_srv_1d_res':f2_fas_srv_1d_res , 'f3_fas_srv_1d_res':f3_fas_srv_1d_res, 'f4_fas_srv_1d_res':f4_fas_srv_1d_res , 
             'f5_fas_srv_1d_res':f5_fas_srv_1d_res , 'f6_fas_srv_1d_res':f6_fas_srv_1d_res , 'f7_fas_srv_1d_res':f7_fas_srv_1d_res ,
             'pgv_srv_1d_res':pgv_srv_1d_res,
             # 'f2_fas_val_srv_1d_res':f2_fas_val_srv_1d_res , 'f3_fas_val_srv_1d_res':f3_fas_val_srv_1d_res, 'f4_fas_val_srv_1d_res':f4_fas_val_srv_1d_res , 
             # 'f5_fas_val_srv_1d_res':f5_fas_val_srv_1d_res , 'f6_fas_val_srv_1d_res':f6_fas_val_srv_1d_res , 'f7_fas_val_srv_1d_res':f7_fas_val_srv_1d_res ,
             
             'f2_fas_srv_1d_ratio':f2_fas_srv_1d_ratio , 'f3_fas_srv_1d_ratio':f3_fas_srv_1d_ratio, 'f4_fas_srv_1d_ratio':f4_fas_srv_1d_ratio , 
             'f5_fas_srv_1d_ratio':f5_fas_srv_1d_ratio , 'f6_fas_srv_1d_ratio':f6_fas_srv_1d_ratio , 'f7_fas_srv_1d_ratio':f7_fas_srv_1d_ratio ,
             # 'f2_fas_srv_1d_ratio_val':f2_fas_srv_1d_ratio_val , 'f3_fas_srv_1d_ratio_val':f3_fas_srv_1d_ratio_val, 'f4_fas_srv_1d_ratio_val':f4_fas_srv_1d_ratio_val , 
             # 'f5_fas_srv_1d_ratio_val':f5_fas_srv_1d_ratio_val , 'f6_fas_srv_1d_ratio_val':f6_fas_srv_1d_ratio_val , 'f7_fas_srv_1d_ratio_val':f7_fas_srv_1d_ratio_val ,
             
             'pgv_srv_1d_ratio':pgv_srv_1d_ratio,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             'pgv_steph_res':pgv_steph_res} #, 'pgv_srv_1d_ratio_val':pgv_srv_1d_ratio_val}

df_res = pd.DataFrame(data=data_res)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res.to_csv(resrat_ff_dir+'srv_1d_res_ff.csv')

###################

data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'lat_region':lat_region,'model':srv_sm_nm_ls,
             'f2_fas_srv_sm_res':f2_fas_srv_sm_res , 'f3_fas_srv_sm_res':f3_fas_srv_sm_res, 'f4_fas_srv_sm_res':f4_fas_srv_sm_res , 
             'f5_fas_srv_sm_res':f5_fas_srv_sm_res , 'f6_fas_srv_sm_res':f6_fas_srv_sm_res , 'f7_fas_srv_sm_res':f7_fas_srv_sm_res ,
             'pgv_srv_sm_res':pgv_srv_sm_res,
             # 'f2_fas_val_srv_sm_res':f2_fas_val_srv_sm_res , 'f3_fas_val_srv_sm_res':f3_fas_val_srv_sm_res, 'f4_fas_val_srv_sm_res':f4_fas_val_srv_sm_res , 
             # 'f5_fas_val_srv_sm_res':f5_fas_val_srv_sm_res , 'f6_fas_val_srv_sm_res':f6_fas_val_srv_sm_res , 'f7_fas_val_srv_sm_res':f7_fas_val_srv_sm_res ,
             
             'f2_fas_srv_sm_ratio':f2_fas_srv_sm_ratio , 'f3_fas_srv_sm_ratio':f3_fas_srv_sm_ratio, 'f4_fas_srv_sm_ratio':f4_fas_srv_sm_ratio , 
             'f5_fas_srv_sm_ratio':f5_fas_srv_sm_ratio , 'f6_fas_srv_sm_ratio':f6_fas_srv_sm_ratio , 'f7_fas_srv_sm_ratio':f7_fas_srv_sm_ratio ,
             # 'f2_fas_srv_sm_ratio_val':f2_fas_srv_sm_ratio_val , 'f3_fas_srv_sm_ratio_val':f3_fas_srv_sm_ratio_val, 'f4_fas_srv_sm_ratio_val':f4_fas_srv_sm_ratio_val , 
             # 'f5_fas_srv_sm_ratio_val':f5_fas_srv_sm_ratio_val , 'f6_fas_srv_sm_ratio_val':f6_fas_srv_sm_ratio_val , 'f7_fas_srv_sm_ratio_val':f7_fas_srv_sm_ratio_val ,
             
             'pgv_srv_sm_ratio':pgv_srv_sm_ratio,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             'pgv_steph_res':pgv_steph_res} #, 'pgv_srv_sm_ratio_val':pgv_srv_sm_ratio_val}

df_res = pd.DataFrame(data=data_res)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res.to_csv(resrat_ff_dir+'srv_sm_res_ff.csv')

#%%

data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] , 'lat_region':lat_region_val,'model':steph_nm_ls_val,
             'f2_fas_val_steph_res':f2_fas_val_steph_res , 'f3_fas_val_steph_res':f3_fas_val_steph_res, 'f4_fas_val_steph_res':f4_fas_val_steph_res , 
             'f5_fas_val_steph_res':f5_fas_val_steph_res , 'f6_fas_val_steph_res':f6_fas_val_steph_res , 'f7_fas_val_steph_res':f7_fas_val_steph_res ,
             'pgv_steph_res_val':pgv_steph_res_val}

df_res_val = pd.DataFrame(data=data_res_val)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res_val.to_csv(resrat_ff_dir+'steph_res_ff_val.csv')
              
###################

data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] , 'lat_region':lat_region_val,'model':eugspe_1d_nm_ls_val,
              'f2_fas_val_eugspe_1d_res':f2_fas_val_eugspe_1d_res , 'f3_fas_val_eugspe_1d_res':f3_fas_val_eugspe_1d_res, 'f4_fas_val_eugspe_1d_res':f4_fas_val_eugspe_1d_res , 
              'f5_fas_val_eugspe_1d_res':f5_fas_val_eugspe_1d_res , 'f6_fas_val_eugspe_1d_res':f6_fas_val_eugspe_1d_res , 'f7_fas_val_eugspe_1d_res':f7_fas_val_eugspe_1d_res ,
              'pgv_eugspe_1d_res_val':pgv_eugspe_1d_res_val,
              'f2_fas_eugspe_1d_ratio_val':f2_fas_eugspe_1d_ratio_val , 'f3_fas_eugspe_1d_ratio_val':f3_fas_eugspe_1d_ratio_val, 'f4_fas_eugspe_1d_ratio_val':f4_fas_eugspe_1d_ratio_val , 
              'f5_fas_eugspe_1d_ratio_val':f5_fas_eugspe_1d_ratio_val , 'f6_fas_eugspe_1d_ratio_val':f6_fas_eugspe_1d_ratio_val , 'f7_fas_eugspe_1d_ratio_val':f7_fas_eugspe_1d_ratio_val ,
              'pgv_eugspe_1d_ratio_val':pgv_eugspe_1d_ratio_val,
              'f2_fas_val_steph_res':f2_fas_val_steph_res , 'f3_fas_val_steph_res':f3_fas_val_steph_res, 'f4_fas_val_steph_res':f4_fas_val_steph_res ,
              'f5_fas_val_steph_res':f5_fas_val_steph_res , 'f6_fas_val_steph_res':f6_fas_val_steph_res , 'f7_fas_val_steph_res':f7_fas_val_steph_res ,
              'pgv_steph_res_val':pgv_steph_res_val}

df_res_val = pd.DataFrame(data=data_res_val)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res_val.to_csv(resrat_ff_dir+'eugspe_1d_res_ff_val.csv')

###################

data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] , 'lat_region':lat_region_val,'model':eugspe_sm_nm_ls_val,
              'f2_fas_val_eugspe_sm_res':f2_fas_val_eugspe_sm_res , 'f3_fas_val_eugspe_sm_res':f3_fas_val_eugspe_sm_res, 'f4_fas_val_eugspe_sm_res':f4_fas_val_eugspe_sm_res , 
              'f5_fas_val_eugspe_sm_res':f5_fas_val_eugspe_sm_res , 'f6_fas_val_eugspe_sm_res':f6_fas_val_eugspe_sm_res , 'f7_fas_val_eugspe_sm_res':f7_fas_val_eugspe_sm_res ,
              'pgv_eugspe_sm_res_val':pgv_eugspe_sm_res_val,
              'f2_fas_eugspe_sm_ratio_val':f2_fas_eugspe_sm_ratio_val , 'f3_fas_eugspe_sm_ratio_val':f3_fas_eugspe_sm_ratio_val, 'f4_fas_eugspe_sm_ratio_val':f4_fas_eugspe_sm_ratio_val , 
              'f5_fas_eugspe_sm_ratio_val':f5_fas_eugspe_sm_ratio_val , 'f6_fas_eugspe_sm_ratio_val':f6_fas_eugspe_sm_ratio_val , 'f7_fas_eugspe_sm_ratio_val':f7_fas_eugspe_sm_ratio_val ,            
             'pgv_eugspe_sm_ratio_val':pgv_eugspe_sm_ratio_val,
             'f2_fas_val_steph_res':f2_fas_val_steph_res , 'f3_fas_val_steph_res':f3_fas_val_steph_res, 'f4_fas_val_steph_res':f4_fas_val_steph_res ,
             'f5_fas_val_steph_res':f5_fas_val_steph_res , 'f6_fas_val_steph_res':f6_fas_val_steph_res , 'f7_fas_val_steph_res':f7_fas_val_steph_res ,
             'pgv_steph_res_val':pgv_steph_res_val}

df_res_val = pd.DataFrame(data=data_res_val)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res_val.to_csv(resrat_ff_dir+'eugspe_sm_res_ff_val.csv')

###################

data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] , 'lat_region':lat_region_val,'model':srv_1d_nm_ls_val,
              'f2_fas_val_srv_1d_res':f2_fas_val_srv_1d_res , 'f3_fas_val_srv_1d_res':f3_fas_val_srv_1d_res, 'f4_fas_val_srv_1d_res':f4_fas_val_srv_1d_res , 
              'f5_fas_val_srv_1d_res':f5_fas_val_srv_1d_res , 'f6_fas_val_srv_1d_res':f6_fas_val_srv_1d_res , 'f7_fas_val_srv_1d_res':f7_fas_val_srv_1d_res ,
              'pgv_srv_1d_res_val':pgv_srv_1d_res_val,
              'f2_fas_srv_1d_ratio_val':f2_fas_srv_1d_ratio_val , 'f3_fas_srv_1d_ratio_val':f3_fas_srv_1d_ratio_val, 'f4_fas_srv_1d_ratio_val':f4_fas_srv_1d_ratio_val , 
              'f5_fas_srv_1d_ratio_val':f5_fas_srv_1d_ratio_val , 'f6_fas_srv_1d_ratio_val':f6_fas_srv_1d_ratio_val , 'f7_fas_srv_1d_ratio_val':f7_fas_srv_1d_ratio_val ,
             'pgv_srv_1d_ratio_val':pgv_srv_1d_ratio_val,
             'f2_fas_val_steph_res':f2_fas_val_steph_res , 'f3_fas_val_steph_res':f3_fas_val_steph_res, 'f4_fas_val_steph_res':f4_fas_val_steph_res ,
             'f5_fas_val_steph_res':f5_fas_val_steph_res , 'f6_fas_val_steph_res':f6_fas_val_steph_res , 'f7_fas_val_steph_res':f7_fas_val_steph_res ,
             'pgv_steph_res_val':pgv_steph_res_val}

df_res_val = pd.DataFrame(data=data_res_val)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res_val.to_csv(resrat_ff_dir+'srv_1d_res_ff_val.csv')

###################

data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] , 'lat_region':lat_region_val,'model':srv_sm_nm_ls_val,
              'f2_fas_val_srv_sm_res':f2_fas_val_srv_sm_res , 'f3_fas_val_srv_sm_res':f3_fas_val_srv_sm_res, 'f4_fas_val_srv_sm_res':f4_fas_val_srv_sm_res , 
              'f5_fas_val_srv_sm_res':f5_fas_val_srv_sm_res , 'f6_fas_val_srv_sm_res':f6_fas_val_srv_sm_res , 'f7_fas_val_srv_sm_res':f7_fas_val_srv_sm_res ,
              'pgv_srv_sm_res_val':pgv_srv_sm_res_val,
              'f2_fas_srv_sm_ratio_val':f2_fas_srv_sm_ratio_val , 'f3_fas_srv_sm_ratio_val':f3_fas_srv_sm_ratio_val, 'f4_fas_srv_sm_ratio_val':f4_fas_srv_sm_ratio_val , 
              'f5_fas_srv_sm_ratio_val':f5_fas_srv_sm_ratio_val , 'f6_fas_srv_sm_ratio_val':f6_fas_srv_sm_ratio_val , 'f7_fas_srv_sm_ratio_val':f7_fas_srv_sm_ratio_val ,
             'pgv_srv_sm_ratio_val':pgv_srv_sm_ratio_val,
             'f2_fas_val_steph_res':f2_fas_val_steph_res , 'f3_fas_val_steph_res':f3_fas_val_steph_res, 'f4_fas_val_steph_res':f4_fas_val_steph_res ,
             'f5_fas_val_steph_res':f5_fas_val_steph_res , 'f6_fas_val_steph_res':f6_fas_val_steph_res , 'f7_fas_val_steph_res':f7_fas_val_steph_res ,
             'pgv_steph_res_val':pgv_steph_res_val}

df_res_val = pd.DataFrame(data=data_res_val)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

df_res_val.to_csv(resrat_ff_dir+'srv_sm_res_ff_val.csv')