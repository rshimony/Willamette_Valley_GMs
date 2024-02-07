#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:57:11 2023

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
        
        fig.plot(x=evlon[3] , y=evlat[3] , style='a0.8c',fill='yellow',pen='black')
        
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

fas_obs_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/fas_ff_obs_springfield.csv')
fas_steph_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/fas_ff_steph_springfield.csv')
psa_obs_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/psa_ff_obs_springfield.csv')
psa_steph_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/psa_ff_steph_springfield.csv')

fas_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_1d/fas_ff_eugspe_1d.csv')
psa_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_1d/psa_ff_eugspe_1d.csv')
fas_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_sm/fas_ff_eugspe_sm.csv')
psa_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_sm/psa_ff_eugspe_sm.csv')

fas_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_1d/fas_ff_srv_1d.csv')
psa_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_1d/psa_ff_srv_1d.csv')
fas_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_sm/fas_ff_srv_sm.csv')
psa_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_sm/psa_ff_srv_sm.csv')

ims_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_1d/ims_ff_eugspe_1d.csv')
ims_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_sm/ims_ff_eugspe_sm.csv')

ims_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_1d/ims_ff_srv_1d.csv')
ims_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_sm/ims_ff_srv_sm.csv')

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
    t10 = np.array(spec_ims_file['T=10'])
    t5 = np.array(spec_ims_file['T=5'])
    t2 = np.array(spec_ims_file['T=2'])
    t1 = np.array(spec_ims_file['T=1'])
    
    t10_3 = []
    t5_3 = []
    t2_3 = []
    t1_3 = []
    st_names_unq = pd.unique(st_name)
    st_lons_unq = pd.unique(st_lons)
    st_lats_unq = pd.unique(st_lats)
    
    for i in range(len(st_names_unq)):
        same_st = st_name == st_names_unq[i]
        same_st_idx = np.where(same_st)[0]
        t10_3.append(t10[same_st_idx])
        t5_3.append(t5[same_st_idx])
        t2_3.append(t2[same_st_idx])
        t1_3.append(t1[same_st_idx])
        
    # return st_names_unq,st_lons_unq,st_lats_unq,t10_3,t5_3,t2_3,t1_3
    
    t10_t = []
    t5_t = []
    t2_t = []
    t1_t = []
    for j in range(len(t10_3)):
        st_t10_t = find_vec_norm(t10_3[j])
        t10_t.append(st_t10_t)
        st_t5_t = find_vec_norm(t5_3[j])
        t5_t.append(st_t5_t)
        st_t2_t = find_vec_norm(t2_3[j])
        t2_t.append(st_t2_t)
        st_t1_t = find_vec_norm(t1_3[j])
        t1_t.append(st_t1_t)
        
    
    return st_names_unq,st_lons_unq,st_lats_unq,t10_t,t5_t,t2_t,t1_t
    
def extract_ims_pgv(ims_file):
    
    st_lons = np.array(ims_file['st_lons'])
    st_lats = np.array(ims_file['st_lats'])
    st_name = np.array(ims_file['st_names'])
    pgv_t_synt = np.array(ims_file['pgv_t_synt'])
    pgv_t_obs = np.array(ims_file['pgv_t_obs'])
    pgv_t_steph = np.array(ims_file['pgv_t_steph'])
    
    return st_name,st_lons,st_lats,pgv_t_synt,pgv_t_obs,pgv_t_steph    


st_name_fas_obs , st_lons_fas_obs , st_lats_fas_obs  , t10_fas_obs , t5_fas_obs , t2_fas_obs , t1_fas_obs = extract_spec_ims_props(fas_obs_file)
st_name_fas_steph , st_lons_fas_steph , st_lats_fas_steph  , t10_fas_steph , t5_fas_steph , t2_fas_steph , t1_fas_steph = extract_spec_ims_props(fas_steph_file)
st_name_psa_obs , st_lons_psa_obs , st_lats_psa_obs  , t10_psa_obs , t5_psa_obs , t2_psa_obs , t1_psa_obs = extract_spec_ims_props(psa_obs_file)
st_name_psa_steph , st_lons_psa_steph , st_lats_psa_steph  , t10_psa_steph , t5_psa_steph , t2_psa_steph , t1_psa_steph = extract_spec_ims_props(psa_steph_file) 

st_name_fas_eugspe_1d , st_lons_fas_eugspe_1d , st_lats_fas_eugspe_1d  , t10_fas_eugspe_1d , t5_fas_eugspe_1d , t2_fas_eugspe_1d , t1_fas_eugspe_1d = extract_spec_ims_props(fas_eugspe_1d_file) 
st_name_psa_eugspe_1d , st_lons_psa_eugspe_1d , st_lats_psa_eugspe_1d  , t10_psa_eugspe_1d , t5_psa_eugspe_1d , t2_psa_eugspe_1d , t1_psa_eugspe_1d = extract_spec_ims_props(psa_eugspe_1d_file)
st_name_fas_eugspe_sm , st_lons_fas_eugspe_sm , st_lats_fas_eugspe_sm  , t10_fas_eugspe_sm , t5_fas_eugspe_sm , t2_fas_eugspe_sm , t1_fas_eugspe_sm = extract_spec_ims_props(fas_eugspe_sm_file)
st_name_psa_eugspe_sm , st_lons_psa_eugspe_sm , st_lats_psa_eugspe_sm  , t10_psa_eugspe_sm , t5_psa_eugspe_sm , t2_psa_eugspe_sm , t1_psa_eugspe_sm = extract_spec_ims_props(psa_eugspe_sm_file)

st_name_fas_srv_1d , st_lons_fas_srv_1d , st_lats_fas_srv_1d  , t10_fas_srv_1d , t5_fas_srv_1d , t2_fas_srv_1d , t1_fas_srv_1d = extract_spec_ims_props(fas_srv_1d_file)
st_name_psa_srv_1d , st_lons_psa_srv_1d , st_lats_psa_srv_1d  , t10_psa_srv_1d , t5_psa_srv_1d , t2_psa_srv_1d , t1_psa_srv_1d = extract_spec_ims_props(psa_srv_1d_file)
st_name_fas_srv_sm , st_lons_fas_srv_sm , st_lats_fas_srv_sm  , t10_fas_srv_sm , t5_fas_srv_sm , t2_fas_srv_sm , t1_fas_srv_sm = extract_spec_ims_props(fas_srv_sm_file)
st_name_psa_srv_sm , st_lons_psa_srv_sm , st_lats_psa_srv_sm  , t10_psa_srv_sm , t5_psa_srv_sm , t2_psa_srv_sm , t1_psa_srv_sm = extract_spec_ims_props(psa_srv_sm_file)

st_name_eugspe_1d , st_lons_eugspe_1d , st_lats_eugspe_1d , pgv_t_synt_eugspe_1d , pgv_t_obs_eugspe_1d , pgv_t_steph_eugspe_1d = extract_ims_pgv(ims_eugspe_1d_file)
st_name_eugspe_sm , st_lons_eugspe_sm , st_lats_eugspe_sm , pgv_t_synt_eugspe_sm , pgv_t_obs_eugspe_sm , pgv_t_steph_eugspe_sm = extract_ims_pgv(ims_eugspe_sm_file)

st_name_srv_1d , st_lons_srv_1d , st_lats_srv_1d , pgv_t_synt_srv_1d , pgv_t_obs_srv_1d , pgv_t_steph_srv_1d = extract_ims_pgv(ims_srv_1d_file)
st_name_srv_sm , st_lons_srv_sm , st_lats_srv_sm , pgv_t_synt_srv_sm , pgv_t_obs_srv_sm , pgv_t_steph_srv_sm = extract_ims_pgv(ims_srv_sm_file)

#%%

# diff_syn = np.array(t2_psa_obs) - np.array(t2_psa_eugspe_1d)
# diff_steph = np.array(t2_psa_obs) - np.array(t2_psa_steph)


#%%


t10_fas_steph_res = np.log(t10_fas_obs) - np.log(t10_fas_steph)
t5_fas_steph_res = np.log(t5_fas_obs) - np.log(t5_fas_steph)
t2_fas_steph_res = np.log(t2_fas_obs) - np.log(t2_fas_steph)
t1_fas_steph_res = np.log(t1_fas_obs) - np.log(t1_fas_steph)

t10_psa_steph_res = np.log(t10_psa_obs) - np.log(t10_psa_steph)
t5_psa_steph_res = np.log(t5_psa_obs) - np.log(t5_psa_steph)
t2_psa_steph_res = np.log(t2_psa_obs) - np.log(t2_psa_steph)
t1_psa_steph_res = np.log(t1_psa_obs) - np.log(t1_psa_steph)


t10_fas_eugspe_1d_res = np.log(t10_fas_obs) - np.log(t10_fas_eugspe_1d)
t5_fas_eugspe_1d_res = np.log(t5_fas_obs) - np.log(t5_fas_eugspe_1d)
t2_fas_eugspe_1d_res = np.log(t2_fas_obs) - np.log(t2_fas_eugspe_1d)
t1_fas_eugspe_1d_res = np.log(t1_fas_obs) - np.log(t1_fas_eugspe_1d)

t10_psa_eugspe_1d_res = np.log(t10_psa_obs) - np.log(t10_psa_eugspe_1d)
t5_psa_eugspe_1d_res = np.log(t5_psa_obs) - np.log(t5_psa_eugspe_1d)
t2_psa_eugspe_1d_res = np.log(t2_psa_obs) - np.log(t2_psa_eugspe_1d)
t1_psa_eugspe_1d_res = np.log(t1_psa_obs) - np.log(t1_psa_eugspe_1d)

t10_fas_eugspe_sm_res = np.log(t10_fas_obs) - np.log(t10_fas_eugspe_sm)
t5_fas_eugspe_sm_res = np.log(t5_fas_obs) - np.log(t5_fas_eugspe_sm)
t2_fas_eugspe_sm_res = np.log(t2_fas_obs) - np.log(t2_fas_eugspe_sm)
t1_fas_eugspe_sm_res = np.log(t1_fas_obs) - np.log(t1_fas_eugspe_sm)

t10_psa_eugspe_sm_res = np.log(t10_psa_obs) - np.log(t10_psa_eugspe_sm)
t5_psa_eugspe_sm_res = np.log(t5_psa_obs) - np.log(t5_psa_eugspe_sm)
t2_psa_eugspe_sm_res = np.log(t2_psa_obs) - np.log(t2_psa_eugspe_sm)
t1_psa_eugspe_sm_res = np.log(t1_psa_obs) - np.log(t1_psa_eugspe_sm)


t10_fas_srv_1d_res = np.log(t10_fas_obs) - np.log(t10_fas_srv_1d)
t5_fas_srv_1d_res = np.log(t5_fas_obs) - np.log(t5_fas_srv_1d)
t2_fas_srv_1d_res = np.log(t2_fas_obs) - np.log(t2_fas_srv_1d)
t1_fas_srv_1d_res = np.log(t1_fas_obs) - np.log(t1_fas_srv_1d)

t10_psa_srv_1d_res = np.log(t10_psa_obs) - np.log(t10_psa_srv_1d)
t5_psa_srv_1d_res = np.log(t5_psa_obs) - np.log(t5_psa_srv_1d)
t2_psa_srv_1d_res = np.log(t2_psa_obs) - np.log(t2_psa_srv_1d)
t1_psa_srv_1d_res = np.log(t1_psa_obs) - np.log(t1_psa_srv_1d)

t10_fas_srv_sm_res = np.log(t10_fas_obs) - np.log(t10_fas_srv_sm)
t5_fas_srv_sm_res = np.log(t5_fas_obs) - np.log(t5_fas_srv_sm)
t2_fas_srv_sm_res = np.log(t2_fas_obs) - np.log(t2_fas_srv_sm)
t1_fas_srv_sm_res = np.log(t1_fas_obs) - np.log(t1_fas_srv_sm)

t10_psa_srv_sm_res = np.log(t10_psa_obs) - np.log(t10_psa_srv_sm)
t5_psa_srv_sm_res = np.log(t5_psa_obs) - np.log(t5_psa_srv_sm)
t2_psa_srv_sm_res = np.log(t2_psa_obs) - np.log(t2_psa_srv_sm)
t1_psa_srv_sm_res = np.log(t1_psa_obs) - np.log(t1_psa_srv_sm)


pgv_t_steph_res = np.log(pgv_t_obs_eugspe_1d) - np.log(pgv_t_steph_eugspe_1d)

pgv_t_eugspe_1d_res = np.log(pgv_t_obs_eugspe_1d) - np.log(pgv_t_synt_eugspe_1d)
pgv_t_eugspe_sm_res = np.log(pgv_t_obs_eugspe_sm) - np.log(pgv_t_synt_eugspe_sm)

pgv_t_srv_1d_res = np.log(pgv_t_obs_srv_1d) - np.log(pgv_t_synt_srv_1d)
pgv_t_srv_sm_res = np.log(pgv_t_obs_srv_sm) - np.log(pgv_t_synt_srv_sm)



t10_fas_eugspe_1d_ratio = np.abs(t10_fas_steph_res)/np.abs(t10_fas_eugspe_1d_res)
t5_fas_eugspe_1d_ratio = np.abs(t5_fas_steph_res)/np.abs(t5_fas_eugspe_1d_res)
t2_fas_eugspe_1d_ratio = np.abs(t2_fas_steph_res)/np.abs(t2_fas_eugspe_1d_res)
t1_fas_eugspe_1d_ratio = np.abs(t1_fas_steph_res)/np.abs(t1_fas_eugspe_1d_res)

t10_psa_eugspe_1d_ratio = np.abs(t10_psa_steph_res)/np.abs(t10_psa_eugspe_1d_res)
t5_psa_eugspe_1d_ratio = np.abs(t5_psa_steph_res)/np.abs(t5_psa_eugspe_1d_res)
t2_psa_eugspe_1d_ratio = np.abs(t2_psa_steph_res)/np.abs(t2_psa_eugspe_1d_res)
t1_psa_eugspe_1d_ratio = np.abs(t1_psa_steph_res)/np.abs(t1_psa_eugspe_1d_res)

t10_fas_eugspe_sm_ratio = np.abs(t10_fas_steph_res)/np.abs(t10_fas_eugspe_sm_res)
t5_fas_eugspe_sm_ratio = np.abs(t5_fas_steph_res)/np.abs(t5_fas_eugspe_sm_res)
t2_fas_eugspe_sm_ratio = np.abs(t2_fas_steph_res)/np.abs(t2_fas_eugspe_sm_res)
t1_fas_eugspe_sm_ratio = np.abs(t1_fas_steph_res)/np.abs(t1_fas_eugspe_sm_res)

t10_psa_eugspe_sm_ratio = np.abs(t10_psa_steph_res)/np.abs(t10_psa_eugspe_sm_res)
t5_psa_eugspe_sm_ratio = np.abs(t5_psa_steph_res)/np.abs(t5_psa_eugspe_sm_res)
t2_psa_eugspe_sm_ratio = np.abs(t2_psa_steph_res)/np.abs(t2_psa_eugspe_sm_res)
t1_psa_eugspe_sm_ratio = np.abs(t1_psa_steph_res)/np.abs(t1_psa_eugspe_sm_res)


t10_fas_srv_1d_ratio = np.abs(t10_fas_steph_res)/np.abs(t10_fas_srv_1d_res)
t5_fas_srv_1d_ratio = np.abs(t5_fas_steph_res)/np.abs(t5_fas_srv_1d_res)
t2_fas_srv_1d_ratio = np.abs(t2_fas_steph_res)/np.abs(t2_fas_srv_1d_res)
t1_fas_srv_1d_ratio = np.abs(t1_fas_steph_res)/np.abs(t1_fas_srv_1d_res)

t10_psa_srv_1d_ratio = np.abs(t10_psa_steph_res)/np.abs(t10_psa_srv_1d_res)
t5_psa_srv_1d_ratio = np.abs(t5_psa_steph_res)/np.abs(t5_psa_srv_1d_res)
t2_psa_srv_1d_ratio = np.abs(t2_psa_steph_res)/np.abs(t2_psa_srv_1d_res)
t1_psa_srv_1d_ratio = np.abs(t1_psa_steph_res)/np.abs(t1_psa_srv_1d_res)

t10_fas_srv_sm_ratio = np.abs(t10_fas_steph_res)/np.abs(t10_fas_srv_sm_res)
t5_fas_srv_sm_ratio = np.abs(t5_fas_steph_res)/np.abs(t5_fas_srv_sm_res)
t2_fas_srv_sm_ratio = np.abs(t2_fas_steph_res)/np.abs(t2_fas_srv_sm_res)
t1_fas_srv_sm_ratio = np.abs(t1_fas_steph_res)/np.abs(t1_fas_srv_sm_res)

t10_psa_srv_sm_ratio = np.abs(t10_psa_steph_res)/np.abs(t10_psa_srv_sm_res)
t5_psa_srv_sm_ratio = np.abs(t5_psa_steph_res)/np.abs(t5_psa_srv_sm_res)
t2_psa_srv_sm_ratio = np.abs(t2_psa_steph_res)/np.abs(t2_psa_srv_sm_res)
t1_psa_srv_sm_ratio = np.abs(t1_psa_steph_res)/np.abs(t1_psa_srv_sm_res)


pgv_t_eugspe_1d_ratio = np.abs(pgv_t_steph_res)/np.abs(pgv_t_eugspe_1d_res)
pgv_t_eugspe_sm_ratio = np.abs(pgv_t_steph_res)/np.abs(pgv_t_eugspe_sm_res)

pgv_t_srv_1d_ratio = np.abs(pgv_t_steph_res)/np.abs(pgv_t_srv_1d_res)
pgv_t_srv_sm_ratio = np.abs(pgv_t_steph_res)/np.abs(pgv_t_srv_sm_res)


grav_lat_idx = np.argwhere(st_lats_fas_obs > 44.98)[:,0]
lat_region = np.zeros_like(st_lats_fas_obs)
lat_region[grav_lat_idx] = 1

lat_region = np.where(lat_region == 1, 'North', 'South')

srv_sm_nm_ls = ["srv_sm"]*len(lat_region)
srv_1d_nm_ls = ["srv_1d"]*len(lat_region)
eugspe_sm_nm_ls = ["eugspe_sm"]*len(lat_region)
eugspe_1d_nm_ls = ["eugspe_1d"]*len(lat_region)
steph_nm_ls = ["steph"]*len(lat_region)

#%%

obs_steph_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/res_maps_figs/'
try:
    os.mkdir(obs_steph_dir)
except FileExistsError:
    print('Directory already exists')
    
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't10_fas_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't5_fas_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't2_fas_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_fas_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't1_fas_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_fas_steph_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't10_psa_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_psa_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't5_psa_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_psa_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't2_psa_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_psa_steph_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 't1_psa_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_psa_steph_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , obs_steph_dir , 'pgv_t_steph_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=pgv_t_steph_res)


eugspe_1d_res_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_1d/res_maps_figs/'
try:
    os.mkdir(eugspe_1d_res_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't10_fas_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't5_fas_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't2_fas_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_fas_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't1_fas_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_fas_eugspe_1d_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't10_psa_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_psa_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't5_psa_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_psa_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't2_psa_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_psa_eugspe_1d_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 't1_psa_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_psa_eugspe_1d_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_res_figs_dir , 'pgv_t_eugspe_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=pgv_t_eugspe_1d_res)


eugspe_sm_res_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_sm/res_maps_figs/'
try:
    os.mkdir(eugspe_sm_res_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't10_fas_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't5_fas_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't2_fas_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_fas_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't1_fas_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_fas_eugspe_sm_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't10_psa_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_psa_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't5_psa_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_psa_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't2_psa_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_psa_eugspe_sm_res)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 't1_psa_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_psa_eugspe_sm_res)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_res_figs_dir , 'pgv_t_eugspe_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=pgv_t_eugspe_sm_res)


srv_1d_res_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_1d/res_maps_figs/'
try:
    os.mkdir(srv_1d_res_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't10_fas_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't5_fas_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't2_fas_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_fas_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't1_fas_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_fas_srv_1d_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't10_psa_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_psa_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't5_psa_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_psa_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't2_psa_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_psa_srv_1d_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 't1_psa_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_psa_srv_1d_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_res_figs_dir , 'pgv_t_srv_1d_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=pgv_t_srv_1d_res)


srv_sm_res_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_sm/res_maps_figs/'
try:
    os.mkdir(srv_sm_res_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't10_fas_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't5_fas_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't2_fas_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_fas_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't1_fas_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_fas_srv_sm_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't10_psa_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t10_psa_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't5_psa_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t5_psa_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't2_psa_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t2_psa_srv_sm_res)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 't1_psa_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=t1_psa_srv_sm_res)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_res_figs_dir , 'pgv_t_srv_sm_res' , st_lons_fas_obs , st_lats_fas_obs , plot_res=pgv_t_srv_sm_res)


################

eugspe_1d_ratio_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_1d/ratio_maps_figs/'
try:
    os.mkdir(eugspe_1d_ratio_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't10_fas_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't5_fas_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't2_fas_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_fas_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't1_fas_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_fas_eugspe_1d_ratio)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't10_psa_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_psa_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't5_psa_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_psa_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't2_psa_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_psa_eugspe_1d_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 't1_psa_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_psa_eugspe_1d_ratio)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_1d_ratio_figs_dir , 'pgv_t_eugspe_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=pgv_t_eugspe_1d_ratio)


eugspe_sm_ratio_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_sm/ratio_maps_figs/'
try:
    os.mkdir(eugspe_sm_ratio_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't10_fas_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't5_fas_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't2_fas_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_fas_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't1_fas_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_fas_eugspe_sm_ratio)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't10_psa_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_psa_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't5_psa_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_psa_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't2_psa_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_psa_eugspe_sm_ratio)
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 't1_psa_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_psa_eugspe_sm_ratio)

plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , eugspe_sm_ratio_figs_dir , 'pgv_t_eugspe_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=pgv_t_eugspe_sm_ratio)


srv_1d_ratio_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_1d/ratio_maps_figs/'
try:
    os.mkdir(srv_1d_ratio_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't10_fas_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't5_fas_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't2_fas_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_fas_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't1_fas_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_fas_srv_1d_ratio)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't10_psa_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_psa_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't5_psa_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_psa_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't2_psa_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_psa_srv_1d_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 't1_psa_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_psa_srv_1d_ratio)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_1d_ratio_figs_dir , 'pgv_t_srv_1d_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=pgv_t_srv_1d_ratio)


srv_sm_ratio_figs_dir = '/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_sm/ratio_maps_figs/'
try:
    os.mkdir(srv_sm_ratio_figs_dir)
except FileExistsError:
    print('Directory already exists')

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't10_fas_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't5_fas_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't2_fas_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_fas_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't1_fas_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_fas_srv_sm_ratio)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't10_psa_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t10_psa_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't5_psa_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t5_psa_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't2_psa_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t2_psa_srv_sm_ratio)
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 't1_psa_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=t1_psa_srv_sm_ratio)

plot_shade_basin_map(srv_lons , srv_lats , srv_depths , srv_sm_ratio_figs_dir , 'pgv_t_srv_sm_ratio' , st_lons_fas_obs , st_lats_fas_obs , plot_res_ratio=pgv_t_srv_sm_ratio)

#%%

data_res = {'st_name':st_name_fas_obs ,'st_lon':st_lons_fas_obs, 'st_lat':st_lats_fas_obs , 'lat_region':lat_region,'model':steph_nm_ls,
             't10_fas_steph_res':t10_fas_steph_res , 't5_fas_steph_res':t5_fas_steph_res, 't2_fas_steph_res':t2_fas_steph_res , 't1_fas_steph_res':t1_fas_steph_res ,
             't10_psa_steph_res':t10_psa_steph_res , 't5_psa_steph_res':t5_psa_steph_res, 't2_psa_steph_res':t2_psa_steph_res , 't1_psa_steph_res':t1_psa_steph_res ,
             'pgv_t_steph_res':pgv_t_steph_res}

df_res = pd.DataFrame(data=data_res)

df_res.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/steph_res_ff.csv')

###################

data_res = {'st_name':st_name_fas_obs ,'st_lon':st_lons_fas_obs, 'st_lat':st_lats_fas_obs , 'lat_region':lat_region,'model':eugspe_1d_nm_ls,
             't10_fas_eugspe_1d_res':t10_fas_eugspe_1d_res , 't5_fas_eugspe_1d_res':t5_fas_eugspe_1d_res, 't2_fas_eugspe_1d_res':t2_fas_eugspe_1d_res , 't1_fas_eugspe_1d_res':t1_fas_eugspe_1d_res ,
             't10_psa_eugspe_1d_res':t10_psa_eugspe_1d_res , 't5_psa_eugspe_1d_res':t5_psa_eugspe_1d_res, 't2_psa_eugspe_1d_res':t2_psa_eugspe_1d_res , 't1_psa_eugspe_1d_res':t1_psa_eugspe_1d_res ,
             'pgv_t_eugspe_1d_res':pgv_t_eugspe_1d_res , 'pgv_t_eugspe_1d_ratio':pgv_t_eugspe_1d_ratio ,
             't10_fas_eugspe_1d_ratio':t10_fas_eugspe_1d_ratio , 't5_fas_eugspe_1d_ratio':t5_fas_eugspe_1d_ratio, 't2_fas_eugspe_1d_ratio':t2_fas_eugspe_1d_ratio , 't1_fas_eugspe_1d_ratio':t1_fas_eugspe_1d_ratio ,
             't10_psa_eugspe_1d_ratio':t10_psa_eugspe_1d_ratio , 't5_psa_eugspe_1d_ratio':t5_psa_eugspe_1d_ratio, 't2_psa_eugspe_1d_ratio':t2_psa_eugspe_1d_ratio , 't1_psa_eugspe_1d_ratio':t1_psa_eugspe_1d_ratio ,
             't10_fas_steph_res':t10_fas_steph_res , 't5_fas_steph_res':t5_fas_steph_res, 't2_fas_steph_res':t2_fas_steph_res , 't1_fas_steph_res':t1_fas_steph_res ,
             't10_psa_steph_res':t10_psa_steph_res , 't5_psa_steph_res':t5_psa_steph_res, 't2_psa_steph_res':t2_psa_steph_res , 't1_psa_steph_res':t1_psa_steph_res ,
             'pgv_t_steph_res':pgv_t_steph_res}

df_res = pd.DataFrame(data=data_res)

df_res.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_1d/eugspe_1d_res_ff.csv')

#################

data_res = {'st_name':st_name_fas_obs ,'st_lon':st_lons_fas_obs, 'st_lat':st_lats_fas_obs , 'lat_region':lat_region,'model':eugspe_sm_nm_ls,
             't10_fas_eugspe_sm_res':t10_fas_eugspe_sm_res , 't5_fas_eugspe_sm_res':t5_fas_eugspe_sm_res, 't2_fas_eugspe_sm_res':t2_fas_eugspe_sm_res , 't1_fas_eugspe_sm_res':t1_fas_eugspe_sm_res ,
             't10_psa_eugspe_sm_res':t10_psa_eugspe_sm_res , 't5_psa_eugspe_sm_res':t5_psa_eugspe_sm_res, 't2_psa_eugspe_sm_res':t2_psa_eugspe_sm_res , 't1_psa_eugspe_sm_res':t1_psa_eugspe_sm_res ,
             'pgv_t_eugspe_sm_res':pgv_t_eugspe_sm_res , 'pgv_t_eugspe_sm_ratio':pgv_t_eugspe_sm_ratio ,
             't10_fas_eugspe_sm_ratio':t10_fas_eugspe_sm_ratio , 't5_fas_eugspe_sm_ratio':t5_fas_eugspe_sm_ratio, 't2_fas_eugspe_sm_ratio':t2_fas_eugspe_sm_ratio , 't1_fas_eugspe_sm_ratio':t1_fas_eugspe_sm_ratio ,
             't10_psa_eugspe_sm_ratio':t10_psa_eugspe_sm_ratio , 't5_psa_eugspe_sm_ratio':t5_psa_eugspe_sm_ratio, 't2_psa_eugspe_sm_ratio':t2_psa_eugspe_sm_ratio , 't1_psa_eugspe_sm_ratio':t1_psa_eugspe_sm_ratio ,
             't10_fas_steph_res':t10_fas_steph_res , 't5_fas_steph_res':t5_fas_steph_res, 't2_fas_steph_res':t2_fas_steph_res , 't1_fas_steph_res':t1_fas_steph_res ,
             't10_psa_steph_res':t10_psa_steph_res , 't5_psa_steph_res':t5_psa_steph_res, 't2_psa_steph_res':t2_psa_steph_res , 't1_psa_steph_res':t1_psa_steph_res ,
             'pgv_t_steph_res':pgv_t_steph_res}

df_res = pd.DataFrame(data=data_res)

df_res.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/eugspe_sm/eugspe_sm_res_ff.csv')

###################

data_res = {'st_name':st_name_fas_obs ,'st_lon':st_lons_fas_obs, 'st_lat':st_lats_fas_obs , 'lat_region':lat_region,'model':srv_1d_nm_ls,
             't10_fas_srv_1d_res':t10_fas_srv_1d_res , 't5_fas_srv_1d_res':t5_fas_srv_1d_res, 't2_fas_srv_1d_res':t2_fas_srv_1d_res , 't1_fas_srv_1d_res':t1_fas_srv_1d_res ,
             't10_psa_srv_1d_res':t10_psa_srv_1d_res , 't5_psa_srv_1d_res':t5_psa_srv_1d_res, 't2_psa_srv_1d_res':t2_psa_srv_1d_res , 't1_psa_srv_1d_res':t1_psa_srv_1d_res ,
             'pgv_t_srv_1d_res':pgv_t_srv_1d_res , 'pgv_t_srv_1d_ratio':pgv_t_srv_1d_ratio ,
             't10_fas_srv_1d_ratio':t10_fas_srv_1d_ratio , 't5_fas_srv_1d_ratio':t5_fas_srv_1d_ratio, 't2_fas_srv_1d_ratio':t2_fas_srv_1d_ratio , 't1_fas_srv_1d_ratio':t1_fas_srv_1d_ratio ,
             't10_psa_srv_1d_ratio':t10_psa_srv_1d_ratio , 't5_psa_srv_1d_ratio':t5_psa_srv_1d_ratio, 't2_psa_srv_1d_ratio':t2_psa_srv_1d_ratio , 't1_psa_srv_1d_ratio':t1_psa_srv_1d_ratio ,
             't10_fas_steph_res':t10_fas_steph_res , 't5_fas_steph_res':t5_fas_steph_res, 't2_fas_steph_res':t2_fas_steph_res , 't1_fas_steph_res':t1_fas_steph_res ,
             't10_psa_steph_res':t10_psa_steph_res , 't5_psa_steph_res':t5_psa_steph_res, 't2_psa_steph_res':t2_psa_steph_res , 't1_psa_steph_res':t1_psa_steph_res ,
             'pgv_t_steph_res':pgv_t_steph_res}

df_res = pd.DataFrame(data=data_res)

df_res.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_1d/srv_1d_res_ff.csv')

###################

data_res = {'st_name':st_name_fas_obs ,'st_lon':st_lons_fas_obs, 'st_lat':st_lats_fas_obs , 'lat_region':lat_region,'model':srv_sm_nm_ls,
             't10_fas_srv_sm_res':t10_fas_srv_sm_res , 't5_fas_srv_sm_res':t5_fas_srv_sm_res, 't2_fas_srv_sm_res':t2_fas_srv_sm_res , 't1_fas_srv_sm_res':t1_fas_srv_sm_res ,
             't10_psa_srv_sm_res':t10_psa_srv_sm_res , 't5_psa_srv_sm_res':t5_psa_srv_sm_res, 't2_psa_srv_sm_res':t2_psa_srv_sm_res , 't1_psa_srv_sm_res':t1_psa_srv_sm_res ,
             'pgv_t_srv_sm_res':pgv_t_srv_sm_res , 'pgv_t_srv_sm_ratio':pgv_t_srv_sm_ratio ,
             't10_fas_srv_sm_ratio':t10_fas_srv_sm_ratio , 't5_fas_srv_sm_ratio':t5_fas_srv_sm_ratio, 't2_fas_srv_sm_ratio':t2_fas_srv_sm_ratio , 't1_fas_srv_sm_ratio':t1_fas_srv_sm_ratio ,
             't10_psa_srv_sm_ratio':t10_psa_srv_sm_ratio , 't5_psa_srv_sm_ratio':t5_psa_srv_sm_ratio, 't2_psa_srv_sm_ratio':t2_psa_srv_sm_ratio , 't1_psa_srv_sm_ratio':t1_psa_srv_sm_ratio ,
             't10_fas_steph_res':t10_fas_steph_res , 't5_fas_steph_res':t5_fas_steph_res, 't2_fas_steph_res':t2_fas_steph_res , 't1_fas_steph_res':t1_fas_steph_res ,
             't10_psa_steph_res':t10_psa_steph_res , 't5_psa_steph_res':t5_psa_steph_res, 't2_psa_steph_res':t2_psa_steph_res , 't1_psa_steph_res':t1_psa_steph_res ,
             'pgv_t_steph_res':pgv_t_steph_res}

df_res = pd.DataFrame(data=data_res)

df_res.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/final_sims/springfield/srv_sm/srv_sm_res_ff.csv')



































