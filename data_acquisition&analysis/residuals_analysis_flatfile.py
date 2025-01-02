#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 19:44:49 2024

@author: rshimony
"""
'''
Creates RESIDUALS values flatfiles.
This script reads IMs flatfiles of a SINGLE EVENT from ALL MODELS (Full and Valley station inventories).
Calculates residual values (observed - synthetics) and residuals ratios (USGS CVM / WV model).
Writing into flatfiles for each model (Full and Valley station inventories).
Writing additional ALL MODEL and MELTED flatfiles:
    ALL MODEL - Stacking all models in a single flatfile (Full and Valley station inventories).
    MELTED - "Melting" all models by IM.
'''

import pandas as pd
import numpy as np
from shapely.geometry.polygon import Polygon
import os

#%%
##Functions and helpers##
# Creating WV polygon
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

def find_vec_norm(Im_value_list):
    '''
    Get the norm of the three components of recording
    '''
    import numpy as np

    comp_1 = Im_value_list[0]**2;
    comp_2 = Im_value_list[1]**2;
    comp_3 = Im_value_list[2]**2;

    vec_norm = np.sqrt(comp_1 + comp_2 + comp_3)

    return(vec_norm)

def extract_spec_ims_props(spec_ims_file):
    '''
    Extract FAS data from a flatfile.
    Create lists of three components for each frequency.
    Create a norm value (t) for eachfrequency bin.
    '''
    
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
    '''
    Extract PGV values from a larger IMs flatfile.
    Extracting norm PGV values for synthetics (both new WV model and USGS CVM) and observed
    '''
    st_lons = np.array(ims_file['st_lons'])
    st_lats = np.array(ims_file['st_lats'])
    st_name = np.array(ims_file['st_names'])
    pgv_t_synt = np.array(ims_file['pgv_t_synt'])
    pgv_t_obs = np.array(ims_file['pgv_t_obs'])
    pgv_t_steph = np.array(ims_file['pgv_t_steph'])
    
    return st_name,st_lons,st_lats,pgv_t_synt,pgv_t_obs,pgv_t_steph 

#%%
# Reading in FAS and IMs flatfiles (both full and Valley stations inventory)
fas_obs_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/fas_ff_obs_springfield.csv')
fas_steph_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/fas_ff_steph_springfield.csv')
fas_val_obs_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/val_fas_ff_obs_springfield.csv')
fas_val_steph_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/val_fas_ff_steph_springfield.csv')

fas_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/fas_ff_eugspe_1d_springfield.csv')
fas_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/fas_ff_eugspe_sm_springfield.csv')
fas_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/fas_ff_srv_1d_springfield.csv')
fas_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/fas_ff_srv_sm_springfield.csv')

fas_val_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/val_fas_ff_eugspe_1d_springfield.csv')
fas_val_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/val_fas_ff_eugspe_sm_springfield.csv')
fas_val_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/val_fas_ff_srv_1d_springfield.csv')
fas_val_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/val_fas_ff_srv_sm_springfield.csv')

ims_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/ims_ff_eugspe_1d.csv')
ims_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/ims_ff_eugspe_sm.csv')
ims_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/ims_ff_srv_1d.csv')
ims_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/ims_ff_srv_sm.csv')

ims_val_eugspe_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/val_ims_ff_eugspe_1d.csv')
ims_val_eugspe_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/val_ims_ff_eugspe_sm.csv')
ims_val_srv_1d_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/val_ims_ff_srv_1d.csv')
ims_val_srv_sm_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/val_ims_ff_srv_sm.csv')

# Output residuals flatfiles directory
resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

#%%
# Getting IMs values from flatfiles
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
## Calculating residuals (observed - synthetics) for all data:
# FAS USGS CVM - full station inventory
f2_fas_steph_res = np.log(fas_obs[3]) - np.log(fas_steph[3])
f3_fas_steph_res = np.log(fas_obs[4]) - np.log(fas_steph[4])
f4_fas_steph_res = np.log(fas_obs[5]) - np.log(fas_steph[5])
f5_fas_steph_res = np.log(fas_obs[6]) - np.log(fas_steph[6])
f6_fas_steph_res = np.log(fas_obs[7]) - np.log(fas_steph[7])
f7_fas_steph_res = np.log(fas_obs[8]) - np.log(fas_steph[8])

# FAS USGS CVM - valley station inventory
f2_fas_val_steph_res = np.log(fas_val_obs[3]) - np.log(fas_val_steph[3])
f3_fas_val_steph_res = np.log(fas_val_obs[4]) - np.log(fas_val_steph[4])
f4_fas_val_steph_res = np.log(fas_val_obs[5]) - np.log(fas_val_steph[5])
f5_fas_val_steph_res = np.log(fas_val_obs[6]) - np.log(fas_val_steph[6])
f6_fas_val_steph_res = np.log(fas_val_obs[7]) - np.log(fas_val_steph[7])
f7_fas_val_steph_res = np.log(fas_val_obs[8]) - np.log(fas_val_steph[8])

# FAS EugSpe 1D - full station inventory
f2_fas_eugspe_1d_res = np.log(fas_obs[3]) - np.log(fas_eugspe_1d[3])
f3_fas_eugspe_1d_res = np.log(fas_obs[4]) - np.log(fas_eugspe_1d[4])
f4_fas_eugspe_1d_res = np.log(fas_obs[5]) - np.log(fas_eugspe_1d[5])
f5_fas_eugspe_1d_res = np.log(fas_obs[6]) - np.log(fas_eugspe_1d[6])
f6_fas_eugspe_1d_res = np.log(fas_obs[7]) - np.log(fas_eugspe_1d[7])
f7_fas_eugspe_1d_res = np.log(fas_obs[8]) - np.log(fas_eugspe_1d[8])

# FAS EugSpe 1D - valley station inventory
f2_fas_val_eugspe_1d_res = np.log(fas_val_obs[3]) - np.log(fas_val_eugspe_1d[3])
f3_fas_val_eugspe_1d_res = np.log(fas_val_obs[4]) - np.log(fas_val_eugspe_1d[4])
f4_fas_val_eugspe_1d_res = np.log(fas_val_obs[5]) - np.log(fas_val_eugspe_1d[5])
f5_fas_val_eugspe_1d_res = np.log(fas_val_obs[6]) - np.log(fas_val_eugspe_1d[6])
f6_fas_val_eugspe_1d_res = np.log(fas_val_obs[7]) - np.log(fas_val_eugspe_1d[7])
f7_fas_val_eugspe_1d_res = np.log(fas_val_obs[8]) - np.log(fas_val_eugspe_1d[8])

# FAS EugSpe SM - full station inventory
f2_fas_eugspe_sm_res = np.log(fas_obs[3]) - np.log(fas_eugspe_sm[3])
f3_fas_eugspe_sm_res = np.log(fas_obs[4]) - np.log(fas_eugspe_sm[4])
f4_fas_eugspe_sm_res = np.log(fas_obs[5]) - np.log(fas_eugspe_sm[5])
f5_fas_eugspe_sm_res = np.log(fas_obs[6]) - np.log(fas_eugspe_sm[6])
f6_fas_eugspe_sm_res = np.log(fas_obs[7]) - np.log(fas_eugspe_sm[7])
f7_fas_eugspe_sm_res = np.log(fas_obs[8]) - np.log(fas_eugspe_sm[8])

# FAS EugSpe SM - valley station inventory
f2_fas_val_eugspe_sm_res = np.log(fas_val_obs[3]) - np.log(fas_val_eugspe_sm[3])
f3_fas_val_eugspe_sm_res = np.log(fas_val_obs[4]) - np.log(fas_val_eugspe_sm[4])
f4_fas_val_eugspe_sm_res = np.log(fas_val_obs[5]) - np.log(fas_val_eugspe_sm[5])
f5_fas_val_eugspe_sm_res = np.log(fas_val_obs[6]) - np.log(fas_val_eugspe_sm[6])
f6_fas_val_eugspe_sm_res = np.log(fas_val_obs[7]) - np.log(fas_val_eugspe_sm[7])
f7_fas_val_eugspe_sm_res = np.log(fas_val_obs[8]) - np.log(fas_val_eugspe_sm[8])

# FAS SRV 1D - full station inventory
f2_fas_srv_1d_res = np.log(fas_obs[3]) - np.log(fas_srv_1d[3])
f3_fas_srv_1d_res = np.log(fas_obs[4]) - np.log(fas_srv_1d[4])
f4_fas_srv_1d_res = np.log(fas_obs[5]) - np.log(fas_srv_1d[5])
f5_fas_srv_1d_res = np.log(fas_obs[6]) - np.log(fas_srv_1d[6])
f6_fas_srv_1d_res = np.log(fas_obs[7]) - np.log(fas_srv_1d[7])
f7_fas_srv_1d_res = np.log(fas_obs[8]) - np.log(fas_srv_1d[8])

# FAS SRV 1D - valley station inventory
f2_fas_val_srv_1d_res = np.log(fas_val_obs[3]) - np.log(fas_val_srv_1d[3])
f3_fas_val_srv_1d_res = np.log(fas_val_obs[4]) - np.log(fas_val_srv_1d[4])
f4_fas_val_srv_1d_res = np.log(fas_val_obs[5]) - np.log(fas_val_srv_1d[5])
f5_fas_val_srv_1d_res = np.log(fas_val_obs[6]) - np.log(fas_val_srv_1d[6])
f6_fas_val_srv_1d_res = np.log(fas_val_obs[7]) - np.log(fas_val_srv_1d[7])
f7_fas_val_srv_1d_res = np.log(fas_val_obs[8]) - np.log(fas_val_srv_1d[8])

# FAS SRV SM - full station inventory
f2_fas_srv_sm_res = np.log(fas_obs[3]) - np.log(fas_srv_sm[3])
f3_fas_srv_sm_res = np.log(fas_obs[4]) - np.log(fas_srv_sm[4])
f4_fas_srv_sm_res = np.log(fas_obs[5]) - np.log(fas_srv_sm[5])
f5_fas_srv_sm_res = np.log(fas_obs[6]) - np.log(fas_srv_sm[6])
f6_fas_srv_sm_res = np.log(fas_obs[7]) - np.log(fas_srv_sm[7])
f7_fas_srv_sm_res = np.log(fas_obs[8]) - np.log(fas_srv_sm[8])

# FAS SRV SM - valley station inventory
f2_fas_val_srv_sm_res = np.log(fas_val_obs[3]) - np.log(fas_val_srv_sm[3])
f3_fas_val_srv_sm_res = np.log(fas_val_obs[4]) - np.log(fas_val_srv_sm[4])
f4_fas_val_srv_sm_res = np.log(fas_val_obs[5]) - np.log(fas_val_srv_sm[5])
f5_fas_val_srv_sm_res = np.log(fas_val_obs[6]) - np.log(fas_val_srv_sm[6])
f6_fas_val_srv_sm_res = np.log(fas_val_obs[7]) - np.log(fas_val_srv_sm[7])
f7_fas_val_srv_sm_res = np.log(fas_val_obs[8]) - np.log(fas_val_srv_sm[8])

# PGV USGS CVM - full station inventory
pgv_steph_res = np.log(pgv_eugspe_1d[4]) - np.log(pgv_eugspe_1d[5])
# PGV USGS CVM - valley station inventory
pgv_steph_res_val = np.log(pgv_val_eugspe_1d[4]) - np.log(pgv_val_eugspe_1d[5])

# PGV EugSpe 1D - full station inventory
pgv_eugspe_1d_res = np.log(pgv_eugspe_1d[4]) - np.log(pgv_eugspe_1d[3])
# PGV EugSpe 1D - valley station inventory
pgv_eugspe_1d_res_val = np.log(pgv_val_eugspe_1d[4]) - np.log(pgv_val_eugspe_1d[3])

# PGV EugSpe SM - full station inventory
pgv_eugspe_sm_res = np.log(pgv_eugspe_sm[4]) - np.log(pgv_eugspe_sm[3])
# PGV EugSpe SM - valley station inventory
pgv_eugspe_sm_res_val = np.log(pgv_val_eugspe_sm[4]) - np.log(pgv_val_eugspe_sm[3])

# PGV SRV 1D - full station inventory
pgv_srv_1d_res = np.log(pgv_srv_1d[4]) - np.log(pgv_srv_1d[3])
# PGV SRV 1D - valley station inventory
pgv_srv_1d_res_val = np.log(pgv_val_srv_1d[4]) - np.log(pgv_val_srv_1d[3])

# PGV SRV SM - full station inventory
pgv_srv_sm_res = np.log(pgv_srv_sm[4]) - np.log(pgv_srv_sm[3])
# PGV SRV SM - valley station inventory
pgv_srv_sm_res_val = np.log(pgv_val_srv_sm[4]) - np.log(pgv_val_srv_sm[3])

#%%
## Calculating residuals ratio (USGS CVM / WV model) for all data:
# FAS EugSpe 1D - full station inventory
f2_fas_eugspe_1d_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_eugspe_1d_res)
f3_fas_eugspe_1d_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_eugspe_1d_res)
f4_fas_eugspe_1d_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_eugspe_1d_res)
f5_fas_eugspe_1d_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_eugspe_1d_res)
f6_fas_eugspe_1d_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_eugspe_1d_res)
f7_fas_eugspe_1d_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_eugspe_1d_res)

# FAS EugSpe 1D - valley station inventory
f2_fas_eugspe_1d_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_eugspe_1d_res)
f3_fas_eugspe_1d_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_eugspe_1d_res)
f4_fas_eugspe_1d_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_eugspe_1d_res)
f5_fas_eugspe_1d_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_eugspe_1d_res)
f6_fas_eugspe_1d_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_eugspe_1d_res)
f7_fas_eugspe_1d_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_eugspe_1d_res)

# FAS EugSpe SM - full station inventory
f2_fas_eugspe_sm_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_eugspe_sm_res)
f3_fas_eugspe_sm_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_eugspe_sm_res)
f4_fas_eugspe_sm_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_eugspe_sm_res)
f5_fas_eugspe_sm_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_eugspe_sm_res)
f6_fas_eugspe_sm_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_eugspe_sm_res)
f7_fas_eugspe_sm_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_eugspe_sm_res)

# FAS EugSpe SM - valley station inventory
f2_fas_eugspe_sm_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_eugspe_sm_res)
f3_fas_eugspe_sm_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_eugspe_sm_res)
f4_fas_eugspe_sm_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_eugspe_sm_res)
f5_fas_eugspe_sm_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_eugspe_sm_res)
f6_fas_eugspe_sm_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_eugspe_sm_res)
f7_fas_eugspe_sm_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_eugspe_sm_res)

# FAS SRV 1D - full station inventory
f2_fas_srv_1d_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_srv_1d_res)
f3_fas_srv_1d_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_srv_1d_res)
f4_fas_srv_1d_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_srv_1d_res)
f5_fas_srv_1d_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_srv_1d_res)
f6_fas_srv_1d_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_srv_1d_res)
f7_fas_srv_1d_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_srv_1d_res)

# FAS SRV 1D - valley station inventory
f2_fas_srv_1d_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_srv_1d_res)
f3_fas_srv_1d_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_srv_1d_res)
f4_fas_srv_1d_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_srv_1d_res)
f5_fas_srv_1d_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_srv_1d_res)
f6_fas_srv_1d_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_srv_1d_res)
f7_fas_srv_1d_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_srv_1d_res)

# FAS SRV SM - full station inventory
f2_fas_srv_sm_ratio = np.abs(f2_fas_steph_res)/np.abs(f2_fas_srv_sm_res)
f3_fas_srv_sm_ratio = np.abs(f3_fas_steph_res)/np.abs(f3_fas_srv_sm_res)
f4_fas_srv_sm_ratio = np.abs(f4_fas_steph_res)/np.abs(f4_fas_srv_sm_res)
f5_fas_srv_sm_ratio = np.abs(f5_fas_steph_res)/np.abs(f5_fas_srv_sm_res)
f6_fas_srv_sm_ratio = np.abs(f6_fas_steph_res)/np.abs(f6_fas_srv_sm_res)
f7_fas_srv_sm_ratio = np.abs(f7_fas_steph_res)/np.abs(f7_fas_srv_sm_res)

# FAS SRV SM - valley station inventory
f2_fas_srv_sm_ratio_val = np.abs(f2_fas_val_steph_res)/np.abs(f2_fas_val_srv_sm_res)
f3_fas_srv_sm_ratio_val = np.abs(f3_fas_val_steph_res)/np.abs(f3_fas_val_srv_sm_res)
f4_fas_srv_sm_ratio_val = np.abs(f4_fas_val_steph_res)/np.abs(f4_fas_val_srv_sm_res)
f5_fas_srv_sm_ratio_val = np.abs(f5_fas_val_steph_res)/np.abs(f5_fas_val_srv_sm_res)
f6_fas_srv_sm_ratio_val = np.abs(f6_fas_val_steph_res)/np.abs(f6_fas_val_srv_sm_res)
f7_fas_srv_sm_ratio_val = np.abs(f7_fas_val_steph_res)/np.abs(f7_fas_val_srv_sm_res)

# PGV WV Models - full station inventory
pgv_eugspe_1d_ratio = np.abs(pgv_steph_res)/np.abs(pgv_eugspe_1d_res)
pgv_eugspe_sm_ratio = np.abs(pgv_steph_res)/np.abs(pgv_eugspe_sm_res)
pgv_srv_1d_ratio = np.abs(pgv_steph_res)/np.abs(pgv_srv_1d_res)
pgv_srv_sm_ratio = np.abs(pgv_steph_res)/np.abs(pgv_srv_sm_res)

# PGV WV Models - valley station inventory
pgv_eugspe_1d_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_eugspe_1d_res_val)
pgv_eugspe_sm_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_eugspe_sm_res_val)
pgv_srv_1d_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_srv_1d_res_val)
pgv_srv_sm_ratio_val = np.abs(pgv_steph_res_val)/np.abs(pgv_srv_sm_res_val)

#%%
##### Creating more data for inserting into residual flatfiles
## Creating arrats of model names
# Full station inventory
srv_sm_nm_ls = ["srv_sm"]*len(pgv_srv_sm_ratio)
srv_1d_nm_ls = ["srv_1d"]*len(pgv_srv_sm_ratio)
eugspe_sm_nm_ls = ["eugspe_sm"]*len(pgv_srv_sm_ratio)
eugspe_1d_nm_ls = ["eugspe_1d"]*len(pgv_srv_sm_ratio)
steph_nm_ls = ["steph"]*len(pgv_srv_sm_ratio)

# Valley station inventory
srv_sm_nm_ls_val = ["srv_sm"]*len(pgv_srv_sm_ratio_val)
srv_1d_nm_ls_val = ["srv_1d"]*len(pgv_srv_sm_ratio_val)
eugspe_sm_nm_ls_val = ["eugspe_sm"]*len(pgv_srv_sm_ratio_val)
eugspe_1d_nm_ls_val = ["eugspe_1d"]*len(pgv_srv_sm_ratio_val)
steph_nm_ls_val = ["steph"]*len(pgv_srv_sm_ratio_val)

## Creating valley in/out array. Added to the full station inventory flatfiles.
full_st = ims_eugspe_1d_file['st_names']
valley_st = ims_val_eugspe_1d_file['st_names']
valley_st_mask = np.isin(full_st,valley_st)

valley_ls = []
for i in valley_st_mask:
    if i == True:
        valley_ls.append('in_valley')
    else:
        valley_ls.append('out_valley')
valley_array = np.array(valley_ls)

#%%
## Writing residuals flatfiles
# USGS CVM Full stations inventory
data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'location':valley_array,'model':steph_nm_ls,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res , 
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             'pgv_steph_res':pgv_steph_res}

df_res = pd.DataFrame(data=data_res)
df_res.to_csv(resrat_ff_dir+'steph_res_ff.csv')
              
# USGS CVM Valley stations inventory
data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] ,'model':steph_nm_ls_val,
             'f2_fas_val_steph_res':f2_fas_val_steph_res , 'f3_fas_val_steph_res':f3_fas_val_steph_res, 'f4_fas_val_steph_res':f4_fas_val_steph_res , 
             'f5_fas_val_steph_res':f5_fas_val_steph_res , 'f6_fas_val_steph_res':f6_fas_val_steph_res , 'f7_fas_val_steph_res':f7_fas_val_steph_res ,
             'pgv_steph_res_val':pgv_steph_res_val}

df_res_val = pd.DataFrame(data=data_res_val)
df_res_val.to_csv(resrat_ff_dir+'steph_res_ff_val.csv')

###################

# EugSpe 1D Full stations inventory
data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'location':valley_array,'model':eugspe_1d_nm_ls,
             'f2_fas_eugspe_1d_res':f2_fas_eugspe_1d_res , 'f3_fas_eugspe_1d_res':f3_fas_eugspe_1d_res, 'f4_fas_eugspe_1d_res':f4_fas_eugspe_1d_res , 
             'f5_fas_eugspe_1d_res':f5_fas_eugspe_1d_res , 'f6_fas_eugspe_1d_res':f6_fas_eugspe_1d_res , 'f7_fas_eugspe_1d_res':f7_fas_eugspe_1d_res ,
             'pgv_eugspe_1d_res':pgv_eugspe_1d_res,
             
             'f2_fas_eugspe_1d_ratio':f2_fas_eugspe_1d_ratio , 'f3_fas_eugspe_1d_ratio':f3_fas_eugspe_1d_ratio, 'f4_fas_eugspe_1d_ratio':f4_fas_eugspe_1d_ratio , 
             'f5_fas_eugspe_1d_ratio':f5_fas_eugspe_1d_ratio , 'f6_fas_eugspe_1d_ratio':f6_fas_eugspe_1d_ratio , 'f7_fas_eugspe_1d_ratio':f7_fas_eugspe_1d_ratio ,
              'pgv_eugspe_1d_ratio':pgv_eugspe_1d_ratio,
              
              'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
              'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
              'pgv_steph_res':pgv_steph_res}

df_res = pd.DataFrame(data=data_res)
df_res.to_csv(resrat_ff_dir+'eugspe_1d_res_ff.csv')

# EugSpe 1D Valley stations inventory
data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] ,'model':eugspe_1d_nm_ls_val,
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
df_res_val.to_csv(resrat_ff_dir+'eugspe_1d_res_ff_val.csv')

###################

# EugSpe SM Full stations inventory
data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'location':valley_array,'model':eugspe_sm_nm_ls,
             'f2_fas_eugspe_sm_res':f2_fas_eugspe_sm_res , 'f3_fas_eugspe_sm_res':f3_fas_eugspe_sm_res, 'f4_fas_eugspe_sm_res':f4_fas_eugspe_sm_res , 
             'f5_fas_eugspe_sm_res':f5_fas_eugspe_sm_res , 'f6_fas_eugspe_sm_res':f6_fas_eugspe_sm_res , 'f7_fas_eugspe_sm_res':f7_fas_eugspe_sm_res ,
             'pgv_eugspe_sm_res':pgv_eugspe_sm_res,
             
             'f2_fas_eugspe_sm_ratio':f2_fas_eugspe_sm_ratio , 'f3_fas_eugspe_sm_ratio':f3_fas_eugspe_sm_ratio, 'f4_fas_eugspe_sm_ratio':f4_fas_eugspe_sm_ratio , 
             'f5_fas_eugspe_sm_ratio':f5_fas_eugspe_sm_ratio , 'f6_fas_eugspe_sm_ratio':f6_fas_eugspe_sm_ratio , 'f7_fas_eugspe_sm_ratio':f7_fas_eugspe_sm_ratio ,
             
             'pgv_eugspe_sm_ratio':pgv_eugspe_sm_ratio,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             'pgv_steph_res':pgv_steph_res}

df_res = pd.DataFrame(data=data_res)
df_res.to_csv(resrat_ff_dir+'eugspe_sm_res_ff.csv')

# EugSpe SM Valley stations inventory
data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] ,'model':eugspe_sm_nm_ls_val,
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
df_res_val.to_csv(resrat_ff_dir+'eugspe_sm_res_ff_val.csv')

###################

# SRV 1D Full stations inventory
data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'location':valley_array,'model':srv_1d_nm_ls,
             'f2_fas_srv_1d_res':f2_fas_srv_1d_res , 'f3_fas_srv_1d_res':f3_fas_srv_1d_res, 'f4_fas_srv_1d_res':f4_fas_srv_1d_res , 
             'f5_fas_srv_1d_res':f5_fas_srv_1d_res , 'f6_fas_srv_1d_res':f6_fas_srv_1d_res , 'f7_fas_srv_1d_res':f7_fas_srv_1d_res ,
             'pgv_srv_1d_res':pgv_srv_1d_res,
             
             'f2_fas_srv_1d_ratio':f2_fas_srv_1d_ratio , 'f3_fas_srv_1d_ratio':f3_fas_srv_1d_ratio, 'f4_fas_srv_1d_ratio':f4_fas_srv_1d_ratio , 
             'f5_fas_srv_1d_ratio':f5_fas_srv_1d_ratio , 'f6_fas_srv_1d_ratio':f6_fas_srv_1d_ratio , 'f7_fas_srv_1d_ratio':f7_fas_srv_1d_ratio ,
             
             'pgv_srv_1d_ratio':pgv_srv_1d_ratio,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             'pgv_steph_res':pgv_steph_res}

df_res = pd.DataFrame(data=data_res)
df_res.to_csv(resrat_ff_dir+'srv_1d_res_ff.csv')

# SRV 1D Valley stations inventory
data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] ,'model':srv_1d_nm_ls_val,
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
df_res_val.to_csv(resrat_ff_dir+'srv_1d_res_ff_val.csv')

###################

# SRV SM Full stations inventory
data_res = {'st_name':fas_obs[0] ,'st_lon':fas_obs[1], 'st_lat':fas_obs[2] , 'location':valley_array,'model':srv_sm_nm_ls,
             'f2_fas_srv_sm_res':f2_fas_srv_sm_res , 'f3_fas_srv_sm_res':f3_fas_srv_sm_res, 'f4_fas_srv_sm_res':f4_fas_srv_sm_res , 
             'f5_fas_srv_sm_res':f5_fas_srv_sm_res , 'f6_fas_srv_sm_res':f6_fas_srv_sm_res , 'f7_fas_srv_sm_res':f7_fas_srv_sm_res ,
             'pgv_srv_sm_res':pgv_srv_sm_res,
             
             'f2_fas_srv_sm_ratio':f2_fas_srv_sm_ratio , 'f3_fas_srv_sm_ratio':f3_fas_srv_sm_ratio, 'f4_fas_srv_sm_ratio':f4_fas_srv_sm_ratio , 
             'f5_fas_srv_sm_ratio':f5_fas_srv_sm_ratio , 'f6_fas_srv_sm_ratio':f6_fas_srv_sm_ratio , 'f7_fas_srv_sm_ratio':f7_fas_srv_sm_ratio ,
             
             'pgv_srv_sm_ratio':pgv_srv_sm_ratio,
             'f2_fas_steph_res':f2_fas_steph_res , 'f3_fas_steph_res':f3_fas_steph_res, 'f4_fas_steph_res':f4_fas_steph_res ,
             'f5_fas_steph_res':f5_fas_steph_res , 'f6_fas_steph_res':f6_fas_steph_res , 'f7_fas_steph_res':f7_fas_steph_res ,
             'pgv_steph_res':pgv_steph_res}

df_res = pd.DataFrame(data=data_res)
df_res.to_csv(resrat_ff_dir+'srv_sm_res_ff.csv')

# SRV SM Valley stations inventory
data_res_val = {'st_name_val':fas_val_obs[0] ,'st_lon_val':fas_val_obs[1], 'st_lat_val':fas_val_obs[2] ,'model':srv_sm_nm_ls_val,
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
df_res_val.to_csv(resrat_ff_dir+'srv_sm_res_ff_val.csv')

#%%
## Reading the newly made residuals flatfiles WITHOUT HEADER
# Full station inventory
steph_res_ff_nohead = pd.read_csv(resrat_ff_dir+'steph_res_ff.csv', header=None)
eugspe_1d_res_ff_nohead = pd.read_csv(resrat_ff_dir+'eugspe_1d_res_ff.csv', header=None)
eugspe_sm_res_ff_nohead = pd.read_csv(resrat_ff_dir+'eugspe_sm_res_ff.csv', header=None)
srv_1d_res_ff_nohead = pd.read_csv(resrat_ff_dir+'srv_1d_res_ff.csv', header=None)
srv_sm_res_ff_nohead = pd.read_csv(resrat_ff_dir+'srv_sm_res_ff.csv', header=None)

# Valley station inventory
valst_steph_res_ff_nohead = pd.read_csv(resrat_ff_dir+'steph_res_ff_val.csv', header=None)
valst_eugspe_1d_res_ff_nohead = pd.read_csv(resrat_ff_dir+'eugspe_1d_res_ff_val.csv', header=None)
valst_eugspe_sm_res_ff_nohead = pd.read_csv(resrat_ff_dir+'eugspe_sm_res_ff_val.csv', header=None)
valst_srv_1d_res_ff_nohead = pd.read_csv(resrat_ff_dir+'srv_1d_res_ff_val.csv', header=None)
valst_srv_sm_res_ff_nohead = pd.read_csv(resrat_ff_dir+'srv_sm_res_ff_val.csv', header=None)

## Writing ALL MODEL flatfiles
# Full station inventory
all_model_ff = pd.concat([eugspe_1d_res_ff_nohead.iloc[1:,1:20],eugspe_sm_res_ff_nohead.iloc[1:,1:20],
                          srv_1d_res_ff_nohead.iloc[1:,1:20],srv_sm_res_ff_nohead.iloc[1:,1:20],steph_res_ff_nohead.iloc[1:,1:20]], ignore_index=True)
all_model_ff.columns=['st_name', 'st_lon', 'st_lat', 'location', 'model', 
                          'f2_fas_res', 'f3_fas_res', 'f4_fas_res', 'f5_fas_res', 'f6_fas_res', 'f7_fas_res', 'pgv_res',
                          'f2_fas_ratio', 'f3_fas_ratio', 'f4_fas_ratio', 'f5_fas_ratio', 'f6_fas_ratio', 'f7_fas_ratio', 'pgv_ratio']
all_model_ff.to_csv(resrat_ff_dir+'all_model_ff.csv')

# Valley station inventory
all_model_val_ff = pd.concat([valst_eugspe_1d_res_ff_nohead.iloc[1:,1:19],valst_eugspe_sm_res_ff_nohead.iloc[1:,1:19],valst_srv_1d_res_ff_nohead.iloc[1:,1:19],
                              valst_srv_sm_res_ff_nohead.iloc[1:,1:19],valst_steph_res_ff_nohead.iloc[1:,1:19]],ignore_index=True)
all_model_val_ff.columns=['st_name', 'st_lon', 'st_lat', 'model', 
                          'f2_fas_res', 'f3_fas_res', 'f4_fas_res', 'f5_fas_res', 'f6_fas_res', 'f7_fas_res', 'pgv_res',
                          'f2_fas_ratio', 'f3_fas_ratio', 'f4_fas_ratio', 'f5_fas_ratio', 'f6_fas_ratio', 'f7_fas_ratio', 'pgv_ratio']
all_model_val_ff.to_csv(resrat_ff_dir+'all_model_val_ff.csv')

## Writing a MELTED flatfile from valley stations
melt_df = all_model_val_ff.melt(id_vars=['st_name','st_lon','st_lat','model'], 
                                value_vars=['f2_fas_res','f3_fas_res','f4_fas_res','f5_fas_res','f6_fas_res','f7_fas_res','pgv_res'],
                                var_name='res_im', value_name='res_val')
melt_df.to_csv(resrat_ff_dir+'melt_val_ff.csv')
















































