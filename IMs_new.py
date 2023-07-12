#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:42:53 2022

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import pandas as pd
import cartopy.crs as ccrs
import glob
import obspy
from geopy import distance

#%%

workingdir = '/Users/rshimony/Desktop/WillametteValley/Models/dist2depth_s_grav_n/'
obs_dir = '/Users/rshimony/Desktop/WillametteValley/val_event_data/observed_data_valevent/'
synthetic_dir = "/Users/rshimony/Desktop/WillametteValley/Models/dist2depth_s_grav_n/results/"
# synthetic_dir = "/Volumes/Roey's_BU/Results_dist2depth_model_novs30_nograv/lower_fmax_model/h100/results/"
stephenson_dir = '/Users/rshimony/Desktop/WillametteValley/4kmdepth_mw4_seismograms/'

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

xpnn = np.array([xpn])
xph1n = np.array([xph1])
xph2n = np.array([xph2])

ypnn = np.array([ypn])
yph1n = np.array([yph1])
yph2n = np.array([yph2])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

#%%

station_file = '/Users/rshimony/Desktop/WillametteValley/recorded_data/metadata/station_inventory_project.csv'

## Directories to be used in pulling data ##
# Need a list of earthquakes and stations
stndata = pd.read_csv(station_file, header = 0)
net = np.array(stndata["network"])[1:]
stnm = np.array(stndata["stName"])[1:]
# chan = np.array(stndata["chan"])
# loc = np.array(['*']* len(net))
slat = np.array(stndata["stlats"])[1:]
slon = np.array(stndata["stlon"])[1:]

#%%

def get_pgv(stream):
    
    '''
    Gives the PGV from a an obspy stream
    '''
    
    import numpy as np
    import obspy
    
    pgv_list = []
    
    for i in range(len(stream)):
        stream_copy = stream.copy()
        data = stream_copy[i].data
        pgv = np.max(np.abs(data))
        pgv_list.append(pgv)
        
    return(pgv_list) 

#%%

def arias_I(stream):
    
    '''
    Calculates Arias intensity from an obspy stream.
    Input is in velocity.
    Output is a list of Ia for all traces in stream and the times for easy plotting.
    '''
    
    import numpy as np
    import obspy
    from scipy import integrate
    
    g = 9.80665
    pi = np.pi
    a = (pi/(2*g))
    a_ints = []
    t_list = []
    max_ais = []
    samprate = stream[0].stats.sampling_rate
    npts = stream[0].stats.npts
    t = np.arange(0, npts / samprate, 1 / samprate)
    
    for i in range(len(stream)):
        stream_copy = stream.copy()
        acc = obspy.Stream.differentiate(stream_copy)
        data = acc[i].data
        acc2 = data**2
        acc2_int = integrate.cumtrapz(acc2)
        ia = a*acc2_int
        samprate = stream[i].stats.sampling_rate
        npts = len(ia)
        t = np.arange(0, npts / samprate, 1 / samprate)
        a_ints.append(ia)
        t_list.append(t)
        max_ai = np.max(np.abs(ia))
        max_ais.append(max_ai)
        
    return(a_ints , t_list , max_ais)

#%%

def make_obspy_from_specfem(specfem_stream_path,example_obspy_stream):
    
    '''
    obspy example should be a stream with 3 traces and the same as the specfem stream
    '''
    
    from obspy.signal.filter import bandpass
    
    specfem_stream = example_obspy_stream.copy()
    
    
    for i in range(3):
        trace2 = np.genfromtxt(specfem_stream_path[i])
        trspecfem = trace2.T[1]
        tspecfem = trace2.T[0]
        stats = obspy.core.trace.Stats()
        stats.sampling_rate = 100
        stats.npts = len(trspecfem)
        trt = obspy.core.trace.Trace(data = trspecfem, header = stats)
        trt = bandpass(trt, 0.1, 0.7*1, trt.stats.sampling_rate, corners=2, zerophase=True)
        
        
        specfem_stream[i].data = trt
        specfem_stream[i].times = tspecfem
        
    return(specfem_stream)

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
#%%

B030_syns = obspy.read(synthetic_dir + 'B030.*')
B032_syns = obspy.read(synthetic_dir + 'B032.*')
LANE_syns = obspy.read(synthetic_dir + 'LANE.*')
BUCK_syns = obspy.read(synthetic_dir + 'BUCK.*')
ALVY_syns = obspy.read(synthetic_dir + 'ALVY.*')
synts_streams = [LANE_syns,B030_syns,B032_syns,BUCK_syns,ALVY_syns]

B030_obs = obspy.read(obs_dir + '*.B030.EH*')
B032_obs = obspy.read(obs_dir + '*.B032.EH*')
LANE_obs = obspy.read(obs_dir + '*.LANE.EN*')
BUCK_obs = obspy.read(obs_dir + '*.BUCK.BH*')
ALVY_obs = obspy.read(obs_dir + '*.ALVY.EN*')
obs_streams = [LANE_obs,B030_obs,B032_obs,BUCK_obs,ALVY_obs]

B030_steph = glob.glob(stephenson_dir + '*.B030.CX*.semv')
B032_steph = glob.glob(stephenson_dir + '*.B032.CX*.semv')
LANE_steph = glob.glob(stephenson_dir + '*.LANE.CX*.semv')
BUCK_steph = glob.glob(stephenson_dir + '*.BUCK.CX*.semv')
ALVY_steph = glob.glob(stephenson_dir + '*.ALVY.CX*.semv')

B030_sf_stream = make_obspy_from_specfem(B030_steph , B030_obs)
B032_sf_stream = make_obspy_from_specfem(B032_steph , B032_obs)
LANE_sf_stream = make_obspy_from_specfem(LANE_steph , LANE_obs)
BUCK_sf_stream = make_obspy_from_specfem(BUCK_steph , BUCK_obs)
ALVY_sf_stream = make_obspy_from_specfem(ALVY_steph , ALVY_obs)
sf_streams = [LANE_sf_stream,B030_sf_stream,B032_sf_stream,BUCK_sf_stream,ALVY_sf_stream]

#%%

pgv_syns_ls = []
pgv_obs_ls = []
pgv_sf_ls = []
ia_syns_ls = []
ia_obs_ls = []
ia_sf_ls = []

for i in range(len(sf_streams)):
    pgv_syns = get_pgv(synts_streams[i])
    pgv_syns_avg = find_vec_norm(pgv_syns)
    pgv_syns_ls.append(pgv_syns_avg)
    
    pgv_obs = get_pgv(obs_streams[i])
    pgv_obs_avg = find_vec_norm(pgv_obs)
    pgv_obs_ls.append(pgv_obs_avg)
    
    pgv_sf = get_pgv(sf_streams[i])
    pgv_sf_avg = find_vec_norm(pgv_sf)
    pgv_sf_ls.append(pgv_sf_avg)

    ia_syns = arias_I(synts_streams[i])[2]
    ia_syns_avg = find_vec_norm(ia_syns)
    ia_syns_ls.append(ia_syns_avg)
    
    ia_obs = arias_I(obs_streams[i])[2]
    ia_obs_avg = find_vec_norm(ia_obs)
    ia_obs_ls.append(ia_obs_avg)
    
    ia_sf = arias_I(sf_streams[i])[2]
    ia_sf_avg = find_vec_norm(ia_sf)
    ia_sf_ls.append(ia_sf_avg)

#%%

res_pgv_syns = []
res_ia_syns = []
res_pgv_ref = []
res_ia_ref = []

for i in range(len(pgv_syns_ls)):
    
    pgv_res_syn = np.log(pgv_obs_ls[i]) - np.log(pgv_syns_ls[i])
    ia_res_syn = np.log(ia_obs_ls[i]) - np.log(pgv_obs_ls[i])
    res_pgv_syns.append(pgv_res_syn)
    res_ia_syns.append(ia_res_syn)
    
    pgv_res_ref = np.log(pgv_obs_ls[i]) - np.log(pgv_sf_ls[i])
    ia_res_ref = np.log(ia_obs_ls[i]) - np.log(ia_sf_ls[i]) 
    res_pgv_ref.append(pgv_res_ref)
    res_ia_ref.append(ia_res_ref)

#%%































