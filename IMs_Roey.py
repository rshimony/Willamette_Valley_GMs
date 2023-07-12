#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 14:39:42 2022

@author: rshimony
"""
#%%
import pandas as pd
import numpy as np
from shapely.geometry.polygon import Polygon

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

def get_pga(stream):
    
    '''
    Gives the PGA from a an obspy stream
    '''
    
    import numpy as np
    import obspy
    
    pga_list = []
    
    for i in range(len(stream)):
        stream_copy = stream.copy()
        acc = obspy.Stream.differentiate(stream_copy)
        data = acc[i].data
        pga = np.max(np.abs(data))
        pga_list.append(pga)
        
    return(pga_list)
        
        
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

def calc_xcorr(stream_1,stream_2):
    
    import numpy as np
    
    #cc = correlate_template(data_1,data_2, normalize='full') 
    
    xcorr_list = []
    
    for i in range(3):
        stream_copy_1 = stream_1.copy()
        stream_copy_2 = stream_2.copy()
        data_1 = stream_copy_1[i].data
        data_2 = stream_copy_2[i].data
        
        data_1 = (data_1 - np.mean(data_1)) / (np.std(data_1))
        data_2 = (data_2 - np.mean(data_2)) / (np.std(data_2))

        cc = np.correlate(data_1,data_2, 'full')/ max(len(data_1), len(data_2))
        indexx = np.argmax(cc)
        xcorr = round(cc[indexx], 4)
        xcorr_list.append(xcorr)

    return xcorr_list


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
# copmute rrup
def compute_rrup(evlon,evlat,evdepth,stlon,stlat,stelv):
    '''
    Compute Rrup given the event and station lon,lat,z positions - ONLY USE ON POINT SOURCES!
    Input:
        evlon:          Array with event longitudes (deg)
        evlat:          Array with event latitudes (deg)
        evdepth:        Array with event depths (km)
        stlon:          Array with station longitudes (deg)
        stlat:          Array with station latitudes (deg)
        stdepth:        Array with station depths (km)
    Output:
        Rrup:           Array with Rrup distances (km)
    '''
    
    from pyproj import Proj
    from numpy import sqrt
    
    #Convert the event and station latitude/longitude to UTM x and y:
    #Make the projection:
    p=Proj(proj='utm',zone='11S',ellps='WGS84',inverse=True)
    
    #Convert the latitude and longitude to UTM X and Y (in meters)    
    evx,evy=p(evlon,evlat)
    stx,sty=p(stlon,stlat)
    
    #Event and station depth is currently in km; convert to m:
    evz=evdepth*1000
    #stations are negative, have positive depth:
    stz=stelv*-1000
    
    #Get distance Rrup - closest distance to rupture.
    #   Since they are almost all point sources, this can just be site/event distance:
    Rrup=sqrt((evx-stx)**2 + (evy-sty)**2 + (evz-stz)**2)    
    
    #Convert back to km:
    Rrup=Rrup/1000
    
    #Return:
    return Rrup

#%%
# # Find arrival times
# def arr_time(stream):
#     '''

#     '''
#     from obspy.core import read
#     from obspy.signal.trigger import plot_trigger
#     from obspy.signal.trigger import classic_sta_lta
    



        
#%%    
import obspy
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd

stream = obspy.read('/Users/rshimony/Desktop/WillametteValley/val_event_data/observed_data_valevent/UW.LANE.EN*.sac')

ai = arias_I(stream)[2]

pga1 = get_pga(stream)

pgv1 = get_pgv(stream)

pga_norm = find_vec_norm(pga1)


#%%

workingdir = '/Users/rshimony/Desktop/WillametteValley/val_event_data/observed_data_valevent/'
synthetic_dir = '/Users/rshimony/Desktop/WillametteValley/Init_Model/talapas_results/no_refin/results/'
stephenson_dir = '/Users/rshimony/Desktop/WillametteValley/4kmdepth_mw4_seismograms/'

B030_syns = obspy.read(synthetic_dir + 'B030_syns.*')
B032_1_syns = obspy.read(synthetic_dir + 'B032_1_syns.*')
LANE_syns = obspy.read(synthetic_dir + 'LANE_syns.*')
BUCK_syns = obspy.read(synthetic_dir + 'BUCK_syns.*')
ALVY_syns = obspy.read(synthetic_dir + 'ALVY_syns.*')
synts_streams = [B030_syns,B032_1_syns,LANE_syns,BUCK_syns,ALVY_syns]

B030_obs = obspy.read(workingdir + '*.B030.EH*')
B032_obs = obspy.read(workingdir + '*.B032.EH*')
LANE_obs = obspy.read(workingdir + '*.LANE.EN*')
BUCK_obs = obspy.read(workingdir + '*.BUCK.BH*')
ALVY_obs = obspy.read(workingdir + '*.ALVY.EN*')
obs_streams = [B030_obs,B032_obs,LANE_obs,BUCK_obs,ALVY_obs]

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
sf_streams = [B030_sf_stream,B032_sf_stream,LANE_sf_stream,BUCK_sf_stream,ALVY_sf_stream]

#%%
for i in range(len(sf_streams)):
    for j in range(len(B030_sf_stream)):
        sf_streams[i][j].write('/Users/rshimony/Desktop/WillametteValley/IMs/' + sf_streams[i][j].id + '_sf.SAC' , format='SAC')


#%%

B027_syns = obspy.read(synthetic_dir + 'B027_8_syns.*')
B027_obs = obspy.read(workingdir + '*.B027.EH*')

ai_27_syn_T = find_vec_norm(arias_I(B027_syns)[2])
ai_27_obs_T = find_vec_norm(arias_I(B027_obs)[2])

pga_27_syn_T = find_vec_norm(get_pga(B027_syns))
pga_27_obs_T = find_vec_norm(get_pga(B027_obs))

pgv_27_syn_T = find_vec_norm(get_pgv(B027_syns))
pgv_27_obs_T = find_vec_norm(get_pgv(B027_obs))

ai_res_syn_obs = np.log(ai_27_obs_T) - np.log(ai_27_syn_T)
pga_res_syn_obs = np.log(pga_27_obs_T) - np.log(pga_27_syn_T)
pgv_res_syn_obs = np.log(pgv_27_obs_T) - np.log(pgv_27_syn_T)

xcorr_syn_obs_T =find_vec_norm(calc_xcorr(B027_syns,B027_obs))

print(ai_res_syn_obs , pga_res_syn_obs , pgv_res_syn_obs , xcorr_syn_obs_T)
data_27 = [ai_res_syn_obs , pga_res_syn_obs , pgv_res_syn_obs , xcorr_syn_obs_T]

df_res_27 = pd.DataFrame(data=data_27)
   
df_res_27.to_csv('/Users/rshimony/Desktop/WillametteValley/IMs/IMs_residuals_27.csv')

#%%

from obspy.signal.trigger import plot_trigger
from obspy.signal.trigger import classic_sta_lta  

trace = B030_syns[0]
df = trace.stats.sampling_rate
cft = classic_sta_lta(trace.data, int(6 * df), int(10 * df))
plot_trigger(trace, cft, 1.5, 0.5)

#%%
arr_file = '/Users/rshimony/Desktop/WillametteValley/IMs/arrivals.csv' 
 
arr_times = pd.read_csv(arr_file)

times_syn = arr_times['Times_synts']
times_obs = arr_times['Time_obs']
times_sf = arr_times['Time_sf']

times_syn_T = [np.average(times_syn[0:3]) , np.average(times_syn[3:6]) , np.average(times_syn[6:9]) , np.average(times_syn[9:12]) , np.average(times_syn[12:15]) , np.average(times_syn[15:18])]
times_obs_T = [np.average(times_obs[0:3]) , np.average(times_obs[3:6]) , np.average(times_obs[6:9]) , np.average(times_obs[9:12]) , np.average(times_obs[12:15]) , np.average(times_obs[15:18])]
times_sf_T = [np.average(times_sf[0:3]) , np.average(times_sf[3:6]) , np.average(times_sf[6:9]) , np.average(times_sf[9:12]) , np.average(times_sf[12:15])]

arrs_res_syn = []
arrs_res_sf = []

for i in range(len(times_syn_T)):
    arrs_res_syn_i = times_obs_T[i] - times_syn_T[i]
    arrs_res_syn.append(arrs_res_syn_i)
    if i < 5:
        arrs_res_sf_i = times_obs_T[i] - times_sf_T[i]
        arrs_res_sf.append(arrs_res_sf_i)
    
arrs_27 =  [arrs_res_syn[-1]]

df_arrs_27 = pd.DataFrame(data=arrs_27)
   
df_arrs_27.to_csv('/Users/rshimony/Desktop/WillametteValley/IMs/arrs_residuals_27.csv')
#%%

synts_ims = np.zeros((len(synts_streams) , 9)) #second number is the amount of Ims x components (3 Ims x 3 comp(E,N,Z))
obs_ims = np.zeros((len(synts_streams) , 9))
sf_ims = np.zeros((len(synts_streams) , 9))

synts_t = np.zeros((len(synts_streams) , 3)) #second number is the amount of Ims (total)
obs_t = np.zeros((len(synts_streams) , 3))
sf_t = np.zeros((len(synts_streams) , 3))


for i in range(len(synts_streams)):
    synts_ims[i][0] = arias_I(synts_streams[i])[2][0]
    synts_ims[i][1] = arias_I(synts_streams[i])[2][1]
    synts_ims[i][2] = arias_I(synts_streams[i])[2][2]
    synts_ims[i][3] = get_pga(synts_streams[i])[0]
    synts_ims[i][4] = get_pga(synts_streams[i])[1]
    synts_ims[i][5] = get_pga(synts_streams[i])[2]
    synts_ims[i][6] = get_pgv(synts_streams[i])[0]
    synts_ims[i][7] = get_pgv(synts_streams[i])[1]
    synts_ims[i][8] = get_pgv(synts_streams[i])[2]
    obs_ims[i][0] = arias_I(obs_streams[i])[2][0]
    obs_ims[i][1] = arias_I(obs_streams[i])[2][1]
    obs_ims[i][2] = arias_I(obs_streams[i])[2][2]
    obs_ims[i][3] = get_pga(obs_streams[i])[0]
    obs_ims[i][4] = get_pga(obs_streams[i])[1]
    obs_ims[i][5] = get_pga(obs_streams[i])[2]
    obs_ims[i][6] = get_pgv(obs_streams[i])[0]
    obs_ims[i][7] = get_pgv(obs_streams[i])[1]
    obs_ims[i][8] = get_pgv(obs_streams[i])[2]
    sf_ims[i][0] = arias_I(sf_streams[i])[2][0]
    sf_ims[i][1] = arias_I(sf_streams[i])[2][1]
    sf_ims[i][2] = arias_I(sf_streams[i])[2][2]
    sf_ims[i][3] = get_pga(sf_streams[i])[0]
    sf_ims[i][4] = get_pga(sf_streams[i])[1]
    sf_ims[i][5] = get_pga(sf_streams[i])[2]
    sf_ims[i][6] = get_pgv(sf_streams[i])[0]
    sf_ims[i][7] = get_pgv(sf_streams[i])[1]
    sf_ims[i][8] = get_pgv(sf_streams[i])[2]
    
    synts_t[i][0] = find_vec_norm(arias_I(synts_streams[i])[2])
    synts_t[i][1] = find_vec_norm(get_pga(synts_streams[i]))
    synts_t[i][2] = find_vec_norm(get_pgv(synts_streams[i]))
    obs_t[i][0] = find_vec_norm(arias_I(obs_streams[i])[2])
    obs_t[i][1] = find_vec_norm(get_pga(obs_streams[i]))
    obs_t[i][2] = find_vec_norm(get_pgv(obs_streams[i]))
    sf_t[i][0] = find_vec_norm(arias_I(sf_streams[i])[2])
    sf_t[i][1] = find_vec_norm(get_pga(sf_streams[i]))
    sf_t[i][2] = find_vec_norm(get_pgv(sf_streams[i]))
    
res_syn_obs = np.zeros((len(synts_streams) , 9))
res_sf_obs = np.zeros((len(synts_streams) , 9))

res_syn_obs_T = np.zeros((len(synts_streams) , 3))
res_sf_obs_T = np.zeros((len(synts_streams) , 3))

xcorr_syn_obs = np.zeros((len(synts_streams) , 3))
xcorr_sf_obs = np.zeros((len(synts_streams) , 3))

xcorr_syn_obs_T = np.zeros(len(synts_streams))
xcorr_sf_obs_T = np.zeros(len(synts_streams))


for i in range(len(synts_streams)):
    for j in range(9):
        res_syn_obs[i][j] = np.log(obs_ims[i][j]) - np.log(synts_ims[i][j])
        res_sf_obs[i][j] = np.log(obs_ims[i][j]) - np.log(sf_ims[i][j])
        
for i in range(len(synts_streams)):
    for j in range(3):
        res_syn_obs_T[i][j] = np.log(obs_t[i][j]) - np.log(synts_t[i][j])
        res_sf_obs_T[i][j] = np.log(obs_t[i][j]) - np.log(sf_t[i][j])
        
for i in range(len(synts_streams)):
    for j in range(3):
        xcorr_syn_obs[i][j] = calc_xcorr(synts_streams[i],obs_streams[i])[j]
        xcorr_sf_obs[i][j] = calc_xcorr(sf_streams[i],obs_streams[i])[j]

for i in range(len(synts_streams)):
    xcorr_syn_obs_T[i] =find_vec_norm(calc_xcorr(synts_streams[i],obs_streams[i]))
    xcorr_sf_obs_T[i] = find_vec_norm(calc_xcorr(sf_streams[i],obs_streams[i]))
       
stationfile = '/Users/rshimony/Desktop/WillametteValley/eve_4_inv/metadata/WV_station_inventory.csv' 
 
stations = pd.read_csv(stationfile)

Net = stations['network']
Name = stations['name']
slat = stations['latitude']
slon = stations['longitude']

sta_lon = []
sta_lat = []
nets = []
sta_name = []

for i in range(len(obs_streams)):
    for j in range(len(Name)):
        if obs_streams[i][0].stats.station == Name[j]:
            sta_lon.append(slon[j])
            sta_lat.append(slat[j])
            nets.append(Net[j])
            sta_name.append(Name[j])
            
            
source_lat = 44.090 #N
source_lon = -122.831 #°W 
source_dep = 4.0 #km

rhyp = []
for i in range(len(sta_name)):
    R = compute_rrup(source_lon,source_lat,source_dep,sta_lon[i],sta_lat[i],0.0)
    rhyp.append(R)
    
    
dict_res = {'net': nets[0], 'station':sta_name[0], 'lat':sta_lat[0],'lon':sta_lon[0],'rhyp':rhyp[0]
          ,'res_syn_Ia_N':res_syn_obs[0][0],'res_syn_Ia_E':res_syn_obs[0][1],'res_syn_Ia_Z':res_syn_obs[0][2],'res_syn_Ia_T':res_syn_obs_T[0][0],
          'res_syn_pga_N':res_syn_obs[0][3],'res_syn_pga_E':res_syn_obs[0][4],'res_syn_pga_Z':res_syn_obs[0][5],'res_syn_pga_T':res_syn_obs_T[0][1],
          'res_syn_pgv_N':res_syn_obs[0][6],'res_syn_pgv_E':res_syn_obs[0][7],'res_syn_pgv_Z':res_syn_obs[0][8],'res_syn_pgv_T':res_syn_obs_T[0][2],
          'res_ref_Ia_N':res_sf_obs[0][0],'res_ref_Ia_E':res_sf_obs[0][1],'res_ref_Ia_Z':res_sf_obs[0][2],'res_ref_Ia_T':res_sf_obs_T[0][0],
          'res_ref_pga_N':res_sf_obs[0][3],'res_ref_pga_E':res_sf_obs[0][4],'res_ref_pga_Z':res_sf_obs[0][5],'res_ref_pga_T':res_sf_obs_T[0][1],
          'res_ref_pgv_N':res_sf_obs[0][6],'res_ref_pgv_E':res_sf_obs[0][7],'res_ref_pgv_Z':res_sf_obs[0][8],'res_ref_pgv_T':res_sf_obs_T[0][2],
          'xcorr_syn_N':xcorr_syn_obs[0][0],'xcorr_syn_E':xcorr_syn_obs[0][1],'xcorr_syn_Z':xcorr_syn_obs[0][2],'xcorr_syn_T':xcorr_syn_obs_T[0],
          'xcorr_ref_N':xcorr_sf_obs[0][0],'xcorr_ref_E':xcorr_sf_obs[0][1],'xcorr_ref_Z':xcorr_sf_obs[0][2],'xcorr_ref_T':xcorr_sf_obs_T[0],
          'arrs_res_syn_T':arrs_res_syn[0],'arrs_res_ref_T':arrs_res_sf[0]}
    
df_res = pd.DataFrame(data=dict_res,index=[0])

for i in range(1,len(sta_name)):

    dict_res_temp = {'net': nets[i], 'station':sta_name[i], 'lat':sta_lat[i],'lon':sta_lon[i],'rhyp':rhyp[i]
              ,'res_syn_Ia_N':res_syn_obs[i][0],'res_syn_Ia_E':res_syn_obs[i][1],'res_syn_Ia_Z':res_syn_obs[i][2],'res_syn_Ia_T':res_syn_obs_T[i][0],
              'res_syn_pga_N':res_syn_obs[i][3],'res_syn_pga_E':res_syn_obs[i][4],'res_syn_pga_Z':res_syn_obs[i][5],'res_syn_pga_T':res_syn_obs_T[i][1],
              'res_syn_pgv_N':res_syn_obs[i][6],'res_syn_pgv_E':res_syn_obs[i][7],'res_syn_pgv_Z':res_syn_obs[i][8],'res_syn_pgv_T':res_syn_obs_T[i][2],
              'res_ref_Ia_N':res_sf_obs[i][0],'res_ref_Ia_E':res_sf_obs[i][1],'res_ref_Ia_Z':res_sf_obs[i][2],'res_ref_Ia_T':res_sf_obs_T[i][0],
              'res_ref_pga_N':res_sf_obs[i][3],'res_ref_pga_E':res_sf_obs[i][4],'res_ref_pga_Z':res_sf_obs[i][5],'res_ref_pga_T':res_sf_obs_T[i][1],
              'res_ref_pgv_N':res_sf_obs[i][6],'res_ref_pgv_E':res_sf_obs[i][7],'res_ref_pgv_Z':res_sf_obs[i][8],'res_ref_pgv_T':res_sf_obs_T[i][2],
              'xcorr_syn_N':xcorr_syn_obs[i][0],'xcorr_syn_E':xcorr_syn_obs[i][1],'xcorr_syn_Z':xcorr_syn_obs[i][2],'xcorr_syn_T':xcorr_syn_obs_T[i],
              'xcorr_ref_N':xcorr_sf_obs[i][0],'xcorr_ref_E':xcorr_sf_obs[i][1],'xcorr_ref_Z':xcorr_sf_obs[i][2],'xcorr_ref_T':xcorr_sf_obs_T[i],
              'arrs_res_syn_T':arrs_res_syn[i],'arrs_res_ref_T':arrs_res_sf[i]}
    
    df_res_temp = pd.DataFrame(data=dict_res_temp , index=[i])

    df_res = pd.concat([df_res, df_res_temp], ignore_index=True)
    
df_res.to_csv('/Users/rshimony/Desktop/WillametteValley/IMs/IMs_residuals_ln_2.csv',index=False)





#%%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import matplotlib.patches as patches
from cartopy.io.img_tiles import GoogleTiles
class ShadedReliefESRI(GoogleTiles):
    # shaded relief
    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
               'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
               z=z, y=y, x=x)
        return url

#%%

flatfile_path = '/Users/rshimony/Desktop/WillametteValley/IMs/IMs_residuals_ln_2.csv'
im_ff = pd.read_csv(flatfile_path)

network = im_ff['net']
station = im_ff['station']
lats = im_ff['lat']
lons = im_ff['lon']
rhypo = im_ff['rhyp']
res_syn_Ia_T = im_ff['res_syn_Ia_T']
res_syn_pga_T = im_ff['res_syn_pga_T']
res_syn_pgv_T = im_ff['res_syn_pgv_T']
res_ref_Ia_T = im_ff['res_ref_Ia_T']
res_ref_pga_T = im_ff['res_ref_pga_T']
res_ref_pgv_T = im_ff['res_ref_pgv_T']
xcorr_syn_T = im_ff['xcorr_syn_T']
xcorr_ref_T = im_ff['xcorr_ref_T']
xcorr_syn_T = im_ff['xcorr_syn_T']
xcorr_ref_T = im_ff['xcorr_ref_T']
arrs_syn_T = im_ff['arrs_res_syn_T']
arrs_ref_T = im_ff['arrs_res_ref_T']

########### Continue from here!!!!
'''
need to:
    1) change the 4 panel res plots to show the arrs instead of xcorr
    2) plot same 4 panel res plots for ref as well
    3) flip the res calc so it will be obs-pred
    4) plot single panel for xcorr
    5) plot ratios (4 panel for all but xcorr, xcorr alone - ???)                  
'''

res_syn_list = [res_syn_pga_T , res_syn_pgv_T , arrs_syn_T , res_syn_Ia_T]
res_ref_list = [res_ref_pga_T , res_ref_pgv_T , arrs_ref_T]
res_names = ['PGA Residuals' , 'PGV Residuals' , 'Arrivals Residuals' , 'Arias_I Residuals' , 'X-Correlation Residuals' , 'X-Correlation Residuals']


valley_poly = np.array([[-123.2641219775928,44.64634967360506],
               [-123.3774486537867,44.5216826327464],
               [-123.1301667224165,43.49355895349841],
               [-123.0758323132938,43.49353022525693],
               [-122.9903835522051,44.45484968462699],
               [-122.7948953662283,44.64386451064287],
               [-123.2641219775928,44.64634967360506]])

# Ia_ratios = np.abs(res_syn_Ia_T[1:] / res_ref_Ia_T[1:])
pga_ratios = np.abs(res_syn_pga_T[1:] / res_ref_pga_T[1:])
pgv_ratios = np.abs(res_syn_pgv_T[1:] / res_ref_pgv_T[1:])
# xcorr_ratios = np.abs(xcorr_syn_T[1:] / xcorr_ref_T[1:])
arrs_ratios = np.abs(arrs_syn_T[1:] / arrs_ref_T[1:])

ratios_list = [pga_ratios , pgv_ratios , arrs_ratios , arrs_ratios]
ratios_names = ['PGA Ratios' , 'PGV Ratios' , 'Arrivals Ratios' , 'Arrivals Ratios']
lon_ratios = lons[1:]
lat_ratios = lats[1:]

#%%

figw, figh = 20.0, 15.0
fig, axs = plt.subplots(2,2,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(figw,figh))
# fig.subplots_adjust(wspace=0.5,hspace=0.1)
plt.subplots_adjust(left=1/figw, right=1-1/figw, bottom=1/figh, top=1-1/figh)

axs=axs.flatten()

for i in range(4):
    if i == 3:
        axs[i].set_axis_off()
        axs[i].patch.set_visible(False)
        axs[i].axis('off')


    else:
        
        axs[i].set_extent([-124.19,-121.51,43.4,46.1])
        axs[i].add_image(ShadedReliefESRI(), 10)
                
        gl_major = axs[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
        # gl_major.xlocator = mticker.FixedLocator(np.arange(-124,-121.5,0.5))
        # gl_major.ylocator = mticker.FixedLocator(np.arange(43.5,46,0.5))
        gl_major.xlabel_style = {'size': 7}
        gl_major.ylabel_style = {'size': 7}
                
        scat=axs[i].scatter(lon_ratios,lat_ratios,s=70,c=ratios_list[i], cmap='seismic',marker='^',alpha=0.7,transform=ccrs.PlateCarree(),edgecolors='k', vmin=0 , vmax = 2)
        
    
        
        
        # rect = patches.Rectangle(xy=[-123.65, 43.49], width=1.15, height=1.16,linewidth=2.5,facecolor='none'
        #                      ,edgecolor='grey',linestyle='--',transform=ccrs.PlateCarree())
        # axs[i].add_patch(rect)
        
        # valley_borders = patches.Polygon(xy=valley_poly,facecolor='none',edgecolor='purple',linewidth=2.5,linestyle='--',
                                         # alpha=0.5,transform=ccrs.PlateCarree(),zorder=4)
        # axs[i].add_patch(valley_borders)
        
        axs[i].plot(xpn,ypn, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
        axs[i].plot(xph1,yph1, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
        axs[i].plot(xph2,yph2, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
                
        axs[i].set_title(ratios_names[i] , pad= 22)
        
        cb=plt.colorbar(scat,fraction=0.022, pad=0.07 , ax=axs[i])
        cb.set_alpha(1)
        cb.set_label(ratios_names[i])
        cb.draw_all()
        fig.suptitle('Residuals Ratios: Model Residuals / Reference Residuals', fontsize=18)

         
plt.savefig('/Users/rshimony/Desktop/WillametteValley/IMs/Ims_ratio_maps.png',dpi=300,bbox_inches='tight') 


plt.close('all')


figw, figh = 20.0, 15.0
fig, axs = plt.subplots(2,2,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(figw,figh), sharey='row' , sharex='col')
# fig.subplots_adjust(wspace=0.5,hspace=0.1)
plt.subplots_adjust(wspace = 0.4 ,hspace = 0.1, right=0.8)

axs=axs.flatten()

for i in range(4):

    axs[i].set_extent([-124.19,-121.51,43.4,46.1])
    axs[i].add_image(ShadedReliefESRI(), 10)
            
    gl_major = axs[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
    # gl_major.xlocator = mticker.FixedLocator(np.arange(-124,-121.5,0.5))
    # gl_major.ylocator = mticker.FixedLocator(np.arange(43.5,46,0.5))
    gl_major.xlabel_style = {'size': 7}
    gl_major.ylabel_style = {'size': 7}
            
    scat=axs[i].scatter(lons,lats,s=70,c=res_syn_list[i], cmap='seismic',marker='^',alpha=0.7,transform=ccrs.PlateCarree(),edgecolors='k' , vmin=-2 , vmax=2)
    

    
    
    # rect = patches.Rectangle(xy=[-123.65, 43.49], width=1.15, height=1.16,linewidth=2.5,facecolor='none'
    #                      ,edgecolor='grey',linestyle='--',transform=ccrs.PlateCarree())
    # axs[i].add_patch(rect)
    
    # valley_borders = patches.Polygon(xy=valley_poly,facecolor='none',edgecolor='purple',linewidth=2.5,linestyle='--',
    #                                  alpha=0.5,transform=ccrs.PlateCarree(),zorder=4)
    # axs[i].add_patch(valley_borders)
    
    axs[i].plot(xpn,ypn, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
    axs[i].plot(xph1,yph1, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
    axs[i].plot(xph2,yph2, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
            
    axs[i].set_title(res_names[i] , pad= 22)
    
    # cb=plt.colorbar(scat,fraction=0.022, pad=0.07 , ax=axs[i])
    # cb.set_alpha(1)
    # cb.set_label(res_names[i])
    # cb.draw_all()
    fig.suptitle('Residual Plots: Observed - Predicted (Model)', fontsize=18)
        
# plt.title('Residual Plots: Observed - Predicted (Initial Model)')
       
cbar_ax = fig.add_axes([0.85, 0.35, 0.02, 0.4])
cbar=fig.colorbar(scat, cax=cbar_ax)
cbar.set_label('ln Residuals',size='13')
cbar.ax.tick_params(labelsize=11)     
plt.savefig('/Users/rshimony/Desktop/WillametteValley/IMs/Ims_maps.png',dpi=300,bbox_inches='tight') 


plt.close('all')





figw, figh = 20.0, 15.0
fig, axs = plt.subplots(2,2,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(figw,figh), sharey='row' , sharex='col')
# fig.subplots_adjust(wspace=0.5,hspace=0.1)
plt.subplots_adjust(wspace = 0.4 ,hspace = 0.1, right=0.8)

axs=axs.flatten()

for i in range(4):
    if i == 3:
        axs[i].set_axis_off()
        axs[i].patch.set_visible(False)
        axs[i].axis('off')


    else:
    
        axs[i].set_extent([-124.19,-121.51,43.4,46.1])
        axs[i].add_image(ShadedReliefESRI(), 10)
                
        gl_major = axs[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
        # gl_major.xlocator = mticker.FixedLocator(np.arange(-124,-121.5,0.5))
        # gl_major.ylocator = mticker.FixedLocator(np.arange(43.5,46,0.5))
        gl_major.xlabel_style = {'size': 7}
        gl_major.ylabel_style = {'size': 7}
                
        scat=axs[i].scatter(lons[1:],lats[1:],s=70,c=res_ref_list[i][1:], cmap='seismic',marker='^',alpha=0.7,transform=ccrs.PlateCarree(),edgecolors='k' , vmin=-2 , vmax=2)
        
    
        
        
        # rect = patches.Rectangle(xy=[-123.65, 43.49], width=1.15, height=1.16,linewidth=2.5,facecolor='none'
        #                      ,edgecolor='grey',linestyle='--',transform=ccrs.PlateCarree())
        # axs[i].add_patch(rect)
        
        # valley_borders = patches.Polygon(xy=valley_poly,facecolor='none',edgecolor='purple',linewidth=2.5,linestyle='--',
        #                                  alpha=0.5,transform=ccrs.PlateCarree(),zorder=4)
        # axs[i].add_patch(valley_borders)
        
        axs[i].plot(xpn,ypn, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1,zorder=4)
        axs[i].plot(xph1,yph1, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1,zorder=4)
        axs[i].plot(xph2,yph2, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1,zorder=4)                

        axs[i].set_title(res_names[i] , pad= 22)
        
        # cb=plt.colorbar(scat,fraction=0.022, pad=0.07 , ax=axs[i])
        # cb.set_alpha(1)
        # cb.set_label(res_names[i])
        # cb.draw_all()
        fig.suptitle('Residual Plots: Observed - Predicted (Reference Model)', fontsize=18)
 
cbar_ax = fig.add_axes([0.85, 0.35, 0.02, 0.4])
cbar=fig.colorbar(scat, cax=cbar_ax)
cbar.set_label('ln Residuals',size='13')
cbar.ax.tick_params(labelsize=11) 
           
plt.savefig('/Users/rshimony/Desktop/WillametteValley/IMs/Ims_maps_ref.png',dpi=300,bbox_inches='tight') 


plt.close('all')









