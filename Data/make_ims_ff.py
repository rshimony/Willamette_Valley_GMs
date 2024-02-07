#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 16:29:52 2023

@author: rshimony
"""

import matplotlib.pyplot as plt
import numpy as np 
import obspy as obs
import pandas as pd
from glob import glob

#%%

st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/snr_stations_fixed.csv')
st_names = np.array(st_file['st_nm'])
st_lons = np.array(st_file['st_lon'])
st_lats = np.array(st_file['st_lat']) 
st_names_unq = pd.unique(st_names)
st_lons_unq = pd.unique(st_lons)
st_lats_unq = pd.unique(st_lats)

valst_file = '/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/valley_snr_st.csv'

obs_dir = '/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/waveforms/observed_data_1hz_snr/'
steph_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/scottsmills/rfile_steph_scottsmills_preempt_freq25/'
synts_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/scottsmills/rfile_steph_basin_scottsmills_eugspe_1d_preempt_freq25/'

outdir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/ims/scottsmills/eugspe_1d/'
ff_name = 'ims_ff_eugspe_1d.csv'

obs_trs = sorted(glob(obs_dir + '*_unfilt.sac'))
# obs_trs_filt = sorted(glob(obs_dir + '*_filt.sac'))
synts_trs = sorted(glob(synts_dir + '*.*v'))
steph_trs = sorted(glob(steph_dir + '*.*v'))

obs_nms = []
for i in range(len(obs_trs)):
    obs_nm = obs_trs[i].split('.')[1]
    obs_nms.append(obs_nm)
obs_nms_unq = np.unique(obs_nms)
    
synt_nms = []
for i in range(len(synts_trs)):   
    synt_nm = (synts_trs[i].split('/')[-1]).split('.')[0]
    synt_nms.append(synt_nm)
synt_nms_unq = np.unique(synt_nms)

steph_nms = []
for i in range(len(steph_trs)):   
    steph_nm = (steph_trs[i].split('/')[-1]).split('.')[0]
    steph_nms.append(steph_nm)
steph_nms_unq = np.unique(steph_nms)

obs_strs = []
obs_strs_filt = []
synts_strs = []
steph_strs = []

for i in range(len(st_names_unq)):
    if st_names_unq[i] in obs_nms_unq:
        obs_str = obs.read(obs_dir + '*' + st_names_unq[i] + '*_unfilt.sac')
        obs_str_cp = obs_str.copy()
        obs_strs.append(obs_str_cp)
        obs_str_filt = obs_str.copy()
        # obs_str_filt.filter('lowpass' ,freq=1.0, corners=4, zerophase=True)
        df_obs = obs_str_filt[0].stats.sampling_rate
        obs_str_filt.filter('lowpass_cheby_2',freq=1.0)
        obs_strs_filt.append(obs_str_filt)
    
    # if st_names_unq[i] in obs_nms_unq:
    #     obs_str_filt = obs.read(obs_dir + '*' + st_names_unq[i] + '*_filt.sac')
    #     obs_str_filt_cp = obs_str_filt.copy()
    #     obs_strs_filt.append(obs_str_filt_cp)
    
    if st_names_unq[i] in synt_nms_unq:
        synts_str = obs.read(synts_dir + st_names_unq[i] + '.*v')
        synts_str_cp = synts_str.copy()
        df_synts = synts_str_cp[0].stats.sampling_rate
        # synts_str_filt = synts_str_cp.filter('lowpass' ,freq=1.0, corners=4, zerophase=True)
        synts_str_filt = synts_str_cp.filter('lowpass_cheby_2',freq=1.0)
        synts_strs.append(synts_str_filt)
        
    if st_names_unq[i] in steph_nms_unq:
        steph_str = obs.read(steph_dir + st_names_unq[i] + '.*v')
        steph_str_cp = steph_str.copy()
        df_steph = steph_str_cp[0].stats.sampling_rate
        # steph_str_filt = steph_str_cp.filter('lowpass', freq=1.0, corners=4, zerophase=True)
        steph_str_filt = steph_str_cp.filter('lowpass_cheby_2',freq=1.0)
        steph_strs.append(steph_str_filt)


obs_strs_3 = []
obs_strs_filt_3 = []
synts_strs_3 = []
steph_strs_3 = []

for i in obs_strs:
    if len(i) >= 3:
       obs_strs_3.append(i[-3:]) 
       
for i in obs_strs_filt:
    if len(i) >= 3:
       obs_strs_filt_3.append(i[-3:]) 
        
for i in synts_strs:
    if len(i) >= 3:
       synts_strs_3.append(i[-3:]) 
       
for i in steph_strs:
    if len(i) >= 3:
       steph_strs_3.append(i[-3:]) 

#%%

# val_st_file = pd.read_csv(valst_file)
# dom_valst_names = np.array(val_st_file['dom_valst_names'])
# dom_valst_lons = np.array(val_st_file['dom_valst_lons'])
# dom_valst_lats = np.array(val_st_file['dom_valst_lats']) 

# valst_obs_strs_3 = []
# valst_obs_strs_filt_3 = []
# valst_synts_strs_3 = []
# valst_steph_strs_3 = []

# for i in range(len(obs_strs_3)):
#     for j in range(len(dom_valst_names)):
        
#         obs_str_name = obs_strs_3[i][0].stats.station
        
#         if obs_str_name == dom_valst_names[j]:
#             valst_obs_strs_3.append(obs_strs_3[i])
#             valst_obs_strs_filt_3.append(obs_strs_filt_3[i])
#             valst_synts_strs_3.append(synts_strs_3[i])
#             valst_steph_strs_3.append(steph_strs_3[i])
            
        

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

def find_vec_norm(Im_value_list):
    '''

    '''
    import numpy as np

    comp_1 = Im_value_list[0]**2;
    comp_2 = Im_value_list[1]**2;
    comp_3 = Im_value_list[2]**2;

    vec_norm = np.sqrt(comp_1 + comp_2 + comp_3)

    return(vec_norm)

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

for i in range(len(synts_strs_3)):
    
    pgv_e_synt = get_pgv(synts_strs_3[i])[0]
    pgv_n_synt = get_pgv(synts_strs_3[i])[1]
    pgv_z_synt = get_pgv(synts_strs_3[i])[2]
    pgv_t_synt = find_vec_norm([pgv_e_synt,pgv_n_synt,pgv_z_synt])

    ai_e_synt = arias_I(synts_strs_3[i])[2][0]
    ai_n_synt = arias_I(synts_strs_3[i])[2][1]
    ai_z_synt = arias_I(synts_strs_3[i])[2][2]
    ai_t_synt = find_vec_norm([ai_e_synt,ai_n_synt,ai_z_synt])

    pgv_e_obs = get_pgv(obs_strs_filt_3[i])[0]
    pgv_n_obs = get_pgv(obs_strs_filt_3[i])[1]
    pgv_z_obs = get_pgv(obs_strs_filt_3[i])[2]
    pgv_t_obs = find_vec_norm([pgv_e_obs,pgv_n_obs,pgv_z_obs])

    ai_e_obs = arias_I(obs_strs_filt_3[i])[2][0]
    ai_n_obs = arias_I(obs_strs_filt_3[i])[2][1]
    ai_z_obs = arias_I(obs_strs_filt_3[i])[2][2]
    ai_t_obs = find_vec_norm([ai_e_obs,ai_n_obs,ai_z_obs])
    
    pgv_e_steph = get_pgv(steph_strs_3[i])[0]
    pgv_n_steph = get_pgv(steph_strs_3[i])[1]
    pgv_z_steph = get_pgv(steph_strs_3[i])[2]
    pgv_t_steph = find_vec_norm([pgv_e_steph,pgv_n_steph,pgv_z_steph])
    
    ai_e_steph = arias_I(steph_strs_3[i])[2][0]
    ai_n_steph = arias_I(steph_strs_3[i])[2][1]
    ai_z_steph = arias_I(steph_strs_3[i])[2][2]
    ai_t_steph = find_vec_norm([ai_e_steph,ai_n_steph,ai_z_steph])
    
    xcorr_e_synt = calc_xcorr((synts_strs_3[i]),(obs_strs_filt_3[i]))[0]
    xcorr_n_synt = calc_xcorr((synts_strs_3[i]),(obs_strs_filt_3[i]))[1]
    xcorr_z_synt = calc_xcorr((synts_strs_3[i]),(obs_strs_filt_3[i]))[2]
    xcorr_t_synt = find_vec_norm([xcorr_e_synt,xcorr_n_synt,xcorr_z_synt])
    
    xcorr_e_steph = calc_xcorr((steph_strs_3[i]),(obs_strs_filt_3[i]))[0]
    xcorr_n_steph = calc_xcorr((steph_strs_3[i]),(obs_strs_filt_3[i]))[1]
    xcorr_z_steph = calc_xcorr((steph_strs_3[i]),(obs_strs_filt_3[i]))[2]
    xcorr_t_steph = find_vec_norm([xcorr_e_steph,xcorr_n_steph,xcorr_z_steph])
    
    dict_ims = {'st_names':st_names_unq[i] , 'st_lons':st_lons_unq[i] , 'st_lats':st_lats_unq[i]
                ,'pgv_e_synt':pgv_e_synt , 'pgv_n_synt':pgv_n_synt , 'pgv_z_synt':pgv_z_synt, 'pgv_t_synt':pgv_t_synt 
                ,'ai_e_synt':ai_e_synt , 'ai_n_synt':ai_n_synt, 'ai_z_synt':ai_z_synt, 'ai_t_synt':ai_t_synt
                ,'pgv_e_obs':pgv_e_obs , 'pgv_n_obs':pgv_n_obs , 'pgv_z_obs':pgv_z_obs, 'pgv_t_obs':pgv_t_obs
                ,'ai_e_obs':ai_e_obs , 'ai_n_obs':ai_n_obs , 'ai_z_obs':ai_z_obs, 'ai_t_obs':ai_t_obs
                ,'pgv_e_steph':pgv_e_steph , 'pgv_n_steph':pgv_n_steph , 'pgv_z_steph':pgv_z_steph, 'pgv_t_steph':pgv_t_steph
                ,'ai_e_steph':ai_e_steph , 'ai_n_steph':ai_n_steph , 'ai_z_steph':ai_z_steph, 'ai_t_steph':ai_t_steph
                ,'xcorr_e_steph':xcorr_e_steph , 'xcorr_n_steph':xcorr_n_steph , 'xcorr_z_steph':xcorr_z_steph, 'xcorr_t_steph':xcorr_t_steph
                ,'xcorr_e_synt':xcorr_e_synt , 'xcorr_n_synt':xcorr_n_synt , 'xcorr_z_synt':xcorr_z_synt, 'xcorr_t_synt':xcorr_t_synt}


    if i == 0:
        df_ims = pd.DataFrame(data=dict_ims,index=[0])
    else:
        df_ims_temp = pd.DataFrame(data=dict_ims , index=[i])
        df_ims = pd.concat([df_ims, df_ims_temp], ignore_index=True)
        
df_ims.to_csv(outdir+ff_name,index=False)

#%%

read_ff = pd.read_csv(outdir+ff_name)
stations_names = np.array(read_ff['st_names'])

val_st_file = pd.read_csv(valst_file)
dom_valst_names = np.array(val_st_file['dom_valst_names'])

valst_ff_mask = np.isin(stations_names,dom_valst_names)
valst_ff = read_ff[valst_ff_mask]

valst_ff.to_csv(outdir+'val_'+ff_name,index=False)

















































