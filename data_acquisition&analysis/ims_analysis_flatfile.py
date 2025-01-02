#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:45:34 2024

@author: rshimony
"""
'''
Calculates intensity measures (IMs - PGV, Arias Intensity, Cross Correlation) from waveform data (observed, new WV model synthetics, USGS CVM synthetics).
Saving the calculated IMs in a flatfile.
This script works for ONE EVENT and ONE MODEL at a time:
    - Observed data
    - USGS CVM synthetics
    - 1 new WV model
The output flatfile holds IMs for the three datasets.
Each IM is calculated for each component seperatly (East(e), North(n), Vertical(z)) and then the norm of the 3 components (t).
Finally, creates a Valley stations only IMs flatfile, out of the newly made full station inventry IMs flatfile.
'''

import numpy as np 
import obspy as obs
import pandas as pd
from glob import glob

#%%
##Functions##
def make_wf_strings(data_dir,data_type,fmax):
    """
    Makes three components, low-pass filtered, waveform strings 
    from a list of single .sac files.

    Inputs:
        Data directory path.
        Data type (Observed "obs" or Synthetic "synts").
        Maxmimum frequency of the low-pass filter.

    Return:
        Strings list
    """
    
    if data_type == 'obs':
        traces = sorted(glob(data_dir + '*_unfilt.sac'))
        trace_nms = []
        for i in range(len(traces)):
            trace_nm = traces[i].split('.')[1]
            trace_nms.append(trace_nm)
        trace_nms_unq = np.unique(trace_nms)
        
        strings = []
        strings_filt = []
        for i in range(len(st_names_unq)):
            if st_names_unq[i] in trace_nms_unq:
                string = obs.read(data_dir + '*' + st_names_unq[i] + '*_unfilt.sac')
                string_cp = string.copy()
                strings.append(string_cp)
                string_filt = string.copy()
                string_filt.filter('lowpass_cheby_2',freq=fmax)
                strings_filt.append(string_filt)
                
        strings_filt_3 = []       
        strings_3 = []
        for i in strings:
            if len(i) >= 3:
               strings_3.append(i[-3:])
        for i in strings_filt:
            if len(i) >= 3:
               strings_filt_3.append(i[-3:])
        
        return strings_3 , strings_filt_3
        
    if data_type == 'synts':
        traces = sorted(glob(data_dir + '*.*v'))  
        trace_nms = []
        for i in range(len(traces)):
            trace_nm = (traces[i].split('/')[-1]).split('.')[0]
            trace_nms.append(trace_nm)
        trace_nms_unq = np.unique(trace_nms)
        
        strings_filt = []
        for i in range(len(st_names_unq)):
            if st_names_unq[i] in trace_nms_unq:
                string = obs.read(data_dir + st_names_unq[i] + '.*v')
                string_cp = string.copy()
                string_filt = string_cp.filter('lowpass_cheby_2',freq=fmax)
                strings_filt.append(string_filt)
            
        strings_3 = []
        for i in strings_filt:
            if len(i) >= 3:
               strings_3.append(i[-3:]) 
               
        return strings_3 
    
def get_pgv(stream):
    
    '''
    Gives the PGV from a an obspy stream
    '''
    
    import numpy as np

    pgv_list = []
    
    for i in range(len(stream)):
        stream_copy = stream.copy()
        data = stream_copy[i].data
        pgv = np.max(np.abs(data))
        pgv_list.append(pgv)
        
    return(pgv_list) 

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
###INPUTS####
# SNR station inventory data
st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/metadata/snr_stations.csv')
st_names = np.array(st_file['st_nm'])
st_lons = np.array(st_file['st_lon'])
st_lats = np.array(st_file['st_lat']) 
st_names_unq = pd.unique(st_names)
st_lons_unq = pd.unique(st_lons)
st_lats_unq = pd.unique(st_lats)

# Valley SNR stations inventory
valst_file = '/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/metadata/valley_snr_st.csv'

# Waveform data
obs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/waveforms/observed_data_1hz_snr/' # Observed
steph_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/steph/' # USGS CVM
synts_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/srv_sm/' # New WV model

# outpur directory
outdir = '/Users/rshimony/Desktop/WillametteValley/wv_project/im_analysis/springfield/'
# output flatfile name
ff_name = 'ims_ff_srv_sm.csv'

#%%
# Getting strings
obs_strs_3,obs_strs_filt_3 = make_wf_strings(obs_dir,'obs',1.0) 
steph_strs_3 =  make_wf_strings(steph_dir,'synts',1.0)
synts_strs_3 = make_wf_strings(synts_dir,'synts',1.0) 

#%%
## Calculating IMS and saving them to a .csv file
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
# Creating a Valley SNR IMs flatfile
read_ff = pd.read_csv(outdir+ff_name)
stations_names = np.array(read_ff['st_names'])

val_st_file = pd.read_csv(valst_file)
dom_valst_names = np.array(val_st_file['dom_valst_names'])

valst_ff_mask = np.isin(stations_names,dom_valst_names)
valst_ff = read_ff[valst_ff_mask]

valst_ff.to_csv(outdir+'val_'+ff_name,index=False)

















































