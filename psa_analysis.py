#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 11:20:19 2023

@author: rshimony
"""

import os
import re

import matplotlib.pyplot as plt
import numpy as np 
import pyrotd
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client
import obspy as obs
import pandas as pd
from glob import glob

#%%

def calc_spectra(trace):
    """
    Calculates average spectra values for displacement, acceleration,
    and velocity waveforms.

    Inputs:
        Obspy trace. 

    Return:
        freq(list): FFT frequencies
        amp(list): FFT amplitudes
    """
    
    import numpy as np
    from mtspec import mtspec

    # Read in file 
    tr = trace
    data = tr.data
    delta = tr.stats.delta
    npts = tr.stats.npts
    

    # Calc spectra amplitudes and frequencies 

    amp_squared, freq =  mtspec(data, delta=delta, time_bandwidth=4, 
                              number_of_tapers=7, nfft=npts, quadratic=True)
    
    # Convert from power spectra to amplitude spectra
    amp = np.sqrt(amp_squared)
    
    return(freq, amp)


st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/domain_stations.csv')
st_names = np.array(st_file['st_name_dom'])
st_lons = np.array(st_file['st_lon_dom'])
st_lats = np.array(st_file['st_lat_dom']) 
st_names_unq = np.unique(st_names)


obs_dir = '/Users/rshimony/Desktop/WillametteValley/salem_eq/observed_data_salem_1hz/'
steph_dir = '/Users/rshimony/Desktop/WillametteValley/Models/rfile_stephenson_salem_topo/'
synts_dir = '/Users/rshimony/Desktop/WillametteValley/Models/eug_spe/eug_spe_singmat_gravfix/'

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

trobs = obs_strs_filt_3[20][0]
t_start = trobs.stats.starttime
dt = trobs.stats.delta
tobs = trobs.times(reftime=UTCDateTime(t_start))
vels = trobs.data
# plt.figure(figsize=(10,2))     
# plt.plot(tobs,trobs, c = 'navy', label = 'Observed_filt',zorder=0)        
      
osc_damping = 0.05
osc_freqs = [0.1,0.2,0.5,1]

osc_resps = pyrotd.calc_spec_accels(
        dt, trobs.data, osc_freqs, osc_damping)
        
#%%

def load_at2(fpath):
    with open(fpath) as fp:
        for _ in range(3):
            next(fp)
        line = next(fp)
        time_step = float(line[17:25])
        accels = np.array([p for l in fp for p in l.split()]).astype(float)
    return time_step, accels

time_step, accels = load_at2('/Users/rshimony/Downloads/RSN8883_14383980_13849090.AT2')     
        
      
osc_damping = 0.05
osc_freqs = [0.1,0.2,0.5,1]

osc_resps = pyrotd.calc_spec_accels(
        time_step, accels, osc_freqs, osc_damping
)

      
        
      
        
      
        
      
        
      
        
      
        