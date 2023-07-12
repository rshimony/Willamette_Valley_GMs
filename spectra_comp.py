#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 15:57:58 2023

@author: rshimony
"""

from obspy.core.utcdatetime import UTCDateTime
import numpy as np
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
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
# freq_vel_filt , amps_vel_filt = calc_spectra(obs_strs_filt_3[20][0])

# freq_steph , amps_steph = calc_spectra(steph_strs_3[20][0])
# freq_synts , amps_synts = calc_spectra(synts_strs_3[20][0])   
    
test_filt_obs = obs_strs_3[20].filter('lowpass', freq=0.9, corners=4, zerophase=True)
test_filt_steph = steph_strs_3[20].filter('lowpass', freq=0.9, corners=4, zerophase=True)
test_filt_synts = synts_strs_3[20].filter('lowpass', freq=0.9, corners=4, zerophase=True)

# test_filt_steph_2 = test_filt_steph.filter('lowpass', freq=1.0, corners=4, zerophase=True)
# test_filt_synts_2 = test_filt_synts.filter('lowpass', freq=1.0, corners=4, zerophase=True)

# test_filt_steph_3 = test_filt_steph_2.filter('lowpass', freq=1.0, corners=4, zerophase=True)
# test_filt_synts_3 = test_filt_synts_2.filter('lowpass', freq=1.0, corners=4, zerophase=True)

freq_obs_filt , amps_obs_filt = calc_spectra(test_filt_obs[0])
freq_steph_filt , amps_steph_filt = calc_spectra(test_filt_steph[0])
freq_synts_filt , amps_synts_filt = calc_spectra(test_filt_synts[0]) 

plt.figure(figsize=(8,10))
plt.loglog(freq_obs_filt,amps_obs_filt, c = 'dimgrey' , label='obs',zorder=0)
plt.loglog(freq_steph_filt,amps_steph_filt, c = 'navy' , label='steph',zorder=1)
plt.loglog(freq_synts_filt,amps_synts_filt,  c = 'maroon' , label='synts',zorder=2)
# plt.loglog(freq_synts_filt,amps_synts_filt,  c = 'green' , label='synts_filt',zorder=3)
# plt.loglog(freq_vel_filt,amps_vel_filt,  c = 'yellow' , label='obs_filt',zorder=3)
plt.xlabel('Freqs [Hz]')
plt.ylabel('Amps')
plt.legend(loc='upper left', bbox_to_anchor=(1, 0.8))
plt.xlim(0.063,(100))
plt.grid(which="both", axis='x') 

#%%      

def subplot_wf_obs(stream , component_idx , subplot_loc , outdir):
    
    freq_vel , amps_vel = calc_spectra(obs_strs_3[i][0])
    freq_vel_filt , amps_vel_filt = calc_spectra(obs_strs_filt_3[i][0])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][0])
    freq_synts , amps_synts = calc_spectra(synts_strs_3[i][0])

    plt.subplot(6,2,1)
    trobs = obs_strs_3[i][0]
    name = st_names_unq[i]
    t_start = trobs.stats.starttime
    samprate_obs = trobs.stats.sampling_rate
    tobs = trobs.times(reftime=UTCDateTime(t_start))
    
    plt.plot(tobs,trobs.data, c = 'dimgrey', label = 'Observed')
    plt.xlim([-5,80])
    plt.title(name + ' ' + trobs.stats.channel ,size=12)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yaxis.offsetText.set_fontsize(10)
    
    limobs = max(abs(trobs.data))
    
    plt.ylim([-1*limobs, limobs])
    plt.ylim([-1*limobs, limobs])
    plt.locator_params(tight=True, nbins=4)
    plt.ylabel('vel (m/s)',fontsize=14)
    plt.tick_params(labelsize=14)
    plt.legend(loc = 'upper right',prop={'size': 10.0})
    
    
    
               
#%%
bot_lim = 0.063

for i in range(len(obs_strs_3)):
        
    # fig = plt.figure(figsize=(10,4))
    plt.figure(figsize=(18,15))
    
    # test_filt_steph = steph_strs_3[i].filter('lowpass', freq=0.8, corners=4, zerophase=True)
    # test_filt_synts = synts_strs_3[i].filter('lowpass', freq=0.8, corners=4, zerophase=True)
    
    freq_vel , amps_vel = calc_spectra(obs_strs_3[i][0])
    freq_vel_filt , amps_vel_filt = calc_spectra(obs_strs_filt_3[i][0])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][0])
    freq_synts , amps_synts = calc_spectra(synts_strs_3[i][0])

    # ax1 = fig.add_subplot(3,2,j+1)
    plt.subplot(6,2,1)
    trobs = obs_strs_3[i][0]
    name = st_names_unq[i]
    t_start = trobs.stats.starttime
    samprate_obs = trobs.stats.sampling_rate
    tobs = trobs.times(reftime=UTCDateTime(t_start))
    
    plt.plot(tobs,trobs.data, c = 'dimgrey', label = 'Observed')
    plt.xlim([-5,80])
    plt.title(name + ' ' + trobs.stats.channel ,size=12)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # ax.yaxis.offsetText.set_fontsize(10)
    
    limobs = max(abs(trobs.data))
    
    plt.ylim([-1*limobs, limobs])
    plt.ylim([-1*limobs, limobs])
    plt.locator_params(tight=True, nbins=4)
    plt.ylabel('vel (m/s)')#,fontsize=14)
    # ax1.tick_params(labelsize=14)
    plt.legend(loc = 'upper right')#,prop={'size': 10.0})
    
    # ax1 = fig.add_subplot(3,2,j+1)
    plt.subplot(6,2,3)
    trobs_filt = obs_strs_filt_3[i][0]
    trsyn = synts_strs_3[i][0]
    trsteph = steph_strs_3[i][0]
    npts_syn = trsyn.stats.npts
    samprate_syn = trsyn.stats.sampling_rate
    tsyn = np.arange(0, npts_syn / samprate_syn, 1 / samprate_syn)
    npts_steph = trsteph.stats.npts
    samprate_steph = trsteph.stats.sampling_rate
    tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)
    
    plt.plot(tobs,trobs_filt.data, c = 'navy', label = 'Observed_filt',zorder=0)
    plt.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
    plt.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',zorder=2)
    plt.xlim([-5,80])
    # plt.title(name + ' ' + trobs.stats.channel)# ,size=12)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # ax1.yaxis.offsetText.set_fontsize(10)
    
    limobs = max(abs(trobs_filt.data))
    limsyn = max(abs(trsyn.data))
    lim = max(limobs, limsyn)
    
    plt.ylim([-1*lim, lim])
    plt.ylim([-1*lim, lim])
    plt.locator_params(tight=True, nbins=4)
    plt.ylabel('vel (m/s)')#,fontsize=14)
    # ax1.tick_params(labelsize=14)
    plt.legend(loc = 'upper right')#,prop={'size': 10.0})
           
    
    # ax2 = fig.add_subplot(3,2,j+2)
    plt.subplot(3,2,2)
    plt.loglog(freq_vel,amps_vel, c = 'dimgrey' , label='Obsereved',zorder=0)
    plt.loglog(freq_vel_filt,amps_vel_filt, c = 'navy' , label='Obsereved filt',zorder=1)
    plt.loglog(freq_synts,amps_synts,  c = 'maroon' , label='Synthetic',zorder=2)
    plt.loglog(freq_steph,amps_steph,  c = 'green' , label='Stephenson',zorder=3)
    plt.xlabel('Freqs [Hz]')
    plt.ylabel('Amps')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.8))
    plt.xlim(bot_lim,(samprate_obs/2))
    plt.grid(which="both", axis='x')   


    freq_vel , amps_vel = calc_spectra(obs_strs_3[i][1])
    freq_vel_filt , amps_vel_filt = calc_spectra(obs_strs_filt_3[i][1])
    # freq_steph , amps_steph = calc_spectra(test_filt_steph[1])
    # freq_synts , amps_synts = calc_spectra(test_filt_synts[1])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][1])
    freq_synts , amps_synts = calc_spectra(synts_strs_3[i][1])

    # ax1 = fig.add_subplot(3,2,j+1)
    plt.subplot(6,2,5)
    trobs = obs_strs_3[i][1]
    name = st_names_unq[i]
    t_start = trobs.stats.starttime
    samprate_obs = trobs.stats.sampling_rate
    tobs = trobs.times(reftime=UTCDateTime(t_start))
    
    plt.plot(tobs,trobs.data, c = 'dimgrey', label = 'Observed')
    plt.xlim([-5,80])
    plt.title(name + ' ' + trobs.stats.channel ,size=12)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # ax1.yaxis.offsetText.set_fontsize(10)
    
    limobs = max(abs(trobs.data))
    
    plt.ylim([-1*limobs, limobs])
    plt.ylim([-1*limobs, limobs])
    plt.locator_params(tight=True, nbins=4)
    plt.ylabel('vel (m/s)')#,fontsize=14)
    # ax1.tick_params(labelsize=14)
    plt.legend(loc = 'upper right')#,prop={'size': 10.0})
    
    # ax1 = fig.add_subplot(3,2,j+1)
    plt.subplot(6,2,7)
    trobs_filt = obs_strs_filt_3[i][1]
    trsyn = synts_strs_3[i][1]
    trsteph = steph_strs_3[i][1]
    npts_syn = trsyn.stats.npts
    samprate_syn = trsyn.stats.sampling_rate
    tsyn = np.arange(0, npts_syn / samprate_syn, 1 / samprate_syn)
    npts_steph = trsteph.stats.npts
    samprate_steph = trsteph.stats.sampling_rate
    tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)

    plt.plot(tobs,trobs_filt.data, c = 'navy', label = 'Observed_filt',zorder=0)
    plt.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
    plt.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',zorder=2)
    plt.xlim([-5,80])
    # plt.title(name + ' ' + trobs.stats.channel)# ,size=12)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # ax1.yaxis.offsetText.set_fontsize(10)
    
    limobs = max(abs(trobs_filt.data))
    limsyn = max(abs(trsyn.data))
    lim = max(limobs, limsyn)
    
    plt.ylim([-1*lim, lim])
    plt.ylim([-1*lim, lim])
    plt.locator_params(tight=True, nbins=4)
    plt.ylabel('vel (m/s)')#,fontsize=14)
    # ax1.tick_params(labelsize=14)
    plt.legend(loc = 'upper right')#,prop={'size': 10.0})
           
    
    # ax2 = fig.add_subplot(3,2,j+2)
    plt.subplot(3,2,4)
    plt.loglog(freq_vel,amps_vel, c = 'dimgrey' , label='Obsereved',zorder=0)
    plt.loglog(freq_vel_filt,amps_vel_filt, c = 'navy' , label='Obsereved filt',zorder=1)
    plt.loglog(freq_synts,amps_synts,  c = 'maroon' , label='Synthetic',zorder=2)
    plt.loglog(freq_steph,amps_steph,  c = 'green' , label='Stephenson',zorder=3)
    plt.xlabel('Freqs [Hz]')
    plt.ylabel('Amps')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.8))
    plt.xlim(bot_lim,(samprate_obs/2))
    plt.grid(which="both", axis='x')        
    
    
    freq_vel , amps_vel = calc_spectra(obs_strs_3[i][2])
    freq_vel_filt , amps_vel_filt = calc_spectra(obs_strs_filt_3[i][2])
    # freq_steph , amps_steph = calc_spectra(test_filt_steph[2])
    # freq_synts , amps_synts = calc_spectra(test_filt_synts[2])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][2])
    freq_synts , amps_synts = calc_spectra(synts_strs_3[i][2])

    # ax1 = fig.add_subplot(3,2,j+1)
    plt.subplot(6,2,9)
    trobs = obs_strs_3[i][2]
    name = st_names_unq[i]
    t_start = trobs.stats.starttime
    samprate_obs = trobs.stats.sampling_rate
    tobs = trobs.times(reftime=UTCDateTime(t_start))
    
    plt.plot(tobs,trobs.data, c = 'dimgrey', label = 'Observed',zorder=0)
    plt.xlim([-5,80])
    plt.title(name + ' ' + trobs.stats.channel ,size=12)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # ax1.yaxis.offsetText.set_fontsize(10)
    
    limobs = max(abs(trobs.data))
    
    plt.ylim([-1*limobs, limobs])
    plt.ylim([-1*limobs, limobs])
    plt.locator_params(tight=True, nbins=4)
    plt.ylabel('vel (m/s)')#,fontsize=14)
    # ax1.tick_params(labelsize=14)
    plt.legend(loc = 'upper right')#,prop={'size': 10.0})
    
    # ax1 = fig.add_subplot(3,2,j+1)
    plt.subplot(6,2,11)
    trobs_filt = obs_strs_filt_3[i][2]
    trsyn = synts_strs_3[i][2]
    trsteph = steph_strs_3[i][2]
    npts_syn = trsyn.stats.npts
    samprate_syn = trsyn.stats.sampling_rate
    tsyn = np.arange(0, npts_syn / samprate_syn, 1 / samprate_syn)
    npts_steph = trsteph.stats.npts
    samprate_steph = trsteph.stats.sampling_rate
    tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)
    
    plt.plot(tobs,trobs_filt.data, c = 'navy', label = 'Observed_filt',zorder=0)
    plt.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
    plt.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',zorder=2)
    plt.xlim([-5,80])
    # plt.title(name + ' ' + trobs.stats.channel)# ,size=12)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # ax1.yaxis.offsetText.set_fontsize(10)
    
    limobs = max(abs(trobs_filt.data))
    limsyn = max(abs(trsyn.data))
    lim = max(limobs, limsyn)
    
    plt.ylim([-1*lim, lim])
    plt.ylim([-1*lim, lim])
    plt.locator_params(tight=True, nbins=4)
    plt.ylabel('vel (m/s)')#,fontsize=14)
    # ax1.tick_params(labelsize=14)
    plt.legend(loc = 'upper right')#,prop={'size': 10.0})
           
    
    # ax2 = fig.add_subplot(3,2,j+2)
    plt.subplot(3,2,6)
    plt.loglog(freq_vel,amps_vel, c = 'dimgrey' , label='Obsereved',zorder=0)
    plt.loglog(freq_vel_filt,amps_vel_filt, c = 'navy' , label='Obsereved filt',zorder=1)
    plt.loglog(freq_synts,amps_synts,  c = 'maroon' , label='Synthetic',zorder=2)
    plt.loglog(freq_steph,amps_steph,  c = 'green' , label='Stephenson',zorder=3)
    plt.xlabel('Freqs [Hz]')
    plt.ylabel('Amps')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.8))
    plt.xlim(bot_lim,(samprate_obs/2))
    plt.grid(which="both", axis='x')  
    
    plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/wf_spec_comp/synts_gp1_eugspe_sm_filt_cheby/'+name+'.png' ,dpi=300,bbox_inches='tight')       
       
       
plt.close('all')       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
      