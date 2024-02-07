#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 19:22:54 2024

@author: rshimony
"""

from obspy.core.utcdatetime import UTCDateTime
import numpy as np
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import obspy as obs
import pandas as pd
from glob import glob
import os

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


st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/salem_eq/metadata/domain_stations.csv')
st_names = np.array(st_file['st_name_dom'])
st_lons = np.array(st_file['st_lon_dom'])
st_lats = np.array(st_file['st_lat_dom']) 
st_names_unq = np.unique(st_names)


obs_dir = '/Users/rshimony/Desktop/WillametteValley/salem_eq/observed_data_salem_1hz/'
steph_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/salem/rfile_stephenson_salem_topo_preempt_freq25/'
eugspe_1d_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/salem/rfile_steph_basin_salem_eugspe_1d_preempt_freq25/'
eugspe_sm_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/salem/rfile_steph_basin_salem_cart_eugspe_sm_preempt_freq25/'
srv_1d_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/salem/rfile_steph_basin_salem_srv_1d_preempt_freq25/'
srv_sm_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/salem/rfile_steph_basin_salem_srv_sm_preempt_freq25/'

working_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/spec_comp_model/salem/'

try:
    os.mkdir(working_dir)
except FileExistsError:
    print('Directory already exists')

#%%

def make_wf_strings(data_dir,data_type,fmax):
    
    if data_type == 'obs':
        traces = sorted(glob(data_dir + '*_unfilt.sac'))
        trace_nms = []
        for i in range(len(traces)):
            trace_nm = traces[i].split('.')[1]
            trace_nms.append(trace_nm)
        trace_nms_unq = np.unique(trace_nms)
        
        strings = []
        if fmax == 1:
            strings_filt = []
            for i in range(len(st_names_unq)):
                if st_names_unq[i] in trace_nms_unq:
                    string = obs.read(data_dir + '*' + st_names_unq[i] + '*_unfilt.sac')
                    string_cp = string.copy()
                    strings.append(string_cp)
                    string_filt = string.copy()
                    string_filt.filter('lowpass_cheby_2',freq=1.0)
                    strings_filt.append(string_filt)
                
        if fmax == 2:
            strings_filt = []
            for i in range(len(st_names_unq)):
                if st_names_unq[i] in trace_nms_unq:
                    string = obs.read(data_dir + '*' + st_names_unq[i] + '*_unfilt.sac')
                    string_cp = string.copy()
                    strings.append(string_cp)
                    string_filt = string.copy()
                    string_filt.filter('lowpass_cheby_2',freq=2.0)
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
        if fmax == 1:
            for i in range(len(st_names_unq)):
                if st_names_unq[i] in trace_nms_unq:
                    string = obs.read(data_dir + st_names_unq[i] + '.*v')
                    string_cp = string.copy()
                    string_filt = string_cp.filter('lowpass_cheby_2',freq=1.0)
                    strings_filt.append(string_filt)
            
        if fmax == 2:
            for i in range(len(st_names_unq)):
                if st_names_unq[i] in trace_nms_unq:
                    string = obs.read(data_dir + st_names_unq[i] + '.*v')
                    string_cp = string.copy()
                    string_filt = string_cp.filter('lowpass_cheby_2',freq=2.0)
                    strings_filt.append(string_filt)
            
        strings_3 = []
        for i in strings_filt:
            if len(i) >= 3:
               strings_3.append(i[-3:]) 
               
        return strings_3 
       
       
#%%

obs_strs_3,obs_strs_filt_3 = make_wf_strings(obs_dir,'obs',1) 
steph_strs_3 =  make_wf_strings(steph_dir,'synts',1)
eugspe_1d_strs_3 = make_wf_strings(eugspe_1d_dir,'synts',1)   
eugspe_sm_strs_3 = make_wf_strings(eugspe_sm_dir,'synts',1)  
srv_1d_strs_3 = make_wf_strings(srv_1d_dir,'synts',1)  
srv_sm_strs_3 = make_wf_strings(srv_sm_dir,'synts',1)
       
#%%

bot_lim = 0.063
linewidth=1.0

for i in range(len(obs_strs_3)):
        
    # plt.figure(figsize=(10,5))
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))
    name = st_names_unq[i]
    fig.suptitle(name,size=24,y=1.05)
    
    freq_obs , amps_obs = calc_spectra(obs_strs_3[i][0])
    freq_obs_filt , amps_obs_filt = calc_spectra(obs_strs_filt_3[i][0])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][0])
    freq_eugspe_1d , amps_eugspe_1d = calc_spectra(eugspe_1d_strs_3[i][0])
    freq_eugspe_sm , amps_eugspe_sm = calc_spectra(eugspe_sm_strs_3[i][0])
    freq_srv_1d , amps_srv_1d = calc_spectra(srv_1d_strs_3[i][0])
    freq_srv_sm , amps_srv_sm = calc_spectra(srv_sm_strs_3[i][0])
    
    trobs = obs_strs_3[i][0]
    samprate_obs = trobs.stats.sampling_rate

    # ax1.subplot(1,2,1)
    # ax3 = plt.subplot2grid((13,5), (2,3), colspan=2, rowspan=4)
    ax1.loglog(freq_obs,amps_obs, c = 'dimgrey' , label='Obsereved',alpha=0.7,linewidth=linewidth)
    ax1.loglog(freq_obs_filt,amps_obs_filt, c = 'navy' , label='Obsereved filt',alpha=0.7,linewidth=linewidth)
    ax1.loglog(freq_eugspe_1d,amps_eugspe_1d,  c = 'darkred' , label='eugspe_1d',alpha=0.7,linewidth=linewidth,ls='--')
    ax1.loglog(freq_eugspe_sm,amps_eugspe_sm,  c = 'chocolate' , label='eugspe_sm',alpha=0.7,linewidth=linewidth,ls='--')
    ax1.loglog(freq_srv_1d,amps_srv_1d,  c = 'gold' , label='srv_1d',alpha=0.7,linewidth=linewidth,ls='-.')
    ax1.loglog(freq_srv_sm,amps_srv_sm,  c = 'darkviolet' , label='srv_sm',alpha=0.7,linewidth=linewidth,ls='-.')
    ax1.loglog(freq_steph,amps_steph,  c = 'green' , label='Stephenson',alpha=0.7,linewidth=linewidth)
    ax1.set_xlabel('Freqs [Hz]',fontsize=14)
    ax1.set_ylabel('Amps',fontsize=14)
    # plt.legend(loc='upper center',bbox_to_anchor=(0.45,  1.6),prop={'size': 11.0,'weight':'bold'},ncol=2)
    ax1.set_xlim(bot_lim,(samprate_obs/2))
    ax1.set_ylim(1e-8,max(amps_obs)+5e-5)
    ax1.grid(which="both", axis='x') 
    ax1.tick_params(axis='both', which='major', labelsize=15)
    # ax1.set_yticks(fontsize=15)
    # ax1.set_xticks(fontsize=15)
    # ax1.xaxis.set_ticklabels([])
    ax1.set_title('East',size=15)


    freq_obs , amps_obs = calc_spectra(obs_strs_3[i][1])
    freq_obs_filt , amps_obs_filt = calc_spectra(obs_strs_filt_3[i][1])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][1])
    freq_eugspe_1d , amps_eugspe_1d = calc_spectra(eugspe_1d_strs_3[i][1])
    freq_eugspe_sm , amps_eugspe_sm = calc_spectra(eugspe_sm_strs_3[i][1])
    freq_srv_1d , amps_srv_1d = calc_spectra(srv_1d_strs_3[i][1])
    freq_srv_sm , amps_srv_sm = calc_spectra(srv_sm_strs_3[i][1])
    
    trobs = obs_strs_3[i][1]
    samprate_obs = trobs.stats.sampling_rate
     
    # ax2.subplot(1,2,2)
    # ax6 = plt.subplot2grid((13,5), (9,3), colspan=2, rowspan=4)
    ax2.loglog(freq_obs,amps_obs, c = 'dimgrey' , label='Obsereved',alpha=0.7,linewidth=linewidth)
    ax2.loglog(freq_obs_filt,amps_obs_filt, c = 'navy' , label='Obsereved filt',alpha=0.7,linewidth=linewidth)
    ax2.loglog(freq_eugspe_1d,amps_eugspe_1d,  c = 'darkred' , label='eugspe_1d',alpha=0.7,linewidth=linewidth,ls='--')
    ax2.loglog(freq_eugspe_sm,amps_eugspe_sm,  c = 'chocolate' , label='eugspe_sm',alpha=0.7,linewidth=linewidth,ls='--')
    ax2.loglog(freq_srv_1d,amps_srv_1d,  c = 'gold' , label='srv_1d',alpha=0.7,linewidth=linewidth,ls='-.')
    ax2.loglog(freq_srv_sm,amps_srv_sm,  c = 'darkviolet' , label='srv_sm',alpha=0.7,linewidth=linewidth,ls='-.')
    ax2.loglog(freq_steph,amps_steph,  c = 'green' , label='Stephenson',alpha=0.7,linewidth=linewidth)
    ax2.set_xlabel('Freqs [Hz]',fontsize=14)
    ax2.set_ylabel('Amps',fontsize=14)
    ax2.legend(loc='upper right',prop={'size': 10.0,'weight':'bold'})
    ax2.set_xlim(bot_lim,(samprate_obs/2))
    ax2.set_ylim(1e-8,max(amps_obs)+5e-5)
    ax2.grid(which="both", axis='x') 
    ax2.tick_params(axis='both', which='major', labelsize=15)
    # ax2.set_yticks(fontsize=15)
    # ax2.set_xticks(fontsize=15)
    # ax2.xaxis.set_ticklabels([])
    ax2.set_title('North',size=15)
    
    # plt.subplots_adjust(hspace=0.6,wspace=1.1)
    fig.tight_layout()
    
    fig.savefig(working_dir+name+'.png' ,dpi=150,bbox_inches='tight') 

plt.close('all') 