#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:09:44 2024

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


st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/snr_stations_fixed.csv')
st_names = np.array(st_file['st_nm'])
st_lons = np.array(st_file['st_lon'])
st_lats = np.array(st_file['st_lat']) 
st_names_unq = np.unique(st_names)


obs_dir = '/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/waveforms/observed_data_1hz_snr/'
steph_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/scottsmills/rfile_steph_scottsmills_preempt_freq25/'
eugspe_1d_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/scottsmills/rfile_steph_basin_scottsmills_eugspe_1d_preempt_freq25/'
eugspe_sm_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/scottsmills/rfile_steph_basin_scottsmills_eugspe_sm_preempt_freq25/'
srv_1d_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/scottsmills/rfile_steph_basin_scottsmills_srv_1d_preempt_freq25/'
srv_sm_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/scottsmills/rfile_steph_basin_scottsmills_srv_sm_preempt_freq25/'

working_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/wf_spec_model_comp_horiz/scottsmills/'

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
        
    plt.figure(figsize=(13,9))
    gridspec.GridSpec(13,5)
    name = st_names_unq[i]
    plt.suptitle(name,size=24,y=0.95)
    plt.figtext(0.05,0.5,'vel (m/s)',fontsize=20,rotation = 90)
    plt.figtext(0.57,0.5,'Amps',fontsize=20,rotation = 90)
    
    freq_obs , amps_obs = calc_spectra(obs_strs_3[i][0])
    freq_obs_filt , amps_obs_filt = calc_spectra(obs_strs_filt_3[i][0])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][0])
    freq_eugspe_1d , amps_eugspe_1d = calc_spectra(eugspe_1d_strs_3[i][0])
    freq_eugspe_sm , amps_eugspe_sm = calc_spectra(eugspe_sm_strs_3[i][0])
    freq_srv_1d , amps_srv_1d = calc_spectra(srv_1d_strs_3[i][0])
    freq_srv_sm , amps_srv_sm = calc_spectra(srv_sm_strs_3[i][0])

    # plt.subplot(6,2,1)
    ax1 = plt.subplot2grid((13,5), (0,0), colspan=3, rowspan=2)
    trobs = obs_strs_3[i][0]
    # name = st_names_unq[i]
    t_start = trobs.stats.starttime
    samprate_obs = trobs.stats.sampling_rate
    tobs = trobs.times(reftime=UTCDateTime(t_start))
    
    plt.plot(tobs,trobs.data, c = 'dimgrey', label = 'Observed',linewidth=linewidth)
    plt.xlim([-5,80])
    plt.title('East' ,size=15)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yticks(fontsize=15)
    # plt.xticks(fontsize=15)
    ax1.xaxis.set_ticklabels([])
    
    limobs = max(abs(trobs.data))
    
    plt.ylim([-1*limobs, limobs])
    plt.ylim([-1*limobs, limobs])
    plt.locator_params(tight=True, nbins=4)
    # plt.ylabel('vel (m/s)',fontsize=14)
    # plt.xlabel('time (sec)',fontsize=14)
    
    # plt.subplot(6,2,3)
    ax2 = plt.subplot2grid((13,5), (2,0), colspan=3, rowspan=4)
    trobs_filt = obs_strs_filt_3[i][0]
    trsyn_eugspe_1d = eugspe_1d_strs_3[i][0]
    trsyn_eugspe_sm = eugspe_sm_strs_3[i][0]
    trsyn_srv_1d = srv_1d_strs_3[i][0]
    trsyn_srv_sm = srv_sm_strs_3[i][0]
    trsteph = steph_strs_3[i][0]
    npts_eugspe_1d = trsyn_eugspe_1d.stats.npts
    npts_eugspe_sm = trsyn_eugspe_sm.stats.npts
    npts_srv_1d = trsyn_srv_1d.stats.npts
    npts_srv_sm = trsyn_srv_sm.stats.npts
    npts_steph = trsteph.stats.npts
    samprate_eugspe_1d = trsyn_eugspe_1d.stats.sampling_rate
    samprate_eugspe_sm = trsyn_eugspe_sm.stats.sampling_rate
    samprate_srv_1d = trsyn_srv_1d.stats.sampling_rate
    samprate_srv_sm = trsyn_srv_sm.stats.sampling_rate
    samprate_steph = trsteph.stats.sampling_rate
    t_eugspe_1d = np.arange(0, npts_eugspe_1d / samprate_eugspe_1d, 1 / samprate_eugspe_1d) 
    t_eugspe_sm = np.arange(0, npts_eugspe_sm / samprate_eugspe_sm, 1 / samprate_eugspe_sm)
    t_srv_1d = np.arange(0, npts_srv_1d / samprate_srv_1d, 1 / samprate_srv_1d)
    t_srv_sm = np.arange(0, npts_srv_sm / samprate_srv_sm, 1 / samprate_srv_sm)
    tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)
    
    plt.plot(tobs,trobs_filt.data, c = 'navy', label = 'Observed_filt',zorder=0,linewidth=linewidth)
    plt.plot(t_eugspe_1d,trsyn_eugspe_1d,  c = 'darkred', label  = 'eugspe_1d',alpha=0.7,linewidth=linewidth,ls='--')
    plt.plot(t_eugspe_sm,trsyn_eugspe_sm,  c = 'chocolate', label  = 'eugspe_sm',alpha=0.7,linewidth=linewidth,ls='--')
    plt.plot(t_srv_1d,trsyn_srv_1d,  c = 'gold', label  = 'srv_1d',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.plot(t_srv_sm,trsyn_srv_sm,  c = 'darkviolet', label  = 'srv_sm',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',alpha=0.7,linewidth=linewidth)
    plt.xlim([-5,80])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yticks(fontsize=15)
    # plt.xticks(fontsize=15)
    ax2.xaxis.set_ticklabels([])
    
    limobs = max(abs(trobs_filt.data))
    limsyn_eugspe_1d = max(abs(trsyn_eugspe_1d.data))
    limsyn_eugspe_sm = max(abs(trsyn_eugspe_sm.data))
    limsyn_srv_1 = max(abs(trsyn_srv_1d.data))
    limsyn_srv_sm = max(abs(trsyn_srv_sm.data))
    lim = max(limobs, limsyn_eugspe_1d , limsyn_eugspe_sm , limsyn_srv_1 , limsyn_srv_sm)
    
    plt.ylim([-1*lim, lim])
    plt.ylim([-1*lim, lim])
    plt.locator_params(tight=True, nbins=4)
    # plt.ylabel('vel (m/s)',fontsize=14)
    # plt.xlabel('time (sec)',fontsize=14)


    # plt.subplot(3,2,2)
    ax3 = plt.subplot2grid((13,5), (2,3), colspan=2, rowspan=4)
    plt.loglog(freq_obs,amps_obs, c = 'dimgrey' , label='Obsereved',alpha=0.7,linewidth=linewidth)
    plt.loglog(freq_obs_filt,amps_obs_filt, c = 'navy' , label='Obsereved filt',alpha=0.7,linewidth=linewidth)
    plt.loglog(freq_eugspe_1d,amps_eugspe_1d,  c = 'darkred' , label='eugspe_1d',alpha=0.7,linewidth=linewidth,ls='--')
    plt.loglog(freq_eugspe_sm,amps_eugspe_sm,  c = 'chocolate' , label='eugspe_sm',alpha=0.7,linewidth=linewidth,ls='--')
    plt.loglog(freq_srv_1d,amps_srv_1d,  c = 'gold' , label='srv_1d',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.loglog(freq_srv_sm,amps_srv_sm,  c = 'darkviolet' , label='srv_sm',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.loglog(freq_steph,amps_steph,  c = 'green' , label='Stephenson',alpha=0.7,linewidth=linewidth)
    # plt.xlabel('Freqs [Hz]',fontsize=14)
    # plt.ylabel('Amps',fontsize=14)
    plt.legend(loc='upper center',bbox_to_anchor=(0.45,  1.6),prop={'size': 11.0,'weight':'bold'},ncol=2)
    plt.xlim(bot_lim,(samprate_obs/2))
    plt.ylim(1e-8,max(amps_obs)+5e-5)
    plt.grid(which="both", axis='x') 
    plt.yticks(fontsize=15)
    # plt.xticks(fontsize=15)
    ax3.xaxis.set_ticklabels([])


    freq_obs , amps_obs = calc_spectra(obs_strs_3[i][1])
    freq_obs_filt , amps_obs_filt = calc_spectra(obs_strs_filt_3[i][1])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][1])
    freq_eugspe_1d , amps_eugspe_1d = calc_spectra(eugspe_1d_strs_3[i][1])
    freq_eugspe_sm , amps_eugspe_sm = calc_spectra(eugspe_sm_strs_3[i][1])
    freq_srv_1d , amps_srv_1d = calc_spectra(srv_1d_strs_3[i][1])
    freq_srv_sm , amps_srv_sm = calc_spectra(srv_sm_strs_3[i][1])

    # plt.subplot(6,2,5)
    ax4= plt.subplot2grid((13,5), (7,0), colspan=3, rowspan=2)
    trobs = obs_strs_3[i][1]
    # name = st_names_unq[i]
    t_start = trobs.stats.starttime
    samprate_obs = trobs.stats.sampling_rate
    tobs = trobs.times(reftime=UTCDateTime(t_start))
    
    plt.plot(tobs,trobs.data, c = 'dimgrey', label = 'Observed',linewidth=linewidth)
    plt.xlim([-5,80])
    plt.title('North' ,size=15)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yticks(fontsize=15)
    # plt.xticks(fontsize=15)
    # ax1.yaxis.offsetText.set_fontsize(10)
    ax4.xaxis.set_ticklabels([])
    
    limobs = max(abs(trobs.data))
    
    plt.ylim([-1*limobs, limobs])
    plt.ylim([-1*limobs, limobs])
    plt.locator_params(tight=True, nbins=4)
    # plt.ylabel('vel (m/s)',fontsize=14)
    # plt.xlabel('time (sec)',fontsize=14)

    
    # plt.subplot(6,2,7)
    ax5 = plt.subplot2grid((13,5), (9,0), colspan=3, rowspan=4)
    trobs_filt = obs_strs_filt_3[i][1]
    trsyn_eugspe_1d = eugspe_1d_strs_3[i][1]
    trsyn_eugspe_sm = eugspe_sm_strs_3[i][1]
    trsyn_srv_1d = srv_1d_strs_3[i][1]
    trsyn_srv_sm = srv_sm_strs_3[i][1]
    trsteph = steph_strs_3[i][1]
    npts_eugspe_1d = trsyn_eugspe_1d.stats.npts
    npts_eugspe_sm = trsyn_eugspe_sm.stats.npts
    npts_srv_1d = trsyn_srv_1d.stats.npts
    npts_srv_sm = trsyn_srv_sm.stats.npts
    npts_steph = trsteph.stats.npts
    samprate_eugspe_1d = trsyn_eugspe_1d.stats.sampling_rate
    samprate_eugspe_sm = trsyn_eugspe_sm.stats.sampling_rate
    samprate_srv_1d = trsyn_srv_1d.stats.sampling_rate
    samprate_srv_sm = trsyn_srv_sm.stats.sampling_rate
    samprate_steph = trsteph.stats.sampling_rate
    t_eugspe_1d = np.arange(0, npts_eugspe_1d / samprate_eugspe_1d, 1 / samprate_eugspe_1d) 
    t_eugspe_sm = np.arange(0, npts_eugspe_sm / samprate_eugspe_sm, 1 / samprate_eugspe_sm)
    t_srv_1d = np.arange(0, npts_srv_1d / samprate_srv_1d, 1 / samprate_srv_1d)
    t_srv_sm = np.arange(0, npts_srv_sm / samprate_srv_sm, 1 / samprate_srv_sm)
    tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)
    
    plt.plot(tobs,trobs_filt.data, c = 'navy', label = 'Observed_filt',zorder=0,linewidth=linewidth)
    plt.plot(t_eugspe_1d,trsyn_eugspe_1d,  c = 'darkred', label  = 'eugspe_1d',alpha=0.7,linewidth=linewidth,ls='--')
    plt.plot(t_eugspe_sm,trsyn_eugspe_sm,  c = 'chocolate', label  = 'eugspe_sm',alpha=0.7,linewidth=linewidth,ls='--')
    plt.plot(t_srv_1d,trsyn_srv_1d,  c = 'gold', label  = 'srv_1d',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.plot(t_srv_sm,trsyn_srv_sm,  c = 'darkviolet', label  = 'srv_sm',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',alpha=0.7,linewidth=linewidth)
    plt.xlim([-5,80])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    # ax5.xaxis.set_ticklabels([])
    
    limobs = max(abs(trobs_filt.data))
    limsyn_eugspe_1d = max(abs(trsyn_eugspe_1d.data))
    limsyn_eugspe_sm = max(abs(trsyn_eugspe_sm.data))
    limsyn_srv_1 = max(abs(trsyn_srv_1d.data))
    limsyn_srv_sm = max(abs(trsyn_srv_sm.data))
    lim = max(limobs, limsyn_eugspe_1d , limsyn_eugspe_sm , limsyn_srv_1 , limsyn_srv_sm)
    
    plt.ylim([-1*lim, lim])
    plt.ylim([-1*lim, lim])
    plt.locator_params(tight=True, nbins=4)
    # plt.ylabel('vel (m/s)',fontsize=14)
    # plt.xlabel('time (sec)',fontsize=14)
     
    
    # plt.subplot(3,2,4)
    ax6 = plt.subplot2grid((13,5), (9,3), colspan=2, rowspan=4)
    plt.loglog(freq_obs,amps_obs, c = 'dimgrey' , label='Obsereved',alpha=0.7,linewidth=linewidth)
    plt.loglog(freq_obs_filt,amps_obs_filt, c = 'navy' , label='Obsereved filt',alpha=0.7,linewidth=linewidth)
    plt.loglog(freq_eugspe_1d,amps_eugspe_1d,  c = 'darkred' , label='eugspe_1d',alpha=0.7,linewidth=linewidth,ls='--')
    plt.loglog(freq_eugspe_sm,amps_eugspe_sm,  c = 'chocolate' , label='eugspe_sm',alpha=0.7,linewidth=linewidth,ls='--')
    plt.loglog(freq_srv_1d,amps_srv_1d,  c = 'gold' , label='srv_1d',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.loglog(freq_srv_sm,amps_srv_sm,  c = 'darkviolet' , label='srv_sm',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.loglog(freq_steph,amps_steph,  c = 'green' , label='Stephenson',alpha=0.7,linewidth=linewidth)
    # plt.xlabel('Freqs [Hz]',fontsize=14)
    # plt.ylabel('Amps',fontsize=14)
    # plt.legend(loc='upper right',prop={'size': 10.0,'weight':'bold'})
    plt.xlim(bot_lim,(samprate_obs/2))
    plt.ylim(1e-8,max(amps_obs)+5e-5)
    plt.grid(which="both", axis='x') 
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    # ax6.xaxis.set_ticklabels([])
    
    plt.subplots_adjust(hspace=0.6,wspace=1.1)
    # plt.tight_layout()
    
    plt.savefig(working_dir+name+'.png' ,dpi=150,bbox_inches='tight') 

plt.close('all')  
















