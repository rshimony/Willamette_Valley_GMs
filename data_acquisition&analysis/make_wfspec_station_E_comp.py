#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:51:58 2024

@author: rshimony
"""
'''
Makes ONE (East) component waveform and spectra comparison figures for ONE SPECIFIC new WV model 
(1 (out of 4) new WV model, USGS CVM and observed data) for each station in a single event.
This scripts takes ONE EVENT and ONE MODEL at a time.
The output figures has 3 panels:
    - the unfiltered waveform 
    - Filtered waveform comparison (synthetic and observed)
    - Spectra comparison
'''

from obspy.core.utcdatetime import UTCDateTime
import numpy as np
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

#####

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
####Inputs#####
# Station inventory 
st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/metadata/snr_stations.csv')
st_names = np.array(st_file['st_nm'])
st_lons = np.array(st_file['st_lon'])
st_lats = np.array(st_file['st_lat']) 
st_names_unq = np.unique(st_names)

# Waveform data directories
obs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/waveforms/observed_data_1hz_snr/'
steph_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/steph/'
synts_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/srv_sm/'

# Output directory
working_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_station_E_comp/springfield/srv_sm/'

# Model and event
model = 'SRV SM'
event = 'Springfield'

try:
    os.mkdir(working_dir)
except FileExistsError:
    print('Directory already exists')

#%%
# Getting strings
obs_strs_3,obs_strs_filt_3 = make_wf_strings(obs_dir,'obs',1) 
steph_strs_3 =  make_wf_strings(steph_dir,'synts',1)
synts_strs_3 = make_wf_strings(synts_dir,'synts',1) 
       
#%%
## Plotting parameters
bot_lim = 0.063 # Bottom limit of spectra plots
linewidth=2.0

for i in range(len(obs_strs_3)):
       
    # Initiating figure
    plt.figure(figsize=(12,4))
    gridspec.GridSpec(6,5)
    name = st_names_unq[i]
    plt.suptitle(name+' - '+event,size=20,y=0.98)
    plt.figtext(0.07,0.46,'vel (m/s)',fontsize=18,rotation = 90)
    plt.figtext(0.57,0.48,'Amps',fontsize=18,rotation = 90)
    
    # Calculating spectra for East component (index = 0)
    freq_obs , amps_obs = calc_spectra(obs_strs_3[i][0])
    freq_obs_filt , amps_obs_filt = calc_spectra(obs_strs_filt_3[i][0])
    freq_steph , amps_steph = calc_spectra(steph_strs_3[i][0])
    freq_synts , amps_synts = calc_spectra(synts_strs_3[i][0])

    # Plotting panel 1 - East unfiltered observed
    ax1 = plt.subplot2grid((6,5), (0,0), colspan=3, rowspan=2)
    trobs = obs_strs_3[i][0]
    t_start = trobs.stats.starttime
    samprate_obs = trobs.stats.sampling_rate
    tobs = trobs.times(reftime=UTCDateTime(t_start))
    plt.plot(tobs,trobs.data, c = 'dimgrey', label = 'Observed',linewidth=linewidth)
    plt.xlim([-5,80])
    plt.title('East' ,size=15)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yticks(fontsize=15)
    ax1.xaxis.set_ticklabels([])
    limobs = max(abs(trobs.data))
    plt.ylim([-1*limobs, limobs])
    plt.ylim([-1*limobs, limobs])
    plt.locator_params(tight=True, nbins=4)
    
    # Plotting panel 2 - East filtered wf comparison
    ax2 = plt.subplot2grid((6,5), (2,0), colspan=3, rowspan=4)
    trobs_filt = obs_strs_filt_3[i][0]
    trsyn = synts_strs_3[i][0]
    trsteph = steph_strs_3[i][0]
    npts_synts = trsyn.stats.npts
    npts_steph = trsteph.stats.npts
    samprate_synts = trsyn.stats.sampling_rate
    samprate_steph = trsteph.stats.sampling_rate
    t_synts = np.arange(0, npts_synts / samprate_synts, 1 / samprate_synts) 
    tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)
    plt.plot(tobs,trobs_filt.data, c = 'navy', label = 'Observed filt',zorder=0,linewidth=linewidth)
    plt.plot(t_synts,trsyn,  c = 'darkred', label  = model,alpha=0.7,linewidth=linewidth,ls='--')
    plt.plot(tsteph,trsteph,  c = 'green', label  = 'USGS CVM',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.xlim([-5,80])
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    limobs = max(abs(trobs_filt.data))
    limsyn_synts = max(abs(trsyn.data))
    lim = max(limobs, limsyn_synts)
    plt.ylim([-1*lim, lim])
    plt.ylim([-1*lim, lim])
    plt.locator_params(tight=True, nbins=4)
    plt.xlabel('time (sec)',fontsize=18)

    # Plotting panel 3 - East spectra comparison
    ax3 = plt.subplot2grid((6,5), (2,3), colspan=2, rowspan=4)
    plt.loglog(freq_obs,amps_obs, c = 'dimgrey' , label='Observed',alpha=0.7,linewidth=linewidth)
    plt.loglog(freq_obs_filt,amps_obs_filt, c = 'navy' , label='Observed filt',alpha=0.7,linewidth=linewidth)
    plt.loglog(freq_synts,amps_synts,  c = 'darkred' , label=model,alpha=0.7,linewidth=linewidth,ls='--')
    plt.loglog(freq_steph,amps_steph,  c = 'green' , label='USGS CVM',alpha=0.7,linewidth=linewidth,ls='-.')
    plt.legend(loc='upper center',bbox_to_anchor=(0.45,  1.6),prop={'size': 13.0,'weight':'bold'},ncol=2)
    plt.xlim(bot_lim,(samprate_obs/2))
    plt.ylim(1e-8,max(amps_obs)+5e-5)
    plt.grid(which="both", axis='x') 
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.xlabel('freq (Hz)',fontsize=18)
    
    # Subplots adjustments
    plt.subplots_adjust(hspace=0.6,wspace=1.1)
    # Saving figure
    plt.savefig(working_dir+name+'.png' ,dpi=150,bbox_inches='tight') 

plt.close('all')  