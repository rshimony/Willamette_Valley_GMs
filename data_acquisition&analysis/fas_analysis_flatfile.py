#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:20:27 2024

@author: rshimony
"""
'''
Calculates Fourier Amplitude Spectra (FAS) binned to specific frequencies (0.2, 0.3, 0.4, 0.5, 0.6, 0.7 Hz)
from waveform data (observed, new WV model synthetics, USGS CVM synthetics).
Saving the calculated binned FAS in a flatfile.
This script works for ONE EVENT and ONE MODEL at a time:
    - Observed data
    - USGS CVM synthetics
    - 1 new WV model
The output flatfile holds binned FAS for the three datasets.
Each IM is calculated for each component seperatly and written to the flatfile.
**Before the calculation part, need to set DATA STRINGS type:
    - Observed filtered
    - USGS CVM (steph)
    - Model synthetics
Finally, creates a Valley stations only FAS flatfile, out of the newly made full station inventry FAS flatfile.
**To get all data, this script need to run 6 times for each event:
    1)Observed filtered
    2)USGS CVM
    3-6)Model synthetics (EugSpe 1D, EugSpe SM, SRV 1D, SRV SM)
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
    
def calc_spectra(trace, data_type):
    """
    **Adapted from Nye et al. (2023)**
    Calculates average spectra values in 20 bins for displacement, acceleration,
    and velocity waveforms.

    Inputs:
        stream: Obspy stream. (change by rshimony to trace)
        data_type(str): Data type used for determining bin edges.
            Options:
                disp
                sm
                vel

    Return:
        bin_means(list): Binned spectra for given station
        freq(list): FFT frequencies
        amp(list): FFT amplitudes
    """
    
    import numpy as np
    from mtspec import mtspec
    from scipy import interpolate
    from scipy.stats import binned_statistic 

    # Read in file 
    tr = trace
    data = tr.data
    delta = tr.stats.delta
    samprate = tr.stats.sampling_rate
    npts = tr.stats.npts
    
    # Determine nyquist frequency
    nyquist = 0.5 * samprate
    

    # Calc spectra amplitudes and frequencies 
        # Switched number of tapers from 7 to 5.  Decreases computation time and
        # results are similar
    amp_squared, freq =  mtspec(data, delta=delta, time_bandwidth=4, 
                              number_of_tapers=5, nfft=npts, quadratic=True)
    
    # Convert from power spectra to amplitude spectra
    amp = np.sqrt(amp_squared)
    
    # Use scipy interpolate function to fill in data in missing bins
    f = interpolate.interp1d(freq, amp)
    freq_new = np.arange(np.min(freq), np.max(freq), 0.0001)
    amp_new = f(freq_new)

    # Remove certain frequencies that are too low or high. 
    indexes = []
    
    for i, val in enumerate(freq_new):
        
        # Remove frequencies below 1/2 length of record
        if val <= 1/(delta*npts*0.5) :
            indexes.append(i)
        
        # Remove frequencies above 10 Hz for sm data because of the way it was processed 
        elif val > 10 and data_type == 'sm':
            indexes.append(i)

        # Remove frequencies above nyquist frequency for disp data
            # (it's already removed in the previous step for sm data)
        elif val > nyquist and data_type == 'disp': 
            indexes.append(i)
    
    # Remove any duplicate indexes
    indexes = np.unique(indexes)
    freq_new = np.delete(freq_new,indexes)
    amp_new = np.delete(amp_new,indexes) 
    
    # Set up bins
    bins = np.logspace(np.log10(0.03), np.log10(1), num=50)
#     if data_type == 'sm':
#         # Starting bins at 0.004 Hz (that is about equal to half the length
#             # of the record for the synthetic and observed data) and ending at
#             # 10 Hz because after that the sm data is unusable due to how it was
#             # processed. 
#         bins = np.logspace(np.log10(0.004), np.log10(10), num=21)
    
#     elif data_type == 'disp':
#         # Starting bins at 0.004 Hz (that is about equal to half the length
#             # of the record for the synthetic and observed data) and ending at
#             # 0.5 Hz because that is the nyquist frequency .
#         bins = np.logspace(np.log10(0.004), np.log10(0.5), num=21)
    
#     else:
#         bins = np.array([0.1,0.2,0.5,1])
    
    bin_means, bin_edges, binnumber = binned_statistic(freq_new,
                                                       amp_new,
                                                       statistic='mean',
                                                       bins=bins)
                
    # for i in range(len(bin_means)):
    #     bin_means[i] = 10**bin_means[i]
        
        
    return(bin_means, freq, amp )
#%%
###INPUTS###
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
eugspe_1d_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/eugspe_1d/' # New WV model
eugspe_sm_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/eugspe_sm/' # New WV model
srv_1d_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/srv_1d/' # New WV model
srv_sm_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/synthetic_data/sac_files/springfield/srv_sm/' # New WV model

# outpur directory
outdir = '/Users/rshimony/Desktop/WillametteValley/wv_project/fas_analysis/springfield/'

#%%
# Getting strings
obs_strs_3,obs_strs_filt_3 = make_wf_strings(obs_dir,'obs',1.0) 
steph_strs_3 =  make_wf_strings(steph_dir,'synts',1.0)
eugspe_1d_strs_3 = make_wf_strings(eugspe_1d_dir,'synts',1.0) 
eugspe_sm_strs_3 = make_wf_strings(eugspe_sm_dir,'synts',1.0) 
srv_1d_strs_3 = make_wf_strings(srv_1d_dir,'synts',1.0) 
srv_sm_strs_3 = make_wf_strings(srv_sm_dir,'synts',1.0) 

#%%
## Calculating FAS bins and saving them to a .csv file

def make_fas_ff(data_strings,ff_name):
    counter_fas = 1
    for i in range(len(obs_strs_filt_3)):
        if counter_fas % 1 == 0:
            print(str(counter_fas) + ' from ' + str(len(obs_strs_filt_3)))
        for j in range(3):
            trobs = data_strings[i][j]
            chan = trobs.stats.channel
            
            bin_m , freqs , amps = calc_spectra(trobs,'vel')
            
            data_fas = [[st_names_unq[i], chan ,st_lons_unq[i], st_lats_unq[i] , 
                     bin_m[27] , bin_m[32] , bin_m[36] , bin_m[39] , bin_m[42] , bin_m[44]] ]
            
            if i == 0 and j == 0:
                df_fas = pd.DataFrame(data=data_fas,columns=['st_name', 'channel', 'st_lon', 'st_lat'
                                                                       , 'f=0.2', 'f=0.3', 'f=0.4', 'f=0.5' , 'f=0.6' , 'f=0.7'])
            else:
                df_fas_temp = pd.DataFrame(data=data_fas,columns=['st_name', 'channel', 'st_lon', 'st_lat'
                                                                       , 'f=0.2', 'f=0.3', 'f=0.4', 'f=0.5' , 'f=0.6' , 'f=0.7'])
                df_fas = pd.concat([df_fas, df_fas_temp], ignore_index=True)
    
    
        counter_fas = counter_fas + 1
    
    df_fas.to_csv(outdir+ff_name,index=False)

make_fas_ff(obs_strs_filt_3,'fas_ff_obs_springfield.csv')
make_fas_ff(steph_strs_3,'fas_ff_steph_springfield.csv')
make_fas_ff(eugspe_1d_strs_3,'fas_ff_eugspe_1d_springfield.csv')
make_fas_ff(eugspe_sm_strs_3,'fas_ff_eugspe_sm_springfield.csv')
make_fas_ff(srv_1d_strs_3,'fas_ff_srv_1d_springfield.csv')
make_fas_ff(srv_sm_strs_3,'fas_ff_srv_sm_springfield.csv')

#%%
# Creating a Valley SNR FAS flatfile
def make_valley_snr_fas_ff(val_ff_name):
    read_ff = pd.read_csv(outdir+val_ff_name)
    stations_names = np.array(read_ff['st_name'])
    
    val_st_file = pd.read_csv(valst_file)
    dom_valst_names = np.array(val_st_file['dom_valst_names'])
    
    valst_ff_mask = np.isin(stations_names,dom_valst_names)
    valst_ff = read_ff[valst_ff_mask]
    
    valst_ff.to_csv(outdir+'val_'+val_ff_name,index=False)

make_valley_snr_fas_ff('fas_ff_obs_springfield.csv')
make_valley_snr_fas_ff('fas_ff_steph_springfield.csv')
make_valley_snr_fas_ff('fas_ff_eugspe_1d_springfield.csv')
make_valley_snr_fas_ff('fas_ff_eugspe_sm_springfield.csv')
make_valley_snr_fas_ff('fas_ff_srv_1d_springfield.csv')
make_valley_snr_fas_ff('fas_ff_srv_sm_springfield.csv')








































