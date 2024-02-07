#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 10:48:42 2023

@author: rshimony
"""

import numpy as np 
import obspy as obs
import pandas as pd
from glob import glob
# import sys
# np.set_printoptions(threshold=sys.maxsize)
#%%

st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/springfield_eq/metadata/snr_stations_fixed.csv')
st_names = np.array(st_file['st_nm'])
st_lons = np.array(st_file['st_lon'])
st_lats = np.array(st_file['st_lat']) 
st_names_unq = pd.unique(st_names)
st_lons_unq = pd.unique(st_lons)
st_lats_unq = pd.unique(st_lats)

valst_file = '/Users/rshimony/Desktop/WillametteValley/springfield_eq/metadata/valley_snr_st.csv'

obs_dir = '/Users/rshimony/Desktop/WillametteValley/springfield_eq/waveforms/observed_data_1hz_snr/'
steph_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/springfield/rfile_steph_springfield_preempt_freq25/'
synts_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/final_sims_sac_files/springfield/rfile_steph_basin_springfield_srv_sm_preempt_freq25/'

outdir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/fas/springfield/'
ff_name = 'fas_ff_steph_springfield.csv'

obs_trs = sorted(glob(obs_dir + '*_unfilt.sac'))
# obs_trs_filt = sorted(glob(obs_dir + '*_filt.sac'))
synts_trs = sorted(glob(synts_dir + '*.*v'))
steph_trs = sorted(glob(steph_dir + '*.*v'))

obs_nms = []
for i in range(len(obs_trs)):
    obs_nm = obs_trs[i].split('.')[1]
    obs_nms.append(obs_nm)
obs_nms_unq = pd.unique(obs_nms)
    
synt_nms = []
for i in range(len(synts_trs)):   
    synt_nm = (synts_trs[i].split('/')[-1]).split('.')[0]
    synt_nms.append(synt_nm)
synt_nms_unq = pd.unique(synt_nms)

steph_nms = []
for i in range(len(steph_trs)):   
    steph_nm = (steph_trs[i].split('/')[-1]).split('.')[0]
    steph_nms.append(steph_nm)
steph_nms_unq = pd.unique(steph_nms)

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

def calc_spectra(trace, data_type):
    """
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

# binms , frequs , ampls = calc_spectra(obs_strs_filt_3[0][0],'vel')


#%%

counter_fas = 1
for i in range(len(obs_strs_filt_3)):
    if counter_fas % 1 == 0:
        print(str(counter_fas) + ' from ' + str(len(obs_strs_filt_3)))
    for j in range(3):
        trobs = steph_strs_3[i][j]
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


#%%

read_ff = pd.read_csv(outdir+ff_name)
stations_names = np.array(read_ff['st_name'])

val_st_file = pd.read_csv(valst_file)
dom_valst_names = np.array(val_st_file['dom_valst_names'])

valst_ff_mask = np.isin(stations_names,dom_valst_names)
valst_ff = read_ff[valst_ff_mask]

valst_ff.to_csv(outdir+'val_'+ff_name,index=False)










































