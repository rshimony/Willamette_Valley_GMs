#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:28:53 2024

@author: rshimony
"""

'''
Download waveforms from a specific earthquake.
Creating station inventory, event metadata, distance btwn event and station files.
Dowload and save raw waveform data, data after instrument correction (unfilt) and low-pass filtered data (filt)
Plot and save all stations' waveforms 
Check SNR of dowloaded waveforms, creates a new SNR station inventory and redownload
Finally, creates a Valley stastions inventory file.
'''

import numpy as np
import matplotlib.pyplot as plt
import obspy as obs
from obspy.clients.fdsn.client import Client
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel
from obspy.core.utcdatetime import UTCDateTime
import pandas as pd
import os
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

#%% Helpers

#### Creating WV polygon from data points stord in csv files.
# Outer outlines + two inner polygons concatenated together.
inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/wv_poly_outlines/basin_depth_poly_outer.csv')
lon_outer = outer_f['wv_poly_lons']
lat_outer = outer_f['wv_poly_lats']    
  
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

####

def comp_SNR(stream, noisedur, sigstart, sigdur):
    
    '''
    Takes an ObsPy stream and computes SNR
    '''
    
    SNR_lst = []
    
    for i in range(len(stream)):
    
        time = stream[i].times()
        data = stream[i].data
        noisestart = 1
        noiseend = noisestart + noisedur
        sigend = sigstart + sigdur
        # find indecies
        noise_ind = np.where((time<=noiseend) & (time>=noisestart))[0]
        noise_amp = data[noise_ind]
        noise_avg = float(np.std(noise_amp))

        sig_ind = np.where((time<=sigend) & (time>=sigstart))[0]
        sig_amp = data[sig_ind]
        sig_avg = float(np.std(sig_amp))

        if noise_avg == 0:
            noise_avg = 1E-10
        if sig_avg == 0:
            sig_avg = 1E-10

        SNR=sig_avg/noise_avg
        SNR_lst.append(SNR)

    
    
    return SNR_lst


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

#%%
############# PARAMETERS for event, station inventory and metadata #############
eq_name = 'salem_eq'

working_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/'+eq_name+'/'
meta_dir = working_dir+'/metadata/'

try:
    os.mkdir(working_dir)
except FileExistsError:
    print('Directory already exists')
try:
    os.mkdir(meta_dir)
except FileExistsError:
    print('Directory already exists')
    

## geographic region:
minlat = 43.4
maxlat = 46.2
minlon = -124.2
maxlon = -121.5

## time constraints for event search:
eventTime = UTCDateTime("2022-10-07T12:52:36.010000Z")
stime = eventTime - 60  # 1 minute before the event
etime = eventTime + 15*60  # 15 minutes after the event

## min depth in km:
min_depth = 1

## min magnitude:
min_mag = 3

## stations of interest:
netcode = '??'
channels = 'HH*,BH*,HN*,EH*,EN*' # high gain accelerometers, could do high gain broadbands ('HH*',)

## CLIENTS
## Set up the client to be IRIS, for downloading:
client = Client('IRIS')
# taup:
model = TauPyModel(model="iasp91")


#%%
############# CATALOG DOWNLOADS ##############

## Find stations metadata for position:
sta_inventory = client.get_stations(network=netcode, channel=channels, starttime=stime, endtime=etime, minlatitude = minlat , maxlatitude = maxlat , 
                                      minlongitude = minlon , maxlongitude = maxlon,level="channel")

## Find earthquakes 
eq_catalog = client.get_events(starttime=stime, endtime=etime,
                        minlatitude=minlat, maxlatitude=maxlat, 
                        minlongitude=minlon, maxlongitude=maxlon,
                        minmagnitude=min_mag,mindepth=min_depth)


#%%
######## GET CATALOG AND METADATA ###########
## Get all station inventory csv:
nt_stas = []
st_lons = []
st_lats = []
st_stas = []
st_elev = []
st_chan = []
ch_samprt = []
for network in sta_inventory:
    for station in network:
        for channel in station:
            nt_stas.append(network.code)
            st_stas.append(station.code)
            st_lons.append(station.longitude)
            st_lats.append(station.latitude)
            st_elev.append(station.elevation)
            st_chan.append(channel.code)
            ch_samprt.append(channel.sample_rate)
stdict = {'network':nt_stas, 'name':st_stas, 'longitude':st_lons,
          'latitude':st_lats, 'elevation':st_elev, 'channel':st_chan, 'samprt':ch_samprt }
stdf = pd.DataFrame(stdict)
stdf_nodup=stdf.drop_duplicates()
stdf_nodup.to_csv(meta_dir + 'station_inventory.csv',index=False)
        
# Extract event information to csv:
ev_lons = []
ev_lats = []
ev_depth = []
ev_origint = []
ev_mag = []
ev_id = []
for event in eq_catalog:
    ev_lons.append(event.origins[0].longitude)
    ev_lats.append(event.origins[0].latitude)
    ev_depth.append(event.origins[0].depth)
    ev_origint.append(event.origins[0].time)
    ev_mag.append(event.magnitudes[0].mag)
    evid_str=str(event.resource_id)
    evid_split = evid_str.split('/')[-1]
    evid = evid_split.split('=')[-1]
    ev_id.append(evid)
evdict = {'catalog #':ev_id, 'longitude':ev_lons, 'latitude':ev_lats,'magnitude':ev_mag,
          'depth':ev_depth,'time':ev_origint}
evdf = pd.DataFrame(evdict)
evdf.to_csv(meta_dir + 'event_catalog.csv',index=False)

#%%
######## get distances btwn stations and events ###########
## make dataframe with station repeated
distances_all = []
sta_df_all = []
net_df_all = []
chan_df_all = []
ev_mag_df_all = []
ev_origintime_df_all = []
ev_collecttime_all = []
ev_endtime_all = []
evdf_lat_all = []
evdf_lon_all = []
evdf_depth_all = []

for i_station in range(len(st_stas)):
    i_stlon = st_lons[i_station]
    i_stlat = st_lats[i_station]
    i_network = nt_stas[i_station]
    i_channel = st_chan[i_station]
    for event in eq_catalog:    
        ## add station name to long array:
        sta_df_all.append(st_stas[i_station])
        net_df_all.append(nt_stas[i_station])
        chan_df_all.append(st_chan[i_station])
        ## append event info:
        evdf_lat_all.append(event.origins[0].latitude)
        evdf_lon_all.append(event.origins[0].longitude)
        evdf_depth_all.append(event.origins[0].depth)
        
        ## get great circle distance between this event, and station:
        ij_grcircle_distance = gps2dist_azimuth(i_stlat, i_stlon, event.origins[0].latitude, event.origins[0].longitude)
        ij_grcircle_distance_km = ij_grcircle_distance[0]/1000
        ij_degrees = kilometers2degrees(ij_grcircle_distance_km)
        ## For each distance, assume 5 km/s average p-wave speed to get p-wave travel time:
        ij_pwv_ttime_allP = model.get_travel_times(source_depth_in_km=event.origins[0].depth/1000,
                                  distance_in_degree=ij_degrees,phase_list=["p","P"])
        ij_pwv_ttime = ij_pwv_ttime_allP[0].time
        
        ## start collecting waveform data 10 seconds before predicted p-wave arrival
        ev_collecttime_all.append(event.origins[0].time + (ij_pwv_ttime - 10))
        ## collect a total of 120 seconds
        ev_endtime_all.append(event.origins[0].time + (ij_pwv_ttime - 10) + 120)
        
        ## append:
        distances_all.append(ij_grcircle_distance_km)
        ev_mag_df_all.append(event.magnitudes[0].mag)
        ev_origintime_df_all.append(event.origins[0].time)
        
    print('st # ' + str(i_station) + ' of ' +str(len(st_stas)) + ' is done')
        
## Save then into a dictionary/dataframe, and sdave to file:
metadata_all = {'distance':distances_all,'network':net_df_all, 'stations':sta_df_all , 'channel': chan_df_all, 'mag':ev_mag_df_all, 
                'origint':ev_origintime_df_all, 'collecttime':ev_collecttime_all, 'endtime':ev_endtime_all,
                'evlat':evdf_lat_all, 'evlon':evdf_lon_all, 'evdepth':evdf_depth_all}

metadata_all_df = pd.DataFrame(metadata_all)
metadata_all_df.to_csv((meta_dir + 'metadata_all.csv'),index=False)

#%%
############# PARAMETERS for collecting, dowload and save waveforms #############
#change these variables for different model runs
examplename = 'observed_data_1hz' # Waveforms data directory
freqmax = 1 # max frequency of low pass filter used in pre-filter
pre_filt = (0.063 , 0.28 , 0.9*freqmax, freqmax)
resp_prefilt_bottom = [0.063 , 0.28] # pre-filter for full waveform (no low-pass filter)

## Initiating waveforms and figures directories
all_wf_dir = working_dir+'/waveforms/'
wf_dir = working_dir+'/waveforms/'+examplename+'/'
all_figs_dir = working_dir+'/figs/'
figs_dir = working_dir+'/figs/'+examplename+'/'

try:
    os.mkdir(all_wf_dir)
except FileExistsError:
    print('Directory already exists')
try:
    os.mkdir(wf_dir)
except FileExistsError:
    print('Directory already exists')
try:
    os.mkdir(all_figs_dir)
except FileExistsError:
    print('Directory already exists')
try:
    os.mkdir(figs_dir)
except FileExistsError:
    print('Directory already exists')

## Reading in the newly made station and event files
station_file = meta_dir+'station_inventory.csv'
event_file = meta_dir+'event_catalog.csv'
event = pd.read_csv(event_file, header = 0)
stndata = pd.read_csv(station_file, header = 0)

## Read station and event metadata
ev_t = event['time'][0]

net = np.array(stndata["network"])
stnm = np.array(stndata["name"])
chnm = np.array(stndata["channel"])
lat = np.array(stndata["latitude"])
lon = np.array(stndata["longitude"])

## Setting up time parameters for waveform fetching
# time for full station inventory
t1 = UTCDateTime(ev_t) - 20
t2 = t1 + (120)#delta in seconds
# time for SNR stations 
t1_val = UTCDateTime(ev_t)
t2_val = t1_val + (75)#delta in seconds

## SNR Parameters

noise_dur = 15 #[sec]
signal_start = 40 #[sec]
signal_dur = 20 #[sec]

snr_threshold = 8
        

#%%
# Fetch waveform from IRIS FDSN web service into a ObsPy stream object
# and automatically attach correct response
# Dowload and save raw waveform data, data after instrument correction (unfilt) and low-pass filtered data (filt)
# Plot and save all stations' waveforms 

dlsuccess = [False] * len(net)

for i in range(len(net)):
    try:
        st = client.get_waveforms(network=net[i], station=stnm[i], location='*',
                                       channel=chnm[i], starttime=t1, endtime=t2,
                                       attach_response=True)
        dlsuccess[i] = True
    except:
        dlsuccess[i] = False
        print(str(i) + ' has nothing')
        
    if dlsuccess[i]==True:

        tr = st[0]
        chan = tr.stats.channel
        t = tr.times(reftime=UTCDateTime(ev_t))
        print(tr)

        tr_raw = tr.copy()
        icorr_sacfile_raw = wf_dir + net[i] + '.' + stnm[i] + '.' + chan + '_raw.sac'
        tr_raw.write(icorr_sacfile_raw, format = 'sac')
        data_raw = tr_raw.data

        tr_unfilt = tr.copy()
        tr_unfilt.detrend(type = 'linear') 
        samprate = tr_unfilt.stats.sampling_rate
        prefilt_unfilt = (resp_prefilt_bottom[0], resp_prefilt_bottom[1], ((samprate/2)-5), (samprate/2))
        tr_unfilt.remove_response(output='VEL' , pre_filt=prefilt_unfilt)
        icorr_sacfile_unfilt = wf_dir + net[i] + '.' + stnm[i] + '.' + chan + '_unfilt.sac'
        tr_unfilt.write(icorr_sacfile_unfilt, format = 'sac')
        data_unfilt = tr_unfilt.data

        tr_filt = tr.copy()
        tr_filt.detrend(type = 'linear')                
        tr_filt.remove_response(output='VEL', pre_filt=pre_filt)
        icorr_sacfile_filt = wf_dir + net[i] + '.' + stnm[i] + '.' + chan + '_filt.sac'
        tr_filt.write(icorr_sacfile_filt, format = 'sac')
        data_filt = tr_filt.data

        fig, axs = plt.subplots(1,3, sharex=True,figsize = (14,2))
        axs[0].plot(t, data_raw)
        axs[1].plot(t, data_unfilt)
        axs[2].plot(t,data_filt)
        axs[0].set_title('raw data')
        axs[1].set_title('data instrument response corrected')
        axs[2].set_title('data filtered')
        plt.tight_layout()
        plt.xlim([-5,80])
        plt.savefig(figs_dir + net[i] + '.' + stnm[i]  + '.' + chan +'.png')
        plt.close()

#%%
## Check SNR of all waveforms dowloaded and make an SNR stations inventory file

# read all unfiltered waveforms
data_strm = obs.read(wf_dir + '*_unfilt.sac')

# Get unique stations' name
data_stations = []
for i in range(len(data_strm)):
    data_stations.append(data_strm[i].stats.station) 
data_stations_unq = pd.unique(data_stations)

# Make 3 components ObsPy stringsfor all stations
strings = []
for i in range(len(data_stations_unq)):
    string = obs.read(wf_dir + '*' + data_stations_unq[i] + '*_unfilt.sac')
    strings.append(string)
strs_3 = []
for i in range(len(strings)):
    if len(strings[i]) >= 3:
        strs_3.append(strings[i][-3:])

# Compute SNR and append SNR strings to a new list
snr_strs = []
snr_strs_nm = []

for i in range(len(strs_3)):
    snr_str = comp_SNR(strs_3[i], noise_dur, signal_start, signal_dur)
    snr_avg_str = find_vec_norm(snr_str) 
    snr_name = strs_3[i][0].stats.station
    if snr_avg_str >= 8:
        snr_strs.append(strs_3[i])
        snr_strs_nm.append(snr_name)

# Get SNR stations location and name
st_lons_uniq = pd.unique(st_lons)
st_lats_uniq = pd.unique(st_lats)
st_name_uniq = pd.unique(st_stas)
snr_lons = []
snr_lats = []
snr_name_unq = []
for i in range(len(snr_strs_nm)):
    for j in range(len(st_lons_uniq)):
        if snr_strs_nm[i] == st_name_uniq[j]:
            snr_lons.append(st_lons_uniq[j])
            snr_lats.append(st_lats_uniq[j])
            snr_name_unq.append(st_name_uniq[j])
            
# Save SNR stations to csv            
ldict = {'st_nm':snr_name_unq, 'st_lon':snr_lons, 'st_lat':snr_lats,}
ldf = pd.DataFrame(ldict)
ldf.to_csv(meta_dir + 'snr_stations.csv',index=False)

#%%
# Fetch ONLY SNR waveform
# Dowload and save raw waveform data, data after instrument correction (unfilt) and low-pass filtered data (filt)
# Plot and save all stations' waveforms 

wf_snr_dir = all_wf_dir+examplename+'_snr/'
figs_snr_dir = all_figs_dir+examplename+'_snr/'

try:
    os.mkdir(wf_snr_dir)
except FileExistsError:
    print('Directory already exists')
try:
    os.mkdir(figs_snr_dir)
except FileExistsError:
    print('Directory already exists')


for i in range(len(snr_name_unq)):
    for j in range(len(st_stas)):
        if snr_name_unq[i] == st_stas[j]:

            st = client.get_waveforms(network=nt_stas[j], station=st_stas[j], location='*',
                                           channel=st_chan[j], starttime=t1_val, endtime=t2_val,
                                           attach_response=True)

            tr = st[0]
            chan = tr.stats.channel
            t = tr.times(reftime=UTCDateTime(ev_t))
            print(tr)

            tr_raw = tr.copy()
            icorr_sacfile_raw = wf_snr_dir + nt_stas[j] + '.' + st_stas[j] + '.' + chan + '_raw.sac'
            tr_raw.write(icorr_sacfile_raw, format = 'sac')
            data_raw = tr_raw.data

            tr_unfilt = tr.copy()
            tr_unfilt.detrend(type = 'linear') 
            samprate = tr_unfilt.stats.sampling_rate
            prefilt_unfilt = (resp_prefilt_bottom[0], resp_prefilt_bottom[1], ((samprate/2)-5), (samprate/2))
            tr_unfilt.remove_response(output='VEL' , pre_filt=prefilt_unfilt)
            icorr_sacfile_unfilt = wf_snr_dir + nt_stas[j] + '.' + st_stas[j] + '.' + chan + '_unfilt.sac'
            tr_unfilt.write(icorr_sacfile_unfilt, format = 'sac')
            data_unfilt = tr_unfilt.data

            tr_filt = tr.copy()
            tr_filt.detrend(type = 'linear')                
            tr_filt.remove_response(output='VEL', pre_filt=pre_filt)
            icorr_sacfile_filt = wf_snr_dir + nt_stas[j] + '.' + st_stas[j] + '.' + chan + '_filt.sac'
            tr_filt.write(icorr_sacfile_filt, format = 'sac')
            data_filt = tr_filt.data

            fig, axs = plt.subplots(1,3, sharex=True,figsize = (14,2))
            axs[0].plot(t, data_raw)
            axs[1].plot(t, data_unfilt)
            axs[2].plot(t,data_filt)
            axs[0].set_title('raw data')
            axs[1].set_title('data instrument response corrected')
            axs[2].set_title('data filtered')
            plt.tight_layout()
            plt.xlim([-5,80])
            plt.savefig(figs_snr_dir + nt_stas[j] + '.' + st_stas[j]  + '.' + chan +'.png')
            plt.close()
#%%
## Create Valley SNR stations inventory
# Read the newly made SNR station inventory file
snr_stations_ff = pd.read_csv(meta_dir+'snr_stations.csv')
snr_st_lat = snr_stations_ff['st_lat']
snr_st_lon = snr_stations_ff['st_lon']
snr_st_name = snr_stations_ff['st_nm']

# Creating a buffered WV polygon to get all valley edge stations
buffered_fullpoly = full_poly.buffer(0.1, join_style=2)

# Getting SNR valley station from the full SNR station inventory
valst_lons = []
valst_lats = []
valst_names = []

for i in range(len(snr_st_lat)):
    grd_point = Point(snr_st_lon[i],snr_st_lat[i])
    if buffered_fullpoly.contains(grd_point) == True:
        
        valst_lons.append(snr_st_lon[i])
        valst_lats.append(snr_st_lat[i])
        valst_names.append(snr_st_name[i])
        
# Saving new Valley SNR station inventory file
dom_valstdict = {'dom_valst_names':valst_names, 'dom_valst_lons':valst_lons, 'dom_valst_lats':valst_lats}
dom_valstdf = pd.DataFrame(dom_valstdict)
dom_valstdf_nodup=dom_valstdf.drop_duplicates()
dom_valstdf_nodup.to_csv(meta_dir + 'valley_snr_st.csv',index=False)

















