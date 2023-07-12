#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:39:39 2023

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import obspy as obs
from obspy.clients.fdsn.client import Client
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel
from obspy.core.utcdatetime import UTCDateTime
import pandas as pd
import os
from scipy import signal
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

#%%
############# PARAMETERS #############
eq_name = 'portorford_eq'

infile = '/Users/rshimony/Desktop/WillametteValley/Models/infile_rfile_steph_portorford.txt'

working_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/'
meta_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/metadata/'

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

## OR central longitude/latitude and radius:
central_lon = -123.05
central_lat = 44.8
search_radius_min = 0 ## in degrees
search_radius_max = 2.5 ## in degrees

## time constraints for event search:
eventTime = UTCDateTime("2019-11-30T01:45:12.0Z")
stime = eventTime - 60  # 1 minute before the event
etime = eventTime + 15*60  # 15 minutes after the event
# stime = UTCDateTime("2000-01-01T000000.000")
# etime = UTCDateTime("2021-9-20T235959.000")

## min depth in km:
min_depth = 1

## stations of interest:
netcode = '??'
# all_stations = 'BUCK'
channels = 'HH*,BH*,HN*,EH*,EN*' # high gain accelerometers, could do high gain broadbands ('HH*',)

## CLIENTS
## Set up the client to be IRIS, for downloading:
client = Client('IRIS')
# taup:
model = TauPyModel(model="iasp91")

#%% Helpers

inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_outer.csv')
lon_outer = outer_f['outer_lon']
lat_outer = outer_f['outer_lat']   
  
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

xpn_m, ypn_m = full_poly.exterior.xy 
xph1_m, yph1_m = full_poly.interiors[0].xy
xph2_m, yph2_m = full_poly.interiors[1].xy

xpnn = np.array([xpn_m])
xph1n = np.array([xph1_m])
xph2n = np.array([xph2_m])

ypnn = np.array([ypn_m])
yph1n = np.array([yph1_m])
yph2n = np.array([yph2_m])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

def comp_SNR(stream, noisedur, sigstart, sigdur):
    
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

    '''
    import numpy as np

    comp_1 = Im_value_list[0]**2;
    comp_2 = Im_value_list[1]**2;
    comp_3 = Im_value_list[2]**2;

    vec_norm = np.sqrt(comp_1 + comp_2 + comp_3)

    return(vec_norm)

#%%
############# CATALOG DOWNLOADS ##############

## Find stations metadata for position:
sta_inventory = client.get_stations(network=netcode, channel=channels, starttime=stime, endtime=etime, minlatitude = minlat , maxlatitude = maxlat , 
                                      minlongitude = minlon , maxlongitude = maxlon,level="channel")

## Find earthquakes 
eq_catalog = client.get_events(starttime=stime, endtime=etime,
                        minlatitude=minlat, maxlatitude=maxlat, 
                        minlongitude=minlon, maxlongitude=maxlon,
                        minmagnitude=3,mindepth=min_depth)

# eq_catalog = client.get_events(starttime=stime, endtime=etime,
#                         latitude=central_lat, longitude=central_lon, 
#                         minradius=search_radius_min, maxradius=search_radius_max,
#                         minmagnitude=3,mindepth=min_depth)


#%%
######## GET CATALOG AND METADATA ###########
## Extract the positions of the infrasound stations:
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
            # print(f'{channel.code[0:-1]}N , {channel.code[0:-1]}2')
stdict = {'network':nt_stas, 'name':st_stas, 'longitude':st_lons,
          'latitude':st_lats, 'elevation':st_elev, 'channel':st_chan, 'samprt':ch_samprt }
stdf = pd.DataFrame(stdict)
stdf_nodup=stdf.drop_duplicates()
stdf_nodup.to_csv(meta_dir + 'station_inventory.csv',index=False)

# nt_stas = np.array(nt_stas)
# st_lons = np.array(st_lons)
# st_lats = np.array(st_lats)
# st_stas = np.array(st_stas)
# st_elev = np.array(st_elev)
# st_chan = np.array(st_chan)
# ch_samprt = np.array(ch_samprt)

## Make a .csv of stations around the basin for the validation process
buffered_fullpoly = full_poly.buffer(0.1, join_style=2)

valst_nt = []
valst_lons = []
valst_lats = []
valst_names = []
valst_elev = []
valst_chan = []
valst_ch_samprt = []

for i in range(len(st_lons)):
    grd_point = Point(st_lons[i],st_lats[i])
    if buffered_fullpoly.contains(grd_point) == True:
        valst_nt.append(nt_stas[i])
        valst_lons.append(st_lons[i])
        valst_lats.append(st_lats[i])
        valst_names.append(st_stas[i])
        valst_elev.append(st_elev[i])
        valst_chan.append(st_chan[i])
        valst_ch_samprt.append(ch_samprt[i])
        

valstdict = {'network':valst_nt, 'name':valst_names, 'longitude':valst_lons,
          'latitude':valst_lats, 'elevation':valst_elev, 'channel':valst_chan, 'samprt':valst_ch_samprt }
valstdf = pd.DataFrame(valstdict)
valstdf_nodup=valstdf.drop_duplicates()
valstdf_nodup.to_csv(meta_dir + 'ploy_st_inventory.csv',index=False)
        
# Extract event information:
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
    
## WRite out event metadata to file to use later:
evdict = {'catalog #':ev_id, 'longitude':ev_lons, 'latitude':ev_lats,'magnitude':ev_mag,
          'depth':ev_depth,'time':ev_origint}
evdf = pd.DataFrame(evdict)
evdf.to_csv(meta_dir + 'event_catalog.csv',index=False)

#%%

plt.figure(figsize=(8,10))
plt.scatter(xpnt,ypnt,c='r',s=2)
plt.scatter(st_lons,st_lats, c='b')
plt.scatter(valst_lons,valst_lats, c='g')
plt.show()



#%%
######## get distances btwn stations and events ###########
## make dataframe with station repeated, use this info for dl-ing waveforms
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
        #ij_ttime = IRISclient.traveltime(model='iasp91', phases=['p'], evdepth=event.origins[0].depth, distkm=[ij_grcircle_distance_km], traveltimeonly=True)
        ## For each distance, assume 5 km/s average p-wave speed to get p-wave travel time:
        #ij_pwv_ttime = ij_grcircle_distance_km / 5.
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

#change these variables for different model runs
examplename = 'observed_data_1hz'
freqmax = 1
pre_filt = (0.063 , 0.28 , 0.9*freqmax, freqmax)
resp_prefilt_bottom = [0.063 , 0.28]

all_wf_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/waveforms/'
wf_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/waveforms/'+examplename+'/'
all_figs_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/figs/'
figs_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/figs/'+examplename+'/'

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

station_file = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/metadata/station_inventory.csv'
event_file = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/metadata/event_catalog.csv'

event = pd.read_csv(event_file, header = 0)
ev_t = event['time'][0]

# Need a list of earthquakes and stations
stndata = pd.read_csv(station_file, header = 0)
net = np.array(stndata["network"])
stnm = np.array(stndata["name"])
chnm = np.array(stndata["channel"])
# loc = np.array(['*']* len(net))
lat = np.array(stndata["latitude"])
lon = np.array(stndata["longitude"])

t1 = UTCDateTime(ev_t) - 20
t2 = t1 + (120)#delta in seconds
t1_val = UTCDateTime(ev_t)
t2_val = t1_val + (75)#delta in seconds

#%%
# Fetch waveform from IRIS FDSN web service into a ObsPy stream object
# and automatically attach correct response

dlsuccess = [False] * len(net)

for i in range(len(net)):
    try:
        st = client.get_waveforms(network=net[i], station=stnm[i], location='*',
                                       channel=chnm[i], starttime=t1, endtime=t2,
                                       attach_response=True)
        dlsuccess[i] = True
    except:
        dlsuccess[i] = False
        print(np.str(i) + ' has nothing')
        
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
#Check SNR and 

data_strm = obs.read(wf_dir + '*_unfilt.sac')

data_stations = []

for i in range(len(data_strm)):
    data_stations.append(data_strm[i].stats.station)
    
data_stations_unq = np.unique(data_stations)

strings = []

for i in range(len(data_stations_unq)):
    string = obs.read(wf_dir + '*' + data_stations_unq[i] + '*_unfilt.sac')
    strings.append(string)

strs_3 = []

for i in range(len(strings)):
    if len(strings[i]) >= 3:
        strs_3.append(strings[i][-3:])
    
snr_strs = []
snr_strs_nm = []

for i in range(len(strs_3)):
    snr_str = comp_SNR(strs_3[i], 15, 40, 20)

    snr_avg_str = find_vec_norm(snr_str)
    
    snr_name = strs_3[i][0].stats.station
    
    if snr_avg_str >= 8:
        snr_strs.append(strs_3[i])
        snr_strs_nm.append(snr_name)

st_lons_uniq = np.unique(st_lons)
st_lats_uniq = np.unique(st_lats)
st_name_uniq = np.unique(st_stas)

snr_lons = []
snr_lats = []
snr_name_unq = []

for i in range(len(snr_strs_nm)):
    for j in range(len(st_lons_uniq)):
        if snr_strs_nm[i] == st_name_uniq[j]:
            snr_lons.append(st_lons_uniq[j])
            snr_lats.append(st_lats_uniq[j])
            snr_name_unq.append(st_name_uniq[j])
            
            
ldict = {'st_nm':snr_name_unq, 'st_lon':snr_lons, 'st_lat':snr_lats,}
ldf = pd.DataFrame(ldict)
ldf.to_csv(meta_dir + 'snr_stations.csv',index=False)

#%%

wf_snr_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/waveforms/'+examplename+'_snr/'
figs_snr_dir = '/Users/rshimony/Desktop/WillametteValley/'+eq_name+'/figs/'+examplename+'_snr/'

try:
    os.mkdir(wf_snr_dir)
except FileExistsError:
    print('Directory already exists')
try:
    os.mkdir(figs_snr_dir)
except FileExistsError:
    print('Directory already exists')


for i in range(86,len(snr_name_unq)):
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
# Attach REC to exsisting SW4 infile

with open(infile, 'a') as fout:
    
    fout.write('\n')
    fout.write('**REC stations**\n')

    for i in range(len(snr_name_unq)):
        name = snr_name_unq[i]
        lat = snr_lats[i]
        lon = snr_lons[i]
    
        fout.write('rec lat=%e lon=%e depth=0 file=%s variables=velocity\n' % (lat,lon,name))
    fout.close()

























