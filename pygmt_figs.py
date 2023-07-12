#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:29:21 2023

@author: rshimony
"""

import obspy as obs
import pandas as pd
import pygmt
import numpy as np
import glob
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point,LineString
from shapely.ops import split
from obspy.core.utcdatetime import UTCDateTime
import matplotlib.pyplot as plt
import geopandas as gpd

#%%

st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/domain_stations.csv')
st_names = np.array(st_file['st_name_dom'])
st_lons = np.array(st_file['st_lon_dom'])
st_lats = np.array(st_file['st_lat_dom']) 

#%%

obs_dir = '/Users/rshimony/Desktop/WillametteValley/salem_eq/observed_data_salem_1hz/'
synts_dir = '/Users/rshimony/Desktop/WillametteValley/Models/gp2_singmat/no_prefilt_wfs/'
steph_dir = '/Users/rshimony/Desktop/WillametteValley/Models/rfile_stephenson_salem_topo/'

obs_trs = sorted(glob.glob(obs_dir + '*_filt.sac'))
synts_trs = sorted(glob.glob(synts_dir + '*.*v'))
steph_trs = sorted(glob.glob(steph_dir + '*.*v'))

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
synts_strs = []
steph_strs = []

for i in range(len(st_names)):
    if st_names[i] in obs_nms_unq:
        obs_str = obs.read(obs_dir + '*' + st_names[i] + '*_filt.sac')
        obs_strs.append(obs_str)
    
    if st_names[i] in synt_nms_unq:
        synts_str = obs.read(synts_dir + st_names[i] + '.*v')
        synts_strs.append(synts_str)
        
    if st_names[i] in steph_nms_unq:
        steph_str = obs.read(steph_dir + st_names[i] + '.*v')
        steph_strs.append(steph_str)

#%%

obs_strs_3 = []
synts_strs_3 = []
steph_strs_3 = []

for i in obs_strs:
    if len(i) >= 3:
       obs_strs_3.append(i[-3:]) 
        
for i in synts_strs:
    if len(i) >= 3:
       synts_strs_3.append(i[-3:]) 
       
for i in steph_strs:
    if len(i) >= 3:
       steph_strs_3.append(i[-3:]) 

    

#%% Functions and data

def get_pgv(stream):
    
    '''
    Gives the PGV from a an obspy stream
    '''
    
    import numpy as np
    import obspy
    
    pgv_list = []
    
    for i in range(len(stream)):
        stream_copy = stream.copy()
        data = stream_copy[i].data
        pgv = np.max(np.abs(data))
        pgv_list.append(pgv)
        
    return(pgv_list) 

def find_vec_norm(Im_value_list):
    '''

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

xpn, ypn = full_poly.exterior.xy
xph1, yph1 = full_poly.interiors[0].xy
xph2, yph2 = full_poly.interiors[1].xy

east_patch_max_lat = 45.0
east_patch_min_lat = 44.6
east_patch_max_lon = -122.1
east_patch_min_lon = -122.65
east_patch_dom = Polygon([[east_patch_min_lon,east_patch_min_lat] , [east_patch_max_lon,east_patch_min_lat] , 
                          [east_patch_max_lon,east_patch_max_lat] , [east_patch_min_lon,east_patch_max_lat]])

xpn_c = []
ypn_c = []

for i in range(len(xpn)):
    grd_pnt = Point(xpn[i],ypn[i])
    if east_patch_dom.contains(grd_pnt) == False:
        xpn_c.append(xpn[i])
        ypn_c.append(ypn[i])

xpnn = np.array([xpn])
xph1n = np.array([xph1])
xph2n = np.array([xph2])

ypnn = np.array([ypn])
yph1n = np.array([yph1])
yph2n = np.array([yph2])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

wv_poly_dict = {'wv_poly_lons':xpn_c,'wv_poly_lats':ypn_c}
wv_poly_df = pd.DataFrame(wv_poly_dict)
wv_poly_df.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/wv_poly.csv',index=False)

### define a function that calculates and outputs the SHIFTED amplitude spectrum

def amp_spec(y,dt):
    '''
    y is the input time series
    dt is the sampling interval
    
    return as outoput the shifted freq vector, and the shifted ampl. spectrum
    '''
    
    y_fft = np.fft.fft(y)
    ampl = np.abs(y_fft)  #this is the ampl. spectrum proper
    f = np.fft.fftfreq(len(y), dt) #this is the un-shifted frequency vector
    
    #now shift because you have OCD
    f = np.fft.fftshift(f)
    ampl = np.fft.fftshift(ampl)*2 #double if you are on;ly gping to plot the positive
    
    return f,ampl

#%%

st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/domain_stations.csv')
st_names = np.array(st_file['st_name_dom'])
st_lons = np.array(st_file['st_lon_dom'])
st_lats = np.array(st_file['st_lat_dom']) 

buffered_fullpoly = full_poly.buffer(0.1, join_style=2)

valst_lons = []
valst_lats = []
valst_names = []

for i in range(len(st_lons)):
    grd_point = Point(st_lons[i],st_lats[i])
    if buffered_fullpoly.contains(grd_point) == True:
        valst_lons.append(st_lons[i])
        valst_lats.append(st_lats[i])
        valst_names.append(st_names[i])

#%%

gdf = gpd.read_file('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/city.zip')

cities = gpd.read_file(gpd.datasets.get_path('naturalearth_cities'))

#%%

depth2units_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/depth2units_wells/depth2units_wells.csv')

lon_units = depth2units_f['Longitude']
lat_units = depth2units_f['Latitude']

#%%

# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

# make color pallets
# pygmt.makecpt(
#     cmap='geo',
#     series='-2000/4000/100',
#     continuous=True)

## topo:
topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=18):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=['WSne', "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(topo_data, cmap='elevation',shading=True,frame=True)
    fig.coast(water="skyblue3")
    
    fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    fig.plot(x=xph1, y=yph1,pen="1p,black")
    fig.plot(x=xph2, y=yph2,pen="1p,black")
    
    fig.plot(x=valst_lons, y=valst_lats, 
              style="t0.42c", 
              fill="navy",
              pen="black",
              transparency=30,
              label = 'Stations') 

    fig.plot(x=lon_units, y=lat_units, 
              style="c0.4c", 
              fill="springgreen3",
              pen="black",
              transparency=30,
              label = 'Wells')

    # fig.text(text=st_names[20], x=st_lons[20]-0.02, y=st_lats[20]+0.08 , font='10p,Helvetica-Bold,magenta4')
    # fig.text(text=st_names[64], x=st_lons[64]-0.06, y=st_lats[64]-0.06, font='10p,Helvetica-Bold,magenta4')
    # fig.text(text=st_names[89], x=st_lons[89]+0.05, y=st_lats[89]+0.06, font='10p,Helvetica-Bold,magenta4')
    fig.legend(position="JTR+jTR+o0.2c", box='+gWhite+p1p',)

fig.show()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wv_map_basin_st_wells_valst.png',dpi=600)


#%%

# trobs = obs_strs_3[89][0]
# trsyn = synts_strs_3[89][0]
# name = st_names[89]
# t_start = trobs.stats.starttime
# tobs = trobs.times(reftime=UTCDateTime(t_start))
# npts = trsyn.stats.npts
# samprate = trsyn.stats.sampling_rate
# tsyn = np.arange(0, npts / samprate, 1 / samprate)

# fig = pygmt.Figure()
# fig.plot(
#     region=[0, 80, -3.973167e-05, 3.973167e-05],
#     projection="X25c/20c",
#     frame="a",
#     x=tobs,
#     y=trobs.data,
#     pen="0.8p,black",
# )
# fig.show()


#%%

'basin_stations_obs = 20-COBRA , 64-NOMA , 89-TAUL3'


fig = plt.figure(figsize=(50, 20))
fig.tight_layout(rect=[0, 0.03, 1, 0.985],h_pad=6.0)
fig.suptitle('GP2 Uniform Filling', fontsize=60 , x=0.5)

ax1 = fig.add_subplot(1,3,1)
# ax1 = plt.subplot(131)
Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/gp2_depth_st.png')
ax1.imshow(Image1)
ax1.axis('off')

ax2 = fig.add_subplot(3,3,2)
# ax2 = plt.subplot(332)
trobs = obs_strs_3[89][0]
trsyn = synts_strs_3[89][0]
name = st_names[89]
t_start = trobs.stats.starttime
tobs = trobs.times(reftime=UTCDateTime(t_start))
npts = trsyn.stats.npts
samprate = trsyn.stats.sampling_rate
tsyn = np.arange(0, npts / samprate, 1 / samprate)

ax2.plot(tobs,trobs.data, c = 'darkslategrey', label = 'Observed',zorder=0)
ax2.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
ax2.set_xlim([-5,80])
ax2.set_title(name + ' ' + trobs.stats.channel ,size=25)

limobs = max(abs(trobs.data))
limsyn = max(abs(trsyn.data))
lim = max(limobs, limsyn)

ax2.set_ylim([-1*lim, lim])
ax2.set_ylim([-1*lim, lim])
ax2.tick_params(labelsize=15)
ax2.locator_params(tight=True, nbins=4)
ax2.set_ylabel('vel (m/s)',fontsize=25)
# ax2.set_xlabel('time (s)',fontsize=25)
ax2.tick_params(labelsize=25)
ax2.legend(loc = 'upper right',prop={'size': 20.0})
        

ax3 = fig.add_subplot(3,3,5)
# ax3 = plt.subplot(335)
trobs = obs_strs_3[64][0]
trsyn = synts_strs_3[64][0]
name = st_names[64]
t_start = trobs.stats.starttime
tobs = trobs.times(reftime=UTCDateTime(t_start))
npts = trsyn.stats.npts
samprate = trsyn.stats.sampling_rate
tsyn = np.arange(0, npts / samprate, 1 / samprate)

ax3.plot(tobs,trobs.data, c = 'darkslategrey', label = 'Observed',zorder=0)
ax3.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
ax3.set_xlim([-5,80])
ax3.set_title(name + ' ' + trobs.stats.channel ,size=25)

limobs = max(abs(trobs.data))
limsyn = max(abs(trsyn.data))
lim = max(limobs, limsyn)

ax3.set_ylim([-1*lim, lim])
ax3.set_ylim([-1*lim, lim])
ax3.tick_params(labelsize=15)
ax3.locator_params(tight=True, nbins=4)
ax3.set_ylabel('vel (m/s)',fontsize=25)
# ax3.set_xlabel('time (s)',fontsize=25)
ax3.tick_params(labelsize=25)
ax3.legend(loc = 'upper right',prop={'size': 20.0})


ax4 = fig.add_subplot(3,3,8)
# ax4 = plt.subplot(338)
trobs = obs_strs_3[20][0]
trsyn = synts_strs_3[20][0]
name = st_names[20]
t_start = trobs.stats.starttime
tobs = trobs.times(reftime=UTCDateTime(t_start))
npts = trsyn.stats.npts
samprate = trsyn.stats.sampling_rate
tsyn = np.arange(0, npts / samprate, 1 / samprate)

ax4.plot(tobs,trobs.data, c = 'darkslategrey', label = 'Observed',zorder=0)
ax4.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
ax4.set_xlim([-5,80])
ax4.set_title(name + ' ' + trobs.stats.channel ,size=25)

limobs = max(abs(trobs.data))
limsyn = max(abs(trsyn.data))
lim = max(limobs, limsyn)

ax4.set_ylim([-1*lim, lim])
ax4.set_ylim([-1*lim, lim])
ax4.tick_params(labelsize=15)
ax4.locator_params(tight=True, nbins=4)
ax4.set_ylabel('vel (m/s)',fontsize=25)
ax4.set_xlabel('time (s)',fontsize=25)
ax4.tick_params(labelsize=25)
ax4.legend(loc = 'upper right',prop={'size': 20.0})


ax5 = fig.add_subplot(1,3,3)
Image5 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/dist2edge_depth_pgv_res.png')
ax5.imshow(Image5)
ax5.axis('off')

# ax5 = fig.add_subplot(3,3,3)
# # ax5 = plt.subplot(333)
# fft_tr_obs = obs_strs_3[89][0]
# fft_tr_data_obs = fft_tr_obs.data
# fft_tr_samprate_obs = fft_tr_obs.stats.delta

# freqs_obs , amps_obs = amp_spec(fft_tr_data_obs,fft_tr_samprate_obs)

# fft_tr_synts = synts_strs_3[89][0]
# fft_tr_data_synts = fft_tr_synts.data
# fft_tr_samprate_synts = fft_tr_synts.stats.delta

# freqs_synts , amps_synts = amp_spec(fft_tr_data_synts,fft_tr_samprate_synts)

# ax5.loglog(freqs_obs,amps_obs, c = 'darkslategrey')
# ax5.loglog(freqs_synts,amps_synts,  c = 'maroon')
# # ax5.set_xlabel('Freqs [Hz]',fontsize=25)
# ax5.set_ylabel('Amplitude',fontsize=25)
# ax5.tick_params(labelsize=25)


# ax6 = fig.add_subplot(3,3,6)
# # ax6 = plt.subplot(336)
# fft_tr_obs = obs_strs_3[64][0]
# fft_tr_data_obs = fft_tr_obs.data
# fft_tr_samprate_obs = fft_tr_obs.stats.delta

# freqs_obs , amps_obs = amp_spec(fft_tr_data_obs,fft_tr_samprate_obs)

# fft_tr_synts = synts_strs_3[64][0]
# fft_tr_data_synts = fft_tr_synts.data
# fft_tr_samprate_synts = fft_tr_synts.stats.delta

# freqs_synts , amps_synts = amp_spec(fft_tr_data_synts,fft_tr_samprate_synts)

# ax6.loglog(freqs_obs,amps_obs, c = 'darkslategrey')
# ax6.loglog(freqs_synts,amps_synts,  c = 'maroon')
# # ax6.set_xlabel('Freqs',fontsize=25)
# ax6.set_ylabel('Amps',fontsize=25)
# ax6.tick_params(labelsize=25)


# ax7 = fig.add_subplot(3,3,9)
# # ax7 = plt.subplot(339)
# fft_tr_obs = obs_strs_3[20][0]
# fft_tr_data_obs = fft_tr_obs.data
# fft_tr_samprate_obs = fft_tr_obs.stats.delta

# freqs_obs , amps_obs = amp_spec(fft_tr_data_obs,fft_tr_samprate_obs)

# fft_tr_synts = synts_strs_3[20][0]
# fft_tr_data_synts = fft_tr_synts.data
# fft_tr_samprate_synts = fft_tr_synts.stats.delta

# freqs_synts , amps_synts = amp_spec(fft_tr_data_synts,fft_tr_samprate_synts)

# ax7.loglog(freqs_obs,amps_obs, c = 'darkslategrey')
# ax7.loglog(freqs_synts,amps_synts,  c = 'maroon')
# ax7.set_xlabel('Freqs',fontsize=25)
# ax7.set_ylabel('Amps',fontsize=25)
# ax7.tick_params(labelsize=25)

# plt.show()
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wf_comp_pygmt_test'
                    ,dpi=100,bbox_inches='tight')
#%%
'basin_stations_obs = 20-COBRA , 64-NOMA , 89-TAUL3'

fig = plt.figure(figsize=(40, 16))
# fig.tight_layout(rect=[0, 0.03, 1, 0.985],h_pad=6.0,w_pad=6.0)
fig.subplots_adjust(wspace=0.05)
fig.suptitle('GP2 Uniform Filling', fontsize=60 , x=0.38,)

ax1 = fig.add_subplot(1,4,1)
# ax1 = plt.subplot(131)
Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/grav_top_depth_st_test.png')
ax1.imshow(Image1)
ax1.axis('off')

ax2 = fig.add_subplot(3,4,3)
# ax2 = plt.subplot(332)
trobs = obs_strs_3[89][0]
trsyn = synts_strs_3[89][0]
trsteph = steph_strs_3[89][0]
name = st_names[89]
t_start = trobs.stats.starttime
tobs = trobs.times(reftime=UTCDateTime(t_start))
npts_syn = trsyn.stats.npts
samprate_syn = trsyn.stats.sampling_rate
tsyn = np.arange(0, npts_syn / samprate_syn, 1 / samprate_syn)
npts_steph = trsteph.stats.npts
samprate_steph = trsteph.stats.sampling_rate
tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)

ax2.plot(tobs,trobs.data, c = 'darkslategrey', label = 'Observed',zorder=0)
ax2.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
ax2.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',zorder=2)
ax2.set_xlim([-5,80])
ax2.set_title(name + ' ' + trobs.stats.channel ,size=25)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax2.yaxis.offsetText.set_fontsize(20)

limobs = max(abs(trobs.data))
limsyn = max(abs(trsyn.data))
lim = max(limobs, limsyn)

ax2.set_ylim([-1*lim, lim])
ax2.set_ylim([-1*lim, lim])
ax2.tick_params(labelsize=15)
ax2.locator_params(tight=True, nbins=4)
ax2.set_ylabel('vel (m/s)',fontsize=28)
# ax2.set_xlabel('time (s)',fontsize=25)
ax2.tick_params(labelsize=28)
ax2.legend(loc = 'upper right',prop={'size': 20.0})
        

ax3 = fig.add_subplot(3,4,7)
# ax3 = plt.subplot(335)
trobs = obs_strs_3[64][0]
trsyn = synts_strs_3[64][0]
trsteph = steph_strs_3[64][0]
name = st_names[64]
t_start = trobs.stats.starttime
tobs = trobs.times(reftime=UTCDateTime(t_start))
npts_syn = trsyn.stats.npts
samprate_syn = trsyn.stats.sampling_rate
tsyn = np.arange(0, npts_syn / samprate_syn, 1 / samprate_syn)
npts_steph = trsteph.stats.npts
samprate_steph = trsteph.stats.sampling_rate
tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)

ax3.plot(tobs,trobs.data, c = 'darkslategrey', label = 'Observed',zorder=0)
ax3.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
ax3.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',zorder=2)
ax3.set_xlim([-5,80])
ax3.set_title(name + ' ' + trobs.stats.channel ,size=25)
ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax3.yaxis.offsetText.set_fontsize(20)

limobs = max(abs(trobs.data))
limsyn = max(abs(trsyn.data))
lim = max(limobs, limsyn)

ax3.set_ylim([-1*lim, lim])
ax3.set_ylim([-1*lim, lim])
ax3.tick_params(labelsize=15)
ax3.locator_params(tight=True, nbins=4)
ax3.set_ylabel('vel (m/s)',fontsize=28)
# ax3.set_xlabel('time (s)',fontsize=25)
ax3.tick_params(labelsize=28)
ax3.legend(loc = 'upper right',prop={'size': 20.0})


ax4 = fig.add_subplot(3,4,11)
# ax4 = plt.subplot(338)
trobs = obs_strs_3[20][0]
trsyn = synts_strs_3[20][0]
trsteph = steph_strs_3[20][0]
name = st_names[20]
t_start = trobs.stats.starttime
tobs = trobs.times(reftime=UTCDateTime(t_start))
npts_syn = trsyn.stats.npts
samprate_syn = trsyn.stats.sampling_rate
tsyn = np.arange(0, npts_syn / samprate_syn, 1 / samprate_syn)
npts_steph = trsteph.stats.npts
samprate_steph = trsteph.stats.sampling_rate
tsteph = np.arange(0, npts_steph / samprate_steph, 1 / samprate_steph)

ax4.plot(tobs,trobs.data, c = 'darkslategrey', label = 'Observed',zorder=0)
ax4.plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
ax4.plot(tsteph,trsteph,  c = 'green', label  = 'Stephenson',zorder=2)
ax4.set_xlim([-5,80])
ax4.set_title(name + ' ' + trobs.stats.channel ,size=25)
ax4.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax4.yaxis.offsetText.set_fontsize(20)

limobs = max(abs(trobs.data))
limsyn = max(abs(trsyn.data))
lim = max(limobs, limsyn)

ax4.set_ylim([-1*lim, lim])
ax4.set_ylim([-1*lim, lim])
ax4.tick_params(labelsize=15)
ax4.locator_params(tight=True, nbins=4)
ax4.set_ylabel('vel (m/s)',fontsize=28)
ax4.set_xlabel('time (s)',fontsize=28)
ax4.tick_params(labelsize=28)
ax4.legend(loc = 'upper right',prop={'size': 20.0})


ax5 = fig.add_subplot(3,4,10)
Image5 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/cs_gp2_sm_cobra.jpg')
ax5.imshow(Image5)
ax5.axis('off')

ax6 = fig.add_subplot(3,4,6)
Image6 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/cs_gp2_sm_noma.jpg')
ax6.imshow(Image6)
ax6.axis('off')

ax7 = fig.add_subplot(3,4,2)
Image7 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/cs_gp1_1d_taul3_test.jpg')
ax7.imshow(Image7)
ax7.axis('off')

# ax8 = fig.add_subplot(2,4,4)
# # ax1 = plt.subplot(131)
# Image8 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/dist2edge_depth_ai_res.png')
# ax8.imshow(Image8)
# ax8.set_position([0.59, 0.52, 0.35, 0.35])
# ax8.axis('off')

# ax9 = fig.add_subplot(2,4,8)
# # ax1 = plt.subplot(131)
# Image9 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/dist2edge_depth_pgv_res.png')
# ax9.imshow(Image9)
# ax9.set_position([0.59, 0.12, 0.35, 0.35])
# ax9.axis('off')

# plt.show()
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/resultsfig_pygmt_test'
                    ,dpi=600,bbox_inches='tight')
#%%

# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

# make color pallets
pygmt.makecpt(
    cmap='geo',
    series='-2000/4000/100',
    continuous=True)

## topo:
topo_data = '@earth_relief_15s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
    fig.coast(shorelines=True, frame="ag", resolution="i")
    
    fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    fig.plot(x=xph1, y=yph1,pen="1p,black")
    fig.plot(x=xph2, y=yph2,pen="1p,black")
    
    fig.coast(water="whitesmoke")

fig.show()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wv_map_shade.png',dpi=600)

#%%

def extract_ifile_basin_data(ifile_path):
    ifile_dists_file = np.genfromtxt(ifile_path)
    
    ifile_dists_lons = []
    ifile_dists_lats = []
    ifile_dists = []
    
    for i in range(1,len(ifile_dists_file)):
    
        ifile_dist = ifile_dists_file[i][2]
        ifile_dist_lons = ifile_dists_file[i][0]
        ifile_dist_lats = ifile_dists_file[i][1]
        
        ifile_dists.append(ifile_dist)
        ifile_dists_lons.append(ifile_dist_lons)
        ifile_dists_lats.append(ifile_dist_lats)
    
    ifile_dists_arr = np.array(ifile_dists)
    ifile_dists_basin = ifile_dists_arr[ifile_dists_arr>0]
    
    ifile_dists_lons_arr = np.array(ifile_dists_lons)
    ifile_dists_lons_basin = ifile_dists_lons_arr[ifile_dists_arr>0]
    
    ifile_dists_lats_arr = np.array(ifile_dists_lats)
    ifile_dists_lats_basin = ifile_dists_lats_arr[ifile_dists_arr>0]
    
    return ifile_dists_lons_basin , ifile_dists_lats_basin , ifile_dists_basin

ifile_grav_top_shallow = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist_grav_top.txt'
ifile_grav_top_eugspe = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/ifile_gp1_eug_spe.txt'
ifile_grav_top = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_gp1_deep_singmat.txt'
ifile_grav_bot = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist_grav_bot.txt'
ifile_gp2 = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_gp2_singmat.txt'
ifile_dist2edge = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist.txt'
ifile_dist2edge_eugspe = '/Users/rshimony/Desktop/WillametteValley/Models/ifiles/ifile_d2e_eug_spe.txt'

# grav_top_lon , grav_top_lat , grav_top_depth = extract_ifile_basin_data(ifile_grav_top)
# grav_top_lon_shallow , grav_top_lat_shallow , grav_top_depth_shallow = extract_ifile_basin_data(ifile_grav_top_shallow)
grav_top_lon_eugspe , grav_top_lat_eugspe , grav_top_depth_eugspe = extract_ifile_basin_data(ifile_grav_top_eugspe)
# grav_bot_lon , grav_bot_lat , grav_bot_depth = extract_ifile_basin_data(ifile_grav_bot)
# gp2_lon , gp2_lat , gp2_depth = extract_ifile_basin_data(ifile_gp2)
# dist2edge_lon , dist2edge_lat , dist2edge_depth = extract_ifile_basin_data(ifile_dist2edge)
# d2e_eugspe_lon , d2e_eugspe_lat , d2e_eugspe_depth = extract_ifile_basin_data(ifile_dist2edge_eugspe)

#%%

depth2units_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/depth2units_wells/depth2units_wells.csv')

lon_units = depth2units_f['Longitude']
lat_units = depth2units_f['Latitude']
srv_depth = depth2units_f['srv_depth']
yamhill_depth = depth2units_f['yamhill_depth']
spencer_depth = depth2units_f['Spencer_depth']
eugene_depth = depth2units_f['eugene_depth']
thickness_srv2yam = depth2units_f['thickness_srv2yam']
thickness_yam2spen = depth2units_f['thickness_yam2spen']
thickness_spen2eug = depth2units_f['thickness_spen2eug']
thickness_yam2eug = depth2units_f['thickness_yam2eug']

unit_depths = [srv_depth , yamhill_depth , spencer_depth , eugene_depth]
unit_thicks = [thickness_srv2yam , thickness_yam2spen , thickness_spen2eug , thickness_yam2eug]

#%%

nonan_depth = eugene_depth[eugene_depth.notna()]
nonan_lons = lon_units[eugene_depth.notna()]
nonan_lats = lat_units[eugene_depth.notna()]

# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

# make color pallets
pygmt.makecpt(
    cmap='geo',
    series='-2000/4000/100',
    continuous=True)

## topo:
topo_data = '@earth_relief_03s' #3 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(topo_data, cmap="geo",shading=True,frame=True)
    fig.coast(shorelines=True, frame="ag", resolution="i")
    
    fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    fig.plot(x=xph1, y=yph1,pen="1p,black")
    fig.plot(x=xph2, y=yph2,pen="1p,black")
    pygmt.makecpt(
        cmap='abyss',reverse=True,
        series=[0,3000],no_bg=True)
    fig.plot(x=nonan_lons, y=nonan_lats, 
              style="c0.25c", 
              fill=nonan_depth,
              pen="black",
              transparency=10,
              cmap=True,
              no_clip='r')
    fig.colorbar(frame="af+lDepth (m)")

fig.show()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wv_eugene_depth_map.png',dpi=600)


#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_dists_lons,ifile_dists_lats,c=ifile_dists,cmap='gist_earth_r')
plt.colorbar()
plt.show()

#%%

plt.figure(figsize=(8,12))
Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wv_map_shade.png')
plt.imshow(Image2)
plt.axis('off')
plt.scatter(ifile_dists_lons_basin,ifile_dists_lats_basin,c=ifile_dists_basin,cmap='jet',alpha=0.5)
plt.colorbar()
plt.show()




#%%

# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

## topo:
topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(topo_data, cmap="white",frame=True)
    fig.coast(water="white",land='white')
    
    # make color pallets
    pygmt.makecpt(
        cmap='abyss',reverse=True,
        series=[0,5000])
    
    # fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    # fig.plot(x=xph1, y=yph1,pen="1p,black")
    # fig.plot(x=xph2, y=yph2,pen="1p,black")
    fig.plot(x=grav_bot_lon, y=grav_bot_lat,
                  style="c0.01c", 
                  fill=grav_bot_depth,
                  cmap=True,
                  transparency=50)
    
    fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")

fig.show()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/grav_bot_depth.png',dpi=600)

#%%

topo_data = '@earth_relief_30s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)
# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]


# Create the first map

fig = pygmt.Figure()

# make color pallets
pygmt.makecpt(
    cmap='abyss',reverse=True,
    series=[0,5000])

with fig.subplot(nrows=1, ncols=4, figsize=("25c", "5c"), autolabel='+JTL', margins="0c"):
    with fig.set_panel(panel=0,fixedlabel='Dist2Edge'):
        with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
            fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
            # fig.coast(shorelines=True, frame="ag", resolution="i")
            
            # make color pallets
            pygmt.makecpt(
                cmap='abyss',reverse=True,
                series=[0,5000])
            
            fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
            fig.plot(x=xph1, y=yph1,pen="1p,black")
            fig.plot(x=xph2, y=yph2,pen="1p,black")
            fig.plot(x=dist2edge_lon, y=dist2edge_lat,
                          style="c0.01c", 
                          fill=dist2edge_depth,
                          cmap=True,
                          transparency=50)
            
            fig.coast(water="whitesmoke")
            fig.colorbar(frame="af+lDepth (m)")


    with fig.set_panel(panel=1,fixedlabel='GP1'):
        with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
            fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
            # fig.coast(shorelines=True, frame="ag", resolution="i")
            
            # make color pallets
            pygmt.makecpt(
                cmap='abyss',reverse=True,
                series=[0,5000])
            
            fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
            fig.plot(x=xph1, y=yph1,pen="1p,black")
            fig.plot(x=xph2, y=yph2,pen="1p,black")
            fig.plot(x=grav_top_lon, y=grav_top_lat,
                          style="c0.01c", 
                          fill=grav_top_depth,
                          cmap=True,
                          transparency=50)
            
            fig.coast(water="whitesmoke")
            fig.colorbar(frame="af+lDepth (m)")

         
    with fig.set_panel(panel=2,fixedlabel='GP2'):
        with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
            fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
            # fig.coast(shorelines=True, frame="ag", resolution="i")
            
            # make color pallets
            pygmt.makecpt(
                cmap='abyss',reverse=True,
                series=[0,5000])
            
            fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
            fig.plot(x=xph1, y=yph1,pen="1p,black")
            fig.plot(x=xph2, y=yph2,pen="1p,black")
            fig.plot(x=gp2_lon, y=gp2_lat,
                          style="c0.01c", 
                          fill=gp2_depth,
                          cmap=True,
                          transparency=50)
            
            fig.coast(water="whitesmoke")
            fig.colorbar(frame="af+lDepth (m)")
            
   
    with fig.set_panel(panel=3,fixedlabel='Eastward Dipping'):
        with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
            fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
            # fig.coast(shorelines=True, frame="ag", resolution="i")
            
            # make color pallets
            pygmt.makecpt(
                cmap='abyss',reverse=True,
                series=[0,5000])
            
            fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
            fig.plot(x=xph1, y=yph1,pen="1p,black")
            fig.plot(x=xph2, y=yph2,pen="1p,black")
            fig.plot(x=grav_bot_lon, y=grav_bot_lat,
                          style="c0.01c", 
                          fill=grav_bot_depth,
                          cmap=True,
                          transparency=50)
            
            fig.coast(water="whitesmoke")
            fig.colorbar(frame="af+lDepth (m)")

fig.show()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wv_models_shade.png',dpi=100)


#%%

fig = plt.figure(figsize=(20, 7))
fig.tight_layout(rect=[0, 0.2, 1, 0.985],w_pad=0.5)
fig.suptitle('WV Basin Models', fontsize=30 , x=0.5)

ax1 = fig.add_subplot(1,4,1)
Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/dist2edge_depth.png')
ax1.imshow(Image1)
ax1.axis('off')
ax1.set_title('Dist2Edge' ,size=16)

ax2 = fig.add_subplot(1,4,2)
Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/grav_top_depth.png')
ax2.imshow(Image2)
ax2.axis('off')
ax2.set_title('Geologic Priors 1' ,size=16)

ax3 = fig.add_subplot(1,4,3)
Image3 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/gp2_depth.png')
ax3.imshow(Image3)
ax3.axis('off')
ax3.set_title('Geologic Priors 2' ,size=16)

ax4 = fig.add_subplot(1,4,4)
Image4 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/grav_bot_depth.png')
ax4.imshow(Image4)
ax4.axis('off')
ax4.set_title('East Dipping Surface' ,size=16)


plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wv_models_gmt.jpg'
                    ,dpi=300,bbox_inches='tight')
# plt.show()

#%%

fig = plt.figure(figsize=(10, 7))
# fig.tight_layout(rect=[0, 0.2, 1, 0.985])
fig.subplots_adjust(wspace=0.005)
# fig.suptitle('Stephenson Model', fontsize=20 , x=0.5,y=0.95)

ax1 = fig.add_subplot(2,3,1)
Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/gp1_1d_depth_pgv_res_ratio.png')
ax1.imshow(Image1)
ax1.axis('off')
ax1.set_title('GP1 1D' ,size=12)

ax2 = fig.add_subplot(2,3,2)
Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/gp2_1d_depth_pgv_res_ratio.png')
ax2.imshow(Image2)
ax2.axis('off')
ax2.set_title('GP2 1D' ,size=12)

ax3 = fig.add_subplot(2,3,3)
Image3 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/d2e_1d_depth_pgv_res_ratio.png')
ax3.imshow(Image3)
ax3.axis('off')
ax3.set_title('D2E 1D' ,size=12)

ax4 = fig.add_subplot(2,3,4)
Image4 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/gp1_sm_depth_pgv_res_ratio.png')
ax4.imshow(Image4)
ax4.axis('off')
ax4.set_title('GP1 Uniform' ,size=12)

ax5 = fig.add_subplot(2,3,5)
Image5 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/gp2_sm_depth_pgv_res_ratio.png')
ax5.imshow(Image5)
ax5.axis('off')
ax5.set_title('GP2 Uniform' ,size=12)

ax6 = fig.add_subplot(2,3,6)
Image6 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/d2e_sm_depth_pgv_res_ratio.png')
ax6.imshow(Image6)
ax6.axis('off')
ax6.set_title('D2E Uniform' ,size=12)

plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/res_ratio_pgv.jpg'
                    ,dpi=300,bbox_inches='tight')

#%%

fig = plt.figure(figsize=(6, 7))
# fig.tight_layout(rect=[0, 0.2, 1, 0.985])
fig.subplots_adjust(wspace=0.005)
fig.suptitle('Stephenson Model', fontsize=20 , x=0.5,y=0.95)

ax1 = fig.add_subplot(2,2,1)
Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/pngs/mag_steph_13.0s.png')
ax1.imshow(Image1)
ax1.axis('off')

ax2 = fig.add_subplot(2,2,2)
Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/pngs/mag_steph_25.0s.png')
ax2.imshow(Image2)
ax2.axis('off')

ax3 = fig.add_subplot(2,2,3)
Image3 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/pngs/mag_steph_48.0s.png')
ax3.imshow(Image3)
ax3.axis('off')

ax4 = fig.add_subplot(2,2,4)
Image4 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/pngs/mag_steph_60.0s.png')
ax4.imshow(Image4)
ax4.axis('off')

plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/waveprop_snaps_steph.jpg'
                    ,dpi=300,bbox_inches='tight')
#%%

ims_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/gp2_singmat/no_prefilt_wfs/ims_ff.csv')

res_st_lons = np.array(ims_file['st_lons'])
res_st_lats = np.array(ims_file['st_lats'])

pgv_t_synt_arr = np.array(ims_file['pgv_t_synt'])
ai_t_synt_arr = np.array(ims_file['ai_t_synt'])
pgv_t_obs_arr = np.array(ims_file['pgv_t_obs'])
ai_t_obs_arr = np.array(ims_file['ai_t_obs'])

xcorr_t_synts = np.array(ims_file['xcorr_t_synt'])

pgv_t_res = np.log(pgv_t_obs_arr) - np.log(pgv_t_synt_arr)
ai_t_res = np.log(ai_t_obs_arr) - np.log(ai_t_synt_arr)


ims_file_steph = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/rfile_stephenson_salem_topo/ims_ff.csv')

pgv_t_steph_arr = np.array(ims_file_steph['pgv_t_steph'])
ai_t_steph_arr = np.array(ims_file_steph['ai_t_steph'])
xcorr_t_steph = np.array(ims_file_steph['xcorr_t_steph'])

pgv_t_res_steph = np.log(pgv_t_obs_arr) - np.log(pgv_t_steph_arr)
ai_t_res_steph = np.log(ai_t_obs_arr) - np.log(ai_t_steph_arr)

pgv_t_res_ratio = np.abs(pgv_t_res_steph)/np.abs(pgv_t_res)
ai_t_res_ratio = np.abs(ai_t_res_steph)/np.abs(ai_t_res)
xcorr_ratio = xcorr_t_steph/xcorr_t_synts


#%%

depth2units_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/depth2units_wells/depth2units_wells.csv')

lon_units = depth2units_f['Longitude']
lat_units = depth2units_f['Latitude']
srv_depth = depth2units_f['srv_depth']
yamhill_depth = depth2units_f['yamhill_depth']
spencer_depth = depth2units_f['Spencer_depth']
eugene_depth = depth2units_f['eugene_depth']
thickness_srv2yam = depth2units_f['thickness_srv2yam']
thickness_yam2spen = depth2units_f['thickness_yam2spen']
thickness_spen2eug = depth2units_f['thickness_spen2eug']
thickness_yam2eug = depth2units_f['thickness_yam2eug']

unit_depths = [srv_depth , yamhill_depth , spencer_depth , eugene_depth]
unit_thicks = [thickness_srv2yam , thickness_yam2spen , thickness_spen2eug , thickness_yam2eug]

#%%

def plot_shade_basin_map(basin_model_lon,basin_model_lat,basin_model_depth,fig_name,
                         plot_station=False , plot_res=False , plot_res_ratio=False , plot_wells=False , plot_all_st=False):
    # # Plot the map using pyGMT
    plot_region = [-124.19,-121.51,43.4,46.1]
    
    # Create the first map
    
    fig = pygmt.Figure()
    
    ## topo:
    topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)
    
    with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
        fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
        fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
        # fig.coast(water='skyblue3')
        
        # make color pallets
        pygmt.makecpt(
            cmap='abyss',reverse=True,
            series=[0,5000])
        
        fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
        fig.plot(x=xph1, y=yph1,pen="1p,black")
        fig.plot(x=xph2, y=yph2,pen="1p,black")
        fig.plot(x=basin_model_lon, y=basin_model_lat,
                      style="c0.01c", 
                      fill=basin_model_depth,
                      cmap=True,
                      transparency=50)
        
        fig.coast(water="whitesmoke")
        fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")
        
        if plot_station == True:
            fig.plot(x=[-124.19, -121.51],y=[st_lats[20]+0.03, st_lats[20]+0.03],pen="1p,red")
            fig.plot(x=[-124.19, -121.51],y=[st_lats[64]-0.09, st_lats[64]-0.09],pen="1p,red")
            fig.plot(x=[-124.19, -121.51],y=[st_lats[89]+0.03, st_lats[89]+0.03],pen="1p,red")
            
            fig.text(text='A', x=-121.6, y=st_lats[89]+0.03 , font='18p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
            fig.text(text='B', x=-121.6, y=st_lats[64]-0.09 , font='18p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
            fig.text(text='C', x=-121.6, y=st_lats[20]+0.03 , font='18p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
            
            fig.plot(x=st_lons[[20,64,89]], y=st_lats[[20,64,89]], 
                      style="t0.45c", 
                      fill="red",
                      pen="black",
                      transparency=30)
            fig.text(text=st_names[20], x=st_lons[20]-0.02, y=st_lats[20]+0.09 , font='15p,Helvetica-Bold,magenta4')
            fig.text(text=st_names[64], x=st_lons[64]-0.28, y=st_lats[64], font='15p,Helvetica-Bold,magenta4')
            fig.text(text=st_names[89], x=st_lons[89]+0.05, y=st_lats[89]+0.08, font='15p,Helvetica-Bold,magenta4')
        
        if plot_all_st == True:
            fig.text(text=st_names, x=np.array(st_lons), y=np.array(st_lats)+0.021 , font='6p,Helvetica-Bold,magenta4')
            fig.plot(x=st_lons, y=st_lats, 
                      style="t0.10c", 
                      fill="red",
                      pen="black",
                      transparency=10)
            fig.plot(x=evlon[4] , y=evlat[4] , style='a0.8c',fill='yellow',pen='black')
        
        if plot_res is not False:
            pygmt.makecpt(
                cmap='bam',
                series = (-2,2))
            fig.plot(x=res_st_lons, y=res_st_lats, 
                      style="t0.65c", 
                      fill=plot_res,
                      pen="black",
                      cmap=True,
                      transparency=0)
            fig.colorbar(position="JMR+o0.5c/0c+w12c/0.8c",frame=["af+lResiduals (Observed - Predicted)"])
            
        if plot_res_ratio is not False:
            pygmt.makecpt(
                cmap='bam',
                series = (0,2))
            fig.plot(x=res_st_lons, y=res_st_lats, 
                      style="t0.55c", 
                      fill=plot_res_ratio,
                      pen="black",
                      cmap=True,
                      transparency=0)
            fig.colorbar(position="JMR+o0.5c/0c+w12c/0.8c",frame=["af+lResiduals Ratio (Stephenson/Model)"])
        
        if plot_wells is not False:
            nonan_depth = plot_wells[plot_wells.notna()]
            nonan_lons = lon_units[plot_wells.notna()]
            nonan_lats = lat_units[plot_wells.notna()]
            pygmt.makecpt(
                cmap='abyss',reverse=True,
                series=[plot_wells.min(),plot_wells.max()],no_bg=True)
            fig.plot(x=nonan_lons, y=nonan_lats, 
                      style="c0.35c", 
                      fill=nonan_depth,
                      pen="black",
                      transparency=10,
                      cmap=True,
                      no_clip='r')
            fig.colorbar(position="JMB+w8c/0.7c",frame=["af+lDepth to Surface (m)"])            
            
    
    fig.show()
    fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/'+str(fig_name)+'.png',dpi=300)
    
# plot_shade_basin_map(grav_top_lon , grav_top_lat , grav_top_depth ,'grav_top_deep_depth')
# plot_shade_basin_map(grav_top_lon_shallow , grav_top_lat_shallow , grav_top_depth_shallow ,'grav_top_shallow_depth')
# plot_shade_basin_map(grav_top_lon_eugspe , grav_top_lat_eugspe , grav_top_depth_eugspe ,'grav_top_eugspe_depth')
# plot_shade_basin_map(grav_bot_lon , grav_bot_lat , grav_bot_depth ,'grav_bot_depth')
# plot_shade_basin_map(gp2_lon , gp2_lat , gp2_depth , 'gp2_deep_depth')
# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'dist2edge_deep_depth')
# plot_shade_basin_map(d2e_eugspe_lon , d2e_eugspe_lat , d2e_eugspe_depth , 'dist2edge_eugspe_depth')

plot_shade_basin_map(grav_top_lon_eugspe , grav_top_lat_eugspe , grav_top_depth_eugspe ,'grav_eugspe_stations' , plot_all_st=True)

# plot_shade_basin_map(grav_top_lon , grav_top_lat , grav_top_depth ,'grav_top_depth_st_test',plot_station=True)
# plot_shade_basin_map(gp2_lon , gp2_lat , gp2_depth , 'gp2_depth_st',plot_station=True)
# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'dist2edge_depth_st',plot_station=True)

# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'dist2edge_depth_pgv_res',plot_res=pgv_t_res)
# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'dist2edge_depth_ai_res',plot_res=ai_t_res)

# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'd2e_1d_depth_pgv_res_ratio',plot_res_ratio=pgv_t_res_ratio)
# plot_shade_basin_map(grav_top_lon , grav_top_lat , grav_top_depth , 'gp1_sm_depth_ai_res_ratio',plot_res_ratio=pgv_t_res_ratio)
# plot_shade_basin_map(gp2_lon , gp2_lat , gp2_depth , 'gp2_sm_depth_ai_res_ratio',plot_res_ratio=pgv_t_res_ratio)

# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'depth_eugene_wells',plot_wells=eugene_depth)
# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'depth_srv_wells',plot_wells=srv_depth)
# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'depth_yamhill_wells',plot_wells=yamhill_depth)
# plot_shade_basin_map(dist2edge_lon , dist2edge_lat , dist2edge_depth , 'depth_spencer_wells',plot_wells=spencer_depth)

#%%
# # Plot the map using pyGMT
plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

## topo:
topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(topo_data, cmap="elevation",shading=True,frame=True)
    fig.coast(water="skyblue1")
    
    fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    fig.plot(x=xph1, y=yph1,pen="1p,black")
    fig.plot(x=xph2, y=yph2,pen="1p,black")
    
    fig.coast(water="whitesmoke")
    fig.colorbar(frame="af+lDepth (m)")
    

    fig.plot(x=st_lons[[20,64,89]], y=st_lats[[20,64,89]], 
              style="t0.35c", 
              fill="red",
              pen="black",
              transparency=30)
    fig.text(text=st_names[20], x=st_lons[20]-0.02, y=st_lats[20]+0.08 , font='10p,Helvetica-Bold,magenta4')
    fig.text(text=st_names[64], x=st_lons[64]-0.06, y=st_lats[64]-0.06, font='10p,Helvetica-Bold,magenta4')
    fig.text(text=st_names[89], x=st_lons[89]+0.05, y=st_lats[89]+0.06, font='10p,Helvetica-Bold,magenta4')
        

    pygmt.makecpt(
        cmap='bam',
        series = (-2,2))
    fig.plot(x=res_st_lons, y=res_st_lats, 
              style="t0.35c", 
              fill=plot_res,
              pen="black",
              cmap=True,
              transparency=0)
    fig.colorbar(position="JMR+o0.5c/0c+w7c/0.5c",frame=["af+lResiduals"])
    
fig.show()
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/'+str(fig_name)+'.png',dpi=600)


#%%

fig,axs = plt.subplots(figsize=(12, 7),nrows=1, ncols=3)
# fig.tight_layout(rect=[0, 0.01, 1, 0.995],w_pad=0.1,pad=0.5)
fig.subplots_adjust(wspace=0.08)
fig.suptitle('WV Depth to Geologic Surfaces', fontsize=22 , x=0.5,y=0.97)

# ax1 = fig.add_subplot(1,3,1)
# Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/dist2edge_depth_srv_wells.png')
# ax1.imshow(Image1)
# ax1.axis('off')
# ax1.set_title('SRV' ,size=16)

Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/depth_srv_wells.png')
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].set_title('SRV' ,size=16)

# ax2 = fig.add_subplot(1,3,2)
# Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/dist2edge_depth_yamhill_wells.png')
# ax2.imshow(Image2)
# ax2.axis('off')
# ax2.set_title('Yamhill' ,size=16)

Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/depth_yamhill_wells.png')
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].set_title('Yamhill' ,size=16)

# ax3 = fig.add_subplot(1,3,3)
# Image3 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/dist2edge_depth_eugene_wells.png')
# ax3.imshow(Image3)
# ax3.axis('off')
# ax3.set_title('Eugene' ,size=16)

Image3 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/depth_eugene_wells.png')
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].set_title('Eugene' ,size=16)

plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/wv_d2surf_sp.jpg'
                    ,dpi=600,bbox_inches='tight')
# plt.show()

#%%

events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/event_catalog_valevents.csv')

evlon = events['longitude']
evlat = events['latitude']
evcode = events['Code']

#%%

# Load the TIGER/Line shapefile data
shapefile = gpd.read_file(
    '/Users/rshimony/Downloads/tl_2019_41_tract.zip'
)

# sf_lons = []
# sf_lats = []
# for i in range(len(shapefile)):
#     sflon = float(shapefile.INTPTLON[i][1:])
#     sflat = float(shapefile.INTPTLAT[i][1:])
#     sf_lons.append(sflon)
#     sf_lats.append(sf_lats)

## topo:
topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=22,FONT_ANNOT_SECONDARY=20,FONT_TITLE=22,FONT_LABEL=23):
# Create a pygmt.Figure object and set the region and projection
    fig = pygmt.Figure()
    fig.basemap(region=[-124.7,-121.51,43.4,46.1], projection='M0/0/10c')
    fig.grdimage(topo_data, cmap="geo",shading=True,frame=True)
    fig.coast(water="skyblue3")
    
    # Plot the shapefile data
    fig.plot(
        data=shapefile,
        style="p0.009p,black",
        pen="black")
    
    fig.plot(x=evlon , y=evlat , style='a0.5c',fill='red',pen='black')
    fig.plot(x=evlon[4] , y=evlat[4] , style='a0.8c',fill='yellow',pen='black')
    
    fig.text(text='Eugene', x=-123.2, y=44.2 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text='Corvalis', x=-123.4, y=44.67 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text='Salem', x=-123.25, y=45.0 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text='Portland', x=-122.65, y=45.7 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text=evcode, x=evlon, y=evlat-0.12 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')

# Show the plot
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/region_map.png',dpi=300)

#%%

# Load the TIGER/Line shapefile data
shapefile = gpd.read_file(
    "https://www2.census.gov/geo/tiger/TIGER2019/TRACT/tl_2019_41_tract.zip"
)

# Load the population data
population = gpd.read_file(
    "https://opendata.arcgis.com/datasets/faf42d3b1a3c40a2a954c75148aefa2f_0.geojson"
)

# Join the shapefile data with the population data
shapefile = shapefile.merge(population[["TRACTCE", "POPULATION"]], on="TRACTCE")

# Create a pygmt.Figure object and set the region and projection
fig = pygmt.Figure()
fig.basemap(region=[-124.19,-121.51,43.4,46.1], projection='M0/0/10c')

# Plot the polygons colored by population
fig.plot(
    data=shapefile.geometry,
    cmap="inferno",
    style="c",
    color=shapefile["POPULATION"],
    pen="black",
)

# Show the plot
fig.show()

#%%

sf = gpd.read_file('/Users/rshimony/Downloads/CensusTracts2020/CensusTracts2020.shp')
sfd = gpd.read_file('/Users/rshimony/Downloads/CensusTracts2020/CensusTracts2020.dbf')

# Create a pygmt.Figure object and set the region and projection
fig = pygmt.Figure()
fig.basemap(region=[-124.19,-121.51,43.4,46.1], projection='M0/0/10c')
# totpop = np.array(sfd.POP20)
# landarea = np.array(sfd.AREALAND20).astype(np.float64)
# sf['pop_den'] = totpop/landarea
# sf.plot(column='pop_den')
# Plot the shapefile data
# linestrings = [geom for geom in sfd.geometry]
# for line in linestrings:
#     x, y = line.exterior.xy
#     fig.plot(x=x, y=y, pen="thin")
    
fig.plot(
    data=sf)

# Show the plot
fig.show()

#%%
totpop = np.array(sfd.POP20)
landarea = np.array(sfd.AREALAND20).astype(np.float64)
sf['pop_den'] = totpop/landarea
sf.plot(column='POP20',cmap='inferno_r')

#%%

import rasterio

tif_file = rasterio.open('/Users/rshimony/Downloads/GHS_POP_P2030_GLOBE_R2022A_54009_100_V1_0_R5_C8/GHS_POP_P2030_GLOBE_R2022A_54009_100_V1_0_R5_C8.tif')
ghs_data = tif_file.read()
#%%
import numpy as np

print("Tiff Boundary", tif_file.bounds)
print("Tiff CRS", tif_file.crs)
print("Data shape", ghs_data.shape)
print("Max value", np.amax(ghs_data))
print("Min value", np.amin(ghs_data))

ghs_data[0][ghs_data[0] < 0.0] = 0.0

#%%

from matplotlib import cm
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap

ourcmap = cm.get_cmap('hot_r', 460)
newcolors = ourcmap(np.linspace(0, 1, 460))
background_colour = np.array([0.9882352941176471, 0.9647058823529412, 0.9607843137254902, 1.0])
newcolors[:1, :] = background_colour
newcmp = ListedColormap(newcolors)


import matplotlib.pyplot as plt
import matplotlib.colors as colors

fig, ax = plt.subplots(facecolor='#FCF6F5FF')
fig.set_size_inches(14, 7)
ax.imshow(ghs_data[0], norm=colors.LogNorm(), cmap=newcmp)
ax.axis('off')
plt.show()

#%%

our_cmap = cm.get_cmap('hot_r', 10)
newcolors = our_cmap(np.linspace(0, 1, 10))
background_colour = np.array([0.9882352941176471, 0.9647058823529412, 0.9607843137254902, 1.0])
newcolors = np.vstack((background_colour, newcolors))
our_cmap = ListedColormap(newcolors)
bounds = [0.0, 1, 5, 10, 20, 50, 100, 200, 1000, 2000, 10000]
norm = colors.BoundaryNorm(bounds, our_cmap.N)

fig, ax = plt.subplots(facecolor='#FCF6F5FF')
fig.set_size_inches(14, 7)
ax.imshow(ghs_data[0], norm=norm, cmap=our_cmap)
ax.axis('off')
plt.show()

#%%
pop_grd_file = '/Users/rshimony/Documents/count.grd'

import xarray as xr

pop_grd = xr.open_dataset(pop_grd_file)

pop_lon = pop_grd

# fig = pygmt.Figure()
# fig.basemap(region=[-124.19,-121.51,43.4,46.1], projection='M0/0/10c')

# fig.plot(data=pop_grd)

# # Show the plot
# fig.show()

#%%

f_path = '/Users/rshimony/Desktop/zzthick4.axyz'

f = np.genfromtxt(f_path)

f_lat = []
f_lon = []
f_z = []

for i in range(len(f)):
    f_lon.append(f[i][0])
    f_lat.append(f[i][1])
    f_z.append(f[i][2])

#%%

plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

## topo:
topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=22,FONT_ANNOT_SECONDARY=20,FONT_TITLE=22,FONT_LABEL=23):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
    # fig.coast(shorelines=True, frame="ag", resolution="i")
    
    # make color pallets
    pygmt.makecpt(
        cmap='abyss',reverse=True,
        series=[0,max(f_z)])

    fig.plot(x=f_lon, y=f_lat,
                  style="c0.01c", 
                  fill=f_z,
                  cmap=True,
                  )
    
    fig.coast(water="whitesmoke")
    fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (km)")

    fig.show()
    fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/grav_data.png',dpi=300)
    
#%%

plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

## topo:
topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
    # fig.basemap(region=plot_region, projection="M15c", frame=None)        
    fig.grdimage(topo_data,region=plot_region, projection="M15c", cmap="grey",shading=True,frame=None)
    # fig.coast(shorelines=True, frame="ag", resolution="i")

    fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    fig.plot(x=xph1, y=yph1,pen="1p,black")
    fig.plot(x=xph2, y=yph2,pen="1p,black")
    
    fig.coast(water="whitesmoke")

    fig.show()
    fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/map4waveprop.png',dpi=600)
#%%
import xarray as xr

plot_region = [-124.19,-121.51,43.4,46.1]

fig = pygmt.Figure()

# image_data = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/gp1_1d/mag_55.0s.csv')

# dataarray = xr.DataArray(data=image_data)

image_path = '/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/image_data.png'

topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=14,FONT_TITLE=14):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"]) 
    fig.grdimage(image_path, cmap="bilbao",shading=True,frame=True)       
    # fig.coast(shorelines=True, frame="ag", resolution="i")
    
    # make color pallets
    # pygmt.makecpt(
    #     cmap='abyss',reverse=True,
    #     series=[0,max(f_z)])
    
    fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    fig.plot(x=xph1, y=yph1,pen="1p,black")
    fig.plot(x=xph2, y=yph2,pen="1p,black")
    
    # fig.plot(x=f_lon, y=f_lat,
    #               style="c0.01c", 
    #               fill=f_z,
    #               cmap=True,
    #               )
    
    # fig.coast(water="whitesmoke")
    # fig.colorbar(frame="af+lDepth (km)")

    fig.show()
    # fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/grav_data.png',dpi=600)

#%%
from glob import glob

image_data_ls = sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/gp1_1d/mag*.csv'))

for i in range(len(image_data_ls)):
    time = image_data_ls[1].split('/')[-1].split('_')[1][:4]
    image_data = pd.read_csv(image_data_ls[i])
    fig, ax = plt.subplots()
    
    ax.imshow(image_data, cmap='Reds')
    ax.axis('off')
    # plt.show()
    plt.close(fig)
    fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/gp1_1d/png/mag_image_'+time+'.png',dpi=100,bbox_inches='tight')

#%%

image_data_ls = sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/gp1_1d/mag*.csv'))

time = image_data_ls[1].split('/')[-1].split('_')[1][:4]
image_data = pd.read_csv(image_data_ls[30])
fig, ax = plt.subplots()

ax.imshow(image_data, cmap='Reds')
# ax.axis('off')
# plt.show()
plt.close(fig)
fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/sw4_image_data/gp1_1d/png/mag_image_'+time+'.png',dpi=300,bbox_inches='tight')
                    
#%%

tiff_f = '/Users/rshimony/Desktop/elev_layer2_hh_200.tif'

# # Plot the map using pyGMT
# plot_region = [-124.19,-121.51,43.4,46.1]

# Create the first map

fig = pygmt.Figure()

## topo:
topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
    # fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(tiff_f, cmap='earth',shading=True,frame=True)
    # fig.coast(water='skyblue3')
    
    # make color pallets
    # pygmt.makecpt(
    #     cmap='abyss',reverse=True,
    #     series=[0,5000])
    
    # fig.plot(x=xpn_c, y=ypn_c,pen="1p,black")
    # fig.plot(x=xph1, y=yph1,pen="1p,black")
    # fig.plot(x=xph2, y=yph2,pen="1p,black")
    # fig.plot(x=basin_model_lon, y=basin_model_lat,
    #               style="c0.01c", 
    #               fill=basin_model_depth,
    #               cmap=True,
    #               transparency=50)
    
    # fig.coast(water="whitesmoke")
    fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")                   
    fig.show()















