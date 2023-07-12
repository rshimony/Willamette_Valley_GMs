#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:30:23 2023

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import obspy as obs
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
# import cartopy.crs as ccrs
from obspy.core.utcdatetime import UTCDateTime

#%%

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
        
valst_lons = np.array(valst_lons)
valst_lats = np.array(valst_lats)
valst_names = np.array(valst_names)

#%%

obs_dir = '/Users/rshimony/Desktop/WillametteValley/salem_eq/observed_data_salem_1hz/'
synts_dir = '/Users/rshimony/Desktop/WillametteValley/Models/d2e_singmat/d2e_singmat_salem/'
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

for i in range(len(valst_names)):
    if valst_names[i] in obs_nms_unq:
        obs_str = obs.read(obs_dir + '*' + valst_names[i] + '*_filt.sac')
        obs_strs.append(obs_str)
    
    if valst_names[i] in synt_nms_unq:
        synts_str = obs.read(synts_dir + valst_names[i] + '.*v')
        synts_strs.append(synts_str)
        
    if valst_names[i] in steph_nms_unq:
        steph_str = obs.read(steph_dir + valst_names[i] + '.*v')
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

#%%

for string in synts_strs_3:
    for trace in string:
        trace.filter('bandpass', freqmin=0.05 , freqmax=1 , corners=4, zerophase=True)
    

#%%

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

xpnn = np.array([xpn])
xph1n = np.array([xph1])
xph2n = np.array([xph2])

ypnn = np.array([ypn])
yph1n = np.array([yph1])
yph2n = np.array([yph2])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

# from cartopy.io.img_tiles import GoogleTiles
# class ShadedReliefESRI(GoogleTiles):
#     # shaded relief
#     def _image_url(self, tile):
#         x, y, z = tile
#         url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
#                'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
#                z=z, y=y, x=x)
#         return url


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

# plt.figure()
# plt.plot(xpnn[0],ypnn[0])
# plt.show()

with open('wv_ol_vs30_400.txt', 'w') as f:
    for i in range(len(ypnn[0])):
        f.write("%s %s\n" % (xpnn[0][i],ypnn[0][i]))
        
ol_file = np.genfromtxt('/Users/rshimony/Documents/My_Notebooks/create_vel_model/wv_ol_vs30_400.txt')
    
#%%
#%%

for i in range(len(synts_strs_3)):
    
    pgv_e_synt = get_pgv(synts_strs_3[i])[0]
    pgv_n_synt = get_pgv(synts_strs_3[i])[1]
    pgv_z_synt = get_pgv(synts_strs_3[i])[2]
    pgv_t_synt = find_vec_norm([pgv_e_synt,pgv_n_synt,pgv_z_synt])

    ai_e_synt = arias_I(synts_strs_3[i])[2][0]
    ai_n_synt = arias_I(synts_strs_3[i])[2][1]
    ai_z_synt = arias_I(synts_strs_3[i])[2][2]
    ai_t_synt = find_vec_norm([ai_e_synt,ai_n_synt,ai_z_synt])

    pgv_e_obs = get_pgv(obs_strs_3[i])[0]
    pgv_n_obs = get_pgv(obs_strs_3[i])[1]
    pgv_z_obs = get_pgv(obs_strs_3[i])[2]
    pgv_t_obs = find_vec_norm([pgv_e_obs,pgv_n_obs,pgv_z_obs])

    ai_e_obs = arias_I(obs_strs_3[i])[2][0]
    ai_n_obs = arias_I(obs_strs_3[i])[2][1]
    ai_z_obs = arias_I(obs_strs_3[i])[2][2]
    ai_t_obs = find_vec_norm([ai_e_obs,ai_n_obs,ai_z_obs])
    
    xcorr_e_synt = calc_xcorr((synts_strs_3[i]),(obs_strs_3[i]))[0]
    xcorr_n_synt = calc_xcorr((synts_strs_3[i]),(obs_strs_3[i]))[1]
    xcorr_z_synt = calc_xcorr((synts_strs_3[i]),(obs_strs_3[i]))[2]
    xcorr_t_synt = find_vec_norm([xcorr_e_synt,xcorr_n_synt,xcorr_z_synt])
    
    dict_ims = {'st_names':valst_names[i] , 'st_lons':valst_lons[i] , 'st_lats':valst_lats[i]
                ,'pgv_e_synt':pgv_e_synt , 'pgv_n_synt':pgv_n_synt , 'pgv_z_synt':pgv_z_synt, 'pgv_t_synt':pgv_t_synt 
                ,'ai_e_synt':ai_e_synt , 'ai_n_synt':ai_n_synt, 'ai_z_synt':ai_z_synt, 'ai_t_synt':ai_t_synt
                ,'pgv_e_obs':pgv_e_obs , 'pgv_n_obs':pgv_n_obs , 'pgv_z_obs':pgv_z_obs, 'pgv_t_obs':pgv_t_obs
                ,'ai_e_obs':ai_e_obs , 'ai_n_obs':ai_n_obs , 'ai_z_obs':ai_z_obs, 'ai_t_obs':ai_t_obs
                ,'xcorr_e_synt':xcorr_e_synt , 'xcorr_n_synt':xcorr_n_synt , 'xcorr_z_synt':xcorr_z_synt, 'xcorr_t_synt':xcorr_t_synt}


    if i == 0:
        df_ims = pd.DataFrame(data=dict_ims,index=[0])
    else:
        df_ims_temp = pd.DataFrame(data=dict_ims , index=[i])
        df_ims = pd.concat([df_ims, df_ims_temp], ignore_index=True)
        
df_ims.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/d2e_singmat/d2e_singmat_salem/ims_ff.csv',index=False)

#%%

for i in range(len(steph_strs_3)):
    
    pgv_e_steph = get_pgv(steph_strs_3[i])[0]
    pgv_n_steph = get_pgv(steph_strs_3[i])[1]
    pgv_z_steph = get_pgv(steph_strs_3[i])[2]
    pgv_t_steph = find_vec_norm([pgv_e_steph,pgv_n_steph,pgv_z_steph])

    ai_e_steph = arias_I(steph_strs_3[i])[2][0]
    ai_n_steph = arias_I(steph_strs_3[i])[2][1]
    ai_z_steph = arias_I(steph_strs_3[i])[2][2]
    ai_t_steph = find_vec_norm([ai_e_steph,ai_n_steph,ai_z_steph])
    
    xcorr_e_steph = calc_xcorr((steph_strs_3[i]),(obs_strs_3[i]))[0]
    xcorr_n_steph = calc_xcorr((steph_strs_3[i]),(obs_strs_3[i]))[1]
    xcorr_z_steph = calc_xcorr((steph_strs_3[i]),(obs_strs_3[i]))[2]
    xcorr_t_steph = find_vec_norm([xcorr_e_steph,xcorr_n_steph,xcorr_z_steph])
    
    dict_ims = {'st_names':valst_names[i] , 'st_lons':valst_lons[i] , 'st_lats':valst_lats[i]
                ,'pgv_e_steph':pgv_e_steph , 'pgv_n_steph':pgv_n_steph , 'pgv_z_steph':pgv_z_steph, 'pgv_t_steph':pgv_t_steph 
                ,'ai_e_steph':ai_e_steph , 'ai_n_steph':ai_n_steph, 'ai_z_steph':ai_z_steph, 'ai_t_steph':ai_t_steph
                ,'xcorr_e_steph':xcorr_e_steph , 'xcorr_n_steph':xcorr_n_steph , 'xcorr_z_steph':xcorr_z_steph, 'xcorr_t_steph':xcorr_t_steph}

    if i == 0:
        df_ims = pd.DataFrame(data=dict_ims,index=[0])
    else:
        df_ims_temp = pd.DataFrame(data=dict_ims , index=[i])
        df_ims = pd.concat([df_ims, df_ims_temp], ignore_index=True)
        
df_ims.to_csv('/Users/rshimony/Desktop/WillametteValley/Models/rfile_stephenson_salem_topo/ims_ff.csv',index=False)
#%%

ims_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/d2e_gp/salem/h50/ims_ff.csv')

pgv_t_synt_arr = np.array(ims_file['pgv_t_synt'])
ai_t_synt_arr = np.array(ims_file['ai_t_synt'])
pgv_t_obs_arr = np.array(ims_file['pgv_t_obs'])
ai_t_obs_arr = np.array(ims_file['ai_t_obs'])

xcorr_t = np.array(ims_file['xcorr_t_synt'])

pgv_t_res = np.log(pgv_t_obs_arr) - np.log(pgv_t_synt_arr)
ai_t_res = np.log(ai_t_obs_arr) - np.log(ai_t_synt_arr)

#%%

plt.figure(figsize=(10, 15))

ax = plt.axes(projection=ShadedReliefESRI().crs)
ax.set_extent([-124.19,-121.51,43.4,46.1])
ax.add_image(ShadedReliefESRI(), 10)

gl_major = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
gl_major.xlabels_top = False
gl_major.ylabels_right = False
gl_major.xlabel_style = {'size': 13}
gl_major.ylabel_style = {'size': 13}

# scat1 = ax.scatter(st_lons,st_lats, s=120 , c=pgv_t_res ,cmap='seismic', transform=ccrs.PlateCarree(),edgecolors='k', marker='^', vmin=-4 , vmax=4)
scat1 = ax.scatter(st_lons,st_lats, s=120 , c=xcorr_t ,cmap='Reds', transform=ccrs.PlateCarree(),edgecolors='k', marker='^', vmin=0 , vmax=1)

cb=plt.colorbar(scat1,fraction=0.035, pad=0.08)
# cb.set_alpha(1)
# cb.set_label('Residuals obs-pred',fontsize=15)
cb.set_label('xcorr',fontsize=15)
# cb.draw_all()

ax.plot(xpn,ypn, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
ax.plot(xph1,yph1, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
ax.plot(xph2,yph2, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)

plt.title('xcorr_t',fontsize=18)
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/d2e_gp/salem/h50/xcorr_t.png',
            dpi=200,bbox_inches='tight')
plt.show()

#%%

n = len(st_names)
fig, (axs) = plt.subplots(n,3,figsize = (18,190))
# fig.subplots_adjust(hspace=0.8)
fig.tight_layout(rect=[0, 0.03, 1, 0.976],h_pad=3.0,w_pad=1.8)
fig.suptitle('Salem_D2E_1D', fontsize=25)

for i in range(len(st_names)):
    name = st_names[i]
    print(name)
    for j in range(len(obs_strs_3[i])):
        trobs = obs_strs_3[i][j]
        trsyn = synts_strs_3[i][j]
        t_start = trobs.stats.starttime
        tobs = trobs.times(reftime=UTCDateTime(t_start))
        npts = trsyn.stats.npts
        samprate = trsyn.stats.sampling_rate
        tsyn = np.arange(0, npts / samprate, 1 / samprate)

        axs[i][j].plot(tobs,trobs.data, c = 'darkslategrey', label = 'Observed',zorder=0)
        axs[i][j].plot(tsyn,trsyn,  c = 'maroon', label  = 'Synthetic',zorder=1)
    
        axs[i][j].set_xlim([-5,80])
    
        axs[i][j].set_title(name + ' ' + trobs.stats.channel ,size=12)
    
        limobs = max(abs(trobs.data))
        limsyn = max(abs(trsyn.data))
        lim = max(limobs, limsyn)
    
        axs[i][j].set_ylim([-1*lim, lim])
        axs[i][j].set_ylim([-1*lim, lim])
        axs[i][j].tick_params(labelsize=8)

        axs[i][j].locator_params(tight=True, nbins=4)
    
        # fig.text(0.5, 0.08, 'time (s)', ha='center', fontsize = 20)
        # fig.text(0.04, 0.5, obs_channel_2 + ' (m/s)', va='center', rotation='vertical', fontsize = 20)
        axs[i][0].set_ylabel('vel (m/s)',fontsize=10)
        axs[i][j].set_xlabel('time (s)',fontsize=10)
        axs[i][0].legend(loc = 'upper right',prop={'size': 6.5})
        
        plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/gp2_1d/Salem/h100/compare_sameax_sharey_comp.pdf'
                    ,dpi=100,bbox_inches='tight')

plt.close('all')


#%%
unfilt_trace_path = '/Users/rshimony/Desktop/WillametteValley/salem_eq/observed_data_salem_1hz/UO.COBRA.HNE_unfilt.sac'
uf_tr = obs.read(unfilt_trace_path)[0]
uf_tr_data = uf_tr.data

fft_tr_obs = obs_strs_3[20][0]
# fft_tr_obs.filter('bandpass', freqmin=0.05 , freqmax=1 , corners=4, zerophase=True)
fft_tr_data_obs = fft_tr_obs.data
fft_tr_samprate_obs = fft_tr_obs.stats.delta

fft_tr_obs_filt = obs_strs_3[20][0]
fft_tr_obs_filt.filter('bandpass', freqmin=0.2 , freqmax=1 , corners=4, zerophase=True)
fft_tr_data_obs_filt = fft_tr_obs_filt.data
fft_tr_samprate_obs = fft_tr_obs.stats.delta

freqs_obs , amps_obs = amp_spec(fft_tr_data_obs,fft_tr_samprate_obs)
freqs_obs_uf , amps_obs_uf = amp_spec(uf_tr_data,fft_tr_samprate_obs)
freqs_obs_filt , amps_obs_filt = amp_spec(fft_tr_data_obs_filt,fft_tr_samprate_obs)

fft_tr_synts = synts_strs_3[20][0]
fft_tr_data_synts = fft_tr_synts.data
fft_tr_samprate_synts = fft_tr_synts.stats.delta

freqs_synts , amps_synts = amp_spec(fft_tr_data_synts,fft_tr_samprate_synts)

plt.figure()
plt.title('Filtered (0.2 Hz) vs Filtered (microseism)')
plt.loglog(freqs_obs,amps_obs, label='iltered (microseism)')
plt.loglog(freqs_obs_filt,amps_obs_filt, label='Filtered (0.2 Hz)')
plt.legend()
plt.xlabel('Freqs [Hz]')
plt.ylabel('Amplitude')

plt.show()

#%%

'basin_stations_obs = 20-COBRA , 64-NOMA , 89-TAUL3'


fig = plt.figure(figsize=(80, 20))
fig.tight_layout(rect=[0, 0.03, 1, 0.985],h_pad=6.0)
fig.suptitle('gp2_singmat_no_prefilt', fontsize=60 , x=0.67)

ax1 = fig.add_subplot(1,3,1)
# ax1 = plt.subplot(131)
ax1 = plt.axes(projection=ShadedReliefESRI().crs)
ax1.set_extent([-124.19,-121.51,43.4,46.1])
ax1.add_image(ShadedReliefESRI(), 10)

gl_major = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
gl_major.xlabels_top = False
gl_major.ylabels_right = False
gl_major.xlabel_style = {'size': 24}
gl_major.ylabel_style = {'size': 28}
   
scat1 = ax1.scatter(st_lons[[20,64,89]],st_lats[[20,64,89]], s=300 , c='orange' , transform=ccrs.PlateCarree(),edgecolors='k', marker='^')
ax1.text(st_lons[20]-0.18, st_lats[20]+0.02,str(st_names[20]), fontsize=22, color='k', transform=ccrs.PlateCarree())
ax1.text(st_lons[64]-0.14, st_lats[64]-0.08,str(st_names[64]), fontsize=22, color='k', transform=ccrs.PlateCarree())
ax1.text(st_lons[89]+0.02, st_lats[89],str(st_names[89]), fontsize=22, color='k', transform=ccrs.PlateCarree())

scat2 = ax1.scatter(-122.551,44.54, s=800 , c='yellow' , transform=ccrs.PlateCarree(),edgecolors='k', marker='*')

ax1.plot(xpn,ypn, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
ax1.plot(xph1,yph1, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)
ax1.plot(xph2,yph2, transform=ccrs.PlateCarree() , c="k" , linewidth=1.5 , alpha=1)


ax2 = fig.add_subplot(3,5,4)
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
        

ax3 = fig.add_subplot(3,5,9)
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


ax4 = fig.add_subplot(3,5,14)
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


ax5 = fig.add_subplot(3,5,5)
# ax5 = plt.subplot(333)
fft_tr_obs = obs_strs_3[89][0]
fft_tr_data_obs = fft_tr_obs.data
fft_tr_samprate_obs = fft_tr_obs.stats.delta

freqs_obs , amps_obs = amp_spec(fft_tr_data_obs,fft_tr_samprate_obs)

fft_tr_synts = synts_strs_3[89][0]
fft_tr_data_synts = fft_tr_synts.data
fft_tr_samprate_synts = fft_tr_synts.stats.delta

freqs_synts , amps_synts = amp_spec(fft_tr_data_synts,fft_tr_samprate_synts)

ax5.loglog(freqs_obs,amps_obs, c = 'darkslategrey')
ax5.loglog(freqs_synts,amps_synts,  c = 'maroon')
# ax5.set_xlabel('Freqs [Hz]',fontsize=25)
ax5.set_ylabel('Amplitude',fontsize=25)
ax5.tick_params(labelsize=25)


ax6 = fig.add_subplot(3,5,10)
# ax6 = plt.subplot(336)
fft_tr_obs = obs_strs_3[64][0]
fft_tr_data_obs = fft_tr_obs.data
fft_tr_samprate_obs = fft_tr_obs.stats.delta

freqs_obs , amps_obs = amp_spec(fft_tr_data_obs,fft_tr_samprate_obs)

fft_tr_synts = synts_strs_3[64][0]
fft_tr_data_synts = fft_tr_synts.data
fft_tr_samprate_synts = fft_tr_synts.stats.delta

freqs_synts , amps_synts = amp_spec(fft_tr_data_synts,fft_tr_samprate_synts)

ax6.loglog(freqs_obs,amps_obs, c = 'darkslategrey')
ax6.loglog(freqs_synts,amps_synts,  c = 'maroon')
# ax6.set_xlabel('Freqs',fontsize=25)
ax6.set_ylabel('Amps',fontsize=25)
ax6.tick_params(labelsize=25)


ax7 = fig.add_subplot(3,5,15)
# ax7 = plt.subplot(339)
fft_tr_obs = obs_strs_3[20][0]
fft_tr_data_obs = fft_tr_obs.data
fft_tr_samprate_obs = fft_tr_obs.stats.delta

freqs_obs , amps_obs = amp_spec(fft_tr_data_obs,fft_tr_samprate_obs)

fft_tr_synts = synts_strs_3[20][0]
fft_tr_data_synts = fft_tr_synts.data
fft_tr_samprate_synts = fft_tr_synts.stats.delta

freqs_synts , amps_synts = amp_spec(fft_tr_data_synts,fft_tr_samprate_synts)

ax7.loglog(freqs_obs,amps_obs, c = 'darkslategrey')
ax7.loglog(freqs_synts,amps_synts,  c = 'maroon')
ax7.set_xlabel('Freqs',fontsize=25)
ax7.set_ylabel('Amps',fontsize=25)
ax7.tick_params(labelsize=25)

# plt.show()
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/gp2_singmat/no_prefilt_wfs/wf_ft_comp.jpg'
                    ,dpi=100,bbox_inches='tight')






