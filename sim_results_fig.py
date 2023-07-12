#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 16:22:07 2023

@author: rshimony
"""

import obspy as obs
import pandas as pd
import pygmt
import numpy as np
import glob
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
from obspy.core.utcdatetime import UTCDateTime
import matplotlib.pyplot as plt

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

# grav_top_lon , grav_top_lat , grav_top_depth = extract_ifile_basin_data(ifile_grav_top)
# grav_top_lon_shallow , grav_top_lat_shallow , grav_top_depth_shallow = extract_ifile_basin_data(ifile_grav_top_shallow)
grav_top_lon_eugspe , grav_top_lat_eugspe , grav_top_depth_eugspe = extract_ifile_basin_data(ifile_grav_top_eugspe)
# grav_bot_lon , grav_bot_lat , grav_bot_depth = extract_ifile_basin_data(ifile_grav_bot)
# gp2_lon , gp2_lat , gp2_depth = extract_ifile_basin_data(ifile_gp2)
# dist2edge_lon , dist2edge_lat , dist2edge_depth = extract_ifile_basin_data(ifile_dist2edge)

#%%

south_depth_gp1deep = grav_top_depth[(grav_top_lat < 44.9795)]
south_depth_gp1shallow = grav_top_depth_shallow[(grav_top_lat_shallow < 44.9795)]
south_depth_gp1eugspe = grav_top_depth_eugspe[(grav_top_lat_eugspe < 44.9795)]
south_depth_gp2 = gp2_depth[(gp2_lat < 44.9795)]
south_depth_d2e = dist2edge_depth[(dist2edge_lat < 44.9795)]

#%%

basin_model_lon = grav_top_lon_eugspe
basin_model_lat = grav_top_lat_eugspe
basin_model_depths = grav_top_depth_eugspe
map_3st_filename = 'gp1_eugspe_st'
# map_ai_res_filename = 'gp1_1d_depth_ai_res'
# map_pgv_res_filename = 'd2e_sm_depth_pgv_res'

synts_dir = '/Users/rshimony/Desktop/WillametteValley/Models/eug_spe_singmat_slu/'

# ims_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/d2e_singmat/d2e_singmat_salem/ims_ff.csv')

outputdir = '/Users/rshimony/Desktop/WillametteValley/Models/pygmt_figs/'

results_fig_suptitle = 'GP1 EugSpe SLU MT Uniform Filling'
results_fig_filename = 'gp1_eugspe_slu_sm'

cs_cobra = '/Users/rshimony/Desktop/WillametteValley/Models/cs_eugspe_gravfix_sm_cobra.jpg'
cs_noma = '/Users/rshimony/Desktop/WillametteValley/Models/cs_eugspe_gravfix_sm_noma.jpg'
cs_taul3 = '/Users/rshimony/Desktop/WillametteValley/Models/cs_eugspe_gravfix_sm_taul3.jpg'


###########################################################################################################
###########################################################################################################

st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/domain_stations.csv')
st_names = np.array(st_file['st_name_dom'])
st_lons = np.array(st_file['st_lon_dom'])
st_lats = np.array(st_file['st_lat_dom']) 


obs_dir = '/Users/rshimony/Desktop/WillametteValley/salem_eq/observed_data_salem_1hz/'
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

    

#Functions and data

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


# ims_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/d2e_gp/salem/h50/ims_ff.csv')

# res_st_lons = np.array(ims_file['st_lons'])
# res_st_lats = np.array(ims_file['st_lats'])

# pgv_t_synt_arr = np.array(ims_file['pgv_t_synt'])
# ai_t_synt_arr = np.array(ims_file['ai_t_synt'])
# pgv_t_obs_arr = np.array(ims_file['pgv_t_obs'])
# ai_t_obs_arr = np.array(ims_file['ai_t_obs'])

# xcorr_t = np.array(ims_file['xcorr_t_synt'])

# pgv_t_res = np.log(pgv_t_obs_arr) - np.log(pgv_t_synt_arr)
# ai_t_res = np.log(ai_t_obs_arr) - np.log(ai_t_synt_arr)


# ims_file_steph = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/rfile_stephenson_salem_topo/ims_ff.csv')

# pgv_t_steph_arr = np.array(ims_file_steph['pgv_t_steph'])
# ai_t_steph_arr = np.array(ims_file_steph['ai_t_steph'])

# pgv_t_res_steph = np.log(pgv_t_obs_arr) - np.log(pgv_t_steph_arr)
# ai_t_res_steph = np.log(ai_t_obs_arr) - np.log(ai_t_steph_arr)

# pgv_t_res_ratio = pgv_t_res_steph/pgv_t_res
# ai_t_res_ratio = ai_t_res_steph/ai_t_res


def plot_shade_basin_map(basin_model_lon,basin_model_lat,basin_model_depth,fig_name,plot_station=False , plot_res=False , plot_wells=False):
    # # Plot the map using pyGMT
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

        if plot_wells is not False:
            nonan_depth = plot_wells[plot_wells.notna()]
            nonan_lons = lon_units[plot_wells.notna()]
            nonan_lats = lat_units[plot_wells.notna()]
            pygmt.makecpt(
                cmap='lajolla',
                series=[plot_wells.min(),plot_wells.max()],no_bg=True)
            fig.plot(x=nonan_lons, y=nonan_lats, 
                      style="c0.25c", 
                      fill=nonan_depth,
                      pen="black",
                      transparency=10,
                      cmap=True,
                      no_clip='r')
            fig.colorbar(position="JMR+o0.5c/0c+w7c/0.5c",frame=["af+lDepth to Surface (m)"])            
            
    
    # fig.show()
    fig.savefig(outputdir+str(fig_name)+'.png',dpi=300)


plot_shade_basin_map(basin_model_lon , basin_model_lat , basin_model_depths ,map_3st_filename ,plot_station=True)

# plot_shade_basin_map(basin_model_lon , basin_model_lat , basin_model_depths , map_pgv_res_filename ,plot_res=pgv_t_res)
# plot_shade_basin_map(basin_model_lon , basin_model_lat , basin_model_depths , map_ai_res_filename ,plot_res=ai_t_res)


'basin_stations_obs = 20-COBRA , 64-NOMA , 89-TAUL3'

fig = plt.figure(figsize=(40, 16))
# fig.tight_layout(rect=[0, 0.03, 1, 0.985],h_pad=6.0,w_pad=6.0)
fig.subplots_adjust(wspace=0.05)
fig.suptitle(results_fig_suptitle, fontsize=60 , x=0.38)

ax1 = fig.add_subplot(1,4,1)
# ax1 = plt.subplot(131)
Image1 = plt.imread(outputdir + map_3st_filename + '.png')
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
ax4.locator_params(tight=True, nbins=4)
ax4.set_ylabel('vel (m/s)',fontsize=28)
ax4.set_xlabel('time (s)',fontsize=28)
ax4.tick_params(labelsize=28)
ax4.legend(loc = 'upper right',prop={'size': 20.0})


ax5 = fig.add_subplot(3,4,10)
Image5 = plt.imread(cs_cobra)
ax5.imshow(Image5)
ax5.axis('off')

ax6 = fig.add_subplot(3,4,6)
Image6 = plt.imread(cs_noma)
ax6.imshow(Image6)
ax6.axis('off')

ax7 = fig.add_subplot(3,4,2)
Image7 = plt.imread(cs_taul3)
ax7.imshow(Image7)
ax7.axis('off')

# ax8 = fig.add_subplot(2,4,4)
# # ax1 = plt.subplot(131)
# Image8 = plt.imread(outputdir + map_ai_res_filename + '.png')
# ax8.imshow(Image8)
# ax8.set_title('Arias Intensity' ,size=25)
# ax8.set_position([0.59, 0.52, 0.35, 0.35])
# ax8.axis('off')

# ax9 = fig.add_subplot(2,4,8)
# # ax1 = plt.subplot(131)
# Image9 = plt.imread(outputdir + map_pgv_res_filename + '.png')
# ax9.imshow(Image9)
# ax9.set_title('PGV' ,size=25)
# ax9.set_position([0.59, 0.12, 0.35, 0.35])
# ax9.axis('off')

# plt.show()
plt.savefig(outputdir+'resultsfig_'+results_fig_filename+'.png' ,dpi=300,bbox_inches='tight')

























