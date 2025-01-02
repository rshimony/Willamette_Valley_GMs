#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 16:24:24 2024

@author: rshimony
"""
'''
Stations' location maps.
This script has THREE PARTS:
    Part 1:
        Plotting stations' location and epicenter from A SINGLE EVENT,
        on top a basin depth regional map.
        This script reads in station inventory file from a specific event, 
        and plots a single map for each station.
    Part 2:
        The last part of this script plots a map of stations used in the WV paper.
    Part 3:
        Check which stations recorded all three events.
        Save these stations' data to a separate 'all_eq' flatfile.
        Plot maps of these stations' location with the three epicenters and basin depth plotted as well.
'''

import pandas as pd
import pygmt
import numpy as np
from shapely.geometry.polygon import Polygon

#%% 
'''
*** PART 1 ***
'''
####INPUTS###
## Station inventory - SINGLE EVENT
st_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/metadata/snr_stations.csv').drop_duplicates()
st_names = np.array(st_file['st_nm'])
st_lons = np.array(st_file['st_lon'])
st_lats = np.array(st_file['st_lat']) 

## Reading in events catalog and extracting event location
events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/event_catalog_valevents.csv')
## idx   EQ ##
##  0    Scotts Mills ##
##  1    Salem ##
##  2    Springfield ##
# This index needs to change accordingly
evlon = events['longitude'][0]
evlat = events['latitude'][0]

## ifiles that contain new WV basin models structure
ifile_eugspe = '/Users/rshimony/Desktop/WillametteValley/wv_project/ifiles/ifile_latlon_eugspe.txt'
ifile_srv = '/Users/rshimony/Desktop/WillametteValley/wv_project/ifiles/ifile_latlon_srv.txt'

## Output directory
outdir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/st_loc_maps/springfield/'

#%%
## Creating WV polygon
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

xpn_m, ypn_m = full_poly.exterior.xy 
xph1_m, yph1_m = full_poly.interiors[0].xy
xph2_m, yph2_m = full_poly.interiors[1].xy

def extract_ifile_basin_data(ifile_path):
    '''
    Extract WV basin depth and grid points location 
    from ifile. 
    Outputs depths, lats and lons arrays
    '''
    
    ifile_f = np.genfromtxt(ifile_path,skip_header=1)
    
    ifile_lons = ifile_f[:,0]
    ifile_lats = ifile_f[:,1]
    ifile_depths = ifile_f[:,2]
    
    ifile_depths_basin = ifile_depths[ifile_depths>0]
    
    ifile_lons_basin = ifile_lons[ifile_depths>0]

    ifile_lats_basin = ifile_lats[ifile_depths>0]
    
    return ifile_lons_basin , ifile_lats_basin , ifile_depths_basin

#%%
## Getting basin depth data for both basin depths (EugSpe and SRV)
eugspe_lons , eugspe_lats , eugspe_depths = extract_ifile_basin_data(ifile_eugspe)
srv_lons , srv_lats , srv_depths = extract_ifile_basin_data(ifile_srv)

#%%
# # Plot the map using pyGMT

def make_st_map(idx):

    plot_region = [-124.19,-121.51,43.4,46.1]
    fig = pygmt.Figure()
    topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)
    
    with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
        fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
        fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
        
        # make color pallets
        pygmt.makecpt(
            cmap='abyss',reverse=True,
            series=[0,5000])
        
        # basin depth
        fig.plot(x=eugspe_lons, y=eugspe_lats,
                      style="c0.01c", 
                      fill=eugspe_depths,
                      cmap=True,
                      transparency=50)
        
        fig.coast(water="whitesmoke")
        fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")
        
        # Event location
        fig.plot(x=evlon , y=evlat , style='a0.8c',fill='yellow',pen='black')
        
        #basin outline
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        
        #stations
        fig.plot(x=st_lons[idx] , y=st_lats[idx] , style='t0.6c',fill='red',pen='black',transparency=30)
    
    # Save the plot
    fig.savefig(outdir+ st_names[idx] +'.png',dpi=200)

for i in range(len(st_names)):
    make_st_map(i)

#%% 
'''
**** PART 2 ****
'''
### Plotting map of stations used in the WV paper ###

## Reading in Salem event stations for plotting
st_file_full = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/salem_eq/metadata/snr_stations.csv').drop_duplicates()
st_names_full = np.array(st_file_full['st_nm'])
st_lons_full = np.array(st_file_full['st_lon'])
st_lats_full = np.array(st_file_full['st_lat']) 

## Reading in events catalog and extracting event location
events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/event_catalog_valevents.csv')
## idx   EQ ##
##  0    Scotts Mills ##
##  1    Salem ##
##  2    Springfield ##
# This index needs to change accordingly
evlon = events['longitude']
evlat = events['latitude']

'''
Station used indeces:
    ALVY - 3
    MONO - 56
    WLOO - 51
    SAIL - 83
    COBRA - 20
    NOMA - 63
    EYES - 29
'''

outdir_paper_st = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/st_loc_maps/'

plot_region = [-124.19,-121.51,43.4,46.1]
fig = pygmt.Figure()
topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)

with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
    fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
    fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
    
    # make color pallets
    pygmt.makecpt(
        cmap='abyss',reverse=True,
        series=[0,5000])
    
    # basin depth
    fig.plot(x=eugspe_lons, y=eugspe_lats,
                  style="c0.01c", 
                  fill=eugspe_depths,
                  cmap=True,
                  transparency=50)
    
    fig.coast(water="whitesmoke")
    fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")
    
    # Event location
    fig.plot(x=evlon , y=evlat , style='a0.8c',fill='yellow',pen='black')
    # Event names
    fig.text(text='Scotts Mills', x=evlon[0]+0.38, y=evlat[0]-0.03 , font='12p,Helvetica-Bold,black')
    fig.text(text='Springfield', x=evlon[2]+0.38, y=evlat[2]-0.03, font='12p,Helvetica-Bold,black')
    fig.text(text='Salem', x=evlon[1]+0.28, y=evlat[1]-0.03, font='12p,Helvetica-Bold,black')
    
    #basin outline
    fig.plot(x=xpn_m, y=ypn_m,pen="0.7p,black")
    fig.plot(x=xph1_m, y=yph1_m,pen="0.7p,black")
    fig.plot(x=xph2_m, y=yph2_m,pen="0.7p,black")
    
    # stations location
    fig.plot(x=st_lons_full[[3,56,51,83,20,29,63]] , y=st_lats_full[[3,56,51,83,20,29,63]] , style='t0.6c',fill='red',pen='black',transparency=30)
    # stations names
    fig.text(text=st_names_full[20], x=st_lons_full[20]-0.18, y=st_lats_full[20]+0.09, font='12p,Helvetica-Bold,firebrick4')
    fig.text(text=st_names_full[29], x=st_lons_full[29]-0.25, y=st_lats_full[29], font='12p,Helvetica-Bold,firebrick4')
    fig.text(text=st_names_full[63], x=st_lons_full[63]-0.06, y=st_lats_full[63]-0.07, font='12p,Helvetica-Bold,firebrick4')
    fig.text(text=st_names_full[3], x=st_lons_full[3], y=st_lats_full[3]-0.09, font='12p,Helvetica-Bold,firebrick4')
    fig.text(text=st_names_full[56], x=st_lons_full[56]-0.25, y=st_lats_full[56]+0.08, font='12p,Helvetica-Bold,firebrick4')
    fig.text(text=st_names_full[51], x=st_lons_full[51]-0.18, y=st_lats_full[51]+0.07, font='12p,Helvetica-Bold,firebrick4')
    fig.text(text=st_names_full[83], x=st_lons_full[83]+0.16, y=st_lats_full[83]-0.07, font='12p,Helvetica-Bold,firebrick4')

# Save the plot  
fig.savefig(outdir_paper_st+'paper_stations.png',dpi=400) 
    
#%% 
'''
**** PART 3 ****
'''
## Output directories
# All EQs Stations ff
output_ff_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/'
# Figures
all_eq_outdir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/st_loc_maps/all_eqs/'

## Reading in all events stations diretories
# Salem
st_file_salem = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/salem_eq/metadata/snr_stations.csv').drop_duplicates()
st_names_salem = np.array(st_file_salem['st_nm'])
st_lons_salem = np.array(st_file_salem['st_lon'])
st_lats_salem = np.array(st_file_salem['st_lat']) 
#Scotts Mills
st_file_scottsmills = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/scottsmills_eq/metadata/snr_stations.csv').drop_duplicates()
st_names_scottsmills = np.array(st_file_scottsmills['st_nm'])
st_lons_scottsmills = np.array(st_file_scottsmills['st_lon'])
st_lats_scottsmills = np.array(st_file_scottsmills['st_lat']) 
# Springfield
st_file_springfield = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/springfield_eq/metadata/snr_stations.csv').drop_duplicates()
st_names_springfield = np.array(st_file_springfield['st_nm'])
st_lons_springfield = np.array(st_file_springfield['st_lon'])
st_lats_springfield = np.array(st_file_springfield['st_lat'])     

## Getting stations that recorded all events
## Saving these into separate flatfile and reading it in for further analysis
salem_scotml_mask = np.isin(st_names_salem,st_names_scottsmills)
salem_scotml_st = st_names_salem[salem_scotml_mask]
sal_scotm_spr_mask = np.isin(salem_scotml_st,st_names_springfield)
sal_scotm_spr_st = salem_scotml_st[sal_scotm_spr_mask]

all_eq_mask = np.isin(st_names_salem,sal_scotm_spr_st)
all_eq_st = st_file_salem[all_eq_mask]

df_all_eq_st = pd.DataFrame(data=all_eq_st)
df_all_eq_st.to_csv(output_ff_dir+'all_eq_st.csv')

st_file_f_all_eq = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/all_eq_st.csv')
st_file_all_eq = st_file_f_all_eq.drop_duplicates()
st_names_all_eq = np.array(st_file_all_eq['st_name_dom'])
st_lons_all_eq = np.array(st_file_all_eq['st_lon_dom'])
st_lats_all_eq = np.array(st_file_all_eq['st_lat_dom']) 
    
## Reading in events catalog and extracting event location
events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/event_catalog_valevents.csv')
## idx   EQ ##
##  0    Scotts Mills ##
##  1    Salem ##
##  2    Springfield ##
evlon = events['longitude']
evlat = events['latitude'] 
    
#%% 
## Plotting all events station maps

plot_region = [-124.19,-121.51,43.4,46.1]
fig = pygmt.Figure()
topo_data = '@earth_relief_03s' #30 arc second global relief (SRTM15+V2.1 @ 1.0 km)

for i in range(len(st_names_all_eq)):
    with pygmt.config(FONT_ANNOT_PRIMARY=18,FONT_ANNOT_SECONDARY=18,FONT_TITLE=18,FONT_LABEL=18):
        fig.basemap(region=plot_region, projection="M0/0/10c", frame=["WSne", "x5+lLatitude(°)", "y5+lLongutide(°)"])        
        fig.grdimage(topo_data, cmap="grey",shading=True,frame=True)
        
        # make color pallets
        pygmt.makecpt(
            cmap='abyss',reverse=True,
            series=[0,5000])
        
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        fig.plot(x=eugspe_lons, y=eugspe_lats,
                      style="c0.01c", 
                      fill=eugspe_depths,
                      cmap=True,
                      transparency=50)
        
        fig.coast(water="whitesmoke")
        fig.colorbar(position="JMB+w8c/0.8c",frame="af+lDepth (m)")
        
        # Epicenters
        fig.plot(x=evlon , y=evlat , style='a0.8c',fill='yellow',pen='black')
        fig.text(text='Scotts Mills', x=evlon[0]+0.38, y=evlat[0]-0.03 , font='12p,Helvetica-Bold,firebrick4')
        fig.text(text='Springfield', x=evlon[2]+0.38, y=evlat[2]-0.03, font='12p,Helvetica-Bold,firebrick4')
        fig.text(text='Salem', x=evlon[1]+0.28, y=evlat[1]-0.03, font='12p,Helvetica-Bold,firebrick4')
        
        #basin outline
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        
        #stations
        fig.plot(x=st_lons_all_eq[i] , y=st_lats_all_eq[i] , style='t0.6c',fill='red',pen='black',transparency=30)
    
    # Save the plot
    fig.savefig(all_eq_outdir+ st_names_all_eq[i] +'.png',dpi=200)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    