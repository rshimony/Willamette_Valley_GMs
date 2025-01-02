#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 19:59:19 2024

@author: rshimony
"""
'''
Reading in IMs residual and residual ratio values from flatfiles,
and plotting on top of regional map.
This script takes in residuals and ratios flatfiles of all models (4 new WV models and USGS CVM) from a SINGLE EVENT.
Outputs IMs residuals and ratio maps of each IM in each model for a specific event.
'''

import pandas as pd
import pygmt
import numpy as np
from shapely.geometry.polygon import Polygon
import os

#%%
####INPUTS####
## Color palettes 
cpt_basin = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/Olive-3D_1.cpt' # Basin depth
cpt_stations = 'vik' # Residuals and ratio values

## Reading in events catalog and extracting event location
events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/event_catalog_valevents.csv')
## idx   EQ ##
##  0    Scotts Mills ##
##  1    Salem ##
##  2    Springfield ##
# This index needs to change accordingly
evlon = events['longitude'][2]
evlat = events['latitude'][2]

## ifile that contains the basin structure
ifile_eugspe = '/Users/rshimony/Desktop/WillametteValley/wv_project/ifiles/ifile_latlon_eugspe.txt'
ifile_srv = '/Users/rshimony/Desktop/WillametteValley/wv_project/ifiles/ifile_latlon_srv.txt'

## Reading in residuals and ratios flatfiles
# Full station inventory
eugspe_1d_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/eugspe_1d_res_ff.csv')
eugspe_sm_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/eugspe_sm_res_ff.csv')
srv_1d_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/srv_1d_res_ff.csv')
srv_sm_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/srv_sm_res_ff.csv')

# Valley station inventory
eugspe_1d_val_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/eugspe_1d_res_ff_val.csv')
eugspe_sm_val_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/eugspe_sm_res_ff_val.csv')
srv_1d_val_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/srv_1d_res_ff_val.csv')
srv_sm_val_df = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/srv_sm_res_ff_val.csv')

## Output directories
# Full station inventory
res_map_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/springfield/'
try:
    os.mkdir(res_map_dir)
except FileExistsError:
    print('Directory already exists')
    
# Valley station inventory
ratio_map_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/ratios_maps/springfield/'
try:
    os.mkdir(ratio_map_dir)
except FileExistsError:
    print('Directory already exists')

#%%
###FUNCTIONS####

## Creating WV basin polygon
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


def plot_shade_basin_map(basin_model_lon,basin_model_lat,basin_model_depth,fig_dir,fig_name,res_st_lons,res_st_lats,
                          plot_res=False ,plot_res_ratio=False):
    '''
    Plot residual or ratio values at all stations on top of a regional map showing WV model basin depth.
    '''
    
    if plot_res is not False:
        
        # # Plot the map using pyGMT
        fig = pygmt.Figure()
        ## topo:
        topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)
        
        plot_region = [-123.65,-122.2,43.8,46.1]
        
        with pygmt.config(FONT_ANNOT_PRIMARY=38,FORMAT_GEO_MAP="ddd"):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["Wsne","ya1f0.25+lLongutide(°)"])        
            fig.grdimage(topo_data, cmap="grey",shading=True)
        with pygmt.config(FONT_ANNOT_PRIMARY=38,FORMAT_GEO_MAP="ddd.xx"):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["wSne", "xa0.8f0.2+lLatitude(°)"]) 

    
        # make color pallets
        pygmt.makecpt(
            cmap=cpt_basin,
            series=[0,5000])
        
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        fig.plot(x=basin_model_lon, y=basin_model_lat,
                      style="c0.01c", 
                      fill=basin_model_depth,
                      cmap=True,
                      transparency=20)
        
        fig.coast(water="whitesmoke")
        with pygmt.config(FONT_ANNOT_PRIMARY=32,FONT_ANNOT_SECONDARY=32,FONT_TITLE=32,FONT_LABEL=32):
            fig.colorbar(position="JMB+w8c/1.1c",frame="af+lDepth (m)")
        
        fig.plot(x=evlon , y=evlat , style='a1.5c',fill='yellow',pen='black')
    
        pygmt.makecpt(
            cmap=cpt_stations,reverse=True,
            series = (-2,2))
        fig.plot(x=res_st_lons, y=res_st_lats, 
                  style="t0.85c", 
                  fill=plot_res,
                  pen="black",
                  cmap=True,
                  transparency=15)
        with pygmt.config(FONT_ANNOT_PRIMARY=32,FONT_ANNOT_SECONDARY=32,FONT_TITLE=32,FONT_LABEL=32):
            fig.colorbar(position="JMR+o0.5c/0c+w18c/0.8c",frame=["af+lResiduals (Observed - Predicted)"])
            
        with pygmt.config(FONT_ANNOT_PRIMARY=32,FONT_ANNOT_SECONDARY=32,FONT_TITLE=32,FONT_LABEL=32,MAP_TICK_PEN_PRIMARY=2.2):
            fig.basemap(map_scale="JBR+o-4.4/-1.3+w50")
        with pygmt.config(FONT_ANNOT_PRIMARY=24,FONT_ANNOT_SECONDARY=22,FONT_TITLE=22,FONT_LABEL=22):
            fig.basemap(rose="JTL+o-2.7c+w1.8c+f1+l")
 
        
    if plot_res_ratio is not False:
        
        # # Plot the map using pyGMT
        fig = pygmt.Figure()
        ## topo:
        topo_data = '@earth_relief_03s' #15 arc second global relief (SRTM15+V2.1 @ 1.0 km)
        
        plot_region = [-124.19,-121.51,43.4,46.1]
        
        with pygmt.config(FONT_ANNOT_PRIMARY=28,FORMAT_GEO_MAP="ddd"):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["Wsne","ya1f0.25+lLongutide(°)"])        
            fig.grdimage(topo_data, cmap="grey",shading=True)
        with pygmt.config(FONT_ANNOT_PRIMARY=28,FORMAT_GEO_MAP="ddd"):
            fig.basemap(region=plot_region, projection="M0/0/10c", frame=["wSne", "xa1f0.2+lLatitude(°)"]) 
            
        # make color pallets
        pygmt.makecpt(
            cmap=cpt_basin,
            series=[0,5000])
        
        fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
        fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
        fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
        fig.plot(x=basin_model_lon, y=basin_model_lat,
                      style="c0.01c", 
                      fill=basin_model_depth,
                      cmap=True,
                      transparency=20)
        
        fig.coast(water="whitesmoke")
        with pygmt.config(FONT_ANNOT_PRIMARY=24,FONT_ANNOT_SECONDARY=24,FONT_TITLE=24,FONT_LABEL=24):
            fig.colorbar(position="JMB+w8c/1.1c",frame="af+lDepth (m)")
        
        fig.plot(x=evlon , y=evlat , style='a1.2c',fill='yellow',pen='black')
        
        pygmt.makecpt(
            cmap=cpt_stations,reverse=True,
            series = (0,2))
        fig.plot(x=res_st_lons, y=res_st_lats, 
                  style="t0.65c", 
                  fill=plot_res_ratio,
                  pen="black",
                  cmap=True,
                  transparency=15)
        with pygmt.config(FONT_ANNOT_PRIMARY=24,FONT_ANNOT_SECONDARY=24,FONT_TITLE=24,FONT_LABEL=24):
            fig.colorbar(position="JMR+o0.5c/0c+w12c/0.8c",frame=["af+lResiduals Ratio (USGS CVM/Model)"])  
    
        with pygmt.config(FONT_ANNOT_PRIMARY=24,FONT_ANNOT_SECONDARY=22,FONT_TITLE=22,FONT_LABEL=20,MAP_TICK_PEN_PRIMARY=1):
            fig.basemap(map_scale="JBR+o-3.4/-1.3+w50")
        with pygmt.config(FONT_ANNOT_PRIMARY=24,FONT_ANNOT_SECONDARY=22,FONT_TITLE=22,FONT_LABEL=22):
            fig.basemap(rose="JTL+o-2.1c+w1.4c+f1+l")
    
    # fig.show()
    fig.savefig(fig_dir+fig_name+'.png',dpi=300)

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

def extract_res_rat_vals(res_rat_df,cvm=False):
    '''
    Extract residuals and ratios values from flatfiles.
    Outputs lists of IM residuals and ratios arrays.
    If input is USGS CVM - only residuals (no ratio)
    '''
    
    f2_res = res_rat_df.iloc[:,5].values
    f3_res = res_rat_df.iloc[:,6].values
    f4_res = res_rat_df.iloc[:,7].values
    f5_res = res_rat_df.iloc[:,8].values
    f6_res = res_rat_df.iloc[:,9].values
    f7_res = res_rat_df.iloc[:,10].values
    pgv_res = res_rat_df.iloc[:,11].values
    
    res_im_ls = [f2_res,f3_res,f4_res,f5_res,f6_res,f7_res,pgv_res]
    
    if cvm == True:
        return res_im_ls
    else:
        f2_rat = res_rat_df.iloc[:,12].values
        f3_rat = res_rat_df.iloc[:,13].values
        f4_rat = res_rat_df.iloc[:,14].values
        f5_rat = res_rat_df.iloc[:,15].values
        f6_rat = res_rat_df.iloc[:,16].values
        f7_rat = res_rat_df.iloc[:,17].values
        pgv_rat = res_rat_df.iloc[:,18].values
        
        
        rat_im_ls = [f2_rat,f3_rat,f4_rat,f5_rat,f6_rat,f7_rat,pgv_rat]
        
        return res_im_ls,rat_im_ls

#%%
## Getting basin depth data for both basin depths (EugSpe and SRV)
eugspe_lons , eugspe_lats , eugspe_depths = extract_ifile_basin_data(ifile_eugspe)
srv_lons , srv_lats , srv_depths = extract_ifile_basin_data(ifile_srv)

## Getting residuals and ratios values
# Full station inventory
usgs_cvm_res = extract_res_rat_vals(eugspe_1d_df,cvm=True)
eugspe_1d_res,eugspe_1d_rat = extract_res_rat_vals(eugspe_1d_df)
eugspe_sm_res,eugspe_sm_rat = extract_res_rat_vals(eugspe_sm_df)
srv_1d_res,srv_1d_rat = extract_res_rat_vals(srv_1d_df)
srv_sm_res,srv_sm_rat = extract_res_rat_vals(srv_sm_df)

# Valley station inventory
usgs_cvm_val_res = extract_res_rat_vals(eugspe_1d_val_df,cvm=True)
eugspe_1d_val_res,eugspe_1d_val_rat = extract_res_rat_vals(eugspe_1d_val_df)
eugspe_sm_val_res,eugspe_sm_val_rat = extract_res_rat_vals(eugspe_sm_val_df)
srv_1d_val_res,srv_1d_val_rat = extract_res_rat_vals(srv_1d_val_df)
srv_sm_val_res,srv_sm_val_rat = extract_res_rat_vals(srv_sm_val_df)

## Getting stations locations
# Full station inventory
st_lons = eugspe_1d_df['st_lon'].values
st_lats = eugspe_1d_df['st_lat'].values

# Valley station inventory
st_lons_val = eugspe_1d_val_df['st_lon_val'].values
st_lats_val = eugspe_1d_val_df['st_lat_val'].values

#%%
## PLOTTING ##

## Residuals
# USGS CVM - Full 
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f2_fas_steph_res' , st_lons , st_lats , plot_res=usgs_cvm_res[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f3_fas_steph_res' , st_lons , st_lats , plot_res=usgs_cvm_res[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f4_fas_steph_res' , st_lons , st_lats , plot_res=usgs_cvm_res[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f5_fas_steph_res' , st_lons , st_lats , plot_res=usgs_cvm_res[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f6_fas_steph_res' , st_lons , st_lats , plot_res=usgs_cvm_res[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f7_fas_steph_res' , st_lons , st_lats , plot_res=usgs_cvm_res[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'pgv_steph_res' , st_lons , st_lats , plot_res=usgs_cvm_res[6])

# USGS CVM - Valley 
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f2_fas_val_steph_res' , st_lons_val , st_lats_val , plot_res=usgs_cvm_val_res[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f3_fas_val_steph_res' , st_lons_val , st_lats_val , plot_res=usgs_cvm_val_res[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f4_fas_val_steph_res' , st_lons_val , st_lats_val , plot_res=usgs_cvm_val_res[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f5_fas_val_steph_res' , st_lons_val , st_lats_val , plot_res=usgs_cvm_val_res[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f6_fas_val_steph_res' , st_lons_val , st_lats_val , plot_res=usgs_cvm_val_res[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'f7_fas_val_steph_res' , st_lons_val , st_lats_val , plot_res=usgs_cvm_val_res[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths*0 , res_map_dir , 'pgv_steph_res_val' , st_lons_val , st_lats_val , plot_res=usgs_cvm_val_res[6])

# EugSpe 1D - Full
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_eugspe_1d_res' , st_lons , st_lats , plot_res=eugspe_1d_res[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_eugspe_1d_res' , st_lons , st_lats , plot_res=eugspe_1d_res[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_eugspe_1d_res' , st_lons , st_lats , plot_res=eugspe_1d_res[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_eugspe_1d_res' , st_lons , st_lats , plot_res=eugspe_1d_res[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_eugspe_1d_res' , st_lons , st_lats , plot_res=eugspe_1d_res[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_eugspe_1d_res' , st_lons , st_lats , plot_res=eugspe_1d_res[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_1d_res' , st_lons , st_lats , plot_res=eugspe_1d_res[6])

# EugSpe 1D - Valley
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_val_eugspe_1d_res' , st_lons_val , st_lats_val , plot_res=eugspe_1d_val_res[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_val_eugspe_1d_res' , st_lons_val , st_lats_val , plot_res=eugspe_1d_val_res[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_val_eugspe_1d_res' , st_lons_val , st_lats_val , plot_res=eugspe_1d_val_res[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_val_eugspe_1d_res' , st_lons_val , st_lats_val , plot_res=eugspe_1d_val_res[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_val_eugspe_1d_res' , st_lons_val , st_lats_val , plot_res=eugspe_1d_val_res[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_val_eugspe_1d_res' , st_lons_val , st_lats_val , plot_res=eugspe_1d_val_res[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_1d_res_val' , st_lons_val , st_lats_val , plot_res=eugspe_1d_val_res[6])

# EugSpe SM - Full
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_eugspe_sm_res' , st_lons , st_lats , plot_res=eugspe_sm_res[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_eugspe_sm_res' , st_lons , st_lats , plot_res=eugspe_sm_res[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_eugspe_sm_res' , st_lons , st_lats , plot_res=eugspe_sm_res[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_eugspe_sm_res' , st_lons , st_lats , plot_res=eugspe_sm_res[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_eugspe_sm_res' , st_lons , st_lats , plot_res=eugspe_sm_res[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_eugspe_sm_res' , st_lons , st_lats , plot_res=eugspe_sm_res[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_sm_res' , st_lons , st_lats , plot_res=eugspe_sm_res[6])

# EugSpe SM - Valley
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f2_fas_val_eugspe_sm_res' , st_lons_val , st_lats_val , plot_res=eugspe_sm_val_res[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f3_fas_val_eugspe_sm_res' , st_lons_val , st_lats_val , plot_res=eugspe_sm_val_res[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f4_fas_val_eugspe_sm_res' , st_lons_val , st_lats_val , plot_res=eugspe_sm_val_res[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f5_fas_val_eugspe_sm_res' , st_lons_val , st_lats_val , plot_res=eugspe_sm_val_res[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f6_fas_val_eugspe_sm_res' , st_lons_val , st_lats_val , plot_res=eugspe_sm_val_res[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'f7_fas_val_eugspe_sm_res' , st_lons_val , st_lats_val , plot_res=eugspe_sm_val_res[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , res_map_dir , 'pgv_eugspe_sm_res_val' , st_lons_val , st_lats_val , plot_res=eugspe_sm_val_res[6])

# SRV 1D - Full
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_srv_1d_res' , st_lons , st_lats , plot_res=srv_1d_res[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_srv_1d_res' , st_lons , st_lats , plot_res=srv_1d_res[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_srv_1d_res' , st_lons , st_lats , plot_res=srv_1d_res[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_srv_1d_res' , st_lons , st_lats , plot_res=srv_1d_res[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_srv_1d_res' , st_lons , st_lats , plot_res=srv_1d_res[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_srv_1d_res' , st_lons , st_lats , plot_res=srv_1d_res[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_1d_res' , st_lons , st_lats , plot_res=srv_1d_res[6])

# SRV 1D - Valley
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_val_srv_1d_res' , st_lons_val , st_lats_val , plot_res=srv_1d_val_res[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_val_srv_1d_res' , st_lons_val , st_lats_val , plot_res=srv_1d_val_res[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_val_srv_1d_res' , st_lons_val , st_lats_val , plot_res=srv_1d_val_res[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_val_srv_1d_res' , st_lons_val , st_lats_val , plot_res=srv_1d_val_res[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_val_srv_1d_res' , st_lons_val , st_lats_val , plot_res=srv_1d_val_res[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_val_srv_1d_res' , st_lons_val , st_lats_val , plot_res=srv_1d_val_res[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_1d_res_val' , st_lons_val , st_lats_val , plot_res=srv_1d_val_res[6])

# SRV SM - Full
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_srv_sm_res' , st_lons , st_lats , plot_res=srv_sm_res[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_srv_sm_res' , st_lons , st_lats , plot_res=srv_sm_res[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_srv_sm_res' , st_lons , st_lats , plot_res=srv_sm_res[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_srv_sm_res' , st_lons , st_lats , plot_res=srv_sm_res[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_srv_sm_res' , st_lons , st_lats , plot_res=srv_sm_res[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_srv_sm_res' , st_lons , st_lats , plot_res=srv_sm_res[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_sm_res' , st_lons , st_lats , plot_res=srv_sm_res[6])

# SRV SM - Valley
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f2_fas_val_srv_sm_res' , st_lons_val , st_lats_val , plot_res=srv_sm_val_res[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f3_fas_val_srv_sm_res' , st_lons_val , st_lats_val , plot_res=srv_sm_val_res[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f4_fas_val_srv_sm_res' , st_lons_val , st_lats_val , plot_res=srv_sm_val_res[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f5_fas_val_srv_sm_res' , st_lons_val , st_lats_val , plot_res=srv_sm_val_res[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f6_fas_val_srv_sm_res' , st_lons_val , st_lats_val , plot_res=srv_sm_val_res[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'f7_fas_val_srv_sm_res' , st_lons_val , st_lats_val , plot_res=srv_sm_val_res[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , res_map_dir , 'pgv_srv_sm_res_val' , st_lons_val , st_lats_val , plot_res=srv_sm_val_res[6])

##################

## Ratios
# EugSpe 1D - Full
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_1d_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_1d_rat[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_1d_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_1d_rat[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_1d_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_1d_rat[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_1d_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_1d_rat[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_1d_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_1d_rat[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_1d_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_1d_rat[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_1d_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_1d_rat[6])

# EugSpe 1D - Valley
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_1d_val_rat[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_1d_val_rat[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_1d_val_rat[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_1d_val_rat[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_1d_val_rat[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_1d_val_rat[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_1d_val_rat[6])

# EugSpe SM - Full
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_sm_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_sm_rat[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_sm_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_sm_rat[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_sm_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_sm_rat[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_sm_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_sm_rat[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_sm_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_sm_rat[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_sm_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_sm_rat[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_sm_ratio' , st_lons , st_lats , plot_res_ratio=eugspe_sm_rat[6])

# EugSpe SM - Valley
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f2_fas_eugspe_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_sm_val_rat[0])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f3_fas_eugspe_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_sm_val_rat[1])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f4_fas_eugspe_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_sm_val_rat[2])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f5_fas_eugspe_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_sm_val_rat[3])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f6_fas_eugspe_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_sm_val_rat[4])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'f7_fas_eugspe_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_sm_val_rat[5])
plot_shade_basin_map(eugspe_lons , eugspe_lats , eugspe_depths , ratio_map_dir , 'pgv_eugspe_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=eugspe_sm_val_rat[6])

# SRV 1D - Full
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_1d_ratio' , st_lons , st_lats , plot_res_ratio=srv_1d_rat[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_1d_ratio' , st_lons , st_lats , plot_res_ratio=srv_1d_rat[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_1d_ratio' , st_lons , st_lats , plot_res_ratio=srv_1d_rat[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_1d_ratio' , st_lons , st_lats , plot_res_ratio=srv_1d_rat[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_1d_ratio' , st_lons , st_lats , plot_res_ratio=srv_1d_rat[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_1d_ratio' , st_lons , st_lats , plot_res_ratio=srv_1d_rat[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_1d_ratio' , st_lons , st_lats , plot_res_ratio=srv_1d_rat[6])

# SRV 1D - Valley
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_1d_val_rat[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_1d_val_rat[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_1d_val_rat[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_1d_val_rat[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_1d_val_rat[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_1d_val_rat[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_1d_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_1d_val_rat[6])

# SRV SM - Full
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_sm_ratio' , st_lons , st_lats , plot_res_ratio=srv_sm_rat[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_sm_ratio' , st_lons , st_lats , plot_res_ratio=srv_sm_rat[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_sm_ratio' , st_lons , st_lats , plot_res_ratio=srv_sm_rat[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_sm_ratio' , st_lons , st_lats , plot_res_ratio=srv_sm_rat[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_sm_ratio' , st_lons , st_lats , plot_res_ratio=srv_sm_rat[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_sm_ratio' , st_lons , st_lats , plot_res_ratio=srv_sm_rat[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_sm_ratio' , st_lons , st_lats , plot_res_ratio=srv_sm_rat[6])

# SRV SM - Valley
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f2_fas_srv_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_sm_val_rat[0])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f3_fas_srv_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_sm_val_rat[1])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f4_fas_srv_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_sm_val_rat[2])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f5_fas_srv_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_sm_val_rat[3])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f6_fas_srv_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_sm_val_rat[4])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'f7_fas_srv_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_sm_val_rat[5])
plot_shade_basin_map(srv_lons , srv_lats , srv_depths , ratio_map_dir , 'pgv_srv_sm_ratio_val' , st_lons_val , st_lats_val , plot_res_ratio=srv_sm_val_rat[6])







































