#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 19:14:56 2024

@author: rshimony
"""

import numpy as np
import pandas as pd
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

'''
OpenQuake shaking calculation script
This script creates OpenQuake input data (including site and basin terms),
computes shaking and saving shaking file
The final part of this script plots the calculated shaking values on top a PyGMT map.

* MIGHT NEED AN EVIRONMENT CHANGE BEFORE RUNNING THE FINAL SECTION
Since PyGMT often can't live in the same environment as OpenQuake
    
Current setting of this script:
    1) Creating OpenQuake input file for mag = 8, Rrup = 100 km, and for full model domain
    2) Extracting ONLT WV datapoints from the full input data file and writing it to a seperate input file
    3) Compute shaking for both input files:
        - Full domain, using Z2.5=nan
        - WV ONLY, using Z2.5 values
    4) Replacing the WV datapoints from the full domain shaking file (that was calculated with Z2.5=nan)
        with the WV datapoints from WV ONLY file (that was calculated using Z2.5 values). *
    5) Saving this to a new shaking file. This new file has shaking values of Z2.5=nan outside of the WV
        and Z2.5 values shaing inside the WV.
    6) Reading in the new shaking files, calculating shaking ratios, extracting ONLY WV data points from the shaking ratio array.
    7) Plotting shakemaps of Basin shakefile (Z2.5 values), No Basin shakefile (Z2.5=nan), Shaking ratio (Basin/No Basin).
        Plotting on top a PyGMT map.
        
    *This setting is set to 'mask' OpenQuake artifacts outside the WV.
'''

#%%
#### Creating OpenQuake Input Data File ####

## Reading in Vs30 file and extracting values and locations
vs30_dom = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/wv_project/vs30/vs30_dom.xygt')
vs30_lon_ls = []
vs30_lat_ls = []
vs30_ls = []
for i in range(len(vs30_dom)):
    vs30_lon = vs30_dom[i][0]
    vs30_lat = vs30_dom[i][1]
    vs30 = vs30_dom[i][2]
    
    vs30_lon_ls.append(vs30_lon)
    vs30_lat_ls.append(vs30_lat)
    vs30_ls.append(vs30) 
    
## Reading in Z-surface file and extracting values and locations (both geographic and cartesian)
z_sur = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z_surface_valley_wv_model.xyz')

z_sur_lon_ls = []
z_sur_lat_ls = []
z_sur_x_ls = []
z_sur_y_ls = []
z_sur_z1_ls = []
z_sur_z25_ls = []

for i in range(len(z_sur)):
    z_sur_lon = z_sur[i][0]
    z_sur_lat = z_sur[i][1]
    z_sur_x = z_sur[i][2]
    z_sur_y = z_sur[i][3]
    z_sur_z1 = z_sur[i][4]
    z_sur_z25 = z_sur[i][5]
    
    z_sur_lon_ls.append(z_sur_lon)
    z_sur_lat_ls.append(z_sur_lat)
    z_sur_x_ls.append(z_sur_x)
    z_sur_y_ls.append(z_sur_y)
    z_sur_z1_ls.append(z_sur_z1)
    z_sur_z25_ls.append(z_sur_z25)

## Initiating OpenQuake input data file and writing to file
#Setting parameters
rrupt = 100 #(km)
mag = 8
oq_data_file = '/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_m8_rrup100.xyz'

lons = np.array(z_sur_lon_ls)
lats = np.array(z_sur_lat_ls)

r_rups = np.full_like(lons,rrupt)
mags = np.full_like(lons,mag)
z_2pt5s = np.array(z_sur_z25_ls)*1000
vs30s = np.array(vs30_ls)

#write to file
out = np.c_[lons,lats,vs30s,z_2pt5s,r_rups,mags]
np.savetxt(oq_data_file,out,fmt='%.5f,%.5f,%d,%d,%.2f,%.2f',header='lon,lat,vs30(m/s),z2pt5(m),Rrupt(km),Mw')


#%%
## Crating WV polygon
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

#%%  
### Extracting ONLY WV datapoints from the openquake data file

def make_valley_only_data_file(full_oq_data_file_path):
    '''
    Reads in OpenQuake data file, extracting its values and locations,
    Using the WV polygon to extract ONLY WV data points from the full data file
    Returns lists of WV ONLY data (values and locations)
    '''
    full_oq_data_file = np.genfromtxt(full_oq_data_file_path, delimiter=",")
    
    oq_data_valley_lon_ls = []
    oq_data_valley_lat_ls = []
    oq_data_valley_vs30_ls = []
    oq_data_valley_z25_ls = []
    oq_data_valley_rrupt_ls = []
    oq_data_valley_mw_ls = []
    
    for i in range(len(full_oq_data_file)):
        oq_data_valley_lon = full_oq_data_file[i][0]
        oq_data_valley_lat = full_oq_data_file[i][1]
        oq_data_valley_vs30 = full_oq_data_file[i][2]
        oq_data_valley_z25 = full_oq_data_file[i][3]
        oq_data_valley_rrupt = full_oq_data_file[i][4]
        oq_data_valley_mw = full_oq_data_file[i][5]
        
        oq_data_valley_lon_ls.append(oq_data_valley_lon)
        oq_data_valley_lat_ls.append(oq_data_valley_lat)
        oq_data_valley_vs30_ls.append(oq_data_valley_vs30)
        oq_data_valley_z25_ls.append(oq_data_valley_z25)
        oq_data_valley_rrupt_ls.append(oq_data_valley_rrupt)
        oq_data_valley_mw_ls.append(oq_data_valley_mw)

    counter = 1
    
    oq_data_valley_only_lon_ls = []
    oq_data_valley_only_lat_ls = []
    oq_data_valley_only_vs30_ls = []
    oq_data_valley_only_z25_ls = []
    oq_data_valley_only_rrupt_ls = []
    oq_data_valley_only_mw_ls = []
    
    for i in range(len(full_oq_data_file)):
        if counter % 100 == 0:
            print(str(counter) + ' from ' + str(len(full_oq_data_file)))
            
        grd_point = Point(oq_data_valley_lon_ls[i],oq_data_valley_lat_ls[i])
    
        if full_poly.contains(grd_point) == True:
            oq_data_valley_only_lon_ls.append(oq_data_valley_lon_ls[i])
            oq_data_valley_only_lat_ls.append(oq_data_valley_lat_ls[i])
            oq_data_valley_only_vs30_ls.append(oq_data_valley_vs30_ls[i])
            oq_data_valley_only_z25_ls.append(oq_data_valley_z25_ls[i])
            oq_data_valley_only_rrupt_ls.append(oq_data_valley_rrupt_ls[i])
            oq_data_valley_only_mw_ls.append(oq_data_valley_mw_ls[i])
                
        counter = counter + 1
    
    return oq_data_valley_only_lon_ls, oq_data_valley_only_lat_ls, oq_data_valley_only_vs30_ls, oq_data_valley_only_z25_ls, oq_data_valley_only_rrupt_ls, oq_data_valley_only_mw_ls

## Initiating the function
# OpenQuake data file path
oq_data_valley_path = '/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_m8_rrup100.xyz'
# Running the function
valley_only_lons, valley_only_lats, valley_only_vs30, valley_only_z25, valley_only_rrupt, valley_only_mw = make_valley_only_data_file(oq_data_valley_path)

## Initiating and writing the ONLY WV data into file
oq_data_valley_only_file = '/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_m8_rrup100_valley.xyz'

lons = np.array(valley_only_lons)
lats = np.array(valley_only_lats)
r_rups = np.array(valley_only_rrupt)
mags = np.array(valley_only_mw)
z_2pt5s = np.array(valley_only_z25)
vs30s = np.array(valley_only_vs30)

#write to file
out = np.c_[lons,lats,vs30s,z_2pt5s,r_rups,mags]
np.savetxt(oq_data_valley_only_file,out,fmt='%.5f,%.5f,%d,%d,%.2f,%.2f',header='lon,lat,vs30(m/s),z2pt5(m),Rrupt(km),Mw')

#%%

import sys
sys.path.append('/Users/rshimony/Documents/open_quake/oq-engine/')

from openquake.hazardlib import imt, const
from openquake.hazardlib.contexts import RuptureContext
from openquake.hazardlib.contexts import DistancesContext
from openquake.hazardlib.contexts import SitesContext
from openquake.hazardlib.gsim.parker_2020 import ParkerEtAl2020SInterB #gmpe for cascadia


def compute_shaking_oq(oq_df,outdir,outfile,basin_term=True):
    '''
    Computes shaking using the Parker 2020 GMPE.
    Takes in an OpenQuake data file and outputs an OpenQuake shaking file (includes the computed shaking columns)
    
    If calculated WITH basin term:
        GMPE is set specifying region = 'Cascadia'
        Z2.5 values are taken from WV Z2.5 surface
        
    If calculated WITHOUT basin term:
        GMPE is set without specifying a region
        Z2.5 values = nan
    '''
    
    ## Define intensity measure and uncertainty
    imt_pgv = imt.PGV()
    uncertaintytype = const.StdDev.TOTAL
    
    if basin_term==True:
        ## Set GMPEs:
        gmpe = ParkerEtAl2020SInterB(region='Cascadia')
    if basin_term==False:
        ## Set GMPEs:
        gmpe = ParkerEtAl2020SInterB()
    
    ## Set the empty arrays:
    median_gmpe = np.array([])
    sd_gmpe = np.array([])
    
    ## Make contexts:
    rctx = RuptureContext()
    dctx = DistancesContext()
    sctx = SitesContext()
    
    
    ## Use OpenQuake dataframe to set parameters as arrays:
    dctx.rrup = oq_df['Rrupt(km)'].values
    rctx.mag = oq_df['Mw'].values
    
    if basin_term==True:
        ## Site context - seems to now need to be from a "site collection", which seems to be a pandas dataframe.
        sitecol_dict = {'sids':np.arange(1,len(oq_df['vs30(m/s)'])+1,1),'vs30':oq_df['vs30(m/s)'].values,
        'vs30measured':np.full_like(oq_df['vs30(m/s)'].values,np.nan),'z1pt0':np.full_like(oq_df['vs30(m/s)'].values,np.nan),
        'z2pt5':oq_df['z2pt5(m)'].values}
        sitecollection = pd.DataFrame(sitecol_dict)
        
    if basin_term==False:
        ## Site context - seems to now need to be from a "site collection", which seems to be a pandas dataframe.
        sitecol_dict = {'sids':np.arange(1,len(oq_df['vs30(m/s)'])+1,1),'vs30':oq_df['vs30(m/s)'].values,
        'vs30measured':np.full_like(oq_df['vs30(m/s)'].values,np.nan),'z1pt0':np.full_like(oq_df['vs30(m/s)'].values,np.nan),
        'z2pt5':np.full_like(oq_df['vs30(m/s)'].values,np.nan)}
        sitecollection = pd.DataFrame(sitecol_dict)
    
    ## Then put into a sites context:
    sctx = SitesContext(sitecol=sitecollection)    
    
    ## Run GMPE:
    ln_median_gmpe,sd_gmpe = gmpe.get_mean_and_stddevs(sctx, rctx, dctx, imt_pgv, [uncertaintytype])
    
    ## Convert median from ln g to g:
    median_gmpe = np.exp(ln_median_gmpe)
    
    ## Add as column in dataframe:
    shaking_df = oq_df.copy()
    shaking_df['PGV_parker_2020_m/s'] = median_gmpe/100
    shaking_df['PGVsd_parker_2020_lng'] = sd_gmpe[0]  ## This is in a list
    
    shaking_df.to_csv(outdir+outfile)

    return shaking_df

#%%
### Creating OpenQuake shaking files:
    # 1) Full domain, Z2.5 = nan
    # 2) WV only, Z2.5 = values
    
## Setting output directory
outputdir = '/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/'

## Reading OpenQuake data files and setting output shaking file name
# 1) Full domain, Z2.5 = nan
oq_df_m8_rrup100_full_domain = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_m8_rrup100.xyz')
oq_df_m8_rrup100_full_domain_z2pt5_nan_file_nm = 'oq_data_shaking_m8_rrup100_full_domain_z2pt5_nan.xyz'
# 2) WV only, Z2.5 = values
oq_df_m8_rrup100_wv_only = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_m8_rrup100_valley.xyz')
oq_df_m8_rrup100_wv_only_z2pt5_vals_file_nm = 'oq_data_shaking_m8_rrup100_wv_only_z2pt5_vals.xyz'

## Running shaking computation function
# 1) Full domain, Z2.5 = nan
compute_shaking_oq(oq_df_m8_rrup100_full_domain , outputdir , oq_df_m8_rrup100_full_domain_z2pt5_nan_file_nm , basin_term=False)
# 2) WV only, Z2.5 = values
compute_shaking_oq(oq_df_m8_rrup100_wv_only , outputdir , oq_df_m8_rrup100_wv_only_z2pt5_vals_file_nm)

#%%
### Taking an OpenQuake shaking file that was calculated without basin terms (Z2.5 = nan), 
### And replacing ONLY THE WV data points with data from an OpenQuake shaking file that was calculated with basin terms.

# Reading the NO BASIN TERMS file
no_basin_term_file = '/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_full_domain_z2pt5_nan.xyz'
no_basin_term_df = pd.read_csv(no_basin_term_file)
# Reading the BASIN file
valley_only_shaking_file = '/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_wv_only_z2pt5_vals.xyz'
valley_only_shaking_df = pd.read_csv(valley_only_shaking_file)

# Extracting NO BASIN TERM locations
no_basin_term_lons = no_basin_term_df['# lon'].values
no_basin_term_lats = no_basin_term_df['lat'].values
# Extracting BASIN locations
valley_only_shaking_lons = valley_only_shaking_df['# lon'].values
valley_only_shaking_lats = valley_only_shaking_df['lat'].values

# Initiating the NEW shaking file (with basin term shaking in the WV grid points)
new_full_valley_shaking_df = no_basin_term_df.copy()

# Creating an array that contains whether a point is in or out of the WV
in_out_valley = []
counter = 1
for i in range(len(no_basin_term_df)):
    if counter % 100 == 0:
        print(str(counter) + ' from ' + str(len(no_basin_term_df)))
        
    grd_point = Point(no_basin_term_lons[i],no_basin_term_lats[i])

    if full_poly.contains(grd_point) == True:
        in_out_valley.append('in_valley')
    else:
        in_out_valley.append('out_valley')
    counter = counter + 1

# Indeces of the in valley points
in_valley_idx = np.where(np.array(in_out_valley)=='in_valley')[0]

# Replacing the WV shaking data points (without basin term) with basin shaking data points
new_full_valley_shaking_df.loc[in_valley_idx,'PGV_parker_2020_m/s'] = valley_only_shaking_df['PGV_parker_2020_m/s'].values

## Writing the new haking file (with basin term shaking in the WV grid points) to file
new_full_valley_shaking_file = '/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_full_domain_z2pt5_vals.xyz'

lons = new_full_valley_shaking_df['# lon'].values
lats = new_full_valley_shaking_df['lat'].values
r_rups = new_full_valley_shaking_df['Rrupt(km)'].values
mags = new_full_valley_shaking_df['Mw'].values
z_2pt5s = new_full_valley_shaking_df['z2pt5(m)'].values
vs30s = new_full_valley_shaking_df['vs30(m/s)'].values
PGV_parker_2020 = new_full_valley_shaking_df['PGV_parker_2020_m/s'].values
PGVsd_parker_2020 = new_full_valley_shaking_df['PGVsd_parker_2020_lng'].values

#write to file
out = np.c_[lons,lats,vs30s,z_2pt5s,r_rups,mags,PGV_parker_2020,PGVsd_parker_2020]
np.savetxt(new_full_valley_shaking_file,out,fmt='%.5f,%.5f,%d,%d,%.2f,%.2f,%.17f,%.16f',header='lon,lat,vs30(m/s),z2pt5(m),Rrupt(km),Mw,PGV_parker_2020_m/s,PGVsd_parker_2020_lng')

#%%
## Reading in full domain shaking files
# Z2.5 values
shakefile_basin = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_full_domain_z2pt5_vals.xyz')
# Z2.5 = nan
shakefile_nobasin = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_full_domain_z2pt5_nan.xyz')
# Shakefile locations
shakefile_lon = shakefile_basin['# lon'].values
shakefile_lat = shakefile_basin['lat'].values

## Calculating shaking ratio (all domain data points)
shakefile_basin_vals = shakefile_basin['PGV_parker_2020_m/s'].values
shakfile_nobasin_vals = shakefile_nobasin['PGV_parker_2020_m/s'].values
shaking_ratio = shakefile_basin['PGV_parker_2020_m/s'].values/shakefile_nobasin['PGV_parker_2020_m/s'].values

#%%
'''
MIGHT NEED AN EVIRONMENT CHANGE BEFORE RUNNING THIS SECTION
Since PyGMT often can't live in the same environment as OpenQuake
'''

## Extracting ONLY WV basin shaking ratio values from full domain shaking file (for BAF calculation)
def extract_basin_values(data_file):
    counter = 1
    
    baf_data_shaking = []
    
    for i in range(len(shakefile_lon)):
        if counter % 100 == 0:
            print(str(counter) + ' from ' + str(len(shakefile_lon)))
            
        grd_point = Point(shakefile_lon[i],shakefile_lat[i])
    
        if full_poly.contains(grd_point) == True:
            baf_data_shaking.append(data_file[i])
                
        counter = counter + 1
    
    return baf_data_shaking

baf_values = extract_basin_values(shaking_ratio)

#%%

import pygmt

## Shakemap plotting parameters
fig_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/baf_shakemaps/openquake/'
cpt_values = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/baf_shakemaps/sw4_synts/Reds_09.cpt'
cpt_series_vels = (0.05,0.18)
cpt_series_ratio = (1.0,1.20)
reverse_val = False
fig_name_no_basin = 'z2pt5_nan_noregion_nobasin'
fig_name_basin = 'z2pt5_full_valley_vals'
fig_name_ratio = 'z2pt5_full_valley_vals_ratio'

## shakemap plotting function
def make_shaking_map(shaking_data, cpt_vals, cpt_sers, reverse, fig_nm, fig_dir, ratio = False):
    '''
    Plotting OpenQuake shaking data on top a PyGMT map.
    Plotting parameters change depending on whether it is ratio value or not (default - not ratio)
    '''
    
    fig = pygmt.Figure()
    ## topo:
    topo_data = '@earth_relief_03s'
    
    plot_region = [-124.19,-121.38,43.4,46.8]
    
    with pygmt.config(FONT_ANNOT_PRIMARY=38,FORMAT_GEO_MAP="ddd"):
        fig.basemap(region=plot_region, projection="M0/0/10c", frame=["Wsne","ya1f0.25+lLongutide(°)"])        
        fig.grdimage(topo_data, cmap="grey",shading=True)
    with pygmt.config(FONT_ANNOT_PRIMARY=38,FORMAT_GEO_MAP="ddd"):
        fig.basemap(region=plot_region, projection="M0/0/10c", frame=["wSne", "xa1f0.25+lLatitude(°)"]) 
    
    
    pygmt.makecpt(
        cmap=cpt_vals,reverse=reverse,
        series = cpt_sers)
    fig.plot(x=shakefile_lon, y=shakefile_lat, 
              style="c0.015c", 
              fill=shaking_data,
              cmap=True,
              transparency=60)
    
    if ratio == False:
        with pygmt.config(FONT_ANNOT_PRIMARY=32,FONT_ANNOT_SECONDARY=32,FONT_TITLE=32,FONT_LABEL=32):
            fig.colorbar(position="JMR+o0.5c/0c+w18c/0.8c",frame=["af+lVelocity [m/s]"])
        
    if ratio == True:
        with pygmt.config(FONT_ANNOT_PRIMARY=32,FONT_ANNOT_SECONDARY=32,FONT_TITLE=32,FONT_LABEL=32):
            fig.colorbar(position="JMR+o0.5c/0c+w18c/0.8c",frame=["af+lRatio [Basin/No Basin]"])
        
    with pygmt.config(FONT_ANNOT_PRIMARY=32,FONT_ANNOT_SECONDARY=32,FONT_TITLE=32,FONT_LABEL=32,MAP_TICK_PEN_PRIMARY=2.2):
        fig.basemap(map_scale="JBR+o-4.4/-1.3+w50")
    with pygmt.config(FONT_ANNOT_PRIMARY=24,FONT_ANNOT_SECONDARY=22,FONT_TITLE=22,FONT_LABEL=22):
        fig.basemap(rose="JTL+o-2.7c+w1.8c+f1+l")
    
    fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
    fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
    fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
    
    if ratio == True: 
        fig.text(x=-122.2,y=43.9,text='BAF = '+str(np.array(baf_values).mean())[:5],font="26p,Helvetica-BoldOblique,black")
    
    # fig.show()
    fig.savefig(fig_dir +fig_nm+ '.png',dpi=300)

#%%
## Running shakemap plotting functions for 
    # 1) shakemap with Z2.5 values
    # 2) shakemap with Z2.5 = nan
    # 3) shaking ratio map (values/nans or basin/no basin)
make_shaking_map(shakefile_basin_vals, cpt_values, cpt_series_vels, reverse_val, fig_name_basin, fig_dir)
make_shaking_map(shakfile_nobasin_vals, cpt_values, cpt_series_vels, reverse_val, fig_name_no_basin, fig_dir)
make_shaking_map(shaking_ratio, cpt_values, cpt_series_ratio, reverse_val, fig_name_ratio, fig_dir, ratio = True)




































