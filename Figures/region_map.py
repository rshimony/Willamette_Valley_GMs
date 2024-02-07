#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 13:06:56 2023

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
from scipy import spatial

#%%

events = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/event_catalog_valevents.csv')

evlon = events['longitude']
evlat = events['latitude']
evcode = events['Code'] #B-ScottstMills, D-Springfield, E-Salem

#%% NEED TO CHANGE TO VALLEY STATIONS

st_file_f_salem = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/domain_stations.csv')
st_file_salem = st_file_f_salem.drop_duplicates()
st_names_salem = np.array(st_file_salem['st_name_dom'])
st_lons_salem = np.array(st_file_salem['st_lon_dom'])
st_lats_salem = np.array(st_file_salem['st_lat_dom']) 

st_file_f_scottsmills = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/snr_stations_fixed.csv')
st_file_scottsmills = st_file_f_scottsmills.drop_duplicates()
st_names_scottsmills = np.array(st_file_scottsmills['st_nm'])
st_lons_scottsmills = np.array(st_file_scottsmills['st_lon'])
st_lats_scottsmills = np.array(st_file_scottsmills['st_lat']) 

st_file_f_springfield = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/springfield_eq/metadata/snr_stations_fixed.csv')
st_file_springfield = st_file_f_springfield.drop_duplicates()
st_names_springfield = np.array(st_file_springfield['st_nm'])
st_lons_springfield = np.array(st_file_springfield['st_lon'])
st_lats_springfield = np.array(st_file_springfield['st_lat']) 

#%%

depth2units_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/depth2units_wells/depth2units_wells.csv')

lon_units = depth2units_f['Longitude']
lat_units = depth2units_f['Latitude']
srv_depth = depth2units_f['srv_depth']
yamhill_depth = depth2units_f['yamhill_depth']
spencer_depth = depth2units_f['Spencer_depth']
eugene_depth = depth2units_f['eugene_depth']

nonan_srv_depth = srv_depth[srv_depth.notna()]
nonan_srv_lons = lon_units[srv_depth.notna()]
nonan_srv_lats = lat_units[srv_depth.notna()]

nonan_eugspe_depth = spencer_depth[spencer_depth.notna()]
nonan_eugspe_lons = lon_units[spencer_depth.notna()]
nonan_eugspe_lats = lat_units[spencer_depth.notna()]

#%%

inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/wv_poly.csv')
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

xpnn = np.array([xpn_m])
xph1n = np.array([xph1_m])
xph2n = np.array([xph2_m])

ypnn = np.array([ypn_m])
yph1n = np.array([yph1_m])
yph2n = np.array([yph2_m])

xpnt = np.concatenate((xpnn, xph1n , xph2n), axis=1)
ypnt = np.concatenate((ypnn, yph1n , yph2n), axis=1)

#%%
f_path = '/Users/rshimony/Desktop/zzthick4.axyz'

# f = np.genfromtxt(f_path)

# f_lat = []
# f_lon = []
# f_z = []

# for i in range(len(f)):
#     f_lon.append(f[i][0])
#     f_lat.append(f[i][1])
#     f_z.append(f[i][2])

xy = np.loadtxt(f_path)
x, y = xy[:, 0].astype('float64'), xy[:, 1].astype('float64')

uniq_f_lat = np.unique(y)
uniq_f_lon = np.unique(x)

f_lat_min = min(uniq_f_lat)
f_lat_max = max(uniq_f_lat)
f_lon_min = min(uniq_f_lon)
f_lon_max = max(uniq_f_lon)

grav_dom = Polygon([[f_lon_min,f_lat_min] , [f_lon_max,f_lat_min] , [f_lon_max,f_lat_max] , [f_lon_min,f_lat_max]])

xpn, ypn = grav_dom.exterior.xy

plt.plot([f_lon_min,f_lon_max,f_lon_min,f_lon_max], [f_lat_min,f_lat_max,f_lat_max,f_lat_min], 'o', color='red', markersize=6)
plt.plot(xpn, ypn, color='blue')
plt.show()


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
    region=[-125,-121.51,43.4,46.1]
    fig.basemap(region=region, projection='M0/0/10c')
    fig.grdimage(topo_data, cmap="geo",shading=True,frame=True)
    fig.coast(water="skyblue3")

    
    #events location
    fig.plot(x=evlon[[1,3,4]] , y=evlat[[1,3,4]] , style='a0.5c',fill='red',pen='black')
    fig.text(text=evcode[[1,3,4]], x=evlon[[1,3,4]], y=evlat[[1,3,4]]-0.12 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    
    #basin outline
    fig.plot(x=xpn_m, y=ypn_m,pen="1p,black")
    fig.plot(x=xph1_m, y=yph1_m,pen="1p,black")
    fig.plot(x=xph2_m, y=yph2_m,pen="1p,black")
    
    # Plot the shapefile data
    fig.plot(
        data=shapefile,
        style="p0.009p,black",
        pen="black")
    
    #McPhee outline
    fig.plot(x=xpn, y=ypn,pen="1p,black,-")
    
    #cities location names
    fig.text(text='Eugene', x=-123.2, y=44.2 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text='Corvalis', x=-123.4, y=44.67 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text='Salem', x=-123.25, y=45.0 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text='Portland', x=-122.65, y=45.7 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    fig.text(text=evcode[[1,3,4]], x=evlon[[1,3,4]], y=evlat[[1,3,4]]-0.12 , font='12p,Helvetica-Bold,royalblue4',pen='black',fill='whitesmoke')
    
    #stations
    fig.plot(x=st_lons_salem , y=st_lats_salem , style='t0.3c',fill='green',pen='black',transparency=50)
    fig.plot(x=st_lons_scottsmills , y=st_lats_scottsmills , style='t0.3c',fill='blue',pen='black',transparency=50)
    fig.plot(x=st_lons_springfield , y=st_lats_springfield , style='t0.3c',fill='yellow',pen='black',transparency=50)
    
    #wells
    fig.plot(x=nonan_srv_lons , y=nonan_srv_lats , style='c0.3c',fill='grey',pen='black',transparency=50)
    fig.plot(x=nonan_eugspe_lons , y=nonan_eugspe_lats , style='c0.3c',fill='purple',pen='black',transparency=50)
    
    # Create an inset map, setting the position to bottom right, and the x- and
    # y-offsets to 0.1 cm, respectively.
    # The inset map contains the Japan main land. "U54S/3c" means UTM projection
    # with a map width of 3 cm. The inset width and height are automatically
    # calculated from the specified ``region`` and ``projection`` parameters.
    # Draws a rectangular box around the inset with a fill color of "white" and
    # a pen of "1p".
    with fig.inset(
        position="jBL+o0.1c",
        box="+gwhite+p1p",
        region=[-127,-118,38,51],
        projection="M3c",
    ):
        # Highlight the Japan area in "lightbrown"
        # and draw its outline with a pen of "0.2p".
        fig.coast(
            land="gray",
            borders=[1, 2],
            shorelines="1/thin",
            water="white",
            # dcw="US+glightbrown+p0.2p",
            area_thresh=10000,
        )
        
        fig.plot(x=xpn_m, y=ypn_m,pen="0.2p,red")
        fig.plot(x=xph1_m, y=yph1_m,pen="0.2,red")
        fig.plot(x=xph2_m, y=yph2_m,pen="0.2,red")
        
        # fig.grdimage(topo_data, cmap="geo",shading=True,frame=True)
        
        # Plot a rectangle ("r") in the inset map to show the area of the main
        # figure. "+s" means that the first two columns are the longitude and
        # latitude of the bottom left corner of the rectangle, and the last two
        # columns the longitude and latitude of the upper right corner.
        rectangle = [[region[0], region[2], region[1], region[3]+0.2]]
        fig.plot(data=rectangle, style="r+s", pen="2p,blue")

# Show the plot
fig.savefig('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/region_data_map.png',dpi=300)


























