#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 13:21:29 2023

@author: rshimony
"""

import fiona
import pandas as pd
import geopandas as gpd

#%%

stations = '/Users/rshimony/Documents/Spatial_Analysis/salem_obs_pgv_GEOG591_fproject.csv'
wells = '/Users/rshimony/Documents/Spatial_Analysis/basin_depth_GEOG591_project.csv'
poly = '/Users/rshimony/Documents/Spatial_Analysis/basin_poly.csv'
stations_clip = '/Users/rshimony/Documents/Spatial_Analysis/station_table.csv'

outdir = '/Users/rshimony/Documents/Spatial_Analysis/'

stationsDf = pd.read_csv(stations,header=0)

stations_c_Df = pd.read_csv(stations_clip,header=0)

wellsDf = pd.read_csv(wells,header=0)

polyDf = pd.read_csv(poly,header=0)

#%%

schema_point = {
    'geometry':'Point',
    'properties':[('Name','str'),
                  ('PGV', 'float')]
}

#open a fiona object
stationsShp = fiona.open(outdir + 'stations.shp', mode='w', driver='ESRI Shapefile',
          schema = schema_point, crs = "EPSG:4326")
#iterate over each row in the dataframe and save record
for index, row in stationsDf.iterrows():
    rowDict = {
        'geometry' : {'type':'Point',
                      'coordinates': (row.st_lons,row.st_lats)},
        'properties': {'Name' : row.st_names,
                        'PGV' : row.pgv_t_obs},
    }
    stationsShp.write(rowDict)
#close fiona object
stationsShp.close()


#%%

schema_point_c = {
    'geometry':'Point',
    'properties':[('Name','str'),
                  ('PGV', 'float'),
                  ('IDW', 'float'),
                  ('krig', 'float')]
}

#open a fiona object
stations_c_Shp = fiona.open(outdir + 'stations_c.shp', mode='w', driver='ESRI Shapefile',
          schema = schema_point_c, crs = "EPSG:4326")
#iterate over each row in the dataframe and save record
for index, row in stations_c_Df.iterrows():
    rowDict = {
        'geometry' : {'type':'Point',
                      'coordinates': (row.st_lons,row.st_lats)},
        'properties': {'Name' : row.Name,
                        'PGV' : row.PGV,
                        'IDW' : row.IDW_raster,
                        'krig' : row.krig_raster},
    }
    stations_c_Shp.write(rowDict)
#close fiona object
stations_c_Shp.close()

#%%

# schema_wells = {
#     'geometry':'Point',
#     'properties':[('Name','str'),
#                   ('Depth', 'int')]
# }

# #open a fiona object
# wellsShp = fiona.open(outdir + 'wells.shp', mode='w', driver='ESRI Shapefile',
#           schema = schema_wells, crs = "EPSG:4326")
# #iterate over each row in the dataframe and save record
# for index, row in wellsDf.iterrows():
#     rowDict = {
#         'geometry' : {'type':'Point',
#                      'coordinates': (row.Longitude,row.Latitude)},
#         'properties': {'Name' : row.SiteName,
#                        'Depth' : row.Depth},
#     }
#     wellsShp.write(rowDict)
# #close fiona object
# wellsShp.close()

#%%

schema_poly = {
    'geometry':'Polygon',
    'properties':[('Name','str')]
}

#open a fiona object
polyShp = fiona.open(outdir + 'poly.shp', mode='w', driver='ESRI Shapefile',
          schema = schema_poly, crs = "EPSG:4326")

#get list of points
xyList = []
rowName = ''
for index, row in polyDf.iterrows():
    xyList.append((row.lon,row.lat))
    rowName = row.name_b
xyList[:5]

#save record and close shapefile
rowDict = {
'geometry' : {'type':'Polygon',
                 'coordinates': [xyList]}, #Here the xyList is in brackets
'properties': {'Name' : rowName},
}
polyShp.write(rowDict)
#close fiona object
polyShp.close()

#%%

data_shp = gpd.read_file('/Users/rshimony/Documents/Spatial_Analysis/poly.shp')


























