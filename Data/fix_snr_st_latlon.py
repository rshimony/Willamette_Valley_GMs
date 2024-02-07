#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 15:05:29 2023

@author: rshimony
"""

import pandas as pd
import numpy as np

#%%

full_st = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/station_inventory.csv')
snr_st = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/snr_stations.csv')

full_st_nm = full_st['name']
full_st_lon = full_st['longitude']
full_st_lat = full_st['latitude']

full_st_nm_unq = pd.unique(full_st_nm)
full_st_lon_unq = pd.unique(full_st_lon)
full_st_lat_unq = pd.unique(full_st_lat)


snr_st_nm = snr_st['st_nm']
snr_st_lon = snr_st['st_lon']
snr_st_lat = snr_st['st_lat']

snr_st_nm_fix = []
snr_st_lon_fix = []
snr_st_lat_fix = []

for i in range(len(snr_st_nm)):
    for j in range(len(full_st_nm_unq)):
        if snr_st_nm[i] == full_st_nm_unq[j]:
            snr_st_nm_fix.append(full_st_nm_unq[j])
            snr_st_lon_fix.append(full_st_lon_unq[j])
            snr_st_lat_fix.append(full_st_lat_unq[j])
            


fix_dict = {'st_nm':snr_st_nm_fix, 'st_lon':snr_st_lon_fix, 'st_lat':snr_st_lat_fix}
fix_df = pd.DataFrame(fix_dict)
fix_df.to_csv('/Users/rshimony/Desktop/WillametteValley/scottsmills_eq/metadata/snr_stations_fixed.csv',index=False)
























































