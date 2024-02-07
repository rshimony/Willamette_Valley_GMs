#!/bin/bash

# To download Topography & bathymetry data (from OpenSWPC)
#if [ ! -e ETOPO1_Bed_g_gmt4.grd.gz ]; then
#   curl -o ETOPO1_Bed_g_gmt4.grd.gz https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd.gz
#fi
#unzip ETOPO1_Bed_g_gmt4.grd.gz

# cut version
#region=-125.9/-120.1/40/48

# Topography grd file
#topo=./ETOPO1_Bed_g_gmt4.grd
#gmt grdcut -R$region $topo -Gtopo.cut_csz.grd

#gmt grd2xyz topo.cut_csz.grd > topo.cut_csz.dat




region=-140.0/-115.0/39/55 
#
# Topography grd file
topo=./ETOPO1_Bed_g_gmt4.grd
gmt grdcut -R$region $topo -Gtopo.csz.grd

gmt grd2xyz topo.csz.grd > topo.csz.dat



#f_lon[f_z>-1000].min()
#Out[203]: -129.80000000000001

#f_lon[f_z>-1000].max()
#Out[204]: -121.0

#f_lat[f_z>-1000].max()
#Out[205]: 51.0

#f_lat[f_z>-1000].min()
#Out[206]: 39.0
