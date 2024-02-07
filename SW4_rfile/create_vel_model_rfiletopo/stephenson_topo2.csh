#!/bin/bash

# To download Topography & bathymetry data (from OpenSWPC)
#if [ ! -e ETOPO1_Bed_g_gmt4.grd.gz ]; then
 #   curl -o ETOPO1_Bed_g_gmt4.grd.gz https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd.gz
#fi
#gunzip ETOPO1_Bed_g_gmt4.grd.gz
region=-124.19/-121.51/43.4/46.1

# Topography grd file
topo=./ETOPO1_Bed_g_gmt4.grd
gmt grdcut -R$region $topo -Gtopo.cut_csz.grd

gmt grd2xyz topo.cut_csz.grd > topo.cut_csz.dat

