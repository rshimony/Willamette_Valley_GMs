#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 16:30:51 2023

@author: rshimony
"""

import matplotlib.pyplot as plt
import sys
import os

flush = sys.stdout.flush()

from pySW4.utils import geo
from pySW4.prep import rfileIO
import numpy as np
import time
import math
from warnings import warn
import xarray as xr
from scipy.interpolate import griddata
import pandas as pd
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point,LineString

#%% Inputs and options
start = time.time()

hhi = (200, 300, 900)
hvi = (100,  300,  900)
block_layers = (0, 1200, 9900, 59400) 

outpath = '/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin_grad'
os.chdir(outpath)

rfileoutput = 'stephenson_rfile_notopo_cartesian_grad_eugspe.rfile'

# rfile input
CVM_L1_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vp_16l1_my_model.nc'    
CVM_L1_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vs_16l1_my_model.nc'
CVM_L2_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vp_16l2_my_model.nc'
CVM_L2_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vs_16l2_my_model.nc'
CVM_L3_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vp_16l3_my_model.nc'
CVM_L3_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/vs_16l3_my_model.nc'

# CVM_L1_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vp_16l1.nc'    
# CVM_L1_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vs_16l1.nc'
# CVM_L2_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vp_16l2.nc'
# CVM_L2_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vs_16l2.nc'
# CVM_L3_P = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vp_16l3.nc'
# CVM_L3_S = '/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/vs_16l3.nc'

'''
coordinated of SW corner - full Stephenson model

lat = 40.2
lon = -129
'''

filenames = [CVM_L1_P, CVM_L1_S, CVM_L1_P, CVM_L1_S, CVM_L2_P, CVM_L2_S, CVM_L3_P, CVM_L3_S]

block1_vp = xr.open_dataset(CVM_L1_P)
block1_vs = xr.open_dataset(CVM_L1_S)

block2_vp = xr.open_dataset(CVM_L2_P)
block2_vs = xr.open_dataset(CVM_L2_S)

block3_vp = xr.open_dataset(CVM_L3_P)
block3_vs = xr.open_dataset(CVM_L3_S)

# Creating an array for the velocity model
ifile_in = '/Users/rshimony/Desktop/WillametteValley/Models/cartesian_ifile/ifile_cartesian_4rfile_cal.txt' #'ALL_japan_path20111110.dat'
n_tif = 2
dx = 100 #(m)
dy = -100 #(m)

vel_model = np.array([\
[0,   -999,  -999,   -999, -999,-999],
[1,  1.987e3, 0.6e3, 1.899e3, 78.34,  39.17],
[2,  1e20, 1e20, 1e20, 1e20,  1e20]])

vel_model_grad = np.array([\
[0,   0,  0,   0, 0,0],
[1,  0.614, 0.487, 0.119, 0,  0],
[2,  0, 0, 0, 0, 0]])
        
# Extent of the resulting rfile (m) [w,s=0; e=x_dim; n=y_dim] 
w = 0;  e = 216e3 #215e3
s = 0;  n = 378e3 #375e3
    
#%% helper functions
def calcQ(vs3d):

    # print('Setting minimum Qs=60\n')

    vs = vs3d.flatten()
    qs = np.zeros(len(vs))

    for i in range(len(vs)):
        if (vs[i] < 1000):
            qs[i] = (0.23*vs[i]) - 80
        elif (vs[i] >= 1000 and vs[i] <= 2500):
            qs[i] = 150 + 0.1*(vs[i] - 1000)
        else:
            qs[i] = 300.

        if qs[i] < 60:
                qs[i] = 60

    qp = 2*qs

    return qp.reshape(vs3d.shape), qs.reshape(vs3d.shape)


def calcRho(vp3d):

    # From Stpehnson et al., 2017 - Using the empirical relations of Brocher (2005)
    vp = vp3d.flatten()/1000. # Put in km/s
    rho = np.zeros(len(vp))

    for i in range(len(vp)):
        rho[i] = 1.6612*vp[i] - 0.4721*vp[i]**2 + 0.0671*vp[i]**3 - 0.0043*vp[i]**4 + 0.000106*vp[i]**5
        if rho[i] < 2.0:
            rho[i] = 2.0
        elif rho[i] > 3.5:
            rho[i] = 3.5

    rho = rho * 1000 # Back in kg/m3

    return rho.reshape(vp3d.shape)

def fix_vp_vs_ratio(min_vp_vs):
    utme = block1_vp.utme.values
    utmn = block1_vp.utmn.values
    
    X,Y = np.meshgrid(utme,utmn)
    X_t = X.T
    Y_t = Y.T
    
    X_t_flat = X_t.flatten()
    Y_t_flat = Y_t.flatten()
    
    # Looping over each layer in block 1
    vp_vs_l = []
    vp_vs_init = []
    for i in range(len(block1_vp.z.values)):
    
        vs_i = block1_vs.vs.values[:,:,i]
        vp_i = block1_vp.vp.values[:,:,i]
        
        vs_flat = vs_i.flatten()
        vp_flat = vp_i.flatten()
        vp_vs = vp_flat/vs_flat
        vp_vs_reshape = vp_vs.reshape(vp_i.shape)
        vp_vs_init.append(vp_vs_reshape)
        
        # training data
        vp_vs_train = vp_vs[vp_vs > min_vp_vs]
        
        utme_train = X_t_flat[vp_vs > min_vp_vs]
        utmn_train = Y_t_flat[vp_vs > min_vp_vs]
        
        # querry at low vp_vs
        utme_low = X_t_flat[(vp_vs < min_vp_vs) & (vp_flat < 1e4)]
        utmn_low = Y_t_flat[(vp_vs < min_vp_vs) & (vp_flat < 1e4)]
        
        # Interpolate vp_vs ratio map
        vp_vs_corr = griddata((utme_train, utmn_train), vp_vs_train, (utme_low, utmn_low), method='linear')
        # print(len(vp_vs_corr))
        
        # determine vp at the corrected locations
        vp_corr = vp_vs_corr * vs_flat[(vp_vs < min_vp_vs) & (vp_flat < 1e4)]
        # print(len(vp_corr))
        
        # print(i)
        # print(len(vp_flat[(vp_vs < min_vp_vs) & (vp_flat < 1e4)]))
        vp_flat[(vp_vs < min_vp_vs) & (vp_flat < 1e4)] = vp_corr
        
        # print(len(vp_vs[(vp_vs < min_vp_vs) & (vp_flat < 1e4)]))
        vp_vs[(vp_vs < min_vp_vs) & (vp_flat < 1e4)] = vp_vs_corr
        
        vp_reshape = vp_flat.reshape(vp_i.shape)
        vp_vs_reshape = vp_vs.reshape(vp_i.shape)
        vp_vs_l.append(vp_vs_reshape)
    
        block1_vp.vp.values[:,:,i] = vp_reshape
        
    return vp_vs_l,vp_vs_init



def plot_CVM_cs(data_list,vmin,vmax,map_label,filename):
    
    import matplotlib.colors as mcolors
    norm = mcolors.TwoSlopeNorm(vcenter = np.sqrt(2))
    
    utme = block1_vp.utme.values
    utmn = block1_vp.utmn.values
    z = block1_vp.z.values
    
    figw, figh = 10.0, 4.0
    fig, axs = plt.subplots(1,4,figsize=(figw,figh), sharey='row')
    # fig.subplots_adjust(wspace=0.5,hspace=0.1)
    # plt.subplots_adjust(wspace = 0.25 ,hspace = 0.15, right=0.8)
    fig.tight_layout(rect=[0, 0.03, 0.85, 0.90])
    
    axs=axs.flatten()
    
    for i in range(4):
        
        if map_label == 'Vp/Vs fixed' or map_label == 'Vp/Vs':
            im = axs[i].imshow(data_list[i].T, cmap ='seismic', vmin = vmin, vmax = vmax, #vp_vs_l
                              extent =[utme.min(), utme.max(), utmn.min(), utmn.max()],
                                interpolation ='nearest', origin ='lower',norm=norm)
            
        else:   
            im = axs[i].imshow(data_list[:,:,i].T, cmap ='viridis', vmin = vmin, vmax = vmax,
                             extent =[utme.min(), utme.max(), utmn.min(), utmn.max()],
                                interpolation ='nearest', origin ='lower',norm=norm)
            
        axs[i].set_xlabel('utme',fontsize=14)
        axs[i].set_ylabel('utmn',fontsize=14)
        axs[i].tick_params(labelsize=12)
        axs[i].set_xticks((450000, 550000))
        axs[i].set_title('Depth = ' + str(z[i]) + 'm',fontsize=10)
    
        fig.suptitle('Stephenson Model horizontal cross sections', fontsize=15)
        # CS = axs[i].contour(X_t,Y_t,vp_vs_l[i],levels=[1,1.4],colors=['lime','k'] ,linewidths=0.7)
    
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar=fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(map_label)
 
    plt.savefig('/Users/rshimony/Documents/Sci_Mat/Stephenson_model/nc_files/my_model_nc/'+filename+'.png',dpi=300,bbox_inches='tight')
    plt.close('all')


def llz2utm(lon,lat,projection_zone='None'):
    '''
    Convert lat,lon to UTM
    '''
    from numpy import zeros,where,chararray
    import utm
    from pyproj import Proj
    from scipy.stats import mode
    
    x=zeros(lon.shape)
    y=zeros(lon.shape)
    zone=zeros(lon.shape)
    b=chararray(lon.shape)
    if projection_zone==None:
        #Determine most suitable UTM zone
        for k in range(len(lon)):
            #x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k]-360)
            x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k])
        zone_mode=mode(zone)
        i=where(zone==zone_mode)[0]
        letter=b[i[0]]
        z=str(int(zone[0]))+letter
    else:
        z=projection_zone
    p = Proj(proj='utm',zone=z,ellps='WGS84')
    x,y=p(lon,lat)
    return x,y

def ifile2tiff(ifile_in):
    print('\nConverting ifile (%s) to tiffs'  %ifile_in)

    data = np.genfromtxt(ifile_in,skip_header=1) 
    nx = len(np.unique(data[:,0]))
    ny = len(np.unique(data[:,1]))
    tlx = 0
    tly = n-s #1000

    proj4 = ('+proj=tmerc +datum=NAD83 +units=m +lon_0=129 +lat_0=30 +no_defs')

    rasterBand = 1
    nodata = np.nan

    # Creating and reading TIFF for each layer
    for i in range(n_tif):
        j = i+1
        
        # Extracting and Reshaping data for each layer
        globals()['elev_layer_tmp%s' % j] = data[:,i+2].reshape(ny,nx)
        globals()['elev_layer%s' % j] = -globals()['elev_layer_tmp%s' % j][::-1]
        
        # Creating TIFF for each layer
        geo.save_GeoTIFF('elev_layer'+str(j)+'.tif', globals()['elev_layer%s' % j], 
                          tlx, tly, dx, dy, epsg=None, proj4=proj4, 
                          dtype=np.float32, nodata=np.nan, rasterBand=1)

        os.system('mkdir -p tiff');
        os.system('mv elev_layer*.tif tiff');
    
    return

def resample_to_hhi(hhi,w, e, s, n):
    fields = '{:>14}|{:>5}|{:>5}|{:>7}|{:>7}'
    data = '{:>14}|{:>5}|{:>5}|{:>7.2f}|{:>7.2f}'
    
    dirs = []
    
    for b, hhi_ in enumerate(hhi):
        
        print('\nResampling grids to ~{} m for block No.{}...'.format(hhi_, b + 1))
    
        path = 'Grids_{}_hh={}m'.format(b + 1, hhi_)
        if not os.path.exists(path):
            os.mkdir(path)
        dirs += [path]
    
        print(fields.format('grid', 'nx', 'ny', 'dx', 'dy'))
        flush
        
        for jj in range(n_tif):
            
            grid = globals()['layer%s' % str(jj+1)]
            outfile = os.path.join(path, grid.name[:-4] + '_hh='+str(hhi_)+'m.tif')
            
            if os.path.exists(outfile):
                print('Skipping grid {}, '
                      'file already exists in {}...'.format(grid.name, path))
                continue
            
            g = grid.copy()
            set_new_extent(g,w, e, s, n)
            g.resample(to=hhi_)
            
            # check for nodata:
            nodata = grid.nodata
            if (grid.elev == nodata).any():
                warn(
                    'Data in grid {} contains '
                    'pixels with `nodata`.'.format(grid.name))
                nodata_grid = grid.copy()
                nodata_grid.elev = (grid.elev == nodata)
                nodata_grid.resample(to=hhi_, order=0)
                idx = nodata_grid.elev
                g.elev[idx] = nodata
    
            print(data.format(g.name.split('.')[0], g.nx, g.ny, g.dx, g.dy))
            flush
    
            g.write_GeoTIFF(outfile)
    
        print
    
    print('Done resampling grids.')
    flush
        
    return dirs  

# def correct_nans_bottom_tiff():
#     for l in range(len(hhi)):
#         path = 'Grids_{}_hh={}m'.format(l+1, hhi[l])
        
#         for i in range(n_tif):
            
#             # Concatenating tif filename
#             tif_file = 'elev_layer'+str(i+1)+'_hh='+str(hhi[l])+'m_old.tif'
            
#             # Reading grid files for dirs[b] or d
#             tiff_data = geo.read_GeoTIFF(os.path.join(path, tif_file))
            
#             tiff_data.elev = tiff_data.elev[::-1]
            
#             tiff_line = tiff_data.elev[1:,:1]
#             nan_arr = np.isnan(tiff_line)
#             nan_true = nan_arr[nan_arr==True]
#             nan_val = len(nan_true)+1
            
#             new_tiff_data = tiff_data
#             new_tiff_data.elev = np.vstack((tiff_data.elev[nan_val:,:], tiff_data.elev[:nan_val,:]))[::-1]
            
#             outfile = os.path.join(path, 'elev_layer'+str(i+1)+'_hh='+str(hhi[l])+'m.tif')
            
#             new_tiff_data.write_GeoTIFF(outfile)
            
#     return


def set_new_extent(self, w, e, s, n, fill_value=None, mask=False):
    """
    Change the extent of a GeoTIFF file, trimming the boundaries or
    padding them as needed.
    .. warning:: **Watch Out!** This operation is performed in place
        on the actual data. The raw data will no longer be
        accessible afterwards. To keep the original data, use the
        :meth:`~pySW4.utils.geo.GeoTIFF.copy` method to create a
        copy of the current object.
    Parameters
    ----------
    w, e, s, n : float
        The west-, east-, south-, and north-most coordinate to keep
        (may also be the xmin, xmax, ymin, ymax coordinate).
    fill_value : float or None
        Value to pad the data with in case extent is expanded beyond
        the gdata. By default this is set to the ``nodata`` value.
    mask : bool
        If set to ``True``, data is masked using the `fill_value`.
    """
    assert w < e, '`w` must be less than `e`'
    assert s < n, '`s` must be less than `n`'

    npw = int((w - self.w) / self.dx)
    npe = int((e - self.w) / self.dx)

    nps = int((s - self.n) / self.dy)
    npn = int((n - self.n) / self.dy)
    
    # trimming all dimensions
    if 0 <= npw < npe <= self.nx and 0 <= npn < nps <= self.ny:
        self.elev = self.elev[npn: nps, npw: npe]

    # expanding or trimming each dimension independently
    else:
        fill_value = self.nodata if fill_value is None else fill_value

        # first handle east and south boundaries
        if npe > self.nx:
            self.elev = np.pad(self.elev, ((0, 0), (0, npe - self.nx)),
                            'constant', constant_values=fill_value)
        else:
            self.elev = self.elev[:, :npe]
        if nps > self.ny:
            self.elev = np.pad(self.elev, ((0, nps - self.ny), (0, 0)),
                            'constant', constant_values=fill_value)
        else:
            self.elev = self.elev[:nps, :]

        # now handle west and north boundaries
        if npw < 0:
            self.elev = np.pad(self.elev, ((0, 0), (-npw, 0)),
                            'constant', constant_values=fill_value)
        else:
            self.elev = self.elev[:, npw:]
        if npn < 0:
            self.elev = np.pad(self.elev, ((-npn, 0), (0, 0)),
                            'constant', constant_values=fill_value)
        else:
            self.elev = self.elev[npn:, :]

        if mask:
            self.elev = np.ma.masked_equal(self.elev, fill_value)

    # update class data
    self.w = self.w + npw * self.dx
    self.n = self.n + npn * self.dy
    
    return


def reading_grids():
    print('\nReading grids...'),
    flush

    for b, d in enumerate(dirs):
        print('{}:'.format(b + 1)),
        print('Printing tifname nx ny min and max topography of each layer (depth is positive downwards)...'),
        flush
        
        lat_dim = []
        lon_dim = []
        
        for i in range(n_tif):
            
            # Concatenating tif filename
            tif_file = 'elev_layer'+str(i+1)+'_hh='+str(hhi[b])+'m.tif'
            
            # Reading grid files for dirs[b] or d
            globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'] = \
                    geo.read_GeoTIFF(os.path.join(d, tif_file))
                
            lat_dim.append(len(globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].elev))
            lon_dim.append(len(globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].elev[0]))
            
            # Converting the elevation to positive downwards
            globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].vele = (
             -1 * globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].elev[::-1])
            
            if min(lat_dim)==max(lat_dim) and min(lon_dim)==max(lon_dim):
                
                print(globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].name,\
                      len(globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].elev),\
                      len(globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].elev[0]),\
                      globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].vele.min(),\
                      globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].vele.max())
            else:
                print('The resampled tiffs do not have the same dimension. Try fixing it before proceeding...')
                exit()
                
        print('Done!')
        print('\n')
           
    # Fixing the topography of all surfaces
    # It is important to make sure surfaces do not cross each other. 
    for b, d in enumerate(dirs):
        print('{}:'.format(b + 1)),
        flush
        
        for i in range(n_tif-1):
            idx = globals()['elev_layer' + str(i + 2)+'_hh='+str(hhi[b])+'_stream'].elev > \
                globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].elev
        
            isFalseIn = True in idx
            
            if isFalseIn == True:
                print('Some layers are crossing each other. Try fixing it before proceeding...')
                exit()
            
        print('\nNo crossing layers...')
            
    return

#%%
def create_k_array_grad(i,j,hhi_,n_layers,zbl,vel_model,vel_model_grad,hvi_,blk):

    # Determing k_array at location index [i,j]
    
    # Get UNIQUE depth interfaces at location index [i,j]
    dep = np.array([])
    sn = np.array([])

    for index in range(n_layers):
        sn = np.append(sn, index+1)
        dep = np.append(dep,globals()['elev_layer'+ str(index+1) +\
                                      '_hh='+str(hhi_)+'_stream'].vele[i,j])

    dep_interf = np.unique(dep)
    
    if np.isnan(np.sum(dep)):
        
        # init the k_array and insert -999's in the entire k_array
        k_array = np.ones((zbl.size, 5))
        k_array[0:, :] = vel_model[0,[3,1,2,4,5]]        
            
    else:
        # determining the layer number that corresponds to the depth interface
        sn_interf = np.array([])
        
        for intf in range(len(dep_interf)):
             sn_interf = np.append(sn_interf,sn[dep == dep_interf[intf]].max())

        # init the k_array
        k_array = np.ones((zbl.size, 5))

        layr_block_index = np.array([])
        for jj in range(zbl.size):

            dd = sn_interf[zbl[jj] > dep_interf]

            if dd.size == 0:
                layr_block_index = np.append(layr_block_index,0)
            else:
                layr_block_index = np.append(layr_block_index,dd.max())    

        for jj in range(zbl.size):
            k_array[jj, :] = vel_model[vel_model[:,0]==layr_block_index[jj],[3,1,2,4,5]]
       
    # calculating k_array gradient correction
    unq_counts = np.unique(layr_block_index,return_counts=True)[1]
    
    idx_ar = np.array([])
    for i in unq_counts:
        id_ar = np.arange(i)
        idx_ar = np.concatenate((idx_ar,id_ar*hvi_))
    
    # Update idx_ar if we are in the second or more blocks for gradient correction
    if blk > 0:
        
        if zbl[0] not in dep_interf:
            
            corr_ar = np.zeros(len(zbl))
            corr_ar[:unq_counts[0]] = zbl[0]
            
            idx_ar = corr_ar+idx_ar
            
            
    # calculate gradient array to correct the k_array
    grad_ar = np.zeros([len(layr_block_index),5])
    j=0
    
    for i in [3,1,2,4,5]:
        grad_ar[:,j] = vel_model_grad[layr_block_index.astype(int),i]
        j=j+1
        
    k_array_grad = grad_ar*idx_ar[:, np.newaxis]   
       
    # Update k_array
    k_array_upd = k_array+k_array_grad
       
    # printing outputs
    np.set_printoptions(suppress=True) 

    #print(sn_interf)
    #print(dep_interf)
    #print(layr_block_index)
    #print(k_array)
    
    
    
    return(k_array_upd)

def poly_latlon2idx(lons,lats,lon0,lat0,diff_lon,diff_lat,diff_m):
    
    lat0 = lat0
    lon0 = lon0

    diff_lon = diff_lon
    diff_lat = diff_lat


        
    delta_lon = lons - lon0
    delta_lat = lats - lat0

    dec_idx_x = delta_lon/diff_lon
    dec_idx_y = delta_lat/diff_lat
    
    x_m = dec_idx_x * diff_m
    y_m = dec_idx_y * diff_m

    idx_poly_points = np.zeros((len(dec_idx_x),2))
    for i in range(len(dec_idx_x)):
        idx_poly_points[i,0] = dec_idx_x[i]
        idx_poly_points[i,1] = dec_idx_y[i]
        
    m_poly_points = np.zeros((len(dec_idx_x),2))
    for i in range(len(dec_idx_x)):
        m_poly_points[i,0] = x_m[i]
        m_poly_points[i,1] = y_m[i]
   
    return(idx_poly_points,m_poly_points)
#%% Fixing the Vp/Vs ratio to be >1.5

min_vp_vs = 1.5
vp_vs_l , vp_vs_init = fix_vp_vs_ratio(min_vp_vs)

# plot_CVM_cs(block1_vp.vp.values,0,8000,'Vp','steph_z_cs_vp_top4')
# plot_CVM_cs(vp_vs_l,0,3,'Vp/Vs fixed','steph_z_cs_vp_vs_fixed_top4')
# plot_CVM_cs(vp_vs_init,0,3,'Vp/Vs','steph_z_cs_vp_vs_top4')
    
#%% Writing file header
rfilename = os.path.join(outpath, rfileoutput)

proj4 = ('+proj=tmerc +datum=NAD83 +units=m +lon_0=-124.19 +lat_0=43.4 +no_defs') ##

magic        = np.int32(1)
precision    = np.int32(4)
attenuation  = np.int32(1)
az           = np.float64(0.0)
lon0         = np.float64(-124.19) 
lat0         = np.float64(43.4) 
nb           = np.int32(len(block_layers)) 

header = (magic, precision, attenuation, az, lon0, lat0, proj4.encode(), nb)

print('\nrfile initialized at {}.\n'
      'Writing file header...\n'
      '======================'.format(rfilename))
flush
with open(rfilename, 'wb') as f:
    rfileIO.write_hdr(f, *header)
    
print('Done!')


#%% Writing block headers
print('Writing block headers')

# Writing block headers
for i in range(nb):
    
    if i == 0:

        # Block No.0 (topo)
        hh0 = np.float64(hhi[i])    # grid size in hori direction
        hv0 = np.float64(hvi[i])    # grid size in vertical direction
        z00 = np.float64(0)         # base z-level for block b.
        nc0 = np.int32(1)           # no of components in block b (in this case 1, just elevation)
        ni0 = np.int32(len(block1_vp.utmn))   # no of grid points along y axis
        nj0 = np.int32(len(block1_vp.utme))   # no of grid points along x axis
        nk0 = np.int32(1)           # no of grid points along the z axis (in this case 1)

        x_extent = (ni0 - 1) * hh0
        y_extent = (nj0 - 1) * hh0

        print ('\nrfile horizontal dimentions are '
               '{:.2f} km x {:.2f} km\n'
               u'with origin at {:.4f}°N, {:.4f}°E\n'
               'ni0: {}, nj0: {}'.format(x_extent * 1e-3, y_extent * 1e-3,
                                         lat0, lon0,
                                         ni0, nj0))

        block0_hdr = [hh0, hv0, z00, nc0, ni0, nj0, nk0]
        print(block0_hdr)

    else:
        # the data block
        if i == 1:
            zmin_blk = block_layers[i-1]
        else:
            zmin_blk = globals()['z'+str(i-1)].max() 
            
        zmax_blk = block_layers[i]

        globals()['z'+str(i)] = np.arange(zmin_blk, zmax_blk + hvi[i-1], hvi[i-1])

        hh_blk = np.float64(hhi[i-1])
        hv_blk = np.float64(hvi[i-1])
        z0_blk = np.float64(zmin_blk)
        nc_blk = np.int32(5)
        ni_blk = np.int32(len(globals()['block'+str(i)+'_vp'].utmn))
        nj_blk = np.int32(len(globals()['block'+str(i)+'_vp'].utme))
        nk_blk = np.int32(globals()['z'+str(i)].size)

        x_extent = (ni_blk - 1) * hh_blk
        y_extent = (nj_blk - 1) * hh_blk

        print ('\nBlock No.{} horizontal dimentions are '
                '{:.2f} km x {:.2f} km\n'
                'with a vertical extent of {} m high to {} m low.\n'
                'ni1: {}, nj1: {}, nk1: {}'.format(i, x_extent * 1e-3, y_extent * 1e-3,
                                          z0_blk, globals()['z'+str(i)].max(),
                                          ni_blk, nj_blk, nk_blk))

        globals()['block'+str(i)+'_hdr'] = [hh_blk, hv_blk, z0_blk, \
                      nc_blk, ni_blk, nj_blk, nk_blk]
        print(globals()['block'+str(i)+'_hdr'])

# Write block headers
print('\nWriting block headers...\n'
        '========================')
flush

with open(rfilename, 'ab') as f:
    for i in range(nb):
        rfileIO.write_block_hdr(f, *globals()['block'+str(i)+'_hdr'])

print('Done!')


#%% Write topo block data section
print('\nWriting topography data...\n'
        '==========================')
flush

with open(rfilename, 'ab') as f:
    (np.zeros([ni0,nj0]).astype(np.float32)).tofile(f)
    
print('Done!')

#%% Make basin polygon

inner1_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner1.csv')
lon_inner1 = inner1_f['inner1_lon']
lat_inner1 = inner1_f['inner1_lat']

inner2_f = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/basin_depth_poly_inner2.csv')
lon_inner2 = inner2_f['inner2_lon']
lat_inner2 = inner2_f['inner2_lat']

outer_f = pd.read_csv('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/wv_poly.csv')
lon_outer = outer_f['wv_poly_lons']
lat_outer = outer_f['wv_poly_lats'] 

# utme_inner1,utmn_inner1 = llz2utm(np.array(lon_inner1),np.array(lat_inner1),projection_zone=10)
# utme_inner2,utmn_inner2 = llz2utm(np.array(lon_inner2),np.array(lat_inner2),projection_zone=10)
# utme_outer,utmn_outer = llz2utm(np.array(lon_outer),np.array(lat_outer),projection_zone=10)

# inner1_poly_utm = np.vstack((utme_inner1, utmn_inner1)).T
# inner2_poly_utm = np.vstack((utme_inner2, utmn_inner2)).T
# outer_poly_utm  = np.vstack((utme_outer, utmn_outer)).T

# full_poly_utm = Polygon(outer_poly_utm, [inner1_poly_utm , inner2_poly_utm])

lat0 = 43.4
lon0 = -124.19

diff_200_lon = 0.0013
diff_200_lat = 0.0009

diff_m = 100

outer_poly_points_idx , outer_poly_points_m = poly_latlon2idx(lon_outer,lat_outer,lon0,lat0,diff_200_lon,diff_200_lat,diff_m)
inner1_poly_points_idx , inner1_poly_points_m = poly_latlon2idx(lon_inner1,lat_inner1,lon0,lat0,diff_200_lon,diff_200_lat,diff_m)
inner2_poly_points_idx , inner2_poly_points_m = poly_latlon2idx(lon_inner2,lat_inner2,lon0,lat0,diff_200_lon,diff_200_lat,diff_m)

full_poly_m = Polygon(outer_poly_points_m, [inner1_poly_points_m , inner2_poly_points_m])


#%%

# Creating tiffs from ifile
ifile2tiff(ifile_in)

#reading the geotiff files
for i in range(n_tif):
    globals()['layer%s' % str(i+1)] = geo.read_GeoTIFF('tiff/elev_layer%s.tif' % str(i+1))

# Resample grids for different block resolutions
dirs = resample_to_hhi(hhi,w, e, s, n)

# Correcting nans at the bottom and placing on top so the tiffs will fit the Stephenson model
# correct_nans_bottom_tiff()

# Reading the grids into a streams such as 'elev_layer1_hh=0.0125_stream' 
reading_grids()

#%% write block No.1 data section
for blk in range(nb-1):
    
    #blk = 1
    vp_blk = globals()['block'+str(blk+1)+'_vp']
    vs_blk = globals()['block'+str(blk+1)+'_vs']
    
    z_blk = globals()['z'+str(blk+1)]
    ni_blk = np.int32(len(vp_blk.utmn))
    nj_blk = np.int32(len(vp_blk.utme))
    hhi_blk = hhi[blk]
    hvi_blk = hvi[blk]
    
    # init the k_array
    k_array = np.ones((z_blk.size, 5))
    
    print('\nWriting data section of Block No.'+str(blk+1)+'\n'
            '==================================')
    print('z_blk =')
    print(z_blk)
    
    
    # Add repeat layers
    if blk == 1 or blk == 2:
        # Open earlier files
        vp_blk_prev = globals()['block'+str(blk)+'_vp']
        vs_blk_prev = globals()['block'+str(blk)+'_vs']

        vp_prev = vp_blk_prev.vp.values
        vs_prev = vs_blk_prev.vs.values

        # Interpolate last layer of previous block onto current spacing
        xi_new, yi_new = np.meshgrid(vp_blk.utme.values, vp_blk.utmn.values)
        xi_prev, yi_prev = np.meshgrid(vp_blk_prev.utme.values, vp_blk_prev.utmn.values)
        
        vp_prev_base = griddata((xi_prev.T.flatten(), yi_prev.T.flatten()), vp_prev[:,:,-1].flatten(), (xi_new, yi_new), method='linear')
        vs_prev_base = griddata((xi_prev.T.flatten(), yi_prev.T.flatten()), vs_prev[:,:,-1].flatten(), (xi_new, yi_new), method='linear')

        if np.isnan(vs_prev_base[0,0]):
            vp_prev_base[0,:] = vp_prev_base[1,:]
        if np.isnan(vs_prev_base[0,0]):
            vs_prev_base[0,:] = vs_prev_base[1,:]
        
        # Add additional layer at the top
        vp_corr = np.dstack((vp_prev_base.T, vp_blk.vp.values))
        vs_corr = np.dstack((vs_prev_base.T, vs_blk.vs.values))

        vp_blk_prev.close()
        vs_blk_prev.close()
        
    lat_in = []
    lon_in = []
    with open(rfilename, 'ab') as f:
        for i in range(ni_blk):
            for j in range(nj_blk):
                
                # Determing k_array at location index [i,j]  #rho,vp,vs,qp,qs
                if blk == 0:
                    vp = vp_blk.vp.values[j,i,:]
                    vs = vs_blk.vs.values[j,i,:]
                else:
                    vp = vp_corr[j,i,:]
                    vs = vs_corr[j,i,:]
                    
                rho = calcRho(vp)
                qp, qs = calcQ(vs)
                
                # remove layers with 1e20 values
                index = np.where(vp>1e7)
                vp[index] = -999
                vs[index] = -999
                rho[index] = -999
                qp[index] = -999
                qs[index] = -999
                
                
                k_array[:,0]= rho
                k_array[:,1]= vp
                k_array[:,2]= vs
                k_array[:,3]= qp
                k_array[:,4]= qs
                
                if k_array.min() < 50 and k_array.min() > -999:
                    print(i,j,k_array.min())
            
                # Check if ni,nj is in the basin
                grd_pnt = Point(j*hhi_blk,i*hhi_blk)
                if full_poly_m.contains(grd_pnt) == True:
                    lat_in.append(i)
                    lon_in.append(j)                
                
                # grd_pnt = Point(vp_blk.utme.values[j],vp_blk.utmn.values[i])
                # if full_poly_utm.contains(grd_pnt) == True:
                #     lat_in.append(vp_blk.utmn.values[i])
                #     lon_in.append(vp_blk.utme.values[j])
                    
                    # Determing k_array at location index [i,j]
                    k_array_basin = create_k_array_grad(i,j,hhi_blk,n_tif,z_blk,vel_model,vel_model_grad,hvi_blk,blk)
                    
                    # delete basement velocity values
                    k_array_basin = np.delete(k_array_basin,(k_array_basin[:,0]>1e7) | (k_array_basin[:,0]<0),0)    

                    # Insert k_array_basin to k_array
                    k_array[0:len(k_array_basin),:] =k_array_basin[:,:]

                # write k_array to file
                k_array.astype(np.float32).tofile(f)
                if k_array.size != (len(z_blk)*5):
                    print(k_array.size)
                    
                
                    
                if j == 0 and i%100 == 0:
                    print('\rWriting row {} of {}, '
                          'column {} of {} of block 1...'.format(i + 1, ni_blk,
                                                                  j + 1, nj_blk)),
                flush
    
    print('\nDone!\n')

# ####################################################################
end = time.time()
time_elaps = end - start
if time_elaps < 60:
    print(f'Duration: {round(time_elaps)} seconds')
else:
    print(f'Duration: {round(time_elaps/60)} minutes')

#%%

plt.figure(figsize=(10,16))
plt.scatter(lon_in,lat_in)
plt.show()

#%% Reading and plotting the rfile
model = rfileIO.read(rfileoutput, 'all', verbose=True)

#%% ---------------------------------------------------------------------------
model.plot_topography(cmap='terrain')
figpath = os.getcwd() +'/fig.rfiletopo.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()

#%%

cs = model.get_cross_section(0, 370, 100, 100) # lat1,lat2,lon1,lon2 in km
fig, ax, cb = cs.plot(property='vp',aspect=10, cmap='jet') #vmin=0, 
figpath = os.getcwd() +'/fig.crossection2.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()


#%%
# ---------------------------------------------------------------------------
# Vs profile at different locations (x,y) within the domain
lat_bin = np.linspace(100.0, 300, 37)
lon_bin = np.linspace(100.0, 200, 21)

fig, ax = plt.subplots()

for i in range(len(lat_bin)):
    for j in range(len(lon_bin)):
        z, properties = model.get_z_profile(lon_bin[j], lat_bin[i])

        p = np.ma.masked_equal(properties.T[2], -999)
        
    
        ax.plot(p, z,linewidth=0.5) 

ax.invert_yaxis()
ax.xaxis.tick_top()
ax.set_xlabel('Vs (m/s)', fontsize=20)
ax.xaxis.set_label_position('top')
ax.set_ylabel('Z, km', fontsize=20)
ax.tick_params(axis='x',labelsize=15)
ax.tick_params(axis='y',labelsize=15)

figpath = os.getcwd() +'/fig.profile_vs.png'
plt.savefig(figpath, bbox_inches='tight', dpi=500) 

#%%
cs = model.get_cross_section(100, 100, 0, 215) # lat1,lat2,lon1,lon2 in km
fig, ax, cb = cs.plot(property='vs',aspect=1, cmap='jet') #vmin=0, 
figpath = os.getcwd() +'/fig.crossection.png'
plt.savefig(figpath, bbox_inches='tight', dpi=500) 
plt.show()

#%%
cs = model.get_cross_section(0, 375, 100, 100) # lat1,lat2,lon1,lon2 in km
fig, ax, cb = cs.plot(property='vs',aspect=2, cmap='jet') #vmin=0, 
figpath = os.getcwd() +'/fig.crossection2.png'
plt.savefig(figpath, bbox_inches='tight', dpi=500) 
plt.show()

















