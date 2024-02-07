#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 14:45:47 2023

@author: rshimony
"""

import numpy as np
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

wk_dir = '/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin'    
os.chdir(wk_dir)

#%% Inputs and options
start = time.time()

hhi = (200, 300, 900)
hvi = (100,  300,  900)
block_layers = (0, 1200, 9900, 59400) 


# Creating an array for the velocity model
ifile_in = '/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/ifile_gp1_deep_singmat_4rfile.txt' #'ALL_japan_path20111110.dat'
n_tif = 2
dx = 100 #(m)
dy = -100 #(m)

vel_model = np.array([\
[0,   -999,  -999,   -999, -999,-999],
[1,  1.7e3, 0.35e3, 1.80e3, 119,  70],
[2,  1e20, 1e20, 1e20, 1e20,  1e20]])

# Extent of the resulting rfile (m) [w,s=0; e=x_dim; n=y_dim] 
w = 0;  e = 216e3 #215e3
s = 0;  n = 378e3 #375e3


#utme: 1080, utmn: 1890
# utm *2 = 216, 378

# i = 1; b = 0; globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].vele
# Out[30]: 
# array([[ -0.,  -0.,  -0., ...,  -0.,  -0.,  -0.],
#        [ nan,  nan,  nan, ...,  nan,  nan,  nan],
#        [ nan,  nan,  nan, ...,  nan,  nan,  nan],
#        ..., 
#        [ -0.,  -0.,  -0., ...,  nan,  nan,  nan],
#        [ -0.,  -0.,  -0., ...,  nan,  nan,  nan],
#        [ -0.,  -0.,  -0., ...,  nan,  nan,  nan]], dtype=float32)

# i = 1; b = 1; globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].vele
# Out[31]: 
# array([[ nan,  nan,  nan, ...,  nan,  nan,  nan],
#        [ nan,  nan,  nan, ...,  nan,  nan,  nan],
#        [ nan,  nan,  nan, ...,  nan,  nan,  nan],
#        ..., 
#        [ -0.,  -0.,  -0., ...,  nan,  nan,  nan],
#        [ -0.,  -0.,  -0., ...,  nan,  nan,  nan],
#        [ -0.,  -0.,  -0., ...,  nan,  nan,  nan]], dtype=float32)

# plt.imshow(globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].vele)
# Out[32]: <matplotlib.image.AxesImage at 0x7faeadd69bb0>

#%% Helper functions    

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

def create_k_array(i,j,hhi_,n_layers,zbl,vel_model):

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
       
    # printing outputs
    np.set_printoptions(suppress=True) 

    #print(sn_interf)
    #print(dep_interf)
    #print(layr_block_index)
    #print(k_array)
    
    return(k_array)


#%%
# Preprocess the tiff files or the ifiles

# Creating tiffs from ifile
ifile2tiff(ifile_in)

#reading the geotiff files
for i in range(n_tif):
    globals()['layer%s' % str(i+1)] = geo.read_GeoTIFF('tiff/elev_layer%s.tif' % str(i+1))

# Resample grids for different block resolutions
dirs = resample_to_hhi(hhi,w, e, s, n)

# Reading the grids into a streams such as 'elev_layer1_hh=0.0125_stream' 
reading_grids()

i = 2
globals()['z'+str(i)] = np.arange(1200, 110000, hvi[i-1])

# write block No.1 data section
for blk in [1]:#range(2):
    
    z_blk = globals()['z'+str(blk+1)]
    ni_blk = np.int32(globals()['elev_layer1_hh='+str(hhi[blk])+'_stream'].ny)
    nj_blk = np.int32(globals()['elev_layer1_hh='+str(hhi[blk])+'_stream'].nx)
    hhi_blk = hhi[blk]
    
    # init the k_array
    k_array = np.ones((z_blk.size, 5))
    
    print('\nWriting data section of Block No.'+str(blk+1)+'\n'
            '==================================')
    print('z_blk =')
    print(z_blk)
    
   
    for i in [700]:#range(ni_blk):
        for j in [350]:#range(nj_blk):
            
            # Determing k_array at location index [i,j]
            k_array = create_k_array(i,j,hhi_blk,n_tif,z_blk,vel_model)
            
            # delete basement velocity values
            k_array = np.delete(k_array,(k_array[:,0]>1e7) | (k_array[:,0]<0),0)










# ####################################################################
end = time.time()
time_elaps = end - start
if time_elaps < 60:
    print(f'Duration: {round(time_elaps)} seconds')
else:
    print(f'Duration: {round(time_elaps/60)} minutes')
