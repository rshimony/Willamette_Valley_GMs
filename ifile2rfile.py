#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 15:46:25 2022

@author: oluwaseunfadugba
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

#%% Inputs and options
start = time.time()

# hhi = (hh, 2 * hh, 4 * hh, 8 * hh)
# hvi = (hv, 2 * hv, 4 * hv, 8 * hv)
hhi = (1000, 2000, 3000, 4000)
hvi = (200,  300,  500, 1000)
block_layers = (-4000, 30000, 70000, 100000, 200000) ##block_layers = (-4000, 30000, 70000, 100000, 200000)

plot_resampled_3D = 1

outpath = '/Users/oluwaseunfadugba/Documents/Projects/TsE_ValerieDiego/TsE_1D_vs_3D/'+\
    '3D_Modeling_using_SW4/2_Creating_Rfile/Creating_rfile_test8b_SF_4/'


plot_orig_tiff = 1
reproject = 0; epsg = 32636

# Extent of the resulting rfile (m) [w,s=0; e=x_dim; n=y_dim] 
w = 0;  e = 1441e3
s = 0;  n = 2040e3

rfileoutput = '3djapan_hh=1000m_hv=2_3_5_1km.rfile'

# Creating an array for the velocity model
ifile_in = 'ALL_japan_path20111110.dat'
n_tif = 23
dx = 1000 #(m)
dy = -1000 #(m)

vel_model = np.array([\
[0,   -999,  -999,   -999, -999,-999],
[1,  1.7e3, 0.35e3, 1.80e3, 119,  70],
[2,  1.8e3, 0.5e3,  1.95e3, 170, 100],
[3,  2.0e3, 0.6e3,  2.00e3, 204, 120],
[4,  2.1e3, 0.7e3,  2.05e3, 238, 140],
[5,  2.2e3, 0.8e3,  2.07e3, 272, 160],
[6,  2.3e3, 0.9e3,  2.10e3, 306, 180],
[7,  2.4e3, 1.0e3,  2.15e3, 340, 200],
[8,  2.7e3, 1.3e3,  2.20e3, 442, 260],
[9,  3.0e3, 1.5e3,  2.25e3, 510, 300],
[10, 3.2e3, 1.7e3,  2.30e3, 578, 340],
[11, 3.5e3, 2.0e3,  2.35e3, 680, 400],
[12, 4.2e3, 2.4e3,  2.45e3, 680, 400],
[13, 5.0e3, 2.9e3,  2.60e3, 680, 400],
[14, 5.5e3, 3.2e3,  2.65e3, 680, 400],
[15, 5.8e3, 3.4e3,  2.70e3, 680, 400],
[16, 6.4e3, 3.8e3,  2.80e3, 680, 400],
[17, 7.5e3, 4.5e3,  3.20e3, 850, 500],
[18, 5.0e3, 2.9e3,  2.40e3, 340, 200],
[19, 6.8e3, 4.0e3,  2.90e3, 510, 300],
[20, 8.0e3, 4.7e3,  3.20e3, 850, 500],
[21, 5.4e3, 2.8e3,  2.60e3, 340, 200],
[22, 6.5e3, 3.5e3,  2.80e3, 510, 300],
[23, 8.1e3, 4.6e3,  3.40e3, 850, 500]])

# printing outputs
np.set_printoptions(suppress=True) 
    
print(vel_model)

os.chdir(outpath)

#%% Checking input parameters
if max(np.remainder(np.abs(w-e),np.abs(hhi)))>0 or \
    max(np.remainder(np.abs(s-n),np.abs(hhi)))>0:

    sys.exit('The domain must be a multiple of the grid spacings')

#%%
def roundup(x):
      return int(math.ceil(x / 1000.0)) * 1000

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
        
def resampl_newextent(stream,w, e, s, n,resam):

    set_new_extent(stream,w, e, s, n)   
    stream.resample(to=resam)   
    print(stream.dx,stream.dy)
    fig = plt.figure(); plt.imshow(stream.elev, cmap='terrain', \
                                   extent=stream.extent) 
    plt.close()
    
def reproj(stream,epsg):
    
    stream.reproject(epsg=32636); 
    stream.name = 'kk'
    fig = plt.figure(); plt.imshow(stream.elev, cmap='terrain', \
                                   extent=stream.extent) 
    plt.close()
           
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

def ifile2tiff(ifile_in):
    print('\nConverting ifile (%s) to tiffs'  %ifile_in)

    data = np.genfromtxt(ifile_in) 
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
        globals()['elev_layer%s' % j] = -globals()['elev_layer_tmp%s' % j]#[::-1]
        
        # Creating TIFF for each layer
        geo.save_GeoTIFF('elev_layer'+str(j)+'.tif', globals()['elev_layer%s' % j], 
                          tlx, tly, dx, dy, epsg=None, proj4=proj4, 
                          dtype=np.float32, nodata=np.nan, rasterBand=1)

        os.system('mkdir -p tiff');
        os.system('mv elev_layer*.tif tiff');
    
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

#%%
min_topo = roundup(-np.nanmin(globals()['elev_layer1_hh='+str(hhi[0])+'_stream'].elev))
max_topo = -roundup(np.nanmax(globals()['elev_layer1_hh='+str(hhi[0])+'_stream'].elev))

if max_topo < block_layers[0]:
    print('Changing the highest topography to', str(max_topo),'m')
    block_layers[0] = max_topo
    
if min_topo > block_layers[1]:
    print('Changing the second block base to the lowest topography (', str(min_topo),'m + 1000 m)')
    block_layers[1] = min_topo + 1000

#%%
rfilename = os.path.join(outpath, rfileoutput)

proj4 = ('+proj=tmerc +datum=NAD83 +units=m +lon_0=129 +lat_0=30 +no_defs') ##

magic        = np.int32(1)
precision    = np.int32(4)
attenuation  = np.int32(1)
az           = np.float64(0.0)
lon0         = np.float64(129) 
lat0         = np.float64(30) 
nb           = np.int32(len(block_layers)) 

header = (magic, precision, attenuation, az, lon0, lat0, proj4.encode(), nb)

print('\nrfile initialized at {}.\n'
      'Writing file header...\n'
      '======================'.format(rfilename))
flush
with open(rfilename, 'wb') as f:
    rfileIO.write_hdr(f, *header)
    
print('Done!')

#%%
print('Writing block headers')

# Writing block headers
for i in range(nb):
    
    if i == 0:

        # Block No.0 (topo)
        hh0 = np.float64(hhi[i])    # grid size in hori direction
        hv0 = np.float64(hvi[i])    # grid size in vertical direction
        z00 = np.float64(0)         # base z-level for block b.
        nc0 = np.int32(1)           # no of components in block b (in this case 1, just elevation)
        ni0 = np.int32(globals()['elev_layer1_hh='+str(hhi[i])+'_stream'].ny)   # no of grid points along y axis
        nj0 = np.int32(globals()['elev_layer1_hh='+str(hhi[i])+'_stream'].nx)   # no of grid points along x axis
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
        ni_blk = np.int32(globals()['elev_layer1_hh='+str(hhi[i-1])+'_stream'].ny)
        nj_blk = np.int32(globals()['elev_layer1_hh='+str(hhi[i-1])+'_stream'].nx)
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

#%%
# Write topo block data section
print('\nWriting topography data...\n'
        '==========================')
flush

with open(rfilename, 'ab') as f:
    (-1 * globals()['elev_layer1_hh='+str(hhi[0])+'_stream'].\
     vele.astype(np.float32)).tofile(f)
    
print('Done!')


#%%
# write block No.1 data section
for blk in range(nb-1):
    
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
    
    with open(rfilename, 'ab') as f:
        for i in range(ni_blk):
            for j in range(nj_blk):
                
                # Determing k_array at location index [i,j]
                k_array = create_k_array(i,j,hhi_blk,n_tif,z_blk,vel_model)
            
                if k_array.min() < 50 and k_array.min() > -999:
                    print(i,j,k_array.min())
            
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
# Reading and plotting the rfile
model = rfileIO.read(rfileoutput, 'all', verbose=True)

#%%
# ---------------------------------------------------------------------------
model.plot_topography(cmap='terrain')
figpath = os.getcwd() +'/fig.rfiletopo.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()

#%%
# ---------------------------------------------------------------------------
cs = model.get_cross_section(0, 2040, 400, 400)
fig, ax, cb = cs.plot(property='vp', vmin=1800,aspect=10, cmap='jet') #vmin=0,
figpath = os.getcwd() +'/fig.crossection1.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()

#%%
cs = model.get_cross_section(0, 2039, 1000, 1000)
fig, ax, cb = cs.plot(property='vp', vmin=1800,aspect=10, cmap='jet') #vmin=0, 
figpath = os.getcwd() +'/fig.crossection2.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()

cs = model.get_cross_section(1500, 1500, 0, 1400)
fig, ax, cb = cs.plot(property='vp', vmin=1800, aspect=5, cmap='jet')
figpath = os.getcwd() +'/fig.crossection3.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()

cs = model.get_cross_section(1000, 1000, 0, 1400)
fig, ax, cb = cs.plot(property='vp', vmin=1800, aspect=5, cmap='jet')
figpath = os.getcwd() +'/fig.crossection4.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()

cs = model.get_cross_section(100, 100, 0, 1400)
fig, ax, cb = cs.plot(property='vp', vmin=1800, aspect=5, cmap='jet')
figpath = os.getcwd() +'/fig.crossection5.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 
plt.show()

#%%
# ---------------------------------------------------------------------------
# z, v = model.get_z_profile(10, 10)
# print(z.shape, v.shape)
# print(z[0], v[1])

# plt.plot(v[:,1], z)

z, properties = model.get_z_profile(10, 10)

labels = ('rho, kg/m^3', 'vp, m/s', 'vs, m/s', 'qp', 'qs')

fig, ax = plt.subplots()

for p, label in zip(properties.T, labels):
    p = np.ma.masked_equal(p, -999)
    ax.plot(p, z, label=label)
ax.invert_yaxis()
ax.xaxis.tick_top()
ax.set_xlabel('Material properties at (10,10)')
ax.xaxis.set_label_position('top')
ax.set_ylabel('Z, km')
ax.legend(loc=10)

figpath = os.getcwd() +'/fig.profile_loc_10_10.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 

#%%
z, properties = model.get_z_profile(1800, 1400)

labels = ('rho, kg/m^3', 'vp, m/s', 'vs, m/s', 'qp', 'qs')

fig, ax = plt.subplots()

for p, label in zip(properties.T, labels):
    p = np.ma.masked_equal(p, -999)
    ax.plot(p, z, label=label)
ax.invert_yaxis()
ax.xaxis.tick_top()
ax.set_xlabel('Material properties at (1849,1440)')
ax.xaxis.set_label_position('top')
ax.set_ylabel('Z, km')
ax.legend(loc=10)

figpath = os.getcwd() +'/fig.profile_loc_1800_1400.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 


 
#%%
# ---------------------------------------------------------------------------
# Vs profile at different locations (x,y) within the domain
lat_bin = np.linspace(0.0, 2040, 601)
lon_bin = np.linspace(0.0, 1400, 401)

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
plt.savefig(figpath, bbox_inches='tight', dpi=200) 

# #%%
# # Vp profile at different locations (x,y) within the domain
# lat_bin = np.linspace(0.0, 2000, 601)
# lon_bin = np.linspace(0.0, 1400, 401)

# fig, ax = plt.subplots()

# for i in range(len(lat_bin)):
#     for j in range(len(lon_bin)):
#         z, properties = model.get_z_profile(lon_bin[j], lat_bin[i])

#         p = np.ma.masked_equal(properties.T[1], -999)
#         ax.plot(p, z,linewidth=0.5) 

# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.set_xlabel('Vp (m/s)', fontsize=20)
# ax.xaxis.set_label_position('top')
# ax.set_ylabel('Z, km', fontsize=20)
# ax.tick_params(axis='x',labelsize=15)
# ax.tick_params(axis='y',labelsize=15)

# figpath = os.getcwd() +'/fig.profile_vp.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 

# #%%
# # Rho profile at different locations (x,y) within the domain
# lat_bin = np.linspace(0.0, 2000, 601)
# lon_bin = np.linspace(0.0, 1400, 401)

# fig, ax = plt.subplots()

# for i in range(len(lat_bin)):
#     for j in range(len(lon_bin)):
#         z, properties = model.get_z_profile(lon_bin[j], lat_bin[i])

#         p = np.ma.masked_equal(properties.T[0], -999)
#         ax.plot(p, z,linewidth=0.5) 

# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.set_xlabel('rho (kg/m^3)', fontsize=20)
# ax.xaxis.set_label_position('top')
# ax.set_ylabel('Z, km', fontsize=20)
# ax.tick_params(axis='x',labelsize=15)
# ax.tick_params(axis='y',labelsize=15)

# figpath = os.getcwd() +'/fig.profile_rho.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 

# #%%
# # ---------------------------------------------------------------------------
# # Qp profile at different locations (x,y) within the domain
# lat_bin = np.linspace(0.0, 2000, 601)
# lon_bin = np.linspace(0.0, 1400, 401)

# fig, ax = plt.subplots()

# for i in range(len(lat_bin)):
#     for j in range(len(lon_bin)):
#         z, properties = model.get_z_profile(lon_bin[j], lat_bin[i])

#         p = np.ma.masked_equal(properties.T[3], -999)
#         ax.plot(p, z,linewidth=0.5) 

# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.set_xlabel('Qp', fontsize=20)
# ax.xaxis.set_label_position('top')
# ax.set_ylabel('Z, km', fontsize=20)
# ax.tick_params(axis='x',labelsize=15)
# ax.tick_params(axis='y',labelsize=15)

# figpath = os.getcwd() +'/fig.profile_qp.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 

# #%%
# # ---------------------------------------------------------------------------
# # Qs profile at different locations (x,y) within the domain
# lat_bin = np.linspace(0.0, 2040, 510)
# lon_bin = np.linspace(0.0, 1400, 360)

# fig, ax = plt.subplots()

# for i in range(len(lon_bin)):
#     for j in range(len(lon_bin)):
#         z, properties = model.get_z_profile(lon_bin[j], lat_bin[i])

#         p = np.ma.masked_equal(properties.T[4], -999)
#         ax.plot(p, z,linewidth=0.5) 

# ax.invert_yaxis()
# ax.xaxis.tick_top()
# ax.set_xlabel('Qs', fontsize=20)
# ax.xaxis.set_label_position('top')
# ax.set_ylabel('Z, km', fontsize=20)
# ax.tick_params(axis='x',labelsize=15)
# ax.tick_params(axis='y',labelsize=15)

# figpath = os.getcwd() +'/fig.profile_qs.png'
# plt.savefig(figpath, bbox_inches='tight', dpi=200) 


# #%%
# if plot_orig_tiff == 1:
    
#     fig = plt.figure()
#     fig.set_size_inches(100, 90)

#     NN = math.ceil(np.sqrt(n_tif))

#     fig.suptitle('Topography of each layer of the 3D velocity model', fontsize=160)


#     for i in range(n_tif):
#         ax = fig.add_subplot(NN, NN, i+1)
#         ax.set_title('Layer '+ str(i+1), fontsize=80)
        
#         plt.imshow(globals()['layer%s' % str(i+1)].elev, cmap='terrain',\
#                    extent=globals()['layer%s' % str(i+1)].extent)
#         plt.colorbar()
#         plt.rcParams.update({'font.size': 50})
    
#     figpath = os.getcwd() +'/fig.topo_for_each_layer.png'
#     plt.savefig(figpath, bbox_inches='tight', dpi=200)

#     plt.show()
    
# #%%
# if plot_resampled_3D == 1:
#     print('Plotting the topography of each resampled version of the 3D velocity model...')
#     # Plotting the topography of each resampled version of the 3D velocity model
#     for b in range(len(hhi)):

#         fig = plt.figure()
#         fig.set_size_inches(100, 90)

#         NN = math.ceil(np.sqrt(n_tif))

#         fig.suptitle('Topography of each resampled version of the 3D velocity model ('+str(hhi[b])+'m)', fontsize=160)

#         for i in range(n_tif):

#             ax = fig.add_subplot(NN, NN, i+1)
#             ax.set_title('Layer '+ str(i+1), fontsize=70)

#             plt.imshow(globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].elev, cmap='terrain', \
#                        extent=globals()['elev_layer' + str(i + 1)+'_hh='+str(hhi[b])+'_stream'].extent)
#             plt.colorbar()
#             plt.rcParams.update({'font.size': 50})

#         figpath = os.getcwd() +'/fig.topo_for_each_layer_in_3Dvel_resampled_'+str(hhi[b])+'m.png'
#         plt.savefig(figpath, bbox_inches='tight', dpi=200)

#         plt.show()
#         plt.close()
        
# #%%
# # Checking the lowest possible average vs for a layer 200 meters thick at the surface
# def Average(lst):
#     return sum(lst) / len(lst)

# # ---------------------------------------------------------------------------
# # Vs profile at different locations (x,y) within the domain
# lat_bin = np.linspace(0.0, 2040, 601)
# lon_bin = np.linspace(0.0, 1400, 401)

# #fig, ax = plt.subplots()
# vs_0_2 = []
# vs_2_4 = []
# vs_4_6 = []
# vs_6_8 = []
# vs_8_10 = []

# for i in range(len(lon_bin)):
#     for j in range(len(lon_bin)):
#         z, properties = model.get_z_profile(lon_bin[j], lat_bin[i])

#         vs =  properties.T[2]
        
#         vs = vs[vs > 0]
#         vs_0_2.append(vs[1])
#         vs_2_4.append(vs[2])
#         vs_4_6.append(vs[3])
#         vs_6_8.append(vs[4])
#         vs_8_10.append(vs[5])
        
# print("Average of vs 0-200m =", round(Average(vs_0_2), 2))
# print("Average of vs 200-400m =", round(Average(vs_2_4), 2))
# print("Average of vs 400-600m =", round(Average(vs_4_6), 2))
# print("Average of vs 600-800m =", round(Average(vs_6_8), 2))
# print("Average of vs 800-1000m =", round(Average(vs_8_10), 2))

#         #p = np.ma.masked_equal(properties.T[2], -999)
#         #ax.plot(p, z,linewidth=0.5) 





# # #%%
# # ####################################################################
# end = time.time()
# time_elaps = end - start
# if time_elaps < 60:
#     print(f'Duration: {round(time_elaps)} seconds')
# else:
#     print(f'Duration: {round(time_elaps/60)} minutes')
    
    
    
#
#%% Reading and plotting the rfile (single-block rfile)
blk1_rfile = '/Users/oluwaseunfadugba/Documents/Projects/'+\
    'TsE_ValerieDiego/TsE_1D_vs_3D/3D_Modeling_using_SW4/'+\
        '2_Creating_Rfile/Creating_rfile_test5_oldworking/'+\
            '3djapan_hh=1000m_hv=200m_blk1.rfile'
model_blk1 = rfileIO.read(blk1_rfile, 'all', verbose=True)

z, properties = model_blk1.get_z_profile(2041,1440)

labels = ('rho, kg/m^3', 'vp, m/s', 'vs, m/s', 'qp', 'qs')

fig, ax = plt.subplots()

for p, label in zip(properties.T, labels):
    p = np.ma.masked_equal(p, -999)
    ax.plot(p, z, label=label)
ax.invert_yaxis()
ax.xaxis.tick_top()
ax.set_xlabel('Material properties for rfile with 1 block at (2041,1440)')
ax.xaxis.set_label_position('top')
ax.set_ylabel('Z, km')
ax.legend(loc=10)

figpath = os.getcwd() +'/fig.profile_rfile_w_blk1_loc_2041_1440.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200) 


z, properties = model_blk1.get_z_profile(1848,1440)

labels = ('rho, kg/m^3', 'vp, m/s', 'vs, m/s', 'qp', 'qs')

fig, ax = plt.subplots()

for p, label in zip(properties.T, labels):
    p = np.ma.masked_equal(p, -999)
    ax.plot(p, z, label=label)
ax.invert_yaxis()
ax.xaxis.tick_top()
ax.set_xlabel('Material properties for rfile with 1 block at (1848,1440)')
ax.xaxis.set_label_position('top')
ax.set_ylabel('Z, km')
ax.legend(loc=10)

figpath = os.getcwd() +'/fig.profile_rfile_w_blk1_loc_1848_1440.png'
plt.savefig(figpath, bbox_inches='tight', dpi=200)    