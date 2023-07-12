#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 13:06:23 2023

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry.polygon import Polygon

#%%

def get_qs(vs):
    """
    Calculate Shear-wave Quality Factor based on Brocher (2008).

    .. note:: If Shear-wave velocity is less-than 0.3 km/s, Shear-wave
              Quality Factor is set to 13.

    Parameters
    ----------
    vs : float or sequence
        Shear-wave velocity in km/s.

    Returns
    -------
    float or sequence
        Shear-wave Quality Factor.
    """
    qs = (-16
          + 104.13 * vs
          - 25.225 * vs**2
          + 8.2184 * vs**3)
    try:
        qs[vs < 0.3] = 13
    except TypeError:
        if vs < 0.3:
            qs = 13

    return qs


def get_qp(qs):
    """Calculate Pressure-wave Quality Factor based on Brocher (2008).

    Parameters
    ----------
    qs : float or sequence
        Shear-wave Quality Factor.

    Returns
    -------
    float or sequence
        Pressure-wave Quality Factor.
    """
    return 2 * qs

#%%

min_vss = [400,800]
min_vps = [1665,2218]
min_rhos = [1735,1996]

vs_grads = [0.633,0.487]
vp_grads = [0.79,0.614]
rho_grads = [0.303,0.119]

min_depth = [0,0]
max_depth = [4553,5500]

min_vs = min_vss[0]
vs_grad = vs_grads[0]

min_vp = min_vps[0]
vp_grad = vp_grads[0]

min_rho = min_rhos[0]
rho_grad = rho_grads[0]

depths = np.linspace(min_depth[0],max_depth[0],10)

vs_depth = (depths*vs_grad) + min_vs 
vp_depth = (depths*vp_grad) + min_vp 
rho_depth = (depths*rho_grad) + min_rho 

qs_depth = get_qs(vs_depth/1000)
qp_depth = get_qp(qs_depth)

#%%

infile = '/Users/rshimony/Desktop/WillametteValley/Models/dist2depth_s_grav_n/salem/infile_dist2depth_s_grav_n_salem_layers_h50.txt'

with open(infile, 'a') as fout:
    
    fout.write('\n')
    fout.write('##materials##\n')

    for i in range(len(vs_depth)):
        vs = vs_depth[i]
        vp = vp_depth[i]
        rho = rho_depth[i]
        qs = qs_depth[i]
        qp = qp_depth[i]
        ids = i+1
    
        fout.write('material id=%.2f vp=%.2f vs=%.2f rho=%.2f qp=%.2f qs=%.2f\n' % (ids,vp,vs,rho,qp,qs))
    fout.close()

#%%

ifile_top = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/ifiles/cal_depth_ifile_dist_grav_top.txt')

ifile_lons = []
ifile_lats = []
ifile_dists_top = []

for i in range(1,len(ifile_top)):
    ifile_lon = ifile_top[i][0]
    ifile_lat = ifile_top[i][1]
    ifile_dist_top = ifile_top[i][2]
    
    ifile_lons.append(ifile_lon)
    ifile_lats.append(ifile_lat)
    ifile_dists_top.append(ifile_dist_top)
    
plt.figure(figsize=(8,10))
plt.scatter(ifile_lons,ifile_lats,c=ifile_dists_top,cmap='gist_earth_r')
plt.colorbar()
plt.show()

#%%

ifile_bot = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/cal_depth_ifile_dist_grav_bot.txt')

ifile_dists_bot = []

for i in range(1,len(ifile_bot)):
    ifile_dist_bot = ifile_bot[i][2]
    
    ifile_dists_bot.append(ifile_dist_bot)

#%%

layers_coef = np.linspace(0.1,0.9,9)
ifile_dists_top_arr = np.array(ifile_dists_top)

layers = []
for i in range(len(layers_coef)):
    layer_i = ifile_dists_top_arr * layers_coef[i]
    layers.append(layer_i)


#%%

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))
nmat = len(layers)+2

with open('cal_depth_ifile_dist_layers.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,nmat)
    f.write(header_line)
    for i in range(len(ifile_dists_top)):
        
        layer_1 = layers[0][i]
        layer_2 = layers[1][i]
        layer_3 = layers[2][i]
        layer_4 = layers[3][i]
        layer_5 = layers[4][i]
        layer_6 = layers[5][i]
        layer_7 = layers[6][i]
        layer_8 = layers[7][i]
        layer_9 = layers[8][i]
        layer_10 = ifile_dists_top[i]
        layer_11 = ifile_dists_bot[i]
        
        if ifile_dists_top[i] > ifile_dists_bot[i]:
            layer_11 = layer_10
        

        f.write("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n" % (ifile_lons[i],ifile_lats[i],layer_1,layer_2,layer_3,layer_4,layer_5,layer_6,layer_7,
                                                              layer_8,layer_9,layer_10,layer_11))


#%%

ifile_dists_file = np.genfromtxt('/Users/rshimony/Desktop/WillametteValley/Models/cal_depth_ifile_dist.txt')

ifile_dists_lons = []
ifile_dists_lats = []
ifile_dists = []

for i in range(1,len(ifile_dists_file)):

    ifile_dist = ifile_dists_file[i][2]
    ifile_dist_lons = ifile_dists_file[i][0]
    ifile_dist_lats = ifile_dists_file[i][1]
    
    ifile_dists.append(ifile_dist)
    ifile_dists_lons.append(ifile_dist_lons)
    ifile_dists_lats.append(ifile_dist_lats)
    
#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_dists_lons,ifile_dists_lats,c=ifile_dists,cmap='gist_earth_r')
plt.colorbar()
plt.show()

#%%

nx=len(np.unique(ifile_lons))
ny=len(np.unique(ifile_lats))
nmat = len(layers)+1

with open('cal_depth_ifile_dist2depth_layers.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,nmat)
    f.write(header_line)
    for i in range(len(ifile_dists_top)):
        
        layer_1 = layers[0][i]
        layer_2 = layers[1][i]
        layer_3 = layers[2][i]
        layer_4 = layers[3][i]
        layer_5 = layers[4][i]
        layer_6 = layers[5][i]
        layer_7 = layers[6][i]
        layer_8 = layers[7][i]
        layer_9 = layers[8][i]
        layer_10 = ifile_dists[i]
        
        f.write("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n" % (ifile_lons[i],ifile_lats[i],layer_1,layer_2,layer_3,layer_4,layer_5,layer_6,layer_7,
                                                              layer_8,layer_9,layer_10))



#%%

inpoly_index_file = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/Models/inpoly_index.csv')
inpoly_index = inpoly_index_file['inpoly_index']
inpoly_lons = inpoly_index_file['inpoly_lons']  
inpoly_lats = inpoly_index_file['inpoly_lats'] 

#%%

inpoly_south_index = []
inpoly_south_lons = []
inpoly_south_lats = []

inpoly_north_index = []
inpoly_north_lons = []
inpoly_north_lats = []

for i in range(len(inpoly_index)):
    if inpoly_lats[i] < 45.28:
        inpoly_south_index.append(inpoly_index[i])
        inpoly_south_lons.append(inpoly_lons[i])
        inpoly_south_lats.append(inpoly_lats[i])
        
    else:
        inpoly_north_index.append(inpoly_index[i])
        inpoly_north_lons.append(inpoly_lons[i])
        inpoly_north_lats.append(inpoly_lats[i])
        

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_dists_lons,ifile_dists_lats,c=ifile_dists)
plt.scatter(inpoly_south_lons,inpoly_south_lats,c='r')
plt.colorbar()
plt.show()

#%%

south_dists = np.array(ifile_dists)[inpoly_south_index]

norht_grav = np.array(ifile_dists_top)[inpoly_north_index]

ifile_south_dists_norht_grav = np.zeros_like(ifile_dists)

ifile_south_dists_norht_grav[inpoly_south_index] = south_dists
ifile_south_dists_norht_grav[inpoly_north_index] = norht_grav

#%%

plt.figure(figsize=(8,10))
plt.scatter(ifile_dists_lons,ifile_dists_lats,c=ifile_south_dists_norht_grav,cmap='gist_earth_r')
plt.colorbar()
plt.show()

#%%

nx=len(np.unique(ifile_dists_lons))
ny=len(np.unique(ifile_dists_lats))



with open('cal_depth_ifile_south_dists_norht_grav.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_south_dists_norht_grav)):
        
        layer_1 = ifile_south_dists_norht_grav[i]
        layer_2 = ifile_dists_bot[i]
        
        if ifile_south_dists_norht_grav[i] > ifile_dists_bot[i]:
            layer_2 = layer_1
        

        f.write("%s %s %s %s\n" % (ifile_dists_lons[i],ifile_dists_lats[i],layer_1,layer_2))


#%%

nx=len(np.unique(ifile_dists_lons))
ny=len(np.unique(ifile_dists_lats))



with open('cal_depth_ifile_gp2_singmat.txt', 'w') as f:
    header_line="%d %d %d\n"%(nx,ny,2)
    f.write(header_line)
    for i in range(len(ifile_south_dists_norht_grav)):
        
        layer_1 = ifile_south_dists_norht_grav[i]
        
        
        f.write("%s %s %s\n" % (ifile_dists_lons[i],ifile_dists_lats[i],layer_1))


#%%

fig = plt.figure(figsize=(35, 15))
fig.tight_layout(rect=[0, 0.2, 1, 0.985],w_pad=4.8)
fig.suptitle('WV Basin Models', fontsize=60 , x=0.5)

ax1 = fig.add_subplot(1,4,1)
scat1 = ax1.scatter(ifile_dists_lons,ifile_dists_lats,c=ifile_dists,cmap='gist_earth_r', vmin=0 , vmax=5000)
ax1.set_xlim(-123.5,-122.2)
ax1.set_ylim(43.8,46.48)
ax1.set_title('Dist2Edge' ,size=25)
ax1.tick_params(labelsize=25)

ax2 = fig.add_subplot(1,4,2)
scat2 = ax2.scatter(ifile_lons,ifile_lats,c=ifile_dists_top,cmap='gist_earth_r', vmin=0 , vmax=5000)
ax2.set_xlim(-123.5,-122.2)
ax2.set_ylim(43.8,46.48)
ax2.set_title('Geologic Priors 1' ,size=25)
ax2.tick_params(labelsize=25)

ax3 = fig.add_subplot(1,4,3)
scat3 = ax3.scatter(ifile_dists_lons,ifile_dists_lats,c=ifile_south_dists_norht_grav,cmap='gist_earth_r', vmin=0 , vmax=5000)
ax3.set_xlim(-123.5,-122.2)
ax3.set_ylim(43.8,46.48)
ax3.set_title('Geologic Priors 2' ,size=25)
ax3.tick_params(labelsize=25)

ax4 = fig.add_subplot(1,4,4)
scat4 = ax4.scatter(ifile_lons,ifile_lats,c=ifile_dists_bot,cmap='gist_earth_r', vmin=0 , vmax=5000)
ax4.set_xlim(-123.5,-122.2)
ax4.set_ylim(43.8,46.48)
ax4.set_title('East Dipping Surface' ,size=25)
ax4.tick_params(labelsize=25)

cbar = fig.colorbar(scat4, ax=ax4)
cbar.ax.set_ylabel('Basin Depth [m]',size='25') 
cbar.ax.tick_params(labelsize=20)

# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# cbar=fig.colorbar(ifile_dists, cax=cbar_ax)
# cbar.set_label('Basin Depth [m]',size='13')
# cbar.ax.tick_params(labelsize=11)

plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/wv_models.jpg'
                    ,dpi=100,bbox_inches='tight')
# plt.show()































