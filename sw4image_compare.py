#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:41:06 2023

@author: rshimony
"""

import numpy as np
import matplotlib.pyplot as plt
import pySW4 as sw4
import pandas as pd
import obspy as obs
from glob import glob

#%%

image_files_rfile=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/no_refinement_figs/rfile_steph_basin_salem_gp1_sm/*.cycle=000*.s.sw4img'))
ims_x_rfile = image_files_rfile[:5]
ims_y_rfile = image_files_rfile[6:-1]
ims_z_rfile = image_files_rfile[14]

image_files_ifile=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/no_refinement_figs/eug_spe_singmat_gravfix/*.cycle=000*.s.sw4img'))
ims_x_ifile = image_files_ifile[:5]
ims_y_ifile = image_files_ifile[6:-1]
ims_z_ifile = image_files_ifile[14]

#%%

im_z_rfile = sw4.read_image(ims_z_rfile)
im_z_rfile_layer = im_z_rfile.patches[0]
im_z_rfile_data = im_z_rfile.patches[0].data

im_z_ifile = sw4.read_image(ims_z_ifile)
im_z_ifile_layer = im_z_ifile.patches[0]
im_z_ifile_data = im_z_ifile.patches[0].data

im_subtract = im_z_rfile_data - im_z_ifile_data
#%%

fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(10,14))
extent=([im_z_rfile_layer.extent[0],im_z_rfile_layer.extent[1],im_z_rfile_layer.extent[2],im_z_rfile_layer.extent[3]])
ax.tick_params(axis='both', which='major', labelsize=25)
ax.set_ylabel('Northing (m)',fontsize=27)
ax.set_xlabel('Easting (m)',fontsize=27)
im=ax.imshow(im_z_rfile_data,extent=extent,cmap='gist_earth')
ax.set_aspect('auto')
cb=fig.colorbar(im)
cb.set_label('Vs [m/sec]',size=23)
cb.ax.tick_params(labelsize=23) 
plt.show()
# fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/cs_eugspe_steph_sm_cobra.jpg',dpi=300,bbox_inches='tight')

#%%

fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(10,14))
extent=([im_z_ifile_layer.extent[0],im_z_ifile_layer.extent[1],im_z_ifile_layer.extent[2],im_z_ifile_layer.extent[3]])
ax.tick_params(axis='both', which='major', labelsize=25)
ax.set_ylabel('Northing (m)',fontsize=27)
ax.set_xlabel('Easting (m)',fontsize=27)
im=ax.imshow(im_z_ifile_data,extent=extent,cmap='gist_earth')
ax.set_aspect('auto')
cb=fig.colorbar(im)
cb.set_label('Vs [m/sec]',size=23)
cb.ax.tick_params(labelsize=23) 
plt.show()


#%%

fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(10,14))
extent=([im_z_ifile_layer.extent[0],im_z_ifile_layer.extent[1],im_z_ifile_layer.extent[2],im_z_ifile_layer.extent[3]])
ax.tick_params(axis='both', which='major', labelsize=25)
ax.set_ylabel('Northing (m)',fontsize=27)
ax.set_xlabel('Easting (m)',fontsize=27)
im=ax.imshow(im_subtract,extent=extent,cmap='seismic',origin='lower')
ax.set_aspect('auto')
cb=fig.colorbar(im)
cb.set_label('Vs [m/sec]',size=23)
cb.ax.tick_params(labelsize=23) 
plt.show()
# fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/cs_eugspe_steph_sm_cobra.jpg',dpi=300,bbox_inches='tight')

#%%

def plot_mapv_sw4image(image_file,fig_name,plot_diff_file = False):
    
    if plot_diff_file is not False:
        image = sw4.read_image(image_file)
        image2 = sw4.read_image(plot_diff_file)
        image_layer = image.patches[0]
        image_data = image.patches[0].data - image2.patches[0].data
        
        cmap = 'seismic'
        cb_label = 'Image Data Diff'        
    
    else:
        image = sw4.read_image(image_file)
        image_layer = image.patches[0]
        image_data = image.patches[0].data
        
        cmap = 'viridis'
        cb_label = 'Vs [m/sec]'
    
    fig,ax=plt.subplots(ncols=1,nrows=1,figsize=(10,14))
    extent=([image_layer.extent[0],image_layer.extent[1],image_layer.extent[2],image_layer.extent[3]])
    ax.tick_params(axis='both', which='major', labelsize=25)
    ax.set_ylabel('Northing (m)',fontsize=27)
    ax.set_xlabel('Easting (m)',fontsize=27)
    im=ax.imshow(image_data,extent=extent,cmap=cmap,origin='lower')
    ax.set_aspect('auto')
    cb=fig.colorbar(im)
    cb.set_label(cb_label,size=23)
    cb.ax.tick_params(labelsize=23) 
    plt.show()
    fig.savefig('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/'+fig_name+'.jpg',dpi=300,bbox_inches='tight')
    
#%%

plot_mapv_sw4image(ims_z_rfile,'rfile_mapv')
plot_mapv_sw4image(ims_z_ifile,'ifile_mapv')
plot_mapv_sw4image(ims_z_rfile,'rfile_ifile_diff',plot_diff_file = ims_z_ifile)

#%%

fig,axs = plt.subplots(figsize=(12, 7),nrows=1, ncols=3)
fig.subplots_adjust(wspace=0.08)

Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/rfile_mapv.jpg')
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].set_title('rfile Vs map view' ,size=16)

Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/ifile_mapv.jpg')
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].set_title('ifile Vs map view' ,size=16)

Image3 = plt.imread('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/rfile_ifile_diff.jpg')
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].set_title('rfile - ifile diff' ,size=16)

plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/rfile_ifile_mapv_comp.jpg'
                    ,dpi=300,bbox_inches='tight')


#%%

rfile_cs_ims_y=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_gp1_sm_steph_y*.jpg'))
rfile_cs_ims_x=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_gp1_sm_stephx*.jpg'))

ifile_cs_ims_y=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_gp1_sm_y*.jpg'))
ifile_cs_ims_x=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_gp1_smx*.jpg'))


rfile_notaper_cs_ims_y=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_rfile_notapery*.jpg'))
rfile_notaper_cs_ims_x=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_rfile_notaperx*.jpg'))

ifile_notaper_cs_ims_y=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_notapery*.jpg'))
ifile_notaper_cs_ims_x=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_notaperx*.jpg'))


rfile_idxpoly_cs_ims_y=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_idxpoly_y*.jpg'))
rfile_idxpoly_cs_ims_x=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_idxpoly_x*.jpg'))

rfile_polyfix_cs_ims_y=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_polyfix_y*.jpg'))
rfile_polyfix_cs_ims_x=sorted(glob('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/cs_polyfix_x*.jpg'))

#%%

fig,axs = plt.subplots(figsize=(32, 42),nrows=len(rfile_notaper_cs_ims_y), ncols=2)
fig.tight_layout()
for i in range(len(rfile_notaper_cs_ims_y)):
    image_r = plt.imread(rfile_polyfix_cs_ims_y[i])
    image_i = plt.imread(ifile_cs_ims_y[i])

    axs[i][0].imshow(image_r)
    axs[i][0].axis('off')

    axs[i][1].imshow(image_i)
    axs[i][1].axis('off')
        
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/rfile_polyfix_ifile_y_cs_comp.jpg'
                    ,dpi=300,bbox_inches='tight')


#%%

fig,axs = plt.subplots(figsize=(32, 42),nrows=len(rfile_cs_ims_x), ncols=2)
fig.tight_layout()
for i in range(len(rfile_cs_ims_x)):
    image_r = plt.imread(rfile_polyfix_cs_ims_x[i])
    image_i = plt.imread(ifile_cs_ims_x[i])

    axs[i][0].imshow(image_r)
    axs[i][0].axis('off')

    axs[i][1].imshow(image_i)
    axs[i][1].axis('off')
        
plt.savefig('/Users/rshimony/Desktop/WillametteValley/Models/rfile_ifile_comp_figs/rfile_polyfix_ifile_x_cs_comp.jpg'
                    ,dpi=300,bbox_inches='tight')

































