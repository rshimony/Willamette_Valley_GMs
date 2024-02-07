#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 12:35:46 2024

@author: rshimony
"""

import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec

#%%

all_eq_st_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/all_eq_st.csv')
all_eq_st_names = np.array(all_eq_st_ff['st_name_dom'])

wf_spec_station_comp_horiz_salem = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/wf_spec_station_comp_horiz/salem/srv_sm/*.png'))
wf_spec_station_comp_horiz_scottsmills = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/wf_spec_station_comp_horiz/scottsmills/srv_sm/*.png'))
wf_spec_station_comp_horiz_springfield = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/wf_spec_station_comp_horiz/springfield/srv_sm/*.png'))

all_eq_single_st_map = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/single_st_map/all_eq/*'))

wf_spec_st_all_eq_salem = []
wf_spec_st_all_eq_scottsmills = []
wf_spec_st_all_eq_springfield = []

for i in range(len(wf_spec_station_comp_horiz_salem)):
    wf_spec_fig_name_salem = (wf_spec_station_comp_horiz_salem[i].split('/')[-1]).split('.')[0]
    if wf_spec_fig_name_salem in all_eq_st_names:
        wf_spec_st_all_eq_salem.append(wf_spec_station_comp_horiz_salem[i])
        
for i in range(len(wf_spec_station_comp_horiz_scottsmills)):
    wf_spec_fig_name_scottsmills = (wf_spec_station_comp_horiz_scottsmills[i].split('/')[-1]).split('.')[0]
    if wf_spec_fig_name_scottsmills in all_eq_st_names:
        wf_spec_st_all_eq_scottsmills.append(wf_spec_station_comp_horiz_scottsmills[i])

for i in range(len(wf_spec_station_comp_horiz_springfield)):
    wf_spec_fig_name_springfield = (wf_spec_station_comp_horiz_springfield[i].split('/')[-1]).split('.')[0]
    if wf_spec_fig_name_springfield in all_eq_st_names:
        wf_spec_st_all_eq_springfield.append(wf_spec_station_comp_horiz_springfield[i])

#%%

all_eq_st_comp_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/subplots/all_eq_st_comp/srv_sm/'

for i in range(len(all_eq_single_st_map)):
    fig,axs = plt.subplots(figsize=(12, 10),nrows=2, ncols=2)
    fig.subplots_adjust(wspace=0.05,hspace=0.02)
    fig.suptitle(all_eq_st_names[i], fontsize=18 ,y=0.93)
    axs = axs.flatten()
    
    Image1 = plt.imread(all_eq_single_st_map[i])
    axs[0].imshow(Image1)
    axs[0].axis('off')
    axs[0].set_title('Station Location',size=15 , pad=0.3)
    
    Image2 = plt.imread(wf_spec_st_all_eq_salem[i])
    axs[1].imshow(Image2)
    axs[1].axis('off')
    axs[1].set_title('Salem' ,size=15 , pad=0.3)
    
    Image3 = plt.imread(wf_spec_st_all_eq_scottsmills[i])
    axs[2].imshow(Image3)
    axs[2].axis('off')
    axs[2].set_title('Scotts Mills' ,size=15 , pad=0.3)
    
    Image4 = plt.imread(wf_spec_st_all_eq_springfield[i])
    axs[3].imshow(Image4)
    axs[3].axis('off')
    axs[3].set_title('Springfield' ,size=15 , pad=0.3)
    
    plt.savefig(all_eq_st_comp_figs_dir + all_eq_st_names[i] + '.jpg'
                        ,dpi=300,bbox_inches='tight')

#%%

eq = 'springfield'
im = 'pgv'

res_maps_val = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/'+eq+'/'+im+'*_val*.png'))
res_maps_val = res_maps_val[:-1]
boxplots = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/models_boxplots/'+eq+'/'+im+'*_val.png'))

ims_maps_boxplots_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/subplots/ims_maps_boxplots/springfield/'

plt.figure(figsize=(6,11))
gridspec.GridSpec(9,8)
# plt.suptitle('FAS Bin = '+im, fontsize=18 ,y=0.93)
plt.suptitle('PGV', fontsize=18 ,y=0.93)
plt.subplots_adjust(hspace=0.4)

ax1 = plt.subplot2grid((9,8), (0,0), colspan=4, rowspan=3)
Image1 = plt.imread(res_maps_val[0])
ax1.imshow(Image1)
ax1.axis('off')
ax1.set_title('eugspe_1d',size=15 , pad=0.3)

ax2 = plt.subplot2grid((9,8), (0,4), colspan=4, rowspan=3)
Image2 = plt.imread(res_maps_val[1])
ax2.imshow(Image2)
ax2.axis('off')
ax2.set_title('eugspe_sm',size=15 , pad=0.3)

ax3 = plt.subplot2grid((9,8), (3,0), colspan=4, rowspan=3)
Image3 = plt.imread(res_maps_val[2])
ax3.imshow(Image3)
ax3.axis('off')
ax3.set_title('srv_1d',size=15 , pad=0.3)

ax4 = plt.subplot2grid((9,8), (3,4), colspan=4, rowspan=3)
Image4 = plt.imread(res_maps_val[3])
ax4.imshow(Image4)
ax4.axis('off')
ax4.set_title('srv_sm',size=15 , pad=0.3)

ax5 = plt.subplot2grid((9,8), (6,0), colspan=8, rowspan=3)
Image5 = plt.imread(boxplots[0])
ax5.imshow(Image5)
ax5.axis('off')

plt.savefig(ims_maps_boxplots_dir + im + '_im_maps_boxplot.jpg' ,dpi=300,bbox_inches='tight')
    
    
#%%

eq = 'springfield'
model = 'srv_sm'

res_maps_fas_val = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/'+eq+'/*fas_val_'+model+'*.png'))
res_maps_pgv_val = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/'+eq+'/pgv_'+model+'*_val*.png'))
res_maps_val_model = res_maps_fas_val + res_maps_pgv_val

all_ims_maps_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/subplots/all_im_maps/'

fig,axs = plt.subplots(figsize=(10, 6),nrows=2, ncols=4)
# fig.subplots_adjust(wspace=0.05)
fig.suptitle(eq +' '+ model, fontsize=18 ,y=0.98)
axs = axs.flatten()

Image1 = plt.imread(res_maps_val_model[0])
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].set_title('FAS Bin = 0.2 Hz',size=12 , pad=0.4)

Image2 = plt.imread(res_maps_val_model[1])
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].set_title('FAS Bin = 0.3 Hz' ,size=12 , pad=0.4)

Image3 = plt.imread(res_maps_val_model[2])
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].set_title('FAS Bin = 0.4 Hz' ,size=12 , pad=0.4)

Image4 = plt.imread(res_maps_val_model[3])
axs[3].imshow(Image4)
axs[3].axis('off')
axs[3].set_title('FAS Bin = 0.5 Hz' ,size=12 , pad=0.4)

Image4 = plt.imread(res_maps_val_model[4])
axs[4].imshow(Image4)
axs[4].axis('off')
axs[4].set_title('FAS Bin = 0.6 Hz' ,size=12 , pad=0.4)

Image4 = plt.imread(res_maps_val_model[5])
axs[5].imshow(Image4)
axs[5].axis('off')
axs[5].set_title('FAS Bin = 0.7 Hz' ,size=12 , pad=0.4)

Image4 = plt.imread(res_maps_val_model[6])
axs[6].imshow(Image4)
axs[6].axis('off')
axs[6].set_title('PGV' ,size=12 , pad=0.4)

axs[7].axis('off')

plt.savefig(all_ims_maps_figs_dir + eq + '_' + model + '.jpg' ,dpi=300,bbox_inches='tight')

#%%

wf_spec_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/wf_spec_model_comp/springfield/*'))
single_st_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/single_st_map/springfield/*'))

wf_spec_comp_st_map_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/subplots/wf_spec_comp_st_map/springfield/'

for i in range(len(wf_spec_figs)):
    plt.figure(figsize=(6,3))
    gridspec.GridSpec(1,6)
    # plt.suptitle('FAS Bin = '+im, fontsize=18 ,y=0.93)
    # plt.suptitle('PGV', fontsize=18 ,y=0.93)
    # plt.subplots_adjust(hspace=0.4)
    
    ax1 = plt.subplot2grid((1,6), (0,0), colspan=4, rowspan=1)
    Image1 = plt.imread(wf_spec_figs[i])
    ax1.imshow(Image1)
    ax1.axis('off')
    # ax1.set_title('eugspe_1d',size=15 , pad=0.3)
    
    ax2 = plt.subplot2grid((1,6), (0,4), colspan=2, rowspan=1)
    Image2 = plt.imread(single_st_figs[i])
    ax2.imshow(Image2)
    ax2.axis('off')
    # ax1.set_title('eugspe_1d',size=15 , pad=0.3)
    
    st_name = (single_st_figs[i].split('/')[-1]).split('.')[0]
    
    plt.savefig(wf_spec_comp_st_map_dir + st_name +'.jpg' ,dpi=800,bbox_inches='tight')

#%%

ratio_map_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/ratio_maps/springfield/*ratio.png'))
ratio_val_boxplot_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ratio_boxplot_valley/springfield/*'))

ratio_maps_val_boxplots_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/subplots/ratio_maps_val_boxplots/springfield/'

for i in range(len(ratio_map_figs)):
    plt.figure(figsize=(6,4))
    gridspec.GridSpec(6,1)
    # plt.suptitle('FAS Bin = '+im, fontsize=18 ,y=0.93)
    # plt.suptitle('PGV', fontsize=18 ,y=0.93)
    # plt.subplots_adjust(hspace=0.4)
    
    ax1 = plt.subplot2grid((6,1), (0,0), colspan=1, rowspan=4)
    Image1 = plt.imread(ratio_map_figs[i])
    ax1.imshow(Image1)
    ax1.axis('off')
    # ax1.set_title('eugspe_1d',size=15 , pad=0.3)
    
    ax2 = plt.subplot2grid((6,1), (4,0), colspan=1, rowspan=2)
    Image2 = plt.imread(ratio_val_boxplot_figs[i])
    ax2.imshow(Image2)
    ax2.axis('off')
    # ax1.set_title('eugspe_1d',size=15 , pad=0.3)
    
    fig_name = (ratio_map_figs[i].split('/')[-1]).split('.')[0]
    
    plt.savefig(ratio_maps_val_boxplots_dir + fig_name +'.jpg' ,dpi=800,bbox_inches='tight')














































