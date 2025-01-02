#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 14:51:58 2024

@author: rshimony
"""
'''
Making final figure subplots combination to use in publication and supplementary material.
Each cell in this script can run separately and creates a set of (or a single) figures.
'''
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec

#%%
'''
All EQs waveform and spectra comparison.
These subplot figures combine waveform and spectra figures from all 3 EQs and a station location map
for all stations that recorded all 3 EQs.
Run this for EACH MODEL SEPARATELY.
Used in supplementary materials.
'''
## Reading all EQs stations flatfile
all_eq_st_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/observed_data/all_eq_st.csv')
all_eq_st_names = np.array(sorted(all_eq_st_ff['st_nm']))

## Reading waveform and spectra comp figures from the 3 EQs
wf_spec_station_comp_horiz_salem = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_station_comp/salem/srv_sm/*.png'))
wf_spec_station_comp_horiz_scottsmills = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_station_comp/scottsmills/srv_sm/*.png'))
wf_spec_station_comp_horiz_springfield = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_station_comp/springfield/srv_sm/*.png'))

## Reading all EQs station location maps
all_eq_single_st_map = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/st_loc_maps/all_eqs/*'))

## Getting waveform and spectra figs of all EQs stations
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

## Plottin all EQs waveform and spectra comparison figs
# Output directory
all_eq_st_comp_figs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/all_eq_st_comp/srv_sm/'
# Plotting
for i in range(len(all_eq_single_st_map)):
    fig,axs = plt.subplots(figsize=(6, 6),nrows=2, ncols=2)
    fig.subplots_adjust(wspace=0.05,hspace=0.0)
    fig.suptitle(all_eq_st_names[i], fontsize=18 ,y=0.99)
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
'''
Combining SINGLE IM maps from ALL EQs to one figure,
along with a boxplot comparison of this IM.
Used in supplementary materials.
'''
## Setting EQs, IMs and fig titles
eqs = ['salem','scottsmills','springfield']
ims = ['f2','f3','f4','f5','f6','f7','pgv']
im_title = ['FAS Bin = 0.2 Hz','FAS Bin = 0.3 Hz','FAS Bin = 0.4 Hz','FAS Bin = 0.5 Hz','FAS Bin = 0.6 Hz','FAS Bin = 0.7 Hz','PGV']

## Plotting
for eq in eqs:
    for im in range(len(ims)):
        
        ims_maps_boxplots_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/ims_maps_boxplots/'+eq+'/'
        
        res_maps_val = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/'+eq+'/'+ims[im]+'*_val*.png'))
        res_maps_val = res_maps_val[:-1]
        boxplots = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/boxplots/res_val_dist/'+eq+'/'+ims[im]+'*_val.png'))
        
        plt.figure(figsize=(3,7))
        gridspec.GridSpec(9,8)
        plt.suptitle(im_title[im], fontsize=18 ,y=1.02)
        plt.subplots_adjust(hspace=0.0)
        
        ax1 = plt.subplot2grid((9,8), (0,0), colspan=4, rowspan=3)
        Image1 = plt.imread(res_maps_val[0])
        ax1.imshow(Image1)
        ax1.axis('off')
        ax1.set_title('EugSpe 1D',size=15 , pad=0.3)
        
        ax2 = plt.subplot2grid((9,8), (0,4), colspan=4, rowspan=3)
        Image2 = plt.imread(res_maps_val[1])
        ax2.imshow(Image2)
        ax2.axis('off')
        ax2.set_title('EugSpe SM',size=15 , pad=0.3)
        
        ax3 = plt.subplot2grid((9,8), (3,0), colspan=4, rowspan=3)
        Image3 = plt.imread(res_maps_val[2])
        ax3.imshow(Image3)
        ax3.axis('off')
        ax3.set_title('SRV 1D',size=15 , pad=0.3)
        
        ax4 = plt.subplot2grid((9,8), (3,4), colspan=4, rowspan=3)
        Image4 = plt.imread(res_maps_val[3])
        ax4.imshow(Image4)
        ax4.axis('off')
        ax4.set_title('SRV SM',size=15 , pad=0.3)
        
        ax5 = plt.subplot2grid((9,8), (6,0), colspan=8, rowspan=3)
        Image5 = plt.imread(boxplots[0])
        ax5.imshow(Image5)
        ax5.axis('off')
        
        plt.savefig(ims_maps_boxplots_dir + ims[im] + '_im_maps_boxplot.jpg' ,dpi=300,bbox_inches='tight')  
    
#%%
'''
Combining all IMs RESIDUALS maps of every EQ and MODEL pair.
Used in supplementary materials.
'''

eqs_models = [['salem','eugspe_1d'],['salem','eugspe_sm'],['salem','srv_1d'],['salem','srv_sm'],
              ['scottsmills','eugspe_1d'],['scottsmills','eugspe_sm'],['scottsmills','srv_1d'],['scottsmills','srv_sm'],
              ['springfield','eugspe_1d'],['springfield','eugspe_sm'],['springfield','srv_1d'],['springfield','srv_sm']]
titles = ['Salem EugSpe 1D','Salem EugSpe SM','Salem SRV 1D','Salem SRV SM',
         'Scotts Mills EugSpe 1D','Scotts Mills EugSpe SM','Scotts Mills SRV 1D','Scotts Mills SRV SM',
         'Springfield EugSpe 1D','Springfield EugSpe SM','Springfield SRV 1D','Springfield SRV SM']

for em in range(len(eqs_models)):
    
    res_maps_fas_val = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/'+eqs_models[em][0]+'/*fas_val_'+eqs_models[em][1]+'*.png'))
    res_maps_pgv_val = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/'+eqs_models[em][0]+'/pgv_'+eqs_models[em][1]+'*_val*.png'))
    res_maps_val_model = res_maps_fas_val + res_maps_pgv_val
    
    all_ims_maps_figs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/all_im_maps/'
    
    fig,axs = plt.subplots(figsize=(8.5, 6),nrows=2, ncols=4)
    fig.subplots_adjust(wspace=0.,hspace=0.1)
    fig.suptitle(titles[em], fontsize=18 ,y=1.02)
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
    
    plt.savefig(all_ims_maps_figs_dir + eqs_models[em][0] + '_' + eqs_models[em][1] + '.jpg' ,dpi=300,bbox_inches='tight')

#%%

'''
Grouping selected IMs maps (f=0.3,0.5,0.7Hz,PGV) from the SRV 1D model simulations for all events. 
FIGURE 6 in publication.
'''

res_maps_salem = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/salem/*_val_srv_1d_res*'))
res_maps_salem.append('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/salem/pgv_srv_1d_res_val.png')
res_maps_scottsmills = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/scottsmills/*_val_srv_1d_res*'))
res_maps_scottsmills.append('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/scottsmills/pgv_srv_1d_res_val.png')
res_maps_springfield = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/springfield/*_val_srv_1d_res*'))
res_maps_springfield.append('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/residuals_maps/springfield/pgv_srv_1d_res_val.png')

grouped_res_maps_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/grouped_res_maps/'

fig,axs = plt.subplots(figsize=(6, 8),nrows=3, ncols=4)
fig.subplots_adjust(wspace=0,hspace=0.17)
axs = axs.flatten()

Image1 = plt.imread(res_maps_salem[1])
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].set_title('f = 0.3Hz',size=9 , pad=0.3)
axs[0].text(0,-50,'a)',size=8,weight='bold')
axs[0].text(1000,2100,'B028',size=5,weight='bold')
axs[0].text(1000,2400,'ALVY',size=5,weight='bold')
axs[0].plot([985, 830], [2050, 1870], color='k', linestyle='-', linewidth=0.5)
axs[0].plot([985, 790], [2350, 2420], color='k', linestyle='-', linewidth=0.5)

Image2 = plt.imread(res_maps_salem[3])
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].set_title('f = 0.5Hz',size=9 , pad=0.3)
axs[1].text(0,-50,'b)',size=8,weight='bold')


Image3 = plt.imread(res_maps_salem[5])
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].set_title('f = 0.7Hz',size=9, pad=0.3)
axs[2].text(0,-50,'c)',size=8,weight='bold')

Image4 = plt.imread(res_maps_salem[6])
axs[3].imshow(Image4)
axs[3].axis('off')
axs[3].set_title('PGV',size=9 , pad=0.3)
axs[3].text(0,-50,'d)',size=8,weight='bold')

Image5 = plt.imread(res_maps_scottsmills[1])
axs[4].imshow(Image5)
axs[4].axis('off')
axs[4].text(0,-50,'e)',size=8,weight='bold')
axs[4].text(1000,2100,'B028',size=5,weight='bold')
axs[4].text(1000,2400,'ALVY',size=5,weight='bold')
axs[4].plot([985, 830], [2050, 1870], color='k', linestyle='-', linewidth=0.5)
axs[4].plot([985, 790], [2350, 2420], color='k', linestyle='-', linewidth=0.5)

Image6 = plt.imread(res_maps_scottsmills[3])
axs[5].imshow(Image6)
axs[5].axis('off')
axs[5].text(0,-50,'f)',size=8,weight='bold')

Image7 = plt.imread(res_maps_scottsmills[5])
axs[6].imshow(Image7)
axs[6].axis('off')
axs[6].text(0,-50,'g)',size=8,weight='bold')

Image8 = plt.imread(res_maps_scottsmills[6])
axs[7].imshow(Image8)
axs[7].axis('off')
axs[7].text(0,-50,'h)',size=8,weight='bold')

Image9 = plt.imread(res_maps_springfield[1])
axs[8].imshow(Image9)
axs[8].axis('off')
axs[8].text(0,-50,'i)',size=8,weight='bold')
axs[8].text(1000,2100,'B028',size=5,weight='bold')
axs[8].text(1000,2400,'ALVY',size=5,weight='bold')
axs[8].plot([985, 830], [2050, 1870], color='k', linestyle='-', linewidth=0.5)
axs[8].plot([985, 790], [2350, 2420], color='k', linestyle='-', linewidth=0.5)

Image10 = plt.imread(res_maps_springfield[3])
axs[9].imshow(Image10)
axs[9].axis('off')
axs[9].text(0,-50,'j)',size=8,weight='bold')

Image11 = plt.imread(res_maps_springfield[5])
axs[10].imshow(Image11)
axs[10].axis('off')
axs[10].text(0,-50,'k)',size=8,weight='bold')

Image12 = plt.imread(res_maps_springfield[6])
axs[11].imshow(Image12)
axs[11].axis('off')
axs[11].text(0,-50,'l)',size=8,weight='bold')

fig.text(0.14,0.90,'Salem',fontsize=12)
fig.text(0.14,0.632,'Scotts Mills',fontsize=12)
fig.text(0.14,0.368,'Springfield',fontsize=12)

plt.savefig(grouped_res_maps_dir + 'grouped_res_maps.png' ,dpi=600,bbox_inches='tight')


#%%
'''
Combining all IMs RATIOS maps of every EQ and MODEL pair.
Used in supplementary materials.
'''
eqs_models = [['salem','eugspe_1d'],['salem','eugspe_sm'],['salem','srv_1d'],['salem','srv_sm'],
              ['scottsmills','eugspe_1d'],['scottsmills','eugspe_sm'],['scottsmills','srv_1d'],['scottsmills','srv_sm'],
              ['springfield','eugspe_1d'],['springfield','eugspe_sm'],['springfield','srv_1d'],['springfield','srv_sm']]
titles = ['Salem EugSpe 1D','Salem EugSpe SM','Salem SRV 1D','Salem SRV SM',
         'Scotts Mills EugSpe 1D','Scotts Mills EugSpe SM','Scotts Mills SRV 1D','Scotts Mills SRV SM',
         'Springfield EugSpe 1D','Springfield EugSpe SM','Springfield SRV 1D','Springfield SRV SM']

for em in range(len(eqs_models)):
    
    rat_maps_fas_val = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/ratios_maps/'+eqs_models[em][0]+'/*'+eqs_models[em][1]+'*ratio.png'))
    
    all_ims_maps_figs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/all_im_maps/'
    
    fig,axs = plt.subplots(figsize=(10, 6),nrows=2, ncols=4)
    fig.subplots_adjust(wspace=0.05)
    fig.suptitle(titles[em], fontsize=18 ,y=0.98)
    axs = axs.flatten()
    
    Image1 = plt.imread(rat_maps_fas_val[0])
    axs[0].imshow(Image1)
    axs[0].axis('off')
    axs[0].set_title('FAS Bin = 0.2 Hz',size=12 , pad=0.4)
    
    Image2 = plt.imread(rat_maps_fas_val[1])
    axs[1].imshow(Image2)
    axs[1].axis('off')
    axs[1].set_title('FAS Bin = 0.3 Hz' ,size=12 , pad=0.4)
    
    Image3 = plt.imread(rat_maps_fas_val[2])
    axs[2].imshow(Image3)
    axs[2].axis('off')
    axs[2].set_title('FAS Bin = 0.4 Hz' ,size=12 , pad=0.4)
    
    Image4 = plt.imread(rat_maps_fas_val[3])
    axs[3].imshow(Image4)
    axs[3].axis('off')
    axs[3].set_title('FAS Bin = 0.5 Hz' ,size=12 , pad=0.4)
    
    Image4 = plt.imread(rat_maps_fas_val[4])
    axs[4].imshow(Image4)
    axs[4].axis('off')
    axs[4].set_title('FAS Bin = 0.6 Hz' ,size=12 , pad=0.4)
    
    Image4 = plt.imread(rat_maps_fas_val[5])
    axs[5].imshow(Image4)
    axs[5].axis('off')
    axs[5].set_title('FAS Bin = 0.7 Hz' ,size=12 , pad=0.4)
    
    Image4 = plt.imread(rat_maps_fas_val[6])
    axs[6].imshow(Image4)
    axs[6].axis('off')
    axs[6].set_title('PGV' ,size=12 , pad=0.4)
    
    axs[7].axis('off')
    
    plt.savefig(all_ims_maps_figs_dir + eqs_models[em][0] + '_' + eqs_models[em][1] + '_rat.jpg' ,dpi=300,bbox_inches='tight')
#%%
'''
For each station in each EQ, combine the 3 component waveform and spectra model comparison with station location map.
Used in supplementary materials.
'''
eqs = ['salem','scottsmills','springfield']

for eq in eqs:
    wf_spec_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_model_comp/'+eq+'/*'))
    single_st_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/st_loc_maps/'+eq+'/*'))
    
    wf_spec_comp_st_map_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/wf_spec_comp_st_map/'+eq+'/'
    
    for i in range(len(wf_spec_figs)):
        plt.figure(figsize=(6,3))
        gridspec.GridSpec(1,6)
        
        ax1 = plt.subplot2grid((1,6), (0,0), colspan=4, rowspan=1)
        Image1 = plt.imread(wf_spec_figs[i])
        ax1.imshow(Image1)
        ax1.axis('off')
        
        ax2 = plt.subplot2grid((1,6), (0,4), colspan=2, rowspan=1)
        Image2 = plt.imread(single_st_figs[i])
        ax2.imshow(Image2)
        ax2.axis('off')
        
        st_name = (single_st_figs[i].split('/')[-1]).split('.')[0]
        
        plt.savefig(wf_spec_comp_st_map_dir + st_name +'.jpg' ,dpi=800,bbox_inches='tight')

#%%
'''
For all IMs in ALL MODELS and ALL EQS, combine RATIO MAPS with Valley or Mountain (out of valley) stations BOXPLOTS.

The second part of this section takes the (previously made, in this section) 
PGV, SRV 1D model combined boxplot ratio maps from the 3 EQs and grouping them together.
This is FIGURE 5 in the publication.
'''
eqs = ['salem','scottsmills','springfield']

for eq in eqs:
    ratio_map_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/ims_maps/ratios_maps/'+eq+'/*ratio.png'))
    ratio_val_boxplot_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/boxplots/ratio_inout_val_dist/'+eq+'/*'))
    
    ratio_maps_val_boxplots_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/'+eq+'/'
    
    for i in range(len(ratio_map_figs)):
        plt.figure(figsize=(6,4))
        gridspec.GridSpec(6,1)
        plt.subplots_adjust(wspace=0,hspace=0)
        
        ax1 = plt.subplot2grid((6,1), (0,0), colspan=1, rowspan=4)
        Image1 = plt.imread(ratio_map_figs[i])
        ax1.imshow(Image1)
        ax1.axis('off')
        
        ax2 = plt.subplot2grid((6,1), (4,0), colspan=1, rowspan=2)
        Image2 = plt.imread(ratio_val_boxplot_figs[i])
        ax2.imshow(Image2)
        ax2.axis('off')
        
        fig_name = (ratio_map_figs[i].split('/')[-1]).split('.')[0]
        
        plt.savefig(ratio_maps_val_boxplots_dir + fig_name +'.jpg' ,dpi=800,bbox_inches='tight')
    
### Plotting FIGURE 5    
salem_ratio_map_boxplot = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/salem/pgv_srv_1d_ratio.jpg'
scottsmills_ratio_map_boxplot = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/scottsmills/pgv_srv_1d_ratio.jpg'
springfield_ratio_map_boxplot = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/springfield/pgv_srv_1d_ratio.jpg'

eq_comp_ratio_maps_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/'

fig,axs = plt.subplots(figsize=(10, 6),nrows=1, ncols=3)
fig.subplots_adjust(wspace=0,hspace=0)
axs = axs.flatten()

Image1 = plt.imread(salem_ratio_map_boxplot)
axs[0].imshow(Image1)
axs[0].set_title('Salem',size=13,pad=0.5)
axs[0].axis('off')
axs[0].text(0,0,'a)',size=15,weight='bold')

Image2 = plt.imread(scottsmills_ratio_map_boxplot)
axs[1].imshow(Image2)
axs[1].set_title('Scotts Mills',size=13,pad=0.5)
axs[1].axis('off')
axs[1].text(0,0,'b)',size=15,weight='bold')

Image3 = plt.imread(springfield_ratio_map_boxplot)
axs[2].imshow(Image3)
axs[2].set_title('Springfield',size=13,pad=0.5)
axs[2].axis('off')
axs[2].text(0,0,'c)',size=15,weight='bold')

plt.savefig(eq_comp_ratio_maps_dir + 'eq_comp_ratio_maps.jpg' ,dpi=300,bbox_inches='tight')

#%%
'''
Making a grouped cross section figure, showing 3 selected W-E cross sections from all 4 new WV models.
This is FIGURE 3 in the publication.
'''

eugspe_1d_cs_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/models_cs/eugspe_1d/*'))
eugspe_sm_cs_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/models_cs/eugspe_sm/*'))
srv_1d_cs_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/models_cs/srv_1d/*'))
srv_sm_cs_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/models_cs/srv_sm/*'))

all_models_cs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/grouped_model_cs/'

plt.figure(figsize=(8,10))
gridspec.GridSpec(13,2)

ax1 = plt.subplot2grid((13,2), (0,0), colspan=1, rowspan=2)
Image1 = plt.imread(eugspe_1d_cs_figs[4])
ax1.imshow(Image1)
ax1.axis('off')
ax1.text(300,75,'A',color='red', weight='bold')
ax1.text(3150,75,"A'",color='red', weight='bold')
ax1.set_title('EugSpe 1D',size=18 , pad=10)

ax2 = plt.subplot2grid((13,2), (2,0), colspan=1, rowspan=2)
Image2 = plt.imread(eugspe_1d_cs_figs[2])
ax2.imshow(Image2)
ax2.axis('off')
ax2.text(300,75,'B',color='red', weight='bold')
ax2.text(3150,75,"B'",color='red', weight='bold')

ax3 = plt.subplot2grid((13,2), (4,0), colspan=1, rowspan=2)
Image3 = plt.imread(eugspe_1d_cs_figs[0])
ax3.imshow(Image3)
ax3.axis('off')
ax3.text(300,75,'C',color='red', weight='bold')
ax3.text(3150,75,"C'",color='red', weight='bold')

ax4 = plt.subplot2grid((13,2), (7,0), colspan=1, rowspan=2)
Image4 = plt.imread(eugspe_sm_cs_figs[4])
ax4.imshow(Image4)
ax4.axis('off')
ax4.text(300,75,'A',color='red', weight='bold')
ax4.text(3150,75,"A'",color='red', weight='bold')
ax4.set_title('EugSpe SM',size=18 , pad=10)

ax5 = plt.subplot2grid((13,2), (9,0), colspan=1, rowspan=2)
Image5 = plt.imread(eugspe_sm_cs_figs[2])
ax5.imshow(Image5)
ax5.axis('off')
ax5.text(300,75,'B',color='red', weight='bold')
ax5.text(3150,75,"B'",color='red', weight='bold')

ax6 = plt.subplot2grid((13,2), (11,0), colspan=1, rowspan=2)
Image6 = plt.imread(eugspe_sm_cs_figs[0])
ax6.imshow(Image6)
ax6.axis('off')
ax6.text(300,75,'C',color='red', weight='bold')
ax6.text(3150,75,"C'",color='red', weight='bold')

ax7 = plt.subplot2grid((13,2), (0,1), colspan=1, rowspan=2)
Image7 = plt.imread(srv_1d_cs_figs[4])
ax7.imshow(Image7)
ax7.axis('off')
ax7.text(300,75,'A',color='red', weight='bold')
ax7.text(3150,75,"A'",color='red', weight='bold')
ax7.set_title('SRV 1D',size=18 , pad=10)

ax8 = plt.subplot2grid((13,2), (2,1), colspan=1, rowspan=2)
Image8 = plt.imread(srv_1d_cs_figs[2])
ax8.imshow(Image8)
ax8.axis('off')
ax8.text(300,75,'B',color='red', weight='bold')
ax8.text(3150,75,"B'",color='red', weight='bold')

ax9 = plt.subplot2grid((13,2), (4,1), colspan=1, rowspan=2)
Image9 = plt.imread(srv_1d_cs_figs[0])
ax9.imshow(Image9)
ax9.axis('off')
ax9.text(300,75,'C',color='red', weight='bold')
ax9.text(3150,75,"C'",color='red', weight='bold')

ax10 = plt.subplot2grid((13,2), (7,1), colspan=1, rowspan=2)
Image10 = plt.imread(srv_sm_cs_figs[4])
ax10.imshow(Image10)
ax10.axis('off')
ax10.set_title('SRV SM',size=18 , pad=10)
ax10.text(300,75,'A',color='red', weight='bold')
ax10.text(3150,75,"A'",color='red', weight='bold')

ax11 = plt.subplot2grid((13,2), (9,1), colspan=1, rowspan=2)
Image11 = plt.imread(srv_sm_cs_figs[2])
ax11.imshow(Image11)
ax11.axis('off')
ax11.text(300,75,'B',color='red', weight='bold')
ax11.text(3150,75,"B'",color='red', weight='bold')

ax12 = plt.subplot2grid((13,2), (11,1), colspan=1, rowspan=2)
Image12 = plt.imread(srv_sm_cs_figs[0])
ax12.imshow(Image12)
ax12.axis('off')
ax12.text(300,75,'C',color='red', weight='bold')
ax12.text(3150,75,"C'",color='red', weight='bold')

plt.savefig(all_models_cs_dir + 'all_models_cs.jpg' ,dpi=800,bbox_inches='tight')

#%%
'''
IM maps comparison figure.
Plotting selected IMs (f=0.3,0.5,0.7Hz, PGV) from the 3 EQs together for comparison.
This is FIGURE 6 in the publication
'''

res_maps_salem = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/salem/*_val_srv_1d_res*'))
res_maps_salem.append('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/salem/pgv_srv_1d_res_val.png')
res_maps_scottsmills = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/scottsmills/*_val_srv_1d_res*'))
res_maps_scottsmills.append('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/scottsmills/pgv_srv_1d_res_val.png')
res_maps_springfield = sorted(glob('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/springfield/*_val_srv_1d_res*'))
res_maps_springfield.append('/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ims_maps/residuals_maps/springfield/pgv_srv_1d_res_val.png')

res_eq_comp_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/subplots/res_eq_comp/'

fig,axs = plt.subplots(figsize=(6, 8),nrows=3, ncols=4)
fig.subplots_adjust(wspace=0,hspace=0.17)
axs = axs.flatten()

Image1 = plt.imread(res_maps_salem[1])
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].set_title('f = 0.3Hz',size=9 , pad=0.3)
axs[0].text(0,-50,'a)',size=8,weight='bold')
axs[0].text(1000,2100,'B028',size=5,weight='bold')
axs[0].text(1000,2400,'ALVY',size=5,weight='bold')
axs[0].plot([985, 830], [2050, 1870], color='k', linestyle='-', linewidth=0.5)
axs[0].plot([985, 790], [2350, 2420], color='k', linestyle='-', linewidth=0.5)

Image2 = plt.imread(res_maps_salem[3])
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].set_title('f = 0.5Hz',size=9 , pad=0.3)
axs[1].text(0,-50,'b)',size=8,weight='bold')


Image3 = plt.imread(res_maps_salem[5])
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].set_title('f = 0.7Hz',size=9, pad=0.3)
axs[2].text(0,-50,'c)',size=8,weight='bold')

Image4 = plt.imread(res_maps_salem[6])
axs[3].imshow(Image4)
axs[3].axis('off')
axs[3].set_title('PGV',size=9 , pad=0.3)
axs[3].text(0,-50,'d)',size=8,weight='bold')

Image5 = plt.imread(res_maps_scottsmills[1])
axs[4].imshow(Image5)
axs[4].axis('off')
axs[4].text(0,-50,'e)',size=8,weight='bold')
axs[4].text(1000,2100,'B028',size=5,weight='bold')
axs[4].text(1000,2400,'ALVY',size=5,weight='bold')
axs[4].plot([985, 830], [2050, 1870], color='k', linestyle='-', linewidth=0.5)
axs[4].plot([985, 790], [2350, 2420], color='k', linestyle='-', linewidth=0.5)

Image6 = plt.imread(res_maps_scottsmills[3])
axs[5].imshow(Image6)
axs[5].axis('off')
axs[5].text(0,-50,'f)',size=8,weight='bold')

Image7 = plt.imread(res_maps_scottsmills[5])
axs[6].imshow(Image7)
axs[6].axis('off')
axs[6].text(0,-50,'g)',size=8,weight='bold')

Image8 = plt.imread(res_maps_scottsmills[6])
axs[7].imshow(Image8)
axs[7].axis('off')
axs[7].text(0,-50,'h)',size=8,weight='bold')

Image9 = plt.imread(res_maps_springfield[1])
axs[8].imshow(Image9)
axs[8].axis('off')
axs[8].text(0,-50,'i)',size=8,weight='bold')
axs[8].text(1000,2100,'B028',size=5,weight='bold')
axs[8].text(1000,2400,'ALVY',size=5,weight='bold')
axs[8].plot([985, 830], [2050, 1870], color='k', linestyle='-', linewidth=0.5)
axs[8].plot([985, 790], [2350, 2420], color='k', linestyle='-', linewidth=0.5)

Image10 = plt.imread(res_maps_springfield[3])
axs[9].imshow(Image10)
axs[9].axis('off')
axs[9].text(0,-50,'j)',size=8,weight='bold')

Image11 = plt.imread(res_maps_springfield[5])
axs[10].imshow(Image11)
axs[10].axis('off')
axs[10].text(0,-50,'k)',size=8,weight='bold')

Image12 = plt.imread(res_maps_springfield[6])
axs[11].imshow(Image12)
axs[11].axis('off')
axs[11].text(0,-50,'l)',size=8,weight='bold')

fig.text(0.14,0.90,'Salem',fontsize=12)
fig.text(0.14,0.632,'Scotts Mills',fontsize=12)
fig.text(0.14,0.368,'Springfield',fontsize=12)

plt.savefig(res_eq_comp_dir + 'res_eq_comp.jpg' ,dpi=300,bbox_inches='tight')

#%%
'''
Wells data figure.
Plotting 2 wells location maps (SRV and EugSpe penetrating wells) 
grouped with digitized well logs and averaged velocity gradient figure.
This is FIGURE 2 in the publication.
'''

fig,axs = plt.subplots(figsize=(10, 12),nrows=1, ncols=3)
axs = axs.flatten()

Image1 = plt.imread('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wells_maps/srv_wells_map.png')
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].set_title('SRV Wells',size=11 , pad=0.3)
axs[0].text(100,-50,'a)',size=15,weight='bold')

Image2 = plt.imread('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wells_maps/eugspe_wells_map.png')
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].set_title('EugSpe Wells',size=11 , pad=0.3)
axs[1].text(100,-50,'b)',size=15,weight='bold')

Image3 = plt.imread('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/digi_wells/digi_wells_stacked.jpg')
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].set_title('Digitized Wells',size=11, pad=1.6)
axs[2].text(100,-50,'c)',size=15,weight='bold')

plt.savefig('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/wells_fig/wells_fig.png' ,dpi=300,bbox_inches='tight')

#%%
'''
Combining wave propagation maps from the selected model (SRV 1D) and the USGS CVM, from the 3 EQs.
Making this combined figure for every imaged time step (0.5s) in the simulation.
This figures are later made into wave propagation video (using FFMPEG) that is used in the supplementary materials.
'''
salem_srv_1d_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/salem_srv_1d_freq8_mags_figs/*.png'),
                           key = lambda x:float(x.split('/')[-1][:-5]))
salem_steph_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/salem_steph_freq8_mags_figs/*.png'),
                           key = lambda x:float(x.split('/')[-1][:-5]))
scottsmills_srv_1d_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/scottsmills_srv_1d_freq8_mags_figs/*.png'),
                           key = lambda x:float(x.split('/')[-1][:-5]))
scottsmills_steph_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/scottsmills_steph_freq8_mags_figs/*.png'),
                           key = lambda x:float(x.split('/')[-1][:-5]))
springfield_srv_1d_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/springfield_srv_1d_freq8_mags_figs/*.png'),
                           key = lambda x:float(x.split('/')[-1][:-5]))
springfield_steph_figs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/wave_propagation/springfield_steph_freq8_mags_figs/*.png'),
                           key = lambda x:float(x.split('/')[-1][:-5]))

wp_prog_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/wave_propagation/'

for i in range(len(salem_srv_1d_figs)):
    fig,axs = plt.subplots(figsize=(10, 8),nrows=2, ncols=3)
    fig.subplots_adjust(wspace=0,hspace=0)
    axs = axs.flatten()
    
    fig_nm = salem_srv_1d_figs[i].split('/')[-1]
    fig_time = salem_srv_1d_figs[i].split('/')[-1][:-4]
    
    Image1 = plt.imread(salem_srv_1d_figs[i])
    axs[0].imshow(Image1[90:,:, :])
    axs[0].axis('off')
    axs[0].set_title('Salem SRV 1D',size=13 , x=0.35)
    
    Image2 = plt.imread(scottsmills_srv_1d_figs[i])
    axs[1].imshow(Image2[90:,:, :])
    axs[1].axis('off')
    axs[1].set_title('Scotts Mills SRV 1D',size=13 , x=0.35)
    
    Image3 = plt.imread(springfield_srv_1d_figs[i])
    axs[2].imshow(Image3[90:,:, :])
    axs[2].axis('off')
    axs[2].set_title('Springfield SRV 1D',size=13, x=0.35)
    
    Image4 = plt.imread(salem_steph_figs[i])
    axs[3].imshow(Image4[90:,:, :])
    axs[3].axis('off')
    axs[3].set_title('Salem CVM',size=13 , x=0.35)
    
    Image5 = plt.imread(scottsmills_steph_figs[i])
    axs[4].imshow(Image5[90:,:, :])
    axs[4].axis('off')
    axs[4].set_title('Scotts Mills CVM',size=13 , x=0.35)
    
    Image6 = plt.imread(springfield_steph_figs[i])
    axs[5].imshow(Image6[90:,:, :])
    axs[5].axis('off')
    axs[5].set_title('Springfield CVM',size=13, x=0.35)
    
    fig.text(0.135,0.82,fig_time,fontsize=22,weight='bold')
    
    plt.savefig(wp_prog_dir + fig_nm ,dpi=300,bbox_inches='tight')

#%%
'''
Grouping two selected timesteps (14s, 25s) wave propagation maps (made in the section above).
This is FIGURE 7 in the publication.
'''

wp_prog_all_eqs = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/wave_propagation/*.png'))

wp_prog_all_eqs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/wave_prop_14_25/'

fig,axs = plt.subplots(figsize=(8, 6),nrows=2, ncols=1)
fig.subplots_adjust(wspace=0,hspace=0.05)
axs = axs.flatten()

Image1 = plt.imread(wp_prog_all_eqs[28])
axs[0].imshow(Image1)
axs[0].axis('off')

Image2 = plt.imread(wp_prog_all_eqs[50])
axs[1].imshow(Image2)
axs[1].axis('off')

fig.text(0.332,0.867,'a)',fontsize=8,weight='bold')
fig.text(0.440,0.870,'b)',fontsize=8,weight='bold')
fig.text(0.559,0.870,'c)',fontsize=8,weight='bold')
fig.text(0.332,0.68,'d)',fontsize=8,weight='bold')
fig.text(0.440,0.68,'e)',fontsize=8,weight='bold')
fig.text(0.563,0.68,'f)',fontsize=8,weight='bold')
fig.text(0.332,0.479,'g)',fontsize=8,weight='bold')
fig.text(0.440,0.481,'h)',fontsize=8,weight='bold')
fig.text(0.562,0.481,'i)',fontsize=8,weight='bold')
fig.text(0.332,0.293,'j)',fontsize=8,weight='bold')
fig.text(0.442,0.293,'k)',fontsize=8,weight='bold')
fig.text(0.564,0.292,'l)',fontsize=8,weight='bold')

plt.savefig(wp_prog_all_eqs_dir + 'wave_prop_14_25.jpg' ,dpi=300,bbox_inches='tight')

#%%
'''
Plotting three figures of selected waveform and spectra comparisons grouped.
This is FIGURES 9, 10, 11 in the publication.
'''

wf_spec_e_salem = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_station_E_comp/salem/srv_1d/*'))
#ALVY - 3, MONO - 55, WLOO - 101, SAIL - 82, NOMA - 62, COBRA - 20, EYES - 29
wf_spec_e_scottsmills = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_station_E_comp/scottsmills/srv_1d/*'))
#ALVY - 5, MONO - 46, WLOO - 86, SAIL - 75
wf_spec_e_springfield = sorted(glob('/Users/rshimony/Desktop/WillametteValley/wv_project/wf_spec_comparisons/wf_spec_station_E_comp/springfield/srv_1d/*'))
#ALVY - 0, MONO - 30

figs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/paper_wf_spec_comp/'

### NOMA, COBRA, EYES - Salem
fig,axs = plt.subplots(figsize=(10, 6),nrows=3, ncols=1)
fig.subplots_adjust(wspace=0,hspace=0)
axs = axs.flatten()

Image1 = plt.imread(wf_spec_e_salem[20])
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].text(-10,90,'a)',size=12,weight='bold')

Image2 = plt.imread(wf_spec_e_salem[29])
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].text(-10,90,'b)',size=12,weight='bold')

Image3 = plt.imread(wf_spec_e_salem[62])
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].text(-10,90,'c)',size=12,weight='bold')

plt.savefig(figs_dir + 'noma_cobra_eyes.jpg' ,dpi=300,bbox_inches='tight')

### WLOO, SAIL - Salem, Scotts Mills
fig,axs = plt.subplots(figsize=(14, 6),nrows=2, ncols=2)
fig.subplots_adjust(wspace=0,hspace=0)
axs = axs.flatten()

Image1 = plt.imread(wf_spec_e_salem[101])
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].text(-10,90,'a)',size=14,weight='bold')

Image2 = plt.imread(wf_spec_e_salem[82])
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].text(-10,90,'b)',size=14,weight='bold')

Image3 = plt.imread(wf_spec_e_scottsmills[84])
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].text(-10,90,'c)',size=14,weight='bold')

Image4 = plt.imread(wf_spec_e_scottsmills[73])
axs[3].imshow(Image4)
axs[3].axis('off')
axs[3].text(-10,90,'d)',size=14,weight='bold')

plt.savefig(figs_dir + 'wloo_sail.jpg' ,dpi=300,bbox_inches='tight')

### ALVY, MONO - Salem, Scotts Mills, Springfield
fig,axs = plt.subplots(figsize=(7.5, 5),nrows=3, ncols=2)
fig.subplots_adjust(wspace=0,hspace=0)
axs = axs.flatten()

Image1 = plt.imread(wf_spec_e_salem[55])
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].text(-10,90,'a)',size=10,weight='bold')

Image2 = plt.imread(wf_spec_e_salem[3])
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].text(-10,90,'b)',size=10,weight='bold')

Image3 = plt.imread(wf_spec_e_scottsmills[44])
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].text(-10,90,'c)',size=10,weight='bold')

Image4 = plt.imread(wf_spec_e_scottsmills[4])
axs[3].imshow(Image4)
axs[3].axis('off')
axs[3].text(-10,90,'d)',size=10,weight='bold')

Image5 = plt.imread(wf_spec_e_springfield[30])
axs[4].imshow(Image5)
axs[4].axis('off')
axs[4].text(-10,90,'e)',size=10,weight='bold')

Image6 = plt.imread(wf_spec_e_springfield[0])
axs[5].imshow(Image6)
axs[5].axis('off')
axs[5].text(-10,90,'f)',size=10,weight='bold')

plt.savefig(figs_dir + 'alvy_mono.jpg' ,dpi=300,bbox_inches='tight')

#%%
'''
Making BAF summary figure.
Grouping data maps used in the OpenQuake BAF process (Vs30 and Z2.5) with the OpenQuake BAF map (top row)
with BAF maps from the simulations of the three EQs (Salme, Scotts Mills, Springfield) (bottom row).
This is FIGURE 12 in the publication.
'''

vs30_map_fig = '/Users/rshimony/Desktop/WillametteValley/wv_project/vs30/vs30_map.png'
z2pt5_map_fig = '/Users/rshimony/Desktop/WillametteValley/wv_project/z_surfaces/z2pt5_valley_map.png'
shaking_ratio_map_fig = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/baf_shakemaps/openquake/z2pt5_full_valley_vals_ratio.png'

baf_fig_salem = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/baf_shakemaps/sw4_synts/Salem_baf.png'
baf_fig_scottsmills = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/baf_shakemaps/sw4_synts/Scotts_Mills_baf.png'
baf_fig_springfield = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/baf_shakemaps/sw4_synts/Springfield_baf.png'


fig,axs = plt.subplots(figsize=(15, 11),nrows=2, ncols=3)
fig.subplots_adjust(wspace=0.,hspace=0.)
axs = axs.flatten()

title_fontsize = 15

Image1 = plt.imread(z2pt5_map_fig)
axs[0].imshow(Image1)
axs[0].axis('off')
axs[0].text(0,150,'a)',size=title_fontsize+10,weight='bold')

Image2 = plt.imread(vs30_map_fig)
axs[1].imshow(Image2)
axs[1].axis('off')
axs[1].text(0,150,'b)',size=title_fontsize+10,weight='bold')

Image3 = plt.imread(shaking_ratio_map_fig)
axs[2].imshow(Image3)
axs[2].axis('off')
axs[2].text(0,150,'c)',size=title_fontsize+10,weight='bold')

Image4 = plt.imread(baf_fig_salem)
axs[3].imshow(Image4)
axs[3].axis('off')
axs[3].text(30,90,'d)',size=title_fontsize+10,weight='bold')

Image5 = plt.imread(baf_fig_scottsmills)
axs[4].imshow(Image5)
axs[4].axis('off')
axs[4].text(30,90,'e)',size=title_fontsize+10,weight='bold')

Image6 = plt.imread(baf_fig_springfield)
axs[5].imshow(Image6)
axs[5].axis('off')
axs[5].text(30,90,'f)',size=title_fontsize+10,weight='bold')

plt.savefig('/Users/rshimony/Desktop/WillametteValley/wv_project/figures/final_subplots/BAF_fig/BAF_fig.png' ,dpi=300,bbox_inches='tight')








