#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:20:45 2024

@author: rshimony
"""
'''
Creates boxplots of residual values distribution in the WV between all models (4 new WV models and USGS CVM) for a SINGLE EVENT (only Valley stations).
Creates boxplots of residual ratio values distribution between Valley and Mountain (outside of valley) stations, for each model in each event.
Creates STACKED boxplot figure, that compares between all models and all events. This is FIGURE 4 in the publication.
Inputs:
    single model residuals and ratios flatfiles
    All models residuals and ratios flatfiles
    Melted residuals and ratios flatfiles of all events
'''

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#%%
###INPUTS###
## Single model residuals and ratios flatfiles
eugspe_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/eugspe_1d_res_ff.csv')
eugspe_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/eugspe_sm_res_ff.csv')
srv_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/srv_1d_res_ff.csv')
srv_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/srv_sm_res_ff.csv')

## All models residuals and ratios flatfiles - Valley station directory
all_model_val_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/all_model_val_ff.csv')

## Melted residuals and ratios flatfiles of all events
salem_melt_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/salem/melt_val_ff.csv')
scottsmills_melt_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/scottsmills/melt_val_ff.csv')
springfield_melt_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/wv_project/residual_ratio_analysis/springfield/melt_val_ff.csv')

## Ouput directories
# Stacked boxplot
stack_boxplot_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/boxplots/'
# Residual values distribution
boxplots_figs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/boxplots/res_val_dist/springfield/'
# Ratio values distribution
val_boxplot_figs_dir = '/Users/rshimony/Desktop/WillametteValley/wv_project/figures/boxplots/ratio_inout_val_dist/springfield/'  
    
#%%

def make_residuals_boxplot(data_ff , im , figs_dir , fig_name , fig_title):
    plt.figure()
    sns.set(font_scale=1)
    boxplot = sns.boxplot(x='model', y=im, data=data_ff, showfliers=False)
    boxplot.set_title(fig_title,size=15)
    boxplot.set_xticklabels(['EugSpe 1D','EugSpe SM','SRV 1D','SRV SM','USGS CVM'])
    boxplot.axhline(y=0, color='magenta', linestyle='--',linewidth=1)
    boxplot.set_xlabel("Model",fontsize=14)
    boxplot.set_ylabel("Residual (Observed - Predicted)",fontsize=14)
    plt.savefig(figs_dir + fig_name + '.png',dpi=150)
    
def make_ratios_boxplot(data_ff , im , figs_dir , fig_name):
    plt.figure(figsize=(6.2,3))
    boxplot = sns.boxplot(x='location', y=im, data=data_ff, showfliers=False)
    boxplot.set_xticklabels(['Valley Stations','Mountain Stations'])
    boxplot.tick_params(labelsize=23)
    boxplot.set_xlabel("Stations Location",fontsize=26)
    boxplot.set_ylabel("Residual Ratio",fontsize=26)
    boxplot.axhline(y=1, color='magenta', linestyle='--',linewidth=2)
    plt.savefig(figs_dir + fig_name + '.png',dpi=150,bbox_inches='tight')
    
    
#%%

make_residuals_boxplot(all_model_val_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res_val' , 'f = 2Hz')
make_residuals_boxplot(all_model_val_ff , 'f3_fas_res' , boxplots_figs_dir , 'f3_fas_res_val' , 'f = 3Hz')
make_residuals_boxplot(all_model_val_ff , 'f4_fas_res' , boxplots_figs_dir , 'f4_fas_res_val' , 'f = 4Hz')
make_residuals_boxplot(all_model_val_ff , 'f5_fas_res' , boxplots_figs_dir , 'f5_fas_res_val' , 'f = 5Hz')
make_residuals_boxplot(all_model_val_ff , 'f6_fas_res' , boxplots_figs_dir , 'f6_fas_res_val' , 'f = 6Hz')
make_residuals_boxplot(all_model_val_ff , 'f7_fas_res' , boxplots_figs_dir , 'f7_fas_res_val' , 'f = 7Hz')
make_residuals_boxplot(all_model_val_ff , 'pgv_res' , boxplots_figs_dir , 'pgv_res_val' , 'PGV')

#%%

make_ratios_boxplot(eugspe_1d_res_ff , 'f2_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f2_fas_eugspe_1d')
make_ratios_boxplot(eugspe_1d_res_ff , 'f3_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f3_fas_eugspe_1d')
make_ratios_boxplot(eugspe_1d_res_ff , 'f4_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f4_fas_eugspe_1d')
make_ratios_boxplot(eugspe_1d_res_ff , 'f5_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f5_fas_eugspe_1d')
make_ratios_boxplot(eugspe_1d_res_ff , 'f6_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f6_fas_eugspe_1d')
make_ratios_boxplot(eugspe_1d_res_ff , 'f7_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f7_fas_eugspe_1d')
make_ratios_boxplot(eugspe_1d_res_ff , 'pgv_eugspe_1d_ratio' , val_boxplot_figs_dir , 'pgv_eugspe_1d')

make_ratios_boxplot(eugspe_sm_res_ff , 'f2_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f2_fas_eugspe_sm')
make_ratios_boxplot(eugspe_sm_res_ff , 'f3_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f3_fas_eugspe_sm')
make_ratios_boxplot(eugspe_sm_res_ff , 'f4_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f4_fas_eugspe_sm')
make_ratios_boxplot(eugspe_sm_res_ff , 'f5_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f5_fas_eugspe_sm')
make_ratios_boxplot(eugspe_sm_res_ff , 'f6_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f6_fas_eugspe_sm')
make_ratios_boxplot(eugspe_sm_res_ff , 'f7_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f7_fas_eugspe_sm')
make_ratios_boxplot(eugspe_sm_res_ff , 'pgv_eugspe_sm_ratio' , val_boxplot_figs_dir , 'pgv_eugspe_sm')

make_ratios_boxplot(srv_1d_res_ff , 'f2_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f2_fas_srv_1d')
make_ratios_boxplot(srv_1d_res_ff , 'f3_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f3_fas_srv_1d')
make_ratios_boxplot(srv_1d_res_ff , 'f4_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f4_fas_srv_1d')
make_ratios_boxplot(srv_1d_res_ff , 'f5_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f5_fas_srv_1d')
make_ratios_boxplot(srv_1d_res_ff , 'f6_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f6_fas_srv_1d')
make_ratios_boxplot(srv_1d_res_ff , 'f7_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f7_fas_srv_1d')
make_ratios_boxplot(srv_1d_res_ff , 'pgv_srv_1d_ratio' , val_boxplot_figs_dir , 'pgv_srv_1d')

make_ratios_boxplot(srv_sm_res_ff , 'f2_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f2_fas_srv_sm')
make_ratios_boxplot(srv_sm_res_ff , 'f3_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f3_fas_srv_sm')
make_ratios_boxplot(srv_sm_res_ff , 'f4_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f4_fas_srv_sm')
make_ratios_boxplot(srv_sm_res_ff , 'f5_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f5_fas_srv_sm')
make_ratios_boxplot(srv_sm_res_ff , 'f6_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f6_fas_srv_sm')
make_ratios_boxplot(srv_sm_res_ff , 'f7_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f7_fas_srv_sm')
make_ratios_boxplot(srv_sm_res_ff , 'pgv_srv_sm_ratio' , val_boxplot_figs_dir , 'pgv_srv_sm')

#%%

salem_single_model = salem_melt_ff[salem_melt_ff['model']=='eugspe_1d']
salem_single_im = salem_single_model[salem_single_model['res_im']=='pgv_res']
n_stations_salem = len(salem_single_im)

scottsmills_single_model = scottsmills_melt_ff[scottsmills_melt_ff['model']=='steph']
scottsmills_single_im = scottsmills_single_model[scottsmills_single_model['res_im']=='f2_fas_res']
n_stations_scottsmills = len(scottsmills_single_im)

springfield_single_model = springfield_melt_ff[springfield_melt_ff['model']=='steph']
springfield_single_im = springfield_single_model[springfield_single_model['res_im']=='f2_fas_res']
n_stations_springfield = len(springfield_single_im)

fig,axs = plt.subplots(figsize=(28, 18),nrows=2, ncols=2)
fig.subplots_adjust(hspace=0.3)
axs = axs.flatten()
sns.set(font_scale=3.2)

sns.boxplot(ax=axs[0],data = salem_melt_ff , x='res_im' , y='res_val' ,hue='model', showfliers=False)
axs[0].set_title('Salem',size=30)
axs[0].set_xticklabels(['f=0.2Hz','f=0.3Hz','f=0.4Hz','f=0.5Hz','f=0.6Hz','f=0.7Hz','PGV'])
axs[0].tick_params(labelsize=26)
axs[0].set_xlabel("Intensity Measure",fontsize=28)
axs[0].set_ylabel("Residual (Observed - Predicted)",fontsize=28)
axs[0].legend_.remove()
axs[0].axhline(y=0, color='magenta', linestyle='--',linewidth=2)
axs[0].text(0,3.1,'N = '+str(n_stations_salem),color='red',weight='bold')
axs[0].text(-0.8,3.7,'a)',size=43,weight='bold')

sns.boxplot(ax=axs[1],data = scottsmills_melt_ff , x='res_im' , y='res_val' ,hue='model', showfliers=False)
axs[1].set_title('Scotts Mills',size=30)
axs[1].set_xticklabels(['f=0.2Hz','f=0.3Hz','f=0.4Hz','f=0.5Hz','f=0.6Hz','f=0.7Hz','PGV'])
axs[1].tick_params(labelsize=26)
axs[1].set_xlabel("Intensity Measure",fontsize=28)
axs[1].set_ylabel("Residual (Observed - Predicted)",fontsize=28)
axs[1].legend_.remove()
axs[1].axhline(y=0, color='magenta', linestyle='--',linewidth=2)
axs[1].text(5,4.6,'N = '+str(n_stations_scottsmills),color='red',weight='bold')
axs[1].text(-0.8,5.5,'b)',size=43,weight='bold')

sns.boxplot(ax=axs[2],data = springfield_melt_ff , x='res_im' , y='res_val' ,hue='model', showfliers=False)
axs[2].set_title('Springfield',size=30)
axs[2].set_xticklabels(['f=0.2Hz','f=0.3Hz','f=0.4Hz','f=0.5Hz','f=0.6Hz','f=0.7Hz','PGV'])
axs[2].legend(handles=axs[2].legend_.legendHandles ,labels=["EugSpe 1D", "EugSpe SM", "SRV 1D", "SRV SM", "USGS CVM"],title = 'Model')
axs[2].tick_params(labelsize=26)
axs[2].set_xlabel("Intensity Measure",fontsize=28)
axs[2].set_ylabel("Residual (Observed - Predicted)",fontsize=28)
sns.move_legend(axs[2], "center right", bbox_to_anchor=(1.8, 0.5))
axs[2].axhline(y=0, color='magenta', linestyle='--',linewidth=2)
axs[2].text(5,3.3,'N = '+str(n_stations_springfield),color='red',weight='bold')
axs[2].text(-0.8,3.9,'c)',size=43,weight='bold')

axs[3].axis('off')

plt.savefig(stack_boxplot_dir + 'all_eq_stack_boxplots.jpg',dpi=300,bbox_inches='tight')
    
#%%

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

























