#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 14:39:45 2024

@author: rshimony
"""

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

#%%

eugspe_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/eugspe_1d_res_ff.csv')
eugspe_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/eugspe_sm_res_ff.csv')
srv_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/srv_1d_res_ff.csv')
srv_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/srv_sm_res_ff.csv')

eugspe_1d_only_rat = eugspe_1d_res_ff.iloc[:,13:20]
eugspe_sm_only_rat = eugspe_sm_res_ff.iloc[:,13:20]
srv_1d_only_rat = srv_1d_res_ff.iloc[:,13:20]
srv_sm_only_rat = srv_sm_res_ff.iloc[:,13:20]

valst_eugspe_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/eugspe_1d_res_ff_val.csv')
valst_eugspe_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/eugspe_sm_res_ff_val.csv')
valst_srv_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/srv_1d_res_ff_val.csv')
valst_srv_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/srv_sm_res_ff_val.csv')


full_st = srv_1d_res_ff['st_name']
valley_st = valst_srv_1d_res_ff['st_name_val']
valley_st_mask = np.isin(full_st,valley_st)

valley_ls = []
for i in valley_st_mask:
    if i == True:
        valley_ls.append('in_valley')
    else:
        valley_ls.append('out_valley')
valley_array = np.array(valley_ls)

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/scottsmills/'

eugspe_1d_only_rat.insert(0, "st_nm", full_st, True)
eugspe_1d_only_rat.insert(1, "location", valley_ls, True)
eugspe_1d_only_rat.to_csv(resrat_ff_dir+'eugspe_1d_valley_loc_ff.csv')

eugspe_sm_only_rat.insert(0, "st_nm", full_st, True)
eugspe_sm_only_rat.insert(1, "location", valley_ls, True)
eugspe_sm_only_rat.to_csv(resrat_ff_dir+'eugspe_sm_valley_loc_ff.csv')

srv_1d_only_rat.insert(0, "st_nm", full_st, True)
srv_1d_only_rat.insert(1, "location", valley_ls, True)
srv_1d_only_rat.to_csv(resrat_ff_dir+'srv_1d_valley_loc_ff.csv')

srv_sm_only_rat.insert(0, "st_nm", full_st, True)
srv_sm_only_rat.insert(1, "location", valley_ls, True)
srv_sm_only_rat.to_csv(resrat_ff_dir+'srv_sm_valley_loc_ff.csv')

#%%
# ,binrange=(0,5)
def make_multi_hist(data_ff , im , hist_type , figs_dir , fig_name , hue):
    if hist_type == 'reg':
        plt.figure()
        sns.set(font_scale=1)
        hist = sns.histplot(data=data_ff, x=im, bins='auto',hue=hue)
        sns.move_legend(hist, 'upper right')
        hist.set_title(fig_name,size=15)
        plt.savefig(figs_dir + fig_name +'_'+ hist_type + '.png',dpi=150)
        
    if hist_type == 'stack':
        plt.figure()
        sns.set(font_scale=1)
        hist = sns.histplot(data=data_ff, x=im, bins='auto',hue=hue,multiple="stack")
        sns.move_legend(hist, 'upper right')
        hist.set_title(fig_name,size=15)
        plt.savefig(figs_dir + fig_name +'_'+ hist_type + '.png',dpi=150)
        
    if hist_type == 'step':
        plt.figure()
        sns.set(font_scale=1)
        hist = sns.histplot(data=data_ff, x=im, bins='auto',hue=hue,element='step')
        sns.move_legend(hist, 'upper right')
        hist.set_title(fig_name,size=15)
        plt.savefig(figs_dir + fig_name +'_'+ hist_type + '.png',dpi=150)
        
val_hist_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ratio_valley_hist/scottsmills/'

make_multi_hist(eugspe_1d_only_rat , 'f2_fas_eugspe_1d_ratio' , 'reg' , val_hist_figs_dir , 'f2_fas_eugspe_1d' ,'location')
make_multi_hist(eugspe_1d_only_rat , 'f3_fas_eugspe_1d_ratio' , 'reg' , val_hist_figs_dir , 'f3_fas_eugspe_1d' ,'location')
make_multi_hist(eugspe_1d_only_rat , 'f4_fas_eugspe_1d_ratio' , 'reg' , val_hist_figs_dir , 'f4_fas_eugspe_1d' ,'location')
make_multi_hist(eugspe_1d_only_rat , 'f5_fas_eugspe_1d_ratio' , 'reg' , val_hist_figs_dir , 'f5_fas_eugspe_1d' ,'location')
make_multi_hist(eugspe_1d_only_rat , 'f6_fas_eugspe_1d_ratio' , 'reg' , val_hist_figs_dir , 'f6_fas_eugspe_1d' ,'location')
make_multi_hist(eugspe_1d_only_rat , 'f7_fas_eugspe_1d_ratio' , 'reg' , val_hist_figs_dir , 'f7_fas_eugspe_1d' ,'location')
make_multi_hist(eugspe_1d_only_rat , 'pgv_eugspe_1d_ratio' , 'reg' , val_hist_figs_dir , 'pgv_eugspe_1d' ,'location')

make_multi_hist(eugspe_sm_only_rat , 'f2_fas_eugspe_sm_ratio' , 'reg' , val_hist_figs_dir , 'f2_fas_eugspe_sm' ,'location')
make_multi_hist(eugspe_sm_only_rat , 'f3_fas_eugspe_sm_ratio' , 'reg' , val_hist_figs_dir , 'f3_fas_eugspe_sm' ,'location')
make_multi_hist(eugspe_sm_only_rat , 'f4_fas_eugspe_sm_ratio' , 'reg' , val_hist_figs_dir , 'f4_fas_eugspe_sm' ,'location')
make_multi_hist(eugspe_sm_only_rat , 'f5_fas_eugspe_sm_ratio' , 'reg' , val_hist_figs_dir , 'f5_fas_eugspe_sm' ,'location')
make_multi_hist(eugspe_sm_only_rat , 'f6_fas_eugspe_sm_ratio' , 'reg' , val_hist_figs_dir , 'f6_fas_eugspe_sm' ,'location')
make_multi_hist(eugspe_sm_only_rat , 'f7_fas_eugspe_sm_ratio' , 'reg' , val_hist_figs_dir , 'f7_fas_eugspe_sm' ,'location')
make_multi_hist(eugspe_sm_only_rat , 'pgv_eugspe_sm_ratio' , 'reg' , val_hist_figs_dir , 'pgv_eugspe_sm' ,'location')

make_multi_hist(srv_1d_only_rat , 'f2_fas_srv_1d_ratio' , 'reg' , val_hist_figs_dir , 'f2_fas_srv_1d' ,'location')
make_multi_hist(srv_1d_only_rat , 'f3_fas_srv_1d_ratio' , 'reg' , val_hist_figs_dir , 'f3_fas_srv_1d' ,'location')
make_multi_hist(srv_1d_only_rat , 'f4_fas_srv_1d_ratio' , 'reg' , val_hist_figs_dir , 'f4_fas_srv_1d' ,'location')
make_multi_hist(srv_1d_only_rat , 'f5_fas_srv_1d_ratio' , 'reg' , val_hist_figs_dir , 'f5_fas_srv_1d' ,'location')
make_multi_hist(srv_1d_only_rat , 'f6_fas_srv_1d_ratio' , 'reg' , val_hist_figs_dir , 'f6_fas_srv_1d' ,'location')
make_multi_hist(srv_1d_only_rat , 'f7_fas_srv_1d_ratio' , 'reg' , val_hist_figs_dir , 'f7_fas_srv_1d' ,'location')
make_multi_hist(srv_1d_only_rat , 'pgv_srv_1d_ratio' , 'reg' , val_hist_figs_dir , 'pgv_srv_1d' ,'location')

make_multi_hist(srv_sm_only_rat , 'f2_fas_srv_sm_ratio' , 'reg' , val_hist_figs_dir , 'f2_fas_srv_sm' ,'location')
make_multi_hist(srv_sm_only_rat , 'f3_fas_srv_sm_ratio' , 'reg' , val_hist_figs_dir , 'f3_fas_srv_sm' ,'location')
make_multi_hist(srv_sm_only_rat , 'f4_fas_srv_sm_ratio' , 'reg' , val_hist_figs_dir , 'f4_fas_srv_sm' ,'location')
make_multi_hist(srv_sm_only_rat , 'f5_fas_srv_sm_ratio' , 'reg' , val_hist_figs_dir , 'f5_fas_srv_sm' ,'location')
make_multi_hist(srv_sm_only_rat , 'f6_fas_srv_sm_ratio' , 'reg' , val_hist_figs_dir , 'f6_fas_srv_sm' ,'location')
make_multi_hist(srv_sm_only_rat , 'f7_fas_srv_sm_ratio' , 'reg' , val_hist_figs_dir , 'f7_fas_srv_sm' ,'location')
make_multi_hist(srv_sm_only_rat , 'pgv_srv_sm_ratio' , 'reg' , val_hist_figs_dir , 'pgv_srv_sm' ,'location')

#%%

def make_boxplot_ims(data_ff , im , figs_dir , fig_name):
    plt.figure()
    sns.set(font_scale=1)
    boxplot = sns.boxplot(x='location', y=im, data=data_ff, showfliers=False)
    boxplot.set_title(fig_name,size=15)
    plt.savefig(figs_dir + fig_name + '.png',dpi=150)
        
val_boxplot_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/ratio_boxplot_valley/scottsmills/'   

make_boxplot_ims(eugspe_1d_only_rat , 'f2_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f2_fas_eugspe_1d')

make_boxplot_ims(eugspe_1d_only_rat , 'f2_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f2_fas_eugspe_1d')
make_boxplot_ims(eugspe_1d_only_rat , 'f3_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f3_fas_eugspe_1d')
make_boxplot_ims(eugspe_1d_only_rat , 'f4_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f4_fas_eugspe_1d')
make_boxplot_ims(eugspe_1d_only_rat , 'f5_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f5_fas_eugspe_1d')
make_boxplot_ims(eugspe_1d_only_rat , 'f6_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f6_fas_eugspe_1d')
make_boxplot_ims(eugspe_1d_only_rat , 'f7_fas_eugspe_1d_ratio' , val_boxplot_figs_dir , 'f7_fas_eugspe_1d')
make_boxplot_ims(eugspe_1d_only_rat , 'pgv_eugspe_1d_ratio' , val_boxplot_figs_dir , 'pgv_eugspe_1d')

make_boxplot_ims(eugspe_sm_only_rat , 'f2_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f2_fas_eugspe_sm')
make_boxplot_ims(eugspe_sm_only_rat , 'f3_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f3_fas_eugspe_sm')
make_boxplot_ims(eugspe_sm_only_rat , 'f4_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f4_fas_eugspe_sm')
make_boxplot_ims(eugspe_sm_only_rat , 'f5_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f5_fas_eugspe_sm')
make_boxplot_ims(eugspe_sm_only_rat , 'f6_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f6_fas_eugspe_sm')
make_boxplot_ims(eugspe_sm_only_rat , 'f7_fas_eugspe_sm_ratio' , val_boxplot_figs_dir , 'f7_fas_eugspe_sm')
make_boxplot_ims(eugspe_sm_only_rat , 'pgv_eugspe_sm_ratio' , val_boxplot_figs_dir , 'pgv_eugspe_sm')

make_boxplot_ims(srv_1d_only_rat , 'f2_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f2_fas_srv_1d')
make_boxplot_ims(srv_1d_only_rat , 'f3_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f3_fas_srv_1d')
make_boxplot_ims(srv_1d_only_rat , 'f4_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f4_fas_srv_1d')
make_boxplot_ims(srv_1d_only_rat , 'f5_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f5_fas_srv_1d')
make_boxplot_ims(srv_1d_only_rat , 'f6_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f6_fas_srv_1d')
make_boxplot_ims(srv_1d_only_rat , 'f7_fas_srv_1d_ratio' , val_boxplot_figs_dir , 'f7_fas_srv_1d')
make_boxplot_ims(srv_1d_only_rat , 'pgv_srv_1d_ratio' , val_boxplot_figs_dir , 'pgv_srv_1d')

make_boxplot_ims(srv_sm_only_rat , 'f2_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f2_fas_srv_sm')
make_boxplot_ims(srv_sm_only_rat , 'f3_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f3_fas_srv_sm')
make_boxplot_ims(srv_sm_only_rat , 'f4_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f4_fas_srv_sm')
make_boxplot_ims(srv_sm_only_rat , 'f5_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f5_fas_srv_sm')
make_boxplot_ims(srv_sm_only_rat , 'f6_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f6_fas_srv_sm')
make_boxplot_ims(srv_sm_only_rat , 'f7_fas_srv_sm_ratio' , val_boxplot_figs_dir , 'f7_fas_srv_sm')
make_boxplot_ims(srv_sm_only_rat , 'pgv_srv_sm_ratio' , val_boxplot_figs_dir , 'pgv_srv_sm')































