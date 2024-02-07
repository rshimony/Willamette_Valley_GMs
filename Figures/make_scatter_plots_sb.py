#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  6 17:42:29 2024

@author: rshimony
"""

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.ticker import MaxNLocator

#%%

steph_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/steph_res_ff.csv')
eugspe_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_1d_res_ff.csv')
eugspe_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_sm_res_ff.csv')
srv_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_1d_res_ff.csv')
srv_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_sm_res_ff.csv')

valst_steph_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/steph_res_ff_val.csv')
valst_eugspe_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_1d_res_ff_val.csv')
valst_eugspe_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_sm_res_ff_val.csv')
valst_srv_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_1d_res_ff_val.csv')
valst_srv_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_sm_res_ff_val.csv')


#%%

def plot_jointplot_res(data_ff , res_syn , res_steph , outdir , fig_name):
    
    plt.figure()
    sns.set(font_scale=2.3)
    
    jp = sns.jointplot(data=data_ff, x=res_syn, y=res_steph , hue='lat_region',height=10, ratio=2)
    
    jp.ax_marg_x.set_xlim(-3, 3)
    jp.ax_marg_y.set_ylim(-3, 3)
    
    jp.ax_marg_y.tick_params(labeltop=True,grid_linewidth=2.5 , labelsize=28, rotation=70)
    jp.ax_marg_y.grid(True, axis='x', ls=':')
    jp.ax_marg_y.xaxis.set_major_locator(MaxNLocator(4))

    jp.ax_marg_x.tick_params(labelleft=True,grid_linewidth=2.5 , labelsize=28)
    jp.ax_marg_x.grid(True, axis='y', ls=':')
    jp.ax_marg_x.yaxis.set_major_locator(MaxNLocator(4))
    
    jp.ax_joint.tick_params(labelsize=32)
    jp.ax_joint.set_ylabel(res_steph,fontsize=32)
    jp.ax_joint.set_xlabel(res_syn,fontsize=32)
    
    jp.ax_joint.plot([-3,3], [-3,3], 'b-', linewidth = 1.5)
    
    jp.refline(x=0, y=0)
    
    plt.savefig(outdir + fig_name + '.png',dpi=150)



jointplot_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/scatter_jointplots/springfield/'
try:
    os.mkdir(jointplot_figs_dir)
except FileExistsError:
    print('Directory already exists')
    
plot_jointplot_res(eugspe_1d_res_ff , 'f2_fas_eugspe_1d_res' , 'f2_fas_steph_res' , jointplot_figs_dir , 'f2_fas_eugspe_1d_res')
plot_jointplot_res(eugspe_1d_res_ff , 'f3_fas_eugspe_1d_res' , 'f3_fas_steph_res' , jointplot_figs_dir , 'f3_fas_eugspe_1d_res')
plot_jointplot_res(eugspe_1d_res_ff , 'f4_fas_eugspe_1d_res' , 'f4_fas_steph_res' , jointplot_figs_dir , 'f4_fas_eugspe_1d_res')
plot_jointplot_res(eugspe_1d_res_ff , 'f5_fas_eugspe_1d_res' , 'f5_fas_steph_res' , jointplot_figs_dir , 'f5_fas_eugspe_1d_res')
plot_jointplot_res(eugspe_1d_res_ff , 'f6_fas_eugspe_1d_res' , 'f6_fas_steph_res' , jointplot_figs_dir , 'f6_fas_eugspe_1d_res')
plot_jointplot_res(eugspe_1d_res_ff , 'f7_fas_eugspe_1d_res' , 'f7_fas_steph_res' , jointplot_figs_dir , 'f7_fas_eugspe_1d_res')
plot_jointplot_res(eugspe_1d_res_ff , 'pgv_eugspe_1d_res' , 'pgv_steph_res' , jointplot_figs_dir , 'pgv_eugspe_1d_res')

plot_jointplot_res(eugspe_sm_res_ff , 'f2_fas_eugspe_sm_res' , 'f2_fas_steph_res' , jointplot_figs_dir , 'f2_fas_eugspe_sm_res')
plot_jointplot_res(eugspe_sm_res_ff , 'f3_fas_eugspe_sm_res' , 'f3_fas_steph_res' , jointplot_figs_dir , 'f3_fas_eugspe_sm_res')
plot_jointplot_res(eugspe_sm_res_ff , 'f4_fas_eugspe_sm_res' , 'f4_fas_steph_res' , jointplot_figs_dir , 'f4_fas_eugspe_sm_res')
plot_jointplot_res(eugspe_sm_res_ff , 'f5_fas_eugspe_sm_res' , 'f5_fas_steph_res' , jointplot_figs_dir , 'f5_fas_eugspe_sm_res')
plot_jointplot_res(eugspe_sm_res_ff , 'f6_fas_eugspe_sm_res' , 'f6_fas_steph_res' , jointplot_figs_dir , 'f6_fas_eugspe_sm_res')
plot_jointplot_res(eugspe_sm_res_ff , 'f7_fas_eugspe_sm_res' , 'f7_fas_steph_res' , jointplot_figs_dir , 'f7_fas_eugspe_sm_res')
plot_jointplot_res(eugspe_sm_res_ff , 'pgv_eugspe_sm_res' , 'pgv_steph_res' , jointplot_figs_dir , 'pgv_eugspe_sm_res')    

plot_jointplot_res(srv_sm_res_ff , 'f2_fas_srv_sm_res' , 'f2_fas_steph_res' , jointplot_figs_dir , 'f2_fas_srv_sm_res')
plot_jointplot_res(srv_sm_res_ff , 'f3_fas_srv_sm_res' , 'f3_fas_steph_res' , jointplot_figs_dir , 'f3_fas_srv_sm_res')
plot_jointplot_res(srv_sm_res_ff , 'f4_fas_srv_sm_res' , 'f4_fas_steph_res' , jointplot_figs_dir , 'f4_fas_srv_sm_res')
plot_jointplot_res(srv_sm_res_ff , 'f5_fas_srv_sm_res' , 'f5_fas_steph_res' , jointplot_figs_dir , 'f5_fas_srv_sm_res')
plot_jointplot_res(srv_sm_res_ff , 'f6_fas_srv_sm_res' , 'f6_fas_steph_res' , jointplot_figs_dir , 'f6_fas_srv_sm_res')
plot_jointplot_res(srv_sm_res_ff , 'f7_fas_srv_sm_res' , 'f7_fas_steph_res' , jointplot_figs_dir , 'f7_fas_srv_sm_res')
plot_jointplot_res(srv_sm_res_ff , 'pgv_srv_sm_res' , 'pgv_steph_res' , jointplot_figs_dir , 'pgv_srv_sm_res')
    
plot_jointplot_res(srv_1d_res_ff , 'f2_fas_srv_1d_res' , 'f2_fas_steph_res' , jointplot_figs_dir , 'f2_fas_srv_1d_res')
plot_jointplot_res(srv_1d_res_ff , 'f3_fas_srv_1d_res' , 'f3_fas_steph_res' , jointplot_figs_dir , 'f3_fas_srv_1d_res')
plot_jointplot_res(srv_1d_res_ff , 'f4_fas_srv_1d_res' , 'f4_fas_steph_res' , jointplot_figs_dir , 'f4_fas_srv_1d_res')
plot_jointplot_res(srv_1d_res_ff , 'f5_fas_srv_1d_res' , 'f5_fas_steph_res' , jointplot_figs_dir , 'f5_fas_srv_1d_res')
plot_jointplot_res(srv_1d_res_ff , 'f6_fas_srv_1d_res' , 'f6_fas_steph_res' , jointplot_figs_dir , 'f6_fas_srv_1d_res')
plot_jointplot_res(srv_1d_res_ff , 'f7_fas_srv_1d_res' , 'f7_fas_steph_res' , jointplot_figs_dir , 'f7_fas_srv_1d_res')
plot_jointplot_res(srv_1d_res_ff , 'pgv_srv_1d_res' , 'pgv_steph_res' , jointplot_figs_dir , 'pgv_srv_1d_res')

####################
 
plot_jointplot_res(valst_eugspe_1d_res_ff , 'f2_fas_val_eugspe_1d_res' , 'f2_fas_val_steph_res' , jointplot_figs_dir , 'f2_fas_val_eugspe_1d_res')
plot_jointplot_res(valst_eugspe_1d_res_ff , 'f3_fas_val_eugspe_1d_res' , 'f3_fas_val_steph_res' , jointplot_figs_dir , 'f3_fas_val_eugspe_1d_res')
plot_jointplot_res(valst_eugspe_1d_res_ff , 'f4_fas_val_eugspe_1d_res' , 'f4_fas_val_steph_res' , jointplot_figs_dir , 'f4_fas_val_eugspe_1d_res')
plot_jointplot_res(valst_eugspe_1d_res_ff , 'f5_fas_val_eugspe_1d_res' , 'f5_fas_val_steph_res' , jointplot_figs_dir , 'f5_fas_val_eugspe_1d_res')
plot_jointplot_res(valst_eugspe_1d_res_ff , 'f6_fas_val_eugspe_1d_res' , 'f6_fas_val_steph_res' , jointplot_figs_dir , 'f6_fas_val_eugspe_1d_res')
plot_jointplot_res(valst_eugspe_1d_res_ff , 'f7_fas_val_eugspe_1d_res' , 'f7_fas_val_steph_res' , jointplot_figs_dir , 'f7_fas_val_eugspe_1d_res')
plot_jointplot_res(valst_eugspe_1d_res_ff , 'pgv_eugspe_1d_res_val' , 'pgv_steph_res_val' , jointplot_figs_dir , 'pgv_eugspe_1d_res_val')

plot_jointplot_res(valst_eugspe_sm_res_ff , 'f2_fas_val_eugspe_sm_res' , 'f2_fas_val_steph_res' , jointplot_figs_dir , 'f2_fas_val_eugspe_sm_res')
plot_jointplot_res(valst_eugspe_sm_res_ff , 'f3_fas_val_eugspe_sm_res' , 'f3_fas_val_steph_res' , jointplot_figs_dir , 'f3_fas_val_eugspe_sm_res')
plot_jointplot_res(valst_eugspe_sm_res_ff , 'f4_fas_val_eugspe_sm_res' , 'f4_fas_val_steph_res' , jointplot_figs_dir , 'f4_fas_val_eugspe_sm_res')
plot_jointplot_res(valst_eugspe_sm_res_ff , 'f5_fas_val_eugspe_sm_res' , 'f5_fas_val_steph_res' , jointplot_figs_dir , 'f5_fas_val_eugspe_sm_res')
plot_jointplot_res(valst_eugspe_sm_res_ff , 'f6_fas_val_eugspe_sm_res' , 'f6_fas_val_steph_res' , jointplot_figs_dir , 'f6_fas_val_eugspe_sm_res')
plot_jointplot_res(valst_eugspe_sm_res_ff , 'f7_fas_val_eugspe_sm_res' , 'f7_fas_val_steph_res' , jointplot_figs_dir , 'f7_fas_val_eugspe_sm_res')
plot_jointplot_res(valst_eugspe_sm_res_ff , 'pgv_eugspe_sm_res_val' , 'pgv_steph_res_val' , jointplot_figs_dir , 'pgv_eugspe_sm_res_val')

plot_jointplot_res(valst_srv_sm_res_ff , 'f2_fas_val_srv_sm_res' , 'f2_fas_val_steph_res' , jointplot_figs_dir , 'f2_fas_val_srv_sm_res')
plot_jointplot_res(valst_srv_sm_res_ff , 'f3_fas_val_srv_sm_res' , 'f3_fas_val_steph_res' , jointplot_figs_dir , 'f3_fas_val_srv_sm_res')
plot_jointplot_res(valst_srv_sm_res_ff , 'f4_fas_val_srv_sm_res' , 'f4_fas_val_steph_res' , jointplot_figs_dir , 'f4_fas_val_srv_sm_res')
plot_jointplot_res(valst_srv_sm_res_ff , 'f5_fas_val_srv_sm_res' , 'f5_fas_val_steph_res' , jointplot_figs_dir , 'f5_fas_val_srv_sm_res')
plot_jointplot_res(valst_srv_sm_res_ff , 'f6_fas_val_srv_sm_res' , 'f6_fas_val_steph_res' , jointplot_figs_dir , 'f6_fas_val_srv_sm_res')
plot_jointplot_res(valst_srv_sm_res_ff , 'f7_fas_val_srv_sm_res' , 'f7_fas_val_steph_res' , jointplot_figs_dir , 'f7_fas_val_srv_sm_res')
plot_jointplot_res(valst_srv_sm_res_ff , 'pgv_srv_sm_res_val' , 'pgv_steph_res_val' , jointplot_figs_dir , 'pgv_srv_sm_res_val')

plot_jointplot_res(valst_srv_1d_res_ff , 'f2_fas_val_srv_1d_res' , 'f2_fas_val_steph_res' , jointplot_figs_dir , 'f2_fas_val_srv_1d_res')
plot_jointplot_res(valst_srv_1d_res_ff , 'f3_fas_val_srv_1d_res' , 'f3_fas_val_steph_res' , jointplot_figs_dir , 'f3_fas_val_srv_1d_res')
plot_jointplot_res(valst_srv_1d_res_ff , 'f4_fas_val_srv_1d_res' , 'f4_fas_val_steph_res' , jointplot_figs_dir , 'f4_fas_val_srv_1d_res')
plot_jointplot_res(valst_srv_1d_res_ff , 'f5_fas_val_srv_1d_res' , 'f5_fas_val_steph_res' , jointplot_figs_dir , 'f5_fas_val_srv_1d_res')
plot_jointplot_res(valst_srv_1d_res_ff , 'f6_fas_val_srv_1d_res' , 'f6_fas_val_steph_res' , jointplot_figs_dir , 'f6_fas_val_srv_1d_res')
plot_jointplot_res(valst_srv_1d_res_ff , 'f7_fas_val_srv_1d_res' , 'f7_fas_val_steph_res' , jointplot_figs_dir , 'f7_fas_val_srv_1d_res')
plot_jointplot_res(valst_srv_1d_res_ff , 'pgv_srv_1d_res_val' , 'pgv_steph_res_val' , jointplot_figs_dir , 'pgv_srv_1d_res_val')


#%%

def make_ims_hist(data_ff,im_str,figs_dir,fig_name):
    plt.figure()
    sns.set(font_scale=1)
    sns.histplot(data=data_ff, x=im_str, bins='auto',binrange=(-3,3))
    plt.savefig(figs_dir + fig_name + '.png',dpi=150)

hist_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/residuals_histograms/springfield/'
try:
    os.mkdir(hist_figs_dir)
except FileExistsError:
    print('Directory already exists')


make_ims_hist(eugspe_1d_res_ff , 'f2_fas_eugspe_1d_res' , hist_figs_dir , 'f2_fas_eugspe_1d_res_hist')
make_ims_hist(eugspe_1d_res_ff , 'f3_fas_eugspe_1d_res' , hist_figs_dir , 'f3_fas_eugspe_1d_res_hist')
make_ims_hist(eugspe_1d_res_ff , 'f4_fas_eugspe_1d_res' , hist_figs_dir , 'f4_fas_eugspe_1d_res_hist')
make_ims_hist(eugspe_1d_res_ff , 'f5_fas_eugspe_1d_res' , hist_figs_dir , 'f5_fas_eugspe_1d_res_hist')
make_ims_hist(eugspe_1d_res_ff , 'f6_fas_eugspe_1d_res' , hist_figs_dir , 'f6_fas_eugspe_1d_res_hist')
make_ims_hist(eugspe_1d_res_ff , 'f7_fas_eugspe_1d_res' , hist_figs_dir , 'f7_fas_eugspe_1d_res_hist')
make_ims_hist(eugspe_1d_res_ff , 'pgv_eugspe_1d_res' , hist_figs_dir , 'pgv_eugspe_1d_res_hist')

make_ims_hist(valst_eugspe_1d_res_ff , 'f5_fas_val_eugspe_1d_res' , hist_figs_dir , 'f5_fas_val_eugspe_1d_res_hist')
make_ims_hist(valst_eugspe_1d_res_ff , 'f5_fas_val_eugspe_1d_res' , hist_figs_dir , 'f5_fas_val_eugspe_1d_res_hist')
make_ims_hist(valst_eugspe_1d_res_ff , 'f5_fas_val_eugspe_1d_res' , hist_figs_dir , 'f5_fas_val_eugspe_1d_res_hist')
make_ims_hist(valst_eugspe_1d_res_ff , 'f5_fas_val_eugspe_1d_res' , hist_figs_dir , 'f5_fas_val_eugspe_1d_res_hist')
make_ims_hist(valst_eugspe_1d_res_ff , 'f5_fas_val_eugspe_1d_res' , hist_figs_dir , 'f5_fas_val_eugspe_1d_res_hist')
make_ims_hist(valst_eugspe_1d_res_ff , 'f5_fas_val_eugspe_1d_res' , hist_figs_dir , 'f5_fas_val_eugspe_1d_res_hist')
make_ims_hist(valst_eugspe_1d_res_ff , 'pgv_eugspe_1d_res_val' , hist_figs_dir , 'pgv_eugspe_1d_res_val_hist')

###

make_ims_hist(eugspe_sm_res_ff , 'f2_fas_eugspe_sm_res' , hist_figs_dir , 'f2_fas_eugspe_sm_res_hist')
make_ims_hist(eugspe_sm_res_ff , 'f3_fas_eugspe_sm_res' , hist_figs_dir , 'f3_fas_eugspe_sm_res_hist')
make_ims_hist(eugspe_sm_res_ff , 'f4_fas_eugspe_sm_res' , hist_figs_dir , 'f4_fas_eugspe_sm_res_hist')
make_ims_hist(eugspe_sm_res_ff , 'f5_fas_eugspe_sm_res' , hist_figs_dir , 'f5_fas_eugspe_sm_res_hist')
make_ims_hist(eugspe_sm_res_ff , 'f6_fas_eugspe_sm_res' , hist_figs_dir , 'f6_fas_eugspe_sm_res_hist')
make_ims_hist(eugspe_sm_res_ff , 'f7_fas_eugspe_sm_res' , hist_figs_dir , 'f7_fas_eugspe_sm_res_hist')
make_ims_hist(eugspe_sm_res_ff , 'pgv_eugspe_sm_res' , hist_figs_dir , 'pgv_eugspe_sm_res_hist')

make_ims_hist(valst_eugspe_sm_res_ff , 'f5_fas_val_eugspe_sm_res' , hist_figs_dir , 'f5_fas_val_eugspe_sm_res_hist')
make_ims_hist(valst_eugspe_sm_res_ff , 'f5_fas_val_eugspe_sm_res' , hist_figs_dir , 'f5_fas_val_eugspe_sm_res_hist')
make_ims_hist(valst_eugspe_sm_res_ff , 'f5_fas_val_eugspe_sm_res' , hist_figs_dir , 'f5_fas_val_eugspe_sm_res_hist')
make_ims_hist(valst_eugspe_sm_res_ff , 'f5_fas_val_eugspe_sm_res' , hist_figs_dir , 'f5_fas_val_eugspe_sm_res_hist')
make_ims_hist(valst_eugspe_sm_res_ff , 'f5_fas_val_eugspe_sm_res' , hist_figs_dir , 'f5_fas_val_eugspe_sm_res_hist')
make_ims_hist(valst_eugspe_sm_res_ff , 'f5_fas_val_eugspe_sm_res' , hist_figs_dir , 'f5_fas_val_eugspe_sm_res_hist')
make_ims_hist(valst_eugspe_sm_res_ff , 'pgv_eugspe_sm_res_val' , hist_figs_dir , 'pgv_eugspe_sm_res_val_hist')

###

make_ims_hist(srv_1d_res_ff , 'f2_fas_srv_1d_res' , hist_figs_dir , 'f2_fas_srv_1d_res_hist')
make_ims_hist(srv_1d_res_ff , 'f3_fas_srv_1d_res' , hist_figs_dir , 'f3_fas_srv_1d_res_hist')
make_ims_hist(srv_1d_res_ff , 'f4_fas_srv_1d_res' , hist_figs_dir , 'f4_fas_srv_1d_res_hist')
make_ims_hist(srv_1d_res_ff , 'f5_fas_srv_1d_res' , hist_figs_dir , 'f5_fas_srv_1d_res_hist')
make_ims_hist(srv_1d_res_ff , 'f6_fas_srv_1d_res' , hist_figs_dir , 'f6_fas_srv_1d_res_hist')
make_ims_hist(srv_1d_res_ff , 'f7_fas_srv_1d_res' , hist_figs_dir , 'f7_fas_srv_1d_res_hist')
make_ims_hist(srv_1d_res_ff , 'pgv_srv_1d_res' , hist_figs_dir , 'pgv_srv_1d_res_hist')

make_ims_hist(valst_srv_1d_res_ff , 'f5_fas_val_srv_1d_res' , hist_figs_dir , 'f5_fas_val_srv_1d_res_hist')
make_ims_hist(valst_srv_1d_res_ff , 'f5_fas_val_srv_1d_res' , hist_figs_dir , 'f5_fas_val_srv_1d_res_hist')
make_ims_hist(valst_srv_1d_res_ff , 'f5_fas_val_srv_1d_res' , hist_figs_dir , 'f5_fas_val_srv_1d_res_hist')
make_ims_hist(valst_srv_1d_res_ff , 'f5_fas_val_srv_1d_res' , hist_figs_dir , 'f5_fas_val_srv_1d_res_hist')
make_ims_hist(valst_srv_1d_res_ff , 'f5_fas_val_srv_1d_res' , hist_figs_dir , 'f5_fas_val_srv_1d_res_hist')
make_ims_hist(valst_srv_1d_res_ff , 'f5_fas_val_srv_1d_res' , hist_figs_dir , 'f5_fas_val_srv_1d_res_hist')
make_ims_hist(valst_srv_1d_res_ff , 'pgv_srv_1d_res_val' , hist_figs_dir , 'pgv_srv_1d_res_val_hist')

###

make_ims_hist(srv_sm_res_ff , 'f2_fas_srv_sm_res' , hist_figs_dir , 'f2_fas_srv_1d_res_hist')
make_ims_hist(srv_sm_res_ff , 'f3_fas_srv_sm_res' , hist_figs_dir , 'f3_fas_srv_1d_res_hist')
make_ims_hist(srv_sm_res_ff , 'f4_fas_srv_sm_res' , hist_figs_dir , 'f4_fas_srv_1d_res_hist')
make_ims_hist(srv_sm_res_ff , 'f5_fas_srv_sm_res' , hist_figs_dir , 'f5_fas_srv_1d_res_hist')
make_ims_hist(srv_sm_res_ff , 'f6_fas_srv_sm_res' , hist_figs_dir , 'f6_fas_srv_1d_res_hist')
make_ims_hist(srv_sm_res_ff , 'f7_fas_srv_sm_res' , hist_figs_dir , 'f7_fas_srv_1d_res_hist')
make_ims_hist(srv_sm_res_ff , 'pgv_srv_sm_res' , hist_figs_dir , 'pgv_srv_1d_res_hist')

make_ims_hist(valst_srv_sm_res_ff , 'f5_fas_val_srv_sm_res' , hist_figs_dir , 'f5_fas_val_srv_sm_res_hist')
make_ims_hist(valst_srv_sm_res_ff , 'f5_fas_val_srv_sm_res' , hist_figs_dir , 'f5_fas_val_srv_sm_res_hist')
make_ims_hist(valst_srv_sm_res_ff , 'f5_fas_val_srv_sm_res' , hist_figs_dir , 'f5_fas_val_srv_sm_res_hist')
make_ims_hist(valst_srv_sm_res_ff , 'f5_fas_val_srv_sm_res' , hist_figs_dir , 'f5_fas_val_srv_sm_res_hist')
make_ims_hist(valst_srv_sm_res_ff , 'f5_fas_val_srv_sm_res' , hist_figs_dir , 'f5_fas_val_srv_sm_res_hist')
make_ims_hist(valst_srv_sm_res_ff , 'f5_fas_val_srv_sm_res' , hist_figs_dir , 'f5_fas_val_srv_sm_res_hist')
make_ims_hist(valst_srv_sm_res_ff , 'pgv_srv_sm_res_val' , hist_figs_dir , 'pgv_srv_sm_res_val_hist')

#%%

def plot_mean_med_scat(data_ff , models_col , value_col , figs_dir , fig_name):
    plt.figure()
    sns.set(font_scale=1)
    vals = data_ff[value_col]
    max_val = max(vals)
    min_val = min(vals)
    
    scat = sns.scatterplot(data=data_ff , x=models_col , y=value_col)
    scat.set_ylim(min_val-0.02, max_val+0.02)
    plt.savefig(figs_dir + fig_name + '.png',dpi=150)

mean_med_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/models_mean_med_scatter/springfield/'
try:
    os.mkdir(mean_med_figs_dir)
except FileExistsError:
    print('Directory already exists')

fas_2_mean_med_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/ims_mean_med_fas2.csv')
fas_3_mean_med_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/ims_mean_med_fas3.csv')
fas_4_mean_med_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/ims_mean_med_fas4.csv')
fas_5_mean_med_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/ims_mean_med_fas5.csv')
fas_6_mean_med_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/ims_mean_med_fas6.csv')
fas_7_mean_med_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/ims_mean_med_fas7.csv')
pgv_mean_med_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/ims_mean_med_pgv.csv')

plot_mean_med_scat(fas_2_mean_med_ff , 'model_names' , 'means' , mean_med_figs_dir , 'fas_2_mean')
plot_mean_med_scat(fas_3_mean_med_ff , 'model_names' , 'means' , mean_med_figs_dir , 'fas_3_mean')
plot_mean_med_scat(fas_4_mean_med_ff , 'model_names' , 'means' , mean_med_figs_dir , 'fas_4_mean')
plot_mean_med_scat(fas_5_mean_med_ff , 'model_names' , 'means' , mean_med_figs_dir , 'fas_5_mean')
plot_mean_med_scat(fas_6_mean_med_ff , 'model_names' , 'means' , mean_med_figs_dir , 'fas_6_mean')
plot_mean_med_scat(fas_7_mean_med_ff , 'model_names' , 'means' , mean_med_figs_dir , 'fas_7_mean')
plot_mean_med_scat(pgv_mean_med_ff , 'model_names' , 'means' , mean_med_figs_dir , 'pgv_mean')

plot_mean_med_scat(fas_2_mean_med_ff , 'model_names' , 'medians' , mean_med_figs_dir , 'fas_2_med')
plot_mean_med_scat(fas_3_mean_med_ff , 'model_names' , 'medians' , mean_med_figs_dir , 'fas_3_med')
plot_mean_med_scat(fas_4_mean_med_ff , 'model_names' , 'medians' , mean_med_figs_dir , 'fas_4_med')
plot_mean_med_scat(fas_5_mean_med_ff , 'model_names' , 'medians' , mean_med_figs_dir , 'fas_5_med')
plot_mean_med_scat(fas_6_mean_med_ff , 'model_names' , 'medians' , mean_med_figs_dir , 'fas_6_med')
plot_mean_med_scat(fas_7_mean_med_ff , 'model_names' , 'medians' , mean_med_figs_dir , 'fas_7_med')
plot_mean_med_scat(pgv_mean_med_ff , 'model_names' , 'medians' , mean_med_figs_dir , 'pgv_med')


#%%

def make_multi_hist(data_ff , im , hist_type , figs_dir , fig_name , hue):
    if hist_type == 'reg':
        plt.figure()
        sns.set(font_scale=1)
        hist = sns.histplot(data=data_ff, x=im, bins='auto',binrange=(-3,3),hue=hue)
        sns.move_legend(hist, 'upper left')
        hist.set_title(fig_name,size=15)
        plt.savefig(figs_dir + fig_name +'_'+ hist_type + '.png',dpi=150)
        
    if hist_type == 'stack':
        plt.figure()
        sns.set(font_scale=1)
        hist = sns.histplot(data=data_ff, x=im, bins='auto',binrange=(-3,3),hue=hue,multiple="stack")
        sns.move_legend(hist, 'upper left')
        hist.set_title(fig_name,size=15)
        plt.savefig(figs_dir + fig_name +'_'+ hist_type + '.png',dpi=150)
        
    if hist_type == 'step':
        plt.figure()
        sns.set(font_scale=1)
        hist = sns.histplot(data=data_ff, x=im, bins='auto',binrange=(-3,3),hue=hue,element='step')
        sns.move_legend(hist, 'upper left')
        hist.set_title(fig_name,size=15)
        plt.savefig(figs_dir + fig_name +'_'+ hist_type + '.png',dpi=150)
        

all_model_val_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/all_model_val_ff.csv')
all_model_nosteph_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/all_model_nosteph_ff.csv')
all_model_nosteph_val_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/all_model_nosteph_val_ff.csv')
all_model_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/all_model_ff.csv')

multi_hist_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/multi_model_hist/springfield/'
try:
    os.mkdir(multi_hist_figs_dir)
except FileExistsError:
    print('Directory already exists')

make_multi_hist(all_model_val_ff , 'f2_fas_res' , 'reg' , multi_hist_figs_dir , 'f2_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f3_fas_res' , 'reg' , multi_hist_figs_dir , 'f3_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f4_fas_res' , 'reg' , multi_hist_figs_dir , 'f4_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f5_fas_res' , 'reg' , multi_hist_figs_dir , 'f5_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f6_fas_res' , 'reg' , multi_hist_figs_dir , 'f6_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f7_fas_res' , 'reg' , multi_hist_figs_dir , 'f7_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'pgv_res' , 'reg' , multi_hist_figs_dir , 'pgv_res_multihist_val' , 'model')

make_multi_hist(all_model_nosteph_ff , 'f2_fas_res' , 'reg' , multi_hist_figs_dir , 'f2_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f3_fas_res' , 'reg' , multi_hist_figs_dir , 'f3_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f4_fas_res' , 'reg' , multi_hist_figs_dir , 'f4_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f5_fas_res' , 'reg' , multi_hist_figs_dir , 'f5_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f6_fas_res' , 'reg' , multi_hist_figs_dir , 'f6_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f7_fas_res' , 'reg' , multi_hist_figs_dir , 'f7_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'pgv_res' , 'reg' , multi_hist_figs_dir , 'pgv_res_multihist_nosteph' , 'model')

make_multi_hist(all_model_nosteph_val_ff , 'f2_fas_res' , 'reg' , multi_hist_figs_dir , 'f2_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f3_fas_res' , 'reg' , multi_hist_figs_dir , 'f3_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f4_fas_res' , 'reg' , multi_hist_figs_dir , 'f4_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f5_fas_res' , 'reg' , multi_hist_figs_dir , 'f5_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f6_fas_res' , 'reg' , multi_hist_figs_dir , 'f6_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f7_fas_res' , 'reg' , multi_hist_figs_dir , 'f7_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'pgv_res' , 'reg' , multi_hist_figs_dir , 'pgv_res_multihist_nosteph_val' , 'model')

make_multi_hist(all_model_ff , 'f2_fas_res' , 'reg' , multi_hist_figs_dir , 'f2_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f3_fas_res' , 'reg' , multi_hist_figs_dir , 'f3_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f4_fas_res' , 'reg' , multi_hist_figs_dir , 'f4_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f5_fas_res' , 'reg' , multi_hist_figs_dir , 'f5_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f6_fas_res' , 'reg' , multi_hist_figs_dir , 'f6_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f7_fas_res' , 'reg' , multi_hist_figs_dir , 'f7_fas_res_multihistl' , 'model')
make_multi_hist(all_model_ff , 'pgv_res' , 'reg' , multi_hist_figs_dir , 'pgv_res_multihist' , 'model')

####

make_multi_hist(all_model_val_ff , 'f2_fas_res' , 'stack' , multi_hist_figs_dir , 'f2_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f3_fas_res' , 'stack' , multi_hist_figs_dir , 'f3_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f4_fas_res' , 'stack' , multi_hist_figs_dir , 'f4_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f5_fas_res' , 'stack' , multi_hist_figs_dir , 'f5_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f6_fas_res' , 'stack' , multi_hist_figs_dir , 'f6_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f7_fas_res' , 'stack' , multi_hist_figs_dir , 'f7_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'pgv_res' , 'stack' , multi_hist_figs_dir , 'pgv_res_multihist_val' , 'model')

make_multi_hist(all_model_nosteph_ff , 'f2_fas_res' , 'stack' , multi_hist_figs_dir , 'f2_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f3_fas_res' , 'stack' , multi_hist_figs_dir , 'f3_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f4_fas_res' , 'stack' , multi_hist_figs_dir , 'f4_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f5_fas_res' , 'stack' , multi_hist_figs_dir , 'f5_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f6_fas_res' , 'stack' , multi_hist_figs_dir , 'f6_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f7_fas_res' , 'stack' , multi_hist_figs_dir , 'f7_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'pgv_res' , 'stack' , multi_hist_figs_dir , 'pgv_res_multihist_nosteph' , 'model')

make_multi_hist(all_model_nosteph_val_ff , 'f2_fas_res' , 'stack' , multi_hist_figs_dir , 'f2_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f3_fas_res' , 'stack' , multi_hist_figs_dir , 'f3_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f4_fas_res' , 'stack' , multi_hist_figs_dir , 'f4_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f5_fas_res' , 'stack' , multi_hist_figs_dir , 'f5_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f6_fas_res' , 'stack' , multi_hist_figs_dir , 'f6_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f7_fas_res' , 'stack' , multi_hist_figs_dir , 'f7_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'pgv_res' , 'stack' , multi_hist_figs_dir , 'pgv_res_multihist_nosteph_val' , 'model')

make_multi_hist(all_model_ff , 'f2_fas_res' , 'stack' , multi_hist_figs_dir , 'f2_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f3_fas_res' , 'stack' , multi_hist_figs_dir , 'f3_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f4_fas_res' , 'stack' , multi_hist_figs_dir , 'f4_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f5_fas_res' , 'stack' , multi_hist_figs_dir , 'f5_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f6_fas_res' , 'stack' , multi_hist_figs_dir , 'f6_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f7_fas_res' , 'stack' , multi_hist_figs_dir , 'f7_fas_res_multihistl' , 'model')
make_multi_hist(all_model_ff , 'pgv_res' , 'stack' , multi_hist_figs_dir , 'pgv_res_multihist' , 'model')

####

make_multi_hist(all_model_val_ff , 'f2_fas_res' , 'step' , multi_hist_figs_dir , 'f2_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f3_fas_res' , 'step' , multi_hist_figs_dir , 'f3_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f4_fas_res' , 'step' , multi_hist_figs_dir , 'f4_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f5_fas_res' , 'step' , multi_hist_figs_dir , 'f5_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f6_fas_res' , 'step' , multi_hist_figs_dir , 'f6_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'f7_fas_res' , 'step' , multi_hist_figs_dir , 'f7_fas_res_multihist_val' , 'model')
make_multi_hist(all_model_val_ff , 'pgv_res' , 'step' , multi_hist_figs_dir , 'pgv_res_multihist_val' , 'model')

make_multi_hist(all_model_nosteph_ff , 'f2_fas_res' , 'step' , multi_hist_figs_dir , 'f2_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f3_fas_res' , 'step' , multi_hist_figs_dir , 'f3_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f4_fas_res' , 'step' , multi_hist_figs_dir , 'f4_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f5_fas_res' , 'step' , multi_hist_figs_dir , 'f5_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f6_fas_res' , 'step' , multi_hist_figs_dir , 'f6_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'f7_fas_res' , 'step' , multi_hist_figs_dir , 'f7_fas_res_multihist_nosteph' , 'model')
make_multi_hist(all_model_nosteph_ff , 'pgv_res' , 'step' , multi_hist_figs_dir , 'pgv_res_multihist_nosteph' , 'model')

make_multi_hist(all_model_nosteph_val_ff , 'f2_fas_res' , 'step' , multi_hist_figs_dir , 'f2_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f3_fas_res' , 'step' , multi_hist_figs_dir , 'f3_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f4_fas_res' , 'step' , multi_hist_figs_dir , 'f4_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f5_fas_res' , 'step' , multi_hist_figs_dir , 'f5_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f6_fas_res' , 'step' , multi_hist_figs_dir , 'f6_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'f7_fas_res' , 'step' , multi_hist_figs_dir , 'f7_fas_res_multihist_nosteph_val' , 'model')
make_multi_hist(all_model_nosteph_val_ff , 'pgv_res' , 'step' , multi_hist_figs_dir , 'pgv_res_multihist_nosteph_val' , 'model')

make_multi_hist(all_model_ff , 'f2_fas_res' , 'step' , multi_hist_figs_dir , 'f2_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f3_fas_res' , 'step' , multi_hist_figs_dir , 'f3_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f4_fas_res' , 'step' , multi_hist_figs_dir , 'f4_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f5_fas_res' , 'step' , multi_hist_figs_dir , 'f5_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f6_fas_res' , 'step' , multi_hist_figs_dir , 'f6_fas_res_multihist' , 'model')
make_multi_hist(all_model_ff , 'f7_fas_res' , 'step' , multi_hist_figs_dir , 'f7_fas_res_multihistl' , 'model')
make_multi_hist(all_model_ff , 'pgv_res' , 'step' , multi_hist_figs_dir , 'pgv_res_multihist' , 'model')

#%%

all_model_val_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/all_model_val_ff.csv')
all_model_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/all_model_ff.csv')

boxplots_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/models_boxplots/springfield/'
try:
    os.mkdir(boxplots_figs_dir)
except FileExistsError:
    print('Directory already exists')

def make_boxplot_ims(data_ff , im , figs_dir , fig_name):
    plt.figure()
    sns.set(font_scale=1)
    boxplot = sns.boxplot(x='model', y=im, data=data_ff, showfliers=False)
    boxplot.set_title(fig_name,size=15)
    plt.savefig(figs_dir + fig_name + '.png',dpi=150)

make_boxplot_ims(all_model_val_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res_val')
make_boxplot_ims(all_model_val_ff , 'f3_fas_res' , boxplots_figs_dir , 'f3_fas_res_val')
make_boxplot_ims(all_model_val_ff , 'f4_fas_res' , boxplots_figs_dir , 'f4_fas_res_val')
make_boxplot_ims(all_model_val_ff , 'f5_fas_res' , boxplots_figs_dir , 'f5_fas_res_val')
make_boxplot_ims(all_model_val_ff , 'f6_fas_res' , boxplots_figs_dir , 'f6_fas_res_val')
make_boxplot_ims(all_model_val_ff , 'f7_fas_res' , boxplots_figs_dir , 'f7_fas_res_val')
make_boxplot_ims(all_model_val_ff , 'pgv_res' , boxplots_figs_dir , 'pgv_res_val')

make_boxplot_ims(all_model_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res')
make_boxplot_ims(all_model_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res')
make_boxplot_ims(all_model_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res')
make_boxplot_ims(all_model_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res')
make_boxplot_ims(all_model_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res')
make_boxplot_ims(all_model_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res')
make_boxplot_ims(all_model_ff , 'f2_fas_res' , boxplots_figs_dir , 'f2_fas_res')

#%%

all_model_val_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/all_model_val_ff.csv')

violin_figs_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/figs/violin_plots/springfield/'

def make_violin_ims(data_ff , im , figs_dir , fig_name , limit=False):
    plt.figure()
    sns.set(font_scale=1)
    violin = sns.violinplot(x='model', y=im, data=data_ff )
    if limit == True:
        violin.set_ylim([-3, 3])
    plt.savefig(figs_dir + fig_name + '.png',dpi=150)


make_violin_ims(all_model_val_ff , 'f2_fas_res' , violin_figs_dir, 'f2_fas_res')
make_violin_ims(all_model_val_ff , 'f3_fas_res' , violin_figs_dir, 'f3_fas_res')
make_violin_ims(all_model_val_ff , 'f4_fas_res' , violin_figs_dir, 'f4_fas_res')
make_violin_ims(all_model_val_ff , 'f5_fas_res' , violin_figs_dir, 'f5_fas_res')
make_violin_ims(all_model_val_ff , 'f6_fas_res' , violin_figs_dir, 'f6_fas_res')
make_violin_ims(all_model_val_ff , 'f7_fas_res' , violin_figs_dir, 'f7_fas_res')
make_violin_ims(all_model_val_ff , 'pgv_res' , violin_figs_dir, 'pgv_res')

make_violin_ims(all_model_val_ff , 'f2_fas_res' , violin_figs_dir, 'f2_fas_res_limit', limit=True)
make_violin_ims(all_model_val_ff , 'f3_fas_res' , violin_figs_dir, 'f3_fas_res_limit', limit=True)
make_violin_ims(all_model_val_ff , 'f4_fas_res' , violin_figs_dir, 'f4_fas_res_limit', limit=True)
make_violin_ims(all_model_val_ff , 'f5_fas_res' , violin_figs_dir, 'f5_fas_res_limit', limit=True)
make_violin_ims(all_model_val_ff , 'f6_fas_res' , violin_figs_dir, 'f6_fas_res_limit', limit=True)
make_violin_ims(all_model_val_ff , 'f7_fas_res' , violin_figs_dir, 'f7_fas_res_limit', limit=True)
make_violin_ims(all_model_val_ff , 'pgv_res' , violin_figs_dir, 'pgv_res_limit', limit=True)












































