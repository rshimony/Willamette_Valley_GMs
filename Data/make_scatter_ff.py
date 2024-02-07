#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 17:56:50 2024

@author: rshimony
"""

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

#%%

steph_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/steph_res_ff.csv')
eugspe_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/eugspe_1d_res_ff.csv')
eugspe_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/eugspe_sm_res_ff.csv')
srv_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/srv_1d_res_ff.csv')
srv_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/srv_sm_res_ff.csv')

valst_steph_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/steph_res_ff_val.csv')
valst_eugspe_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/eugspe_1d_res_ff_val.csv')
valst_eugspe_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/eugspe_sm_res_ff_val.csv')
valst_srv_1d_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/srv_1d_res_ff_val.csv')
valst_srv_sm_res_ff = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/srv_sm_res_ff_val.csv')

#%% make arrays of mean and median residuals for each model (valley stations only)

def mean_median_arrays(flatfile):
    im_means = []
    im_medians = []
    model_name = []
    for i in range(6,13):
        im_col = flatfile.iloc[:,i]
        im_mean = np.mean(im_col)
        im_median = np.median(im_col)
        im_means.append(im_mean)
        im_medians.append(im_median)
        model_name.append(flatfile.iloc[0,5])
    im_means_arr = np.array(im_means)
    im_medians_arr = np.array(im_medians)
    model_name_arr = np.array(model_name)
    
    return(im_means_arr,im_medians_arr,model_name_arr)

steph_mean_med = mean_median_arrays(valst_steph_res_ff)
eugspe_1d_mean_med = mean_median_arrays(valst_eugspe_1d_res_ff)
eugspe_sm_mean_med = mean_median_arrays(valst_eugspe_sm_res_ff)
srv_1d_mean_med = mean_median_arrays(valst_srv_1d_res_ff)
srv_sm_mean_med = mean_median_arrays(valst_srv_sm_res_ff)


fas_2_mean = np.array([eugspe_1d_mean_med[0][0] , eugspe_sm_mean_med[0][0] , srv_1d_mean_med[0][0] , srv_sm_mean_med[0][0] , steph_mean_med[0][0]])
fas_3_mean = np.array([eugspe_1d_mean_med[0][1] , eugspe_sm_mean_med[0][1] , srv_1d_mean_med[0][1] , srv_sm_mean_med[0][1] , steph_mean_med[0][1]])
fas_4_mean = np.array([eugspe_1d_mean_med[0][2] , eugspe_sm_mean_med[0][2] , srv_1d_mean_med[0][2] , srv_sm_mean_med[0][2] , steph_mean_med[0][2]])
fas_5_mean = np.array([eugspe_1d_mean_med[0][3] , eugspe_sm_mean_med[0][3] , srv_1d_mean_med[0][3] , srv_sm_mean_med[0][3] , steph_mean_med[0][3]])
fas_6_mean = np.array([eugspe_1d_mean_med[0][4] , eugspe_sm_mean_med[0][4] , srv_1d_mean_med[0][4] , srv_sm_mean_med[0][4] , steph_mean_med[0][4]])
fas_7_mean = np.array([eugspe_1d_mean_med[0][5] , eugspe_sm_mean_med[0][5] , srv_1d_mean_med[0][5] , srv_sm_mean_med[0][5] , steph_mean_med[0][5]])
pgv_mean = np.array([eugspe_1d_mean_med[0][6] , eugspe_sm_mean_med[0][6] , srv_1d_mean_med[0][6] , srv_sm_mean_med[0][6] , steph_mean_med[0][6]])

fas_2_med = np.array([eugspe_1d_mean_med[1][0] , eugspe_sm_mean_med[1][0] , srv_1d_mean_med[1][0] , srv_sm_mean_med[1][0] , steph_mean_med[1][0]])
fas_3_med = np.array([eugspe_1d_mean_med[1][1] , eugspe_sm_mean_med[1][1] , srv_1d_mean_med[1][1] , srv_sm_mean_med[1][1] , steph_mean_med[1][1]])
fas_4_med = np.array([eugspe_1d_mean_med[1][2] , eugspe_sm_mean_med[1][2] , srv_1d_mean_med[1][2] , srv_sm_mean_med[1][2] , steph_mean_med[1][2]])
fas_5_med = np.array([eugspe_1d_mean_med[1][3] , eugspe_sm_mean_med[1][3] , srv_1d_mean_med[1][3] , srv_sm_mean_med[1][3] , steph_mean_med[1][3]])
fas_6_med = np.array([eugspe_1d_mean_med[1][4] , eugspe_sm_mean_med[1][4] , srv_1d_mean_med[1][4] , srv_sm_mean_med[1][4] , steph_mean_med[1][4]])
fas_7_med = np.array([eugspe_1d_mean_med[1][5] , eugspe_sm_mean_med[1][5] , srv_1d_mean_med[1][5] , srv_sm_mean_med[1][5] , steph_mean_med[1][5]])
pgv_med = np.array([eugspe_1d_mean_med[1][6] , eugspe_sm_mean_med[1][6] , srv_1d_mean_med[1][6] , srv_sm_mean_med[1][6] , steph_mean_med[1][6]])

model_names = np.array([eugspe_1d_mean_med[2][0] , eugspe_sm_mean_med[2][0] , srv_1d_mean_med[2][0] , srv_sm_mean_med[2][0] , steph_mean_med[2][0]])

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/salem/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')
    

data_mean_med_fas2 = {'model_names':model_names , 'means':fas_2_mean , 'medians':fas_2_med}
df_mean_med_fas2 = pd.DataFrame(data=data_mean_med_fas2)
df_mean_med_fas2.to_csv(resrat_ff_dir+'ims_mean_med_fas2.csv')

data_mean_med_fas3 = {'model_names':model_names , 'means':fas_3_mean , 'medians':fas_3_med}
df_mean_med_fas3 = pd.DataFrame(data=data_mean_med_fas3)
df_mean_med_fas3.to_csv(resrat_ff_dir+'ims_mean_med_fas3.csv')

data_mean_med_fas4 = {'model_names':model_names , 'means':fas_4_mean , 'medians':fas_4_med}
df_mean_med_fas4 = pd.DataFrame(data=data_mean_med_fas4)
df_mean_med_fas4.to_csv(resrat_ff_dir+'ims_mean_med_fas4.csv')

data_mean_med_fas5 = {'model_names':model_names , 'means':fas_5_mean , 'medians':fas_5_med}
df_mean_med_fas5 = pd.DataFrame(data=data_mean_med_fas5)
df_mean_med_fas5.to_csv(resrat_ff_dir+'ims_mean_med_fas5.csv')

data_mean_med_fas6 = {'model_names':model_names , 'means':fas_6_mean , 'medians':fas_6_med}
df_mean_med_fas6 = pd.DataFrame(data=data_mean_med_fas6)
df_mean_med_fas6.to_csv(resrat_ff_dir+'ims_mean_med_fas6.csv')

data_mean_med_fas7 = {'model_names':model_names , 'means':fas_7_mean , 'medians':fas_7_med}
df_mean_med_fas7 = pd.DataFrame(data=data_mean_med_fas7)
df_mean_med_fas7.to_csv(resrat_ff_dir+'ims_mean_med_fas7.csv')

data_mean_med_pgv = {'model_names':model_names , 'means':pgv_mean , 'medians':pgv_med}
df_mean_med_pgv = pd.DataFrame(data=data_mean_med_pgv)
df_mean_med_pgv.to_csv(resrat_ff_dir+'ims_mean_med_pgv.csv')

#%%

resrat_ff_dir = '/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/'
try:
    os.mkdir(resrat_ff_dir)
except FileExistsError:
    print('Directory already exists')

steph_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/steph_res_ff.csv', header=None)
eugspe_1d_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_1d_res_ff.csv', header=None)
eugspe_sm_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_sm_res_ff.csv', header=None)
srv_1d_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_1d_res_ff.csv', header=None)
srv_sm_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_sm_res_ff.csv', header=None)

valst_steph_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/steph_res_ff_val.csv', header=None)
valst_eugspe_1d_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_1d_res_ff_val.csv', header=None)
valst_eugspe_sm_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/eugspe_sm_res_ff_val.csv', header=None)
valst_srv_1d_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_1d_res_ff_val.csv', header=None)
valst_srv_sm_res_ff_nohead = pd.read_csv('/Users/rshimony/Desktop/WillametteValley/final_outputs/flatfiles/res_rat/springfield/srv_sm_res_ff_val.csv', header=None)


all_model_val_ff = pd.concat([valst_eugspe_1d_res_ff_nohead.iloc[1:,1:20],valst_eugspe_sm_res_ff_nohead.iloc[1:,1:20],valst_srv_1d_res_ff_nohead.iloc[1:,1:20],
                              valst_srv_sm_res_ff_nohead.iloc[1:,1:20],valst_steph_res_ff_nohead.iloc[1:,1:20]],ignore_index=True)
all_model_val_ff.columns=['st_name', 'st_lon', 'st_lat', 'lat_region', 'model', 
                          'f2_fas_res', 'f3_fas_res', 'f4_fas_res', 'f5_fas_res', 'f6_fas_res', 'f7_fas_res', 'pgv_res',
                          'f2_fas_ratio', 'f3_fas_ratio', 'f4_fas_ratio', 'f5_fas_ratio', 'f6_fas_ratio', 'f7_fas_ratio', 'pgv_ratio']
all_model_val_ff.to_csv(resrat_ff_dir+'all_model_val_ff.csv')


all_model_ff = pd.concat([eugspe_1d_res_ff_nohead.iloc[1:,1:20],eugspe_sm_res_ff_nohead.iloc[1:,1:20],
                          srv_1d_res_ff_nohead.iloc[1:,1:20],srv_sm_res_ff_nohead.iloc[1:,1:20],steph_res_ff_nohead.iloc[1:,1:20]], ignore_index=True)
all_model_ff.columns=['st_name', 'st_lon', 'st_lat', 'lat_region', 'model', 
                          'f2_fas_res', 'f3_fas_res', 'f4_fas_res', 'f5_fas_res', 'f6_fas_res', 'f7_fas_res', 'pgv_res',
                          'f2_fas_ratio', 'f3_fas_ratio', 'f4_fas_ratio', 'f5_fas_ratio', 'f6_fas_ratio', 'f7_fas_ratio', 'pgv_ratio']
all_model_ff.to_csv(resrat_ff_dir+'all_model_ff.csv')


all_model_nosteph_val_ff = pd.concat([valst_eugspe_1d_res_ff_nohead.iloc[1:,1:20],valst_eugspe_sm_res_ff_nohead.iloc[1:,1:20],valst_srv_1d_res_ff_nohead.iloc[1:,1:20]
                                      ,valst_srv_sm_res_ff_nohead.iloc[1:,1:20]], ignore_index=True)
all_model_nosteph_val_ff.columns=['st_name', 'st_lon', 'st_lat', 'lat_region', 'model', 
                          'f2_fas_res', 'f3_fas_res', 'f4_fas_res', 'f5_fas_res', 'f6_fas_res', 'f7_fas_res', 'pgv_res',
                          'f2_fas_ratio', 'f3_fas_ratio', 'f4_fas_ratio', 'f5_fas_ratio', 'f6_fas_ratio', 'f7_fas_ratio', 'pgv_ratio']
all_model_nosteph_val_ff.to_csv(resrat_ff_dir+'all_model_nosteph_val_ff.csv')


all_model_nosteph_ff = pd.concat([eugspe_1d_res_ff_nohead.iloc[1:,1:20],eugspe_sm_res_ff_nohead.iloc[1:,1:20],srv_1d_res_ff_nohead.iloc[1:,1:20]
                                  ,srv_sm_res_ff_nohead.iloc[1:,1:20]], ignore_index=True)
all_model_nosteph_ff.columns=['st_name', 'st_lon', 'st_lat', 'lat_region', 'model', 
                          'f2_fas_res', 'f3_fas_res', 'f4_fas_res', 'f5_fas_res', 'f6_fas_res', 'f7_fas_res', 'pgv_res',
                          'f2_fas_ratio', 'f3_fas_ratio', 'f4_fas_ratio', 'f5_fas_ratio', 'f6_fas_ratio', 'f7_fas_ratio', 'pgv_ratio']
all_model_nosteph_ff.to_csv(resrat_ff_dir+'all_model_nosteph_ff.csv')




























