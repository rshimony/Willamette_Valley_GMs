#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 13:37:26 2023

@author: oluwaseunfadugba
"""

import os
import numpy as np
import glob 
import time       
import pandas as pd

start = time.time()
current_dir = os.getcwd()
working_dir = '/home/ofadugba/waves'
#working_dir = '/Users/oluwaseunfadugba/Documents/Projects/TsE_ValerieDiego/TsE_1D_vs_3D/3D_Modeling_using_SW4/Running_Simulations/4_Model_results/'
username = 'ofadugba'
starttime = '2020-01-01'

#%%
os.chdir(working_dir)

if os.path.exists('jobs.txt')==False:
    os.system('sacct -S '+starttime+' -u '+username+' --format=JobID,JobName,ReqMem,MaxRSS,Elapsed,Partition,Account,AllocCPUS,State,ExitCode > jobs.txt')

job_id_list = np.array([])
n_grdpts_list = np.array([])
total_seismic_moment_list = np.array([])
mag_list = np.array([])
n_sources_list = np.array([])
hmin_list = np.array([])
PPW_list = np.array([])
goal_time_list = np.array([])
n_timesteps_list = np.array([])
fmax_list = np.array([])
n_cores_list = np.array([])
max_mem_req_list = np.array([])
max_mem_used_GB_list = np.array([])
run_time_list = np.array([])
description_list = np.array([])
username_list = np.array([])

# loop over slurm outputs
for filename in glob.iglob(working_dir+'/**/slurm*',recursive = True):
    
    # find only simulations that finished
    if os.system('grep -q "program sw4 finished" '+ filename) == 0:
        
        # extracting needed parameters
        with open(filename, 'r') as fp:
            # read all lines using readline()
            lines = fp.readlines()
            
            job_id = filename.split('/')[-1].split('.')[0].split('-')[-1]
            username = filename.split('/')[2]
            
            for row in lines:
                # check if string present on a current line
                if row[0:25] ==  ' Making Output Directory:':
                    description = str(row.split(' ')[-1])
                elif row[0:20] == 'Total number of grid':
                    n_grdpts = row.split(' ')[-1].split('\n')[0]
                elif row[29:35] == '<=  Vs':
                    minvs = row.split(' ')[7]
                elif row[21:29] == 'minVs/h=':
                    hmin = row[7:19]
                    minvs_h = float(row[29:].split(' ')[0])
                elif row[6:15] == 'max freq=':
                    fmax = float(row[15:27])
                elif row[0:14] == 'Running sw4 on':
                    n_cores = row.split(' ')[-2]
                elif row[1:13] == 'Start Time =':
                    goal_time = float(row.split(' ')[-1])
                elif row[1:23] == 'Number of time steps =':
                    n_timesteps = row.split(' ')[-3]
                elif row[2:22] ==  'Total seismic moment':
                    total_seismic_moment = row.split(' ')[-3]
                elif row[2:18] == 'Moment magnitude':
                    mag = float(row.split(' ')[-1])
                elif row[2:26] == 'Number of moment sources':
                    n_sources = float(row.split(' ')[-1])
                elif row[3:31] == 'Execution time, solver phase':
                    splt_line = row[31:].split(' ')
                    hours = 0; minutes=0; seconds=0
                    
                    for i in range(len(splt_line)):
                        if splt_line[i] == 'hours':
                            hours = float(splt_line[i-1])
                    
                        if splt_line[i] == 'minutes':
                            minutes = float(splt_line[i-1])
                            
                        if splt_line[i] == 'seconds':
                            seconds = float(splt_line[i-1])
                    run_time = hours+(minutes/60)+(seconds/3600)
                    
            PPW = minvs_h/fmax
            
        # open and extract jog history on talapas
        with open('jobs.txt', 'r') as fp_job:
            # read all lines using readline()
            lines = fp_job.readlines()
            
            for row in lines:
                # check if string present on a current line
                # extract maximum memory used
                if row[0:10] ==  job_id+'.0':
                    sw4_line = row.split(' ')
                        
                    for i in range(len(sw4_line)):
                        
                        if len(sw4_line[i]) >0 and sw4_line[i][-1] == 'K':
                            max_mem_used_GB = float(sw4_line[i][:-1])* float(n_cores)/1e-6                  
                                 
                # extract maximum memory requested
                if row[0:10] ==  job_id+'  ':
                    sw4_lin = row.split(' ')
                    for i in range(len(sw4_lin)):
                        if len(sw4_lin[i]) >0:
                            if sw4_lin[i][-1] == 'M' or sw4_lin[i][-1] == 'G':
                                max_mem_req = sw4_lin[i]
      
        # concatenate variables
        job_id_list = np.append(job_id_list,job_id)
        n_grdpts_list = np.append(n_grdpts_list,n_grdpts)
        total_seismic_moment_list = np.append(total_seismic_moment_list,total_seismic_moment)
        mag_list = np.append(mag_list,mag)
        n_sources_list = np.append(n_sources_list,n_sources)
        hmin_list = np.append(hmin_list,hmin)
        PPW_list = np.append(PPW_list,PPW)
        goal_time_list = np.append(goal_time_list,goal_time)
        n_timesteps_list = np.append(n_timesteps_list,n_timesteps)
        fmax_list = np.append(fmax_list,fmax)
        n_cores_list = np.append(n_cores_list,n_cores)
        max_mem_req_list = np.append(max_mem_req_list,max_mem_req)
        max_mem_used_GB_list = np.append(max_mem_used_GB_list,max_mem_used_GB)
        run_time_list = np.append(run_time_list,run_time)
        description_list = np.append(description_list,description)
        username_list = np.append(username_list,username)
    
######################### Put together dataframe ############################  
# First, make a dictionary for main part of dataframe:
dataset_dict = {'job_id':job_id_list,     
                'n_grdpts':n_grdpts_list, 
                'total_seismic_moment':total_seismic_moment_list,    
                'mag':mag_list,  
                'n_sources':n_sources_list,           
                'hmin':hmin_list,
                'PPW':PPW_list,  
                'goal_time':goal_time_list,               
                'n_timesteps':n_timesteps_list,             
                'fmax':fmax_list,            
                'n_cores':n_cores_list,           
                'max_mem_req':max_mem_req_list,
                'max_mem_used_GB':max_mem_used_GB_list,     
                'run_time':run_time_list, 
                'description':description_list,     
                'username':username_list}    
    
# Make main dataframe
flatfile_res_df = pd.DataFrame(data=dataset_dict)

# Save df to file
flatfile_res_df.to_csv(current_dir+'/SW4_scaling.csv',index=False)    
    
#os.system('rm -rf jobs.txt')
        
# ####################################################################
end = time.time()
time_elaps = end - start
if time_elaps < 60:
    print('Duration: '+str(round(time_elaps))+' seconds')
else:
    print('Duration: '+str(round(time_elaps/60))+' minutes')

# instead of using glob. It is slower though!
# n_gridpoints = subprocess.check_output('grep "(without ghost points):" '+ filename, shell=True).decode('ascii').split(' ')[-1]
# dirname = subprocess.check_output('grep "Making Output Directory" '+ filename, shell=True).decode('ascii').split(' ')[-1]
# minvs = subprocess.check_output('grep "<=  Vs" '+ filename, shell=True).decode('ascii').split(' ')[7]
