#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 12:11:38 2023

@author: rshimony
"""

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(11, 3.5))
# fig.tight_layout(rect=[0, 0.2, 1, 0.985])
fig.subplots_adjust(wspace=0.005)
# fig.suptitle('rfile Figures', fontsize=20 , x=0.5,y=0.95)

ax1 = fig.add_subplot(2,3,1)
Image1 = plt.imread('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/fig.profile_vs_basin.png')
ax1.imshow(Image1)
ax1.axis('off')
ax1.set_title('Vs profile - basin' ,size=8)

ax2 = fig.add_subplot(2,3,2)
Image2 = plt.imread('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/fig.crossection2_basin.png')
ax2.imshow(Image2)
ax2.axis('off')
ax2.set_title('S-N cross section - basin' ,size=8)

ax3 = fig.add_subplot(2,3,3)
Image3 = plt.imread('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/fig.crossection_basin.png')
ax3.imshow(Image3)
ax3.axis('off')
ax3.set_title('W-E cross section - basin' ,size=8)

ax4 = fig.add_subplot(2,3,4)
Image4 = plt.imread('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfiletopo/fig.profile_vs_steph.png')
ax4.imshow(Image4)
ax4.axis('off')
ax4.set_title('Vs profile - stephenson' ,size=8)

ax5 = fig.add_subplot(2,3,5)
Image5 = plt.imread('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfiletopo/fig.crossection2_steph.png')
ax5.imshow(Image5)
ax5.axis('off')
ax5.set_title('S-N cross section - stephenson' ,size=8)

ax6 = fig.add_subplot(2,3,6)
Image6 = plt.imread('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfiletopo/fig.crossection_steph.png')
ax6.imshow(Image6)
ax6.axis('off')
ax6.set_title('W-E cross section - stephenson' ,size=8)

plt.savefig('/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/rfile_comp_fig.jpg'
                    ,dpi=500,bbox_inches='tight')

#%%

plt.figure(figsize=(8,6))
tiff = '/Users/rshimony/Documents/My_Notebooks/create_vel_model_rfilebasin/tiff/elev_layer2.tif'
tiff_data = plt.imread(tiff)
plt.imshow(tiff_data,cmap='viridis')
plt.show()

































