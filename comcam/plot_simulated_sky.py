#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: plot images directly from the butler
##########

# system imports
import time
import sys
from sys import exit
import glob
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table
from scipy.stats import skew
import matplotlib
import os
import shutil
from astropy.visualization import (ImageNormalize,PercentileInterval,HistEqStretch)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc 
import lsst.daf.butler as dafButler
#from conversion import *
from lsst.ip.isr import AssembleCcdTask
assembleTask = AssembleCcdTask()

#output_data='/sdf/home/t/tguillem/public_html/Operations_Rehearsal/OR4/simulation/04000039/'
output_data='/sdf/home/t/tguillem/public_html/Operations_Rehearsal/OR4/simulation/00013404/'
os.makedirs(output_data,exist_ok=True)

rafts=['R22']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
det_ccds=['000' ,'001' ,'002' ,'003' ,'004' ,'005' ,'006' ,'007' ,'008']
arr_sky_ccd=[]
mean_sky_central_ccd=[]

os.makedirs(output_data+'R22',exist_ok=True)
shutil.copy2('index_files/index_raft_plots.html',output_data+'R22/index.html')
os.makedirs(output_data+'raft_level',exist_ok=True)

#inpath = '/sdf/data/rubin/user/jchiang/ops_rehearsal_3/vignetted_flats/flats/'
inpath = '/sdf/data/rubin/shared/ops-rehearsals/ops-rehearsal-4/image_sims/rehearsal_nights/pass_1/'

for i_raft in range(len(rafts)):
    for i_ccd in range(len(ccds)):
        #path_ccd = inpath + '04000039/eimage_04000039-0-r-' + rafts[i_raft] + '_' +ccds[i_ccd] + '-det' + det_ccds[i_ccd] + '.fits'
        path_ccd = inpath + '00013404/eimage_00013404-0-i-' + rafts[i_raft] + '_' +ccds[i_ccd] + '-det' + det_ccds[i_ccd] + '.fits'
        file_ccd=pyfits.open(path_ccd)
        arr_sky = file_ccd[0].data
        arr_sky_ccd.append(arr_sky)
        mean_sky=np.mean(arr_sky_ccd)
        std_sky=np.std(arr_sky_ccd)
        max_sky=mean_sky+1*std_sky
        min_sky=mean_sky-1*std_sky
        plt.figure(figsize=(20, 18))
        plt.imshow(arr_sky, origin='lower', vmin=min_sky, vmax=max_sky, cmap='gist_gray')
        plt.colorbar()
        #plt.title(title)
        plt.savefig(output_data+'R22/image_'+ccds[i_ccd]+'.png')
        plt.close()    

        mean_sky_central=np.mean(arr_sky_ccd[1000:1100,1000:1100])
        mean_sky_central_ccd.append(mean_sky_central)
        
sys.exit()
        
#raft-level
arr_sky_raft_S0=np.concatenate((arr_sky_ccd[0],arr_sky_ccd[1],arr_sky_ccd[2]),axis=1)
arr_sky_raft_S1=np.concatenate((arr_sky_ccd[3],arr_sky_ccd[4],arr_sky_ccd[5]),axis=1)
arr_sky_raft_S2=np.concatenate((arr_sky_ccd[6],arr_sky_ccd[7],arr_sky_ccd[8]),axis=1)
arr_sky_raft=np.concatenate((arr_sky_raft_S0,arr_sky_raft_S1,arr_sky_raft_S2),axis=0)
#print(arr_sky_raft)
plt.figure(figsize=(20, 18))
mean_sky=np.mean(arr_sky_raft)
max_sky=np.max(arr_sky_raft)+100
min_sky=np.min(arr_sky_raft)-100
max_sky=3300
min_sky=3100
plt.imshow(arr_sky_raft, origin='lower', vmin=min_sky, vmax=max_sky, cmap='gist_gray')
plt.colorbar()
#plt.title(title)
plt.savefig(output_data+'/raft_level/'+'image_raft.png')
plt.close()

#1D histogram
histo_arr_sky_raft = arr_sky_raft.flat
plt.figure()
plt.hist(histo_arr_sky_raft, range=[min_sky,max_sky], bins=1000, label='ADU', histtype='step', color = 'black')
#plt.ylim([0.8, 1.1])
plt.xlabel('ADU')
plt.ylabel('pixels')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.savefig(output_data+'/raft_level/'+'histo_raft.png')        
        
sys.exit()
