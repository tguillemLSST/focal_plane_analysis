#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: extract gradient+vignetting from screen flats
##########

# system imports
import time
import sys
from sys import exit
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table
from scipy import stats
import matplotlib
import os
import shutil
from astropy.visualization import (ImageNormalize,PercentileInterval,HistEqStretch)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc 
import lsst.daf.butler as dafButler
from conversion import *
from geom import *
import pickle
#to fit
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.integrate import simps
from scipy.stats import norm
from scipy.stats import lognorm
from scipy import optimize
from scipy.optimize import minimize
from numpy import exp, linspace, random

print('Configuration arguments: ', str(sys.argv))
str_raft = str(sys.argv[1])
str_all_sensors = str(sys.argv[2])

output_data='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/flat_r/vignetting_illumination_combined/BL/itl/'
os.makedirs(output_data,exist_ok=True)

#butler access
#repo='embargo_new'
repo='/repo/main'
print('=========butler: ' + repo)

butler = dafButler.Butler(repo,writeable=False)
registry = butler.registry

#for a set of batch jobs
rafts=[str_raft]
ccds=[str_all_sensors]

#4 regions to run separately (actually 5 because BL is splitted between e2v and ITL)
#TR
#rafts=['R22', 'R23', 'R24', 'R32', 'R33', 'R34']
#BR
#rafts=['R22', 'R23', 'R24', 'R12', 'R13', 'R14']
#BL
#e2v
#rafts=['R21', 'R22', 'R11', 'R12']
#ITL
rafts=['R20', 'R10', 'R01', 'R02']
#TL
#rafts=['R21', 'R22', 'R30', 'R31', 'R32']
ccds=['S00', 'S01', 'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
#ccds=['S01']

#camera geometry
camera = LsstCam.getCamera()

os.makedirs(output_data+'/'+collection_flats[i_filter]+'/results',exist_ok=True)

vignetting_radius=[]
vignetting_value=[]
#loop over all CCDs
for i_raft in range(len(rafts)):
    outpath_final_flat = output_data+'/'+rafts[i_raft]
    os.makedirs(outpath_final_flat,exist_ok=True)
        
    for i_ccd in range(len(ccds)):
        detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
        detector_full_name = rafts[i_raft] + '_' +ccds[i_ccd]
        if(detector_full_name=='R30_S12'):
            continue
        title = rafts[i_raft]+' '+ccds[i_ccd]
        ##############for now : test directly over a flat exposure
        detector_full = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=2025072700161'
        datasetRefs=list(registry.queryDatasets(datasetType='post_isr_image', instrument='LSSTCam', collections='LSSTCam/runs/nightlyValidation', where=f"detector.full_name={detector_full}"))
        postISRCCD = butler.get(datasetRefs[0])
        arr_flat = postISRCCD.getImage().getArray()
            
        #gradient+vignetting profile versus radius
        #print(arr_flat.shape)
        n_bins_x=len(arr_flat[0,:])
        n_bins_y=len(arr_flat[:,0])
        #to check secondary bands
        arr_2nd = np.copy(arr_flat)
        arr_2nd[:,:]=0
        n_bins_x=len(arr_flat[0,:])
        n_bins_y=len(arr_flat[:,0])
        #can use only 1% of the pixels to speed up things
        for i_x in range(n_bins_x):
            if(i_x%100!=5):
                continue
            for i_y in range(n_bins_y):
                if(i_y%100!=5):
                        continue
                position_mm = pixel_to_focal(float(i_x),float(i_y),camera[detector_full_name])
                x_mm=position_mm[0][0]
                y_mm=position_mm[1][0]
                radius=np.sqrt(x_mm**2+y_mm**2)
                value=arr_flat[i_y,i_x]
                vignetting_radius.append(radius)
                vignetting_value.append(value)
                if(radius<350 and value<0.45):
                    arr_2nd[i_y,i_x]=1
                #print(vignetting_radius)
                #print(vignetting_value)

        #plot scatter and median versus radius
        bin_means, bin_edges, binnumber = stats.binned_statistic(vignetting_radius, vignetting_value, statistic='median', range=[(0, 370)], bins=740)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width/2
        #renormalize bin_means (to start at 1 at r=0)
        average_vignetting = []
        for i in range(len(bin_means)):
            if(bin_means[i]>0.8):
                average_vignetting.append(bin_means[i])
            if(len(average_vignetting)>10):
                break
        #print(average_vignetting)
        average_vignetting_v = np.mean(average_vignetting[1:10])
        arr_vignetting_value = np.array(vignetting_value)
        arr_vignetting_value = 1/average_vignetting_v*arr_vignetting_value
        bin_means = 1/average_vignetting_v*bin_means

        plt.figure()
        plt.scatter(vignetting_radius,arr_vignetting_value, marker='.',color = 'blue', s=0.3, alpha=1, label='pixels')
        plt.scatter(bin_centers[1:-1], bin_means[1:-1], marker='.', color = 'red', s=100, label='Model')
        plt.plot(bin_centers, bin_means, color = 'red')
        plt.xlim([0,370])
        plt.ylim([0,1.2])
        plt.xlabel('radius')
        plt.ylabel('vignetting')
        plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
        plt.legend()
        plt.savefig(outpath_final_flat+'/vignetting_'+ccds[i_ccd]+'.png', bbox_inches='tight')

        #save results to parametrize the model later
        t_variance = Table([bin_means,bin_centers],names=('bin_means','bin_centers'),meta={'name': 'Variances'})
        t_variance.write(outpath_final_flat+'/Variance_'+ccds[i_ccd]+'.fits', overwrite=True)
            
print('DONE')
sys.exit()
