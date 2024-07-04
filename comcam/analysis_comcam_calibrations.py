#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: plot calibrations from the butler
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

output_data='/sdf/home/t/tguillem/public_html/Operations_Rehearsal/OR4/calibrations/studies/'
os.makedirs(output_data,exist_ok=True)
os.makedirs(output_data+'/results',exist_ok=True)

#butler access
#repo='/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/'
repo='/sdf/group/rubin/repo/embargo/butler.yaml'
print('=========butler: ' + repo)
butler = dafButler.Butler(repo)
registry = butler.registry

#OR3
collection_bias='u/jchiang/bias_70240217_w_2024_07/20240218T190659Z'
#collection_flat='u/jchiang/flat_70240217_w_2024_07/20240218T195122Z'
#OR4
collection_dark='u/jchiang/dark_70240217_w_2024_07/20240218T191310Z'
#collection_flat='u/jchiang/flat_70240417_w_2024_15/20240418T050546Z'
collection_flat='LSSTComCamSim/calib/DM-44910/or4Flats/flat.20240620a'

rafts=['R22']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
ccds=['S11']

#loop over all CCDs
for i_raft in range(len(rafts)):
    #create raft directory
    outpath_final = output_data+rafts[i_raft]+'/'
    os.makedirs(outpath_final,exist_ok=True)
    
    for i_ccd in range(len(ccds)):
        #detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
        detector_id=i_ccd
        #bias
        calib_bias = butler.get('bias', instrument='LSSTComCamSim', detector=detector_id, collections=collection_bias)
        title = rafts[i_raft]+' '+ccds[i_ccd]
        # get full array
        arr_bias = calib_bias.getImage().getArray()

        #dark
        calib_dark = butler.get('dark', instrument='LSSTComCamSim', detector=detector_id, collections=collection_dark)
        title = rafts[i_raft]+' '+ccds[i_ccd]
        # get full array
        arr_dark = calib_dark.getImage().getArray()

        #per-pixel correlation plot ==> too noisy
        #plt.figure()
        #plt.scatter(arr_bias.flat[0:100],arr_dark.flat[0:100], marker='.',color = 'blue', s=10, alpha=0.3, label='pixels')
        #plt.xlim([-5, 5])
        #plt.ylim([0, 0.05])
        #plt.xlabel('bias')
        #plt.ylabel('dark')
        #plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
        #plt.legend()
        #plt.savefig(outpath_final+'bias_dark_correlation.png', bbox_inches='tight')

        #get per amplifier results
        detector = calib_bias.getDetector()
        #get size
        amplifier = detector['C00']
        sub_im0 = calib_bias.getMaskedImage()[amplifier.getBBox()]
        arr_amp0 = sub_im0.getImage().getArray()
        im_x_size = arr_amp0.shape[1]
        im_y_size = arr_amp0.shape[0]
        var_total =  np.zeros(16)
        mean_total = np.zeros(16)
        mean_corner = np.zeros(16)
        var_corner = np.zeros(16)
        var_line = np.zeros((16,im_y_size))
        mean_line = np.zeros((16,im_y_size))
        var_column =  np.zeros((16,im_x_size))
        mean_column = np.zeros((16,im_x_size))
        mean_line_bias = np.zeros((16,im_y_size))
        mean_line_dark = np.zeros((16,im_y_size))
        amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
        amps=['C03']
        detector = calib_bias.getDetector()
        up_line = im_y_size - 100
        up_column = im_x_size - 100
        #print(up_line)
        #print(up_column)
        for i_amp in range(len(amps)):
                amplifier = detector[amps[i_amp]]
                amplifier_ref = detector['C16']
                sub_im_bias = calib_bias.getMaskedImage()[amplifier_ref.getBBox()]
                arr_amp_bias = sub_im_bias.getImage().getArray()
                sub_im_dark = calib_dark.getMaskedImage()[amplifier.getBBox()]
                arr_amp_dark = sub_im_dark.getImage().getArray()
                mean_line_bias[i_amp,:] = np.mean(arr_amp_bias,axis=1)
                mean_line_dark[i_amp,:] = np.mean(arr_amp_dark,axis=1)
                plt.figure()
                plt.scatter(mean_line_bias[i_amp],mean_line_dark[i_amp], marker='.',color = 'blue', s=10, alpha=0.3, label='lines')
                plt.xlim([-2, 2])
                plt.ylim([-0.02, 0.04])
                plt.xlabel('bias')
                plt.ylabel('dark')
                plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
                plt.legend()
                plt.savefig(outpath_final+'bias_dark_correlation.png', bbox_inches='tight')
        
print('DONE')
sys.exit()
