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

output_data='/sdf/home/t/tguillem/public_html/Operations_Rehearsal/OR4/calibrations/test/'
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
#collection_flat='u/jchiang/flat_70240417_w_2024_15/20240418T050546Z'
collection_flat='LSSTComCamSim/calib/DM-44910/or4Flats/flat.20240620a'

rafts=['R22']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']

#loop over all CCDs
for i_raft in range(len(rafts)):
    #create raft directory
    outpath_final_bias = output_data+'/bias/'+rafts[i_raft]
    os.makedirs(outpath_final_bias,exist_ok=True)
    shutil.copy2('index_files/index_raft_plots.html',outpath_final_bias+'/index.html')
    outpath_final_flat = output_data+'/flat/'+rafts[i_raft]
    os.makedirs(outpath_final_flat,exist_ok=True)
    shutil.copy2('index_files/index_raft_plots.html',outpath_final_flat+'/index.html')

    for i_ccd in range(len(ccds)):
        #detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
        detector_id=i_ccd
        #bias
        calib_1 = butler.get('bias', instrument='LSSTComCamSim', detector=detector_id, collections=collection_bias)
        title = rafts[i_raft]+' '+ccds[i_ccd]
        # get full array
        arr_1 = calib_1.getImage().getArray()
        plt.figure(figsize=[25,20])
        #plt.gcf().subplots_adjust(left = 0.05, bottom = 0.05,right = 0.95, top = 0.95)
        plt.imshow(arr_1, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
        #plt.imshow(arr_1, origin='lower', vmin=-1, vmax=1, cmap = 'hot')
        #plt.imshow(arr_1, origin='lower', vmin=0.005, vmax=0.020, cmap = 'hot')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_final_bias+'/image_'+ccds[i_ccd]+'.png')
        plt.close()

        #flat
        calib_flat = butler.get('flat', instrument='LSSTComCamSim', physical_filter='i_06', detector=detector_id, collections=collection_flat)
        #calib_flat_md = butler_1.get('flat.metadata', instrument='LSSTComCamSim', physical_filter='i_06', detector=detector_id)
        #print(calib_falt_md)
        # get full array
        arr_flat = calib_flat.getImage().getArray()
        plt.figure(figsize=[25,20])
        plt.imshow(arr_flat, origin='lower', vmin=0.95, vmax=1.05, cmap = 'hot')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_final_flat+'/image_'+ccds[i_ccd]+'.png')
        plt.close()
        
print('DONE')
sys.exit()
