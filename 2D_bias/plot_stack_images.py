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
import matplotlib
import os
import shutil
from astropy.visualization import (ImageNormalize,PercentileInterval,HistEqStretch)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc 
import lsst.daf.butler as dafButler

print('Configuration arguments: ', str(sys.argv))
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
repo_path=str(sys.argv[4])
output_data=str(sys.argv[5])

os.makedirs(output_data,exist_ok=True)

run=str_run_all
raft=str_raft
ccd=str_all_sensors

#butler access
repo='/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/'
print('=========butler: ' + repo)
###calibrations
plot_calibrations = False
#collection = "LSSTCam/calib"
collection = "u/tguillem/DM-37455/master_bias_1D_20230116a"#DM-37455/bias.20230116a"
calibType='bias'
###user collections
plot_images = True
#collection ="u/tguillem/DM-37455/biasCorr.20230117a"
collection ="u/tguillem/DM-37455/biasCorr_1D.20230118a"
datasetType='cpBiasProc'
butler = dafButler.Butler(repo,collections=collection)
registry = butler.registry

#runs=('13162')
rafts=['R01','R02','R03','R10','R11','R12','R13','R14','R20','R21','R22','R23','R24','R30','R31','R32','R33','R34','R41','R42','R43']
#rafts=['R14']
ccds=['S00','S01','S02','S10','S11','S12','S20','S21','S22']
#ccds=['S22']

if plot_calibrations:
    detectorId = 1
    calib = butler.get(calibType, instrument='LSSTCam', detector=detectorId)
    #title = rafts[i_raft]+' '+ccds[i_ccd]+' exp='+str(exposure)
    plt.figure()
    # get full array
    arr = calib.getImage().getArray()
    plt.imshow(arr, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
    plt.colorbar()
    #plt.title(title)
    plt.savefig(output_data+'/image.png')
    plt.close()

if plot_images:
    #loop over all CCDs
    for i_raft in range(len(rafts)):
        #create raft directory
        outpath_final = output_data+rafts[i_raft]
        os.makedirs(outpath_final,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final+'/index.html')
        
        for i_ccd in range(len(ccds)):
            detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\''
            #print(detector)
            #datasetRefs=list(registry.queryDatasets(datasetType=datasetType, instrument='LSSTCam', collections=collection, where=f"detector.full_name={detector}"))
            datasetRefs=list(registry.queryDatasets(datasetType=datasetType, instrument='LSSTCam', collections=collection))
            print(datasetRefs[0])

            # get image
            print(datasetRefs[0].dataId.full)
            exposure = datasetRefs[0].dataId['exposure']
            raw = butler.get(datasetRefs[0])
            #title = 'Raft '+rafts[i_raft]+' CCD '+ccds[i_ccd]+' exposure '+str(exposure)
            title = rafts[i_raft]+' '+ccds[i_ccd]+' exp='+str(exposure)
            plt.figure()
            # get full array
            arr = raw.getImage().getArray()
            plt.imshow(arr, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final+'/image_'+ccds[i_ccd]+'.png')
            plt.close()

print('DONE')
sys.exit()
