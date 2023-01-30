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
from conversion import *

print('Configuration arguments: ', str(sys.argv))
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
repo_path=str(sys.argv[4])
output_data=str(sys.argv[5])
###HACK###
collection_selected=repo_path

os.makedirs(output_data,exist_ok=True)
os.makedirs(output_data+'/results',exist_ok=True)

run=str_run_all
raft=str_raft
ccd=str_all_sensors

#butler access
repo='/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/'
print('=========butler: ' + repo)

#select collection
collection=collection_selected
###calibrations
plot_calibrations = False
calibType='bias'
#collection = "LSSTCam/calib"
#collection = "u/tguillem/DM-37455/master_bias_1D_20230123c"#DM-37455/bias.20230116a"
###user collections
plot_images = True
datasetType='cpBiasProc'
#collection ="u/tguillem/DM-37455/biasCorr.20230117a"
#collection ="u/tguillem/DM-37455/biasCorr_1D.20230118a"
#collection="u/tguillem/DM-37455/run_13162_2D_MB.20230119a"
#collection="u/tguillem/DM-37455/run_13162_2D.20230120b"
print('=========collection: ' + collection_selected)
butler = dafButler.Butler(repo,collections=collection)
registry = butler.registry

#runs=('13162')
rafts=['R01','R02','R03','R10','R11','R12','R13','R14','R20','R21','R22','R23','R24','R30','R31','R32','R33','R34','R41','R42','R43']
rafts=['R14']
ccds=['S00','S01','S02','S10','S11','S12','S20','S21','S22']
#ccds=['S22']

if plot_calibrations:
    #loop over all CCDs
    for i_raft in range(len(rafts)):
        #create raft directory
        outpath_final = output_data+rafts[i_raft]
        os.makedirs(outpath_final,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final+'/index.html')

        for i_ccd in range(len(ccds)):
            detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
            calib = butler.get(calibType, instrument='LSSTCam', detector=detector_id)
            title = rafts[i_raft]+' '+ccds[i_ccd]
            plt.figure()
            # get full array
            arr = calib.getImage().getArray()
            plt.imshow(arr, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final+'/image_'+ccds[i_ccd]+'.png')
            plt.close()

if plot_images:
    #loop over all CCDs
    for i_raft in range(len(rafts)):
        #create raft directory
        outpath_final = output_data+rafts[i_raft]+''
        os.makedirs(outpath_final,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final+'/index.html')
        
        for i_ccd in range(len(ccds)):
            #HACK: pass the seq_num here
            detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=3021121200160'
            #without seq_num
            #detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\''
            print(detector)
            datasetRefs=list(registry.queryDatasets(datasetType=datasetType, instrument='LSSTCam', collections=collection, where=f"detector.full_name={detector}"))
            print(datasetRefs)
            
            # get image
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

            #loop over amplifiers
            
            #get per amplifier results
            detector = raw.getDetector()
            #get size
            amplifier = detector['C00']
            sub_im0 = raw.getMaskedImage()[amplifier.getBBox()]
            arr_amp = sub_im0.getImage().getArray()
            im_x_size = arr_amp.shape[1]
            im_y_size = arr_amp.shape[0]
            var_total =  np.zeros(16)
            mean_total = np.zeros(16)
            var_line = np.zeros((16,im_y_size))
            mean_line = np.zeros((16,im_y_size))
            var_column =  np.zeros((16,im_x_size))
            mean_column = np.zeros((16,im_x_size))
            amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
            detector = raw.getDetector()
            for i_amp in range(len(amps)):
                amplifier = detector[i_amp]
                sub_im0 = raw.getMaskedImage()[amplifier.getBBox()]
                arr_amp = sub_im0.getImage().getArray()
                #plt.figure()
                #plt.imshow(arr_amp, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
                #plt.colorbar()
                #plt.savefig(outpath_final+'/raw_amp.png')
                var_total[i_amp] = np.var(arr_amp)
                mean_total[i_amp] = np.mean(arr_amp)
                var_line[i_amp,:] = np.var(arr_amp,axis=1)
                mean_line[i_amp,:] = np.mean(arr_amp,axis=1)
                var_column[i_amp,:] = np.var(arr_amp,axis=0)
                mean_column[i_amp,:] = np.mean(arr_amp,axis=0)
                
            #write results    
            t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column'), meta={'name': 'Variances'})
            #print(t_variance)
            t_variance.write(output_data + '/results/Variances_' + rafts[i_raft] + '_'+ ccds[i_ccd] + '.fits', overwrite=True)
            
print('DONE')
sys.exit()
