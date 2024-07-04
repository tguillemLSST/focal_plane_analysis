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

output_data='/sdf/home/t/tguillem/public_html/Operations_Rehearsal/OR4/debug/'
os.makedirs(output_data,exist_ok=True)

#butler access
#repo_1='/sdf/data/rubin/repo/ops-rehearsal-3-prep/butler.yaml'
#print('=========butler 1: ' + repo_1)
#repo='embargo_or4.yaml'
#print('=========butler : ' + repo)

#collections
#collection_1='LSSTComCamSim/raw/all'#runs/nightlyvalidation/20240403/d_2024_03_29/DM-43612'
collection_1='LSSTComCamSim/raw/all'
datasetType_1_1='raw'
#OR3
#collection_2='LSSTComCamSim/runs/nightlyvalidation/20240403/d_2024_03_29/DM-43612'
#OR4
collection_2='LSSTComCamSim/nightlyValidation'
#collection_2='LSSTComCamSim/prompt/output-2024-06-25'

datasetType_2_1='postISRCCD'
datasetType_2_2='calexp'
datasetType_2_3='calexpBackground'

#print('=========collections: ' + collection_1 + ' and ' + collection_2)
#butler_1 = dafButler.Butler(repo_1)
#registry_1 = butler_1.registry
repo = 'embargo_or4'
butler = dafButler.Butler(repo,writeable=False)
registry = butler.registry
#old
#butler = dafButler.Butler(repo)
#registry = butler.registry

#CCD
rafts=['R22']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
ccds=['S01', 'S02']
title='R22 S01-S02'

#get the list of exposures
#exposures=[]
#datasetRefs=list(registry_1.queryDatasets(datasetType=datasetType, instrument='LSSTComCamSim', collections=collection_1))
#raw_exps = sorted(registry_1.queryDatasets('raw', instrument='LSSTComCamSim',where="exposure.day_obs=20240117",collections=collection_1))
#print(raw_exps)

#by hand
#exposures = ['7024062500131']
exposures = ['7024062500500']
arr_raw_ccd=[]
arr_postISRCCD_ccd=[]
arr_calexp_ccd=[]
arr_calexpBackground_ccd=[]
for i_exp in range(len(exposures)):
    #create exposure directory
    str_exp = str(exposures[i_exp])
    #loop over all CCDs
    for i_raft in range(len(rafts)):

        #CCD-level analysis
        for i_ccd in range(len(ccds)):
            detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND visit.id=' + str_exp

            #calexpBackground
            datasetRefs=list(registry.queryDatasets(datasetType=datasetType_2_3, instrument='LSSTComCamSim', collections=collection_2, where=f"detector.full_name={detector}"))
            calexpBackground = butler.get(datasetRefs[0])
            arr_calexpBackground = calexpBackground.getImage().getArray()
            arr_calexpBackground_ccd.append(arr_calexpBackground)

        #continuity studies
        arr_calexpBackground_S01S02=np.concatenate((arr_calexpBackground_ccd[0],arr_calexpBackground_ccd[1]),axis=1)
        print(arr_calexpBackground_S01S02.shape)
        plt.figure(figsize=(20, 18))
        mean_calexpBackground=np.mean(arr_calexpBackground_S01S02)
        max_calexpBackground=np.max(arr_calexpBackground_S01S02)
        min_calexpBackground=np.min(arr_calexpBackground_S01S02)
        max_calexpBackground=4005
        min_calexpBackground=3995
        plt.imshow(arr_calexpBackground_S01S02, origin='lower', vmin=min_calexpBackground, vmax=max_calexpBackground, cmap='gist_gray')
        plt.ylim((0,2000))
        plt.xlim((3000,5000))
        plt.colorbar()
        plt.title(title)
        plt.savefig(output_data+'image_raft.png')
        plt.close()

        
print('DONE')
sys.exit()
