#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: plot defects
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
import lsst.afw.image as afwImage
import numpy as np
from conversion import *
from matplotlib.colors import ListedColormap

#butler access
repo = '/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/'
print('=========butler: ' + repo)

outpath = "/sps/lsst/users/tguillem/web/debug/defects_stack/flat/"
os.makedirs(outpath,exist_ok=True)

#rafts = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43', 'R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
#ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
rafts = ['R14']
ccds=['S21']
#collection = ['u/lsstccs/defects_13391_w_2023_24']
#collection = ['u/lsstccs/eo_dark_defects_13391_w_2023_24']
collection = ['u/tguillem/run_6_validation/defects_13391_flat_20230906a']

butler = dafButler.Butler(repo, collections=collection)
#3 objects: defects/eoBrightDefects/eoDarkDefects

for i_raft in range(len(rafts)):
    outpath_raft = outpath+rafts[i_raft]+'/'
    os.makedirs(outpath_raft,exist_ok=True)
    for i_ccd in range(len(ccds)):
        detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
        detector_name = dict_detector[str(detector_id)]
        butler = dafButler.Butler(repo, collections=collection)
        defects = butler.get('defects',instrument='LSSTCam',detector=detector_id)
        
        #code snippet from C. Waters
        #e2v
        #ncol=4096
        #nrow=4004
        #ITL
        ncol=4072
        nrow=4000
        realization = afwImage.MaskedImageI(ncol, nrow)  # instantiate an image of the proper size
        defects.maskPixels(realization)   # apply mask to that image.
        defectArray = realization.getMask().getArray()  # get array for plotting/etc.
        print(defectArray.shape)
        np.set_printoptions(threshold=sys.maxsize)
        #print(defectArray)

        for i in range(nrow):
            print('========row ' + str(i))
            for j in range(ncol):
                if(defectArray[i][j]==1):
                    print('row ' + str(i) + ' column ' + str(j) + ' :' + str(defectArray[i][j])) 
                    
        #plot
        defectArray_copy = defectArray[:,:]#5:10,100:120]
        cmapmine = ListedColormap(['w', 'r'], N=2)
        plt.figure()
        plt.imshow(defectArray_copy, origin='lower', vmin=0, vmax=1, cmap='tab10') #cmap='tab10' / cmap=cmapmine
        #bug
        #plt.imshow(defectArray_copy, origin='lower', vmin=0, vmax=1, cmap=cmapmine)# interpolation='none') #extent=[0,4072,0,100])
        plt.colorbar()
        plt.title(collection[0]+'\n'+ detector_name)
        plt.savefig(outpath_raft+'image_'+ccds[i_ccd]+'.png')
        plt.close()

sys.exit()
