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
import pickle
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
from geom import *
from scipy.interpolate import make_interp_spline
from scipy.interpolate import splrep, splev
from scipy import interpolate

print('Configuration arguments: ', str(sys.argv))
str_raft = str(sys.argv[1])
str_all_sensors = str(sys.argv[2])

output_data='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/full_flat_correction/test/'
os.makedirs(output_data,exist_ok=True)

#butler access
repo='embargo_new'
#repo='/repo/main'
print('=========butler: ' + repo)

butler = dafButler.Butler(repo,writeable=False)
registry = butler.registry

#for a set of batch jobs
rafts=[str_raft]
ccds=[str_all_sensors]

#CCD
#rafts=['R01','R02','R03','R10','R11','R12','R13','R14','R20','R21','R22','R23','R24','R30','R31','R32','R33','R34','R41','R42','R43']
#rafts=['R01']
#ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
#ccds=['S01']

camera = LsstCam.getCamera()

#combined correction
inpath='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/vignetting_illumination_combined/model/'
file = open(inpath+'vignetting_splines.pkl', 'rb')
pclk_TR = pickle.load(file)
pclk_BR = pickle.load(file)
pclk_TL = pickle.load(file)
pclk_BL = pickle.load(file)

def bspl_TR(x):
    return interpolate.splev(x, pclk_TR)
def bspl_BR(x):
    return interpolate.splev(x, pclk_BR)
def bspl_BL(x):
    return interpolate.splev(x, pclk_BL)
def bspl_TL(x):
    return interpolate.splev(x, pclk_TL)

#########
def illumination_vignetting_correction(x, y):
    """
    Parameters
    ----------
    x, y : pixel position in mm

    Returns
    -------
    c : full correction
    """
    theta_deg = np.arctan2(y, x) * 180 / np.pi
    radius=np.sqrt(x**2+y**2)
    c=1.0
    #stop at some radius
    if(radius>364.75):
        radius=364.75
    #at a given radius, linear interpolation versus the angle between the 4 references
    if(theta_deg>=-45 and theta_deg<45):
        c = np.abs(-45-theta_deg)/90*bspl_TR(radius) + np.abs(45-theta_deg)/90*bspl_BR(radius)
    elif(theta_deg>=45 and theta_deg<135):
        c = np.abs(45-theta_deg)/90*bspl_TL(radius) + np.abs(135-theta_deg)/90*bspl_TR(radius)
    elif(theta_deg>=135 or theta_deg<-135):
        if(theta_deg>0):
            c = np.abs(45+180-theta_deg)/90*bspl_TL(radius) + np.abs(theta_deg-135)/90*bspl_BL(radius)
        if(theta_deg<0):
            c = np.abs(-135-theta_deg)/90*bspl_TL(radius) + np.abs(90-(-135-theta_deg))/90*bspl_BL(radius)
    elif(theta_deg>=-135 or theta_deg<-45):
         c = np.abs(-45-theta_deg)/90*bspl_BL(radius) + np.abs(-135-theta_deg)/90*bspl_BR(radius)
    return c
#########

#get illumination center
position_center_mm = pixel_to_focal(float(2000),float(2000),camera['R22_S11'])
x_center_mm=position_center_mm[0][0]
y_center_mm=position_center_mm[1][0]

#loop over all CCDs
for i_raft in range(len(rafts)):
    outpath_final_flat = output_data+'/'+rafts[i_raft]
    os.makedirs(outpath_final_flat,exist_ok=True)
    shutil.copy2('index_files/index_raft_plots.html',outpath_final_flat+'/index.html')
    
    for i_ccd in range(len(ccds)):
        detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
        detector_full_name = rafts[i_raft] + '_' +ccds[i_ccd]
            
        title = rafts[i_raft]+' '+ccds[i_ccd]
        ##############HACK : test directly over a flat exposure
        detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=2025071100888'
        datasetRefs=list(registry.queryDatasets(datasetType='post_isr_image', instrument='LSSTCam', collections='LSSTCam/runs/nightlyValidation', where=f"detector.full_name={detector}"))
        postISRCCD = butler.get(datasetRefs[0])
        arr_flat = postISRCCD.getImage().getArray()
        v_vmin=0.8
        v_vmax=1.05
        n_bins_x=len(arr_flat[0,:])
        n_bins_y=len(arr_flat[:,0])
        arr_2nd = np.copy(arr_flat)
        arr_2nd[:,:]=0
        for i_x in range(n_bins_x):
            if(i_x%10!=8):
                continue
            for i_y in range(n_bins_y):
                if(i_y%10!=8):
                    continue
                position_mm = pixel_to_focal(float(i_x),float(i_y),camera[detector_full_name])
                x_mm=position_mm[0][0]
                y_mm=position_mm[1][0]
                distance=np.sqrt((x_mm-x_center_mm)**2+(y_mm-y_center_mm)**2)
                vignetting_corr = illumination_vignetting_correction(x_mm,y_mm)
                arr_2nd[i_y,i_x]=vignetting_corr

        plt.figure(figsize=[25,20])
        #plt.imshow(arr_flat, origin='lower', vmin=v_vmin, vmax=v_vmax, cmap = 'hot')
        plt.imshow(arr_2nd, origin='lower', vmin=v_vmin, vmax=v_vmax, cmap = 'hot')
        #plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_final_flat+'/image_'+ccds[i_ccd]+'.png', bbox_inches='tight')
        plt.close()
        
print('DONE')
sys.exit()
