#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: build a vignetting model
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

inpath='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/vignetting/20250618/u/erykoff/LSSTCam/calib/DM-51432/skyflat-griz/'
file_S22_g=inpath+'flatGen-z.20250617a/R34/Variance_S22.fits'
file_S22_r=inpath+'flatGen-r.20250617a/R34/Variance_S22.fits'
file_S22_i=inpath+'flatGen-g.20250617a/R34/Variance_S22.fits'
file_S22_z=inpath+'flatGen-i.20250617a/R34/Variance_S22.fits'
variance_S22_g = Table.read(file_S22_g)
variance_S22_r = Table.read(file_S22_r)
variance_S22_i = Table.read(file_S22_i)
variance_S22_z = Table.read(file_S22_z)

outpath='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/vignetting/20250618/model/'
os.makedirs(outpath,exist_ok=True)

#check points
#for i in range(len(variance_S22_g['bin_centers'])):
#    print(str(i) + ' ' + str(variance_S22_g['bin_centers'][i]) + ' | ' + str(variance_S22_g['bin_means'][i]))
#sys.exit()

#correction factors (old method)
#corr_S21 = variance_S11['bin_means'][290]/variance_S21['bin_means'][290]
#corr_S22 = corr_S21 * variance_S21['bin_means'][330]/variance_S22['bin_means'][330]
##combine in a single array    
#variance_combined = np.copy(variance_S11['bin_means'])
#for i in range(len(variance_S11['bin_means'])):
#    if(i>290 and i<330):
#        variance_combined[i]=corr_S21*variance_S21['bin_means'][i]
#    if(i>329):
#        variance_combined[i]=corr_S22*variance_S22['bin_means'][i]
#    #if(i==369):
#    #    variance_combined[i]=0
#for i in range(len(variance_S11['bin_centers'])):
#    print(str(i) + ' | ' + str(variance_S11['bin_centers'][i]) + ' | ' + str(variance_combined[i]))
        
#plot with only data points
#plt.figure()
#plt.scatter(variance_S11['bin_centers'][259:308], variance_S11['bin_means'][259:308], marker='.', color = 'black', s=40, label='Model from R34_S11')
#plt.scatter(variance_S21['bin_centers'][279:331], corr_S21*variance_S21['bin_means'][279:331], marker='.', color = 'blue', s=40, label='Model from R34_S21')
#plt.scatter(variance_S22['bin_centers'][315:], corr_S22*variance_S22['bin_means'][315:], marker='.', color = 'red', s=40, label='Model from R34_S22')
#plt.xlim([250,370])
#plt.ylim([0,1.2])
#plt.xlabel('Radius [mm]')
#plt.ylabel('Vignetting model')
#plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
#plt.legend()
#plt.title('Vignetting model from an i-band sky exposure')
#plt.savefig(outpath+'vignetting_model.png', bbox_inches='tight')

#batoid model
vignetting_coeffs = np.asarray(
    [
        -2.04748298e+12,
        4.62036195e+12,
        -4.55318392e+12,
        2.55519946e+12,
        -8.86592878e+11,
        1.89254514e+11,
        -2.11087631e+10,
        -2.68228152e+08,
        4.87993883e+08,
        -8.03764403e+07,
        6.99808127e+06,
        -3.58577957e+05,
        1.05491604e+04,
        -1.60565953e+02,
        9.96009337e-01,
        9.98941038e-01,
    ],)

#single spline
from scipy.interpolate import make_interp_spline
bspl_g = make_interp_spline(variance_S22_g['bin_centers'][0:729], variance_S22_g['bin_means'][0:729], k=3)
bspl_r = make_interp_spline(variance_S22_r['bin_centers'][0:729], variance_S22_r['bin_means'][0:729], k=3)
bspl_i = make_interp_spline(variance_S22_i['bin_centers'][0:729], variance_S22_i['bin_means'][0:729], k=3)
bspl_z = make_interp_spline(variance_S22_z['bin_centers'][0:729], variance_S22_z['bin_means'][0:729], k=3)
plt.figure()
#plt.scatter(variance_S22_g['bin_centers'][0:729], variance_S22_g['bin_means'][0:729], marker='.', color = 'black', s=40, label='g band')
#plt.scatter(variance_S22_r['bin_centers'][0:729], variance_S22_r['bin_means'][0:729], marker='.', color = 'red', s=40, label='r band')
#plt.scatter(variance_S22_i['bin_centers'][0:729], variance_S22_i['bin_means'][0:729], marker='.', color = 'blue', s=40, label='i band')
#plt.scatter(variance_S22_z['bin_centers'][0:729], variance_S22_z['bin_means'][0:729], marker='.', color = 'brown', s=40, label='z band')
plt.xlim([230,370])
xx = np.linspace(0, 364.25, 1000)
plt.plot(xx, bspl_g(xx),color = 'black', label='g band')
plt.plot(xx, bspl_r(xx),color = 'red', label='r band')
plt.plot(xx, bspl_i(xx),color = 'blue', label='i band')
plt.plot(xx, bspl_z(xx),color = 'brown', label='z band')
plt.plot(xx, np.polyval(vignetting_coeffs, xx/1000.), color = 'purple', linestyle='dashed', label='Batoid')
plt.ylim([0,1.2])
plt.xlabel('Radius [mm]')
plt.ylabel('Vignetting model')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.title('Vignetting model from sky flats')
plt.savefig(outpath+'vignetting_model_combined.png', bbox_inches='tight')

#save in a pickle file
file = open(outpath+'vignetting_splines.pkl', 'wb')
pickle.dump(bspl_g, file)
pickle.dump(bspl_r, file)
pickle.dump(bspl_i, file)
pickle.dump(bspl_z, file)
file.close()

#test to read back the values
file = open(outpath+'vignetting_splines.pkl', 'rb')
bspl_g=pickle.load(file)
bspl_r=pickle.load(file)
bspl_i=pickle.load(file)
bspl_z=pickle.load(file)
print(bspl_g(350))
print(bspl_r(350))
print(bspl_i(350))
print(bspl_z(350))

sys.exit()
