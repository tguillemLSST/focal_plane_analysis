#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: build the gradient+vignetting model
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

outpath='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/vignetting_illumination_combined/model/'
os.makedirs(outpath,exist_ok=True)

#read the 5 inputs 
inpath='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/vignetting_illumination_combined/'
variance_TR = Table.read(inpath+'TR/e2v/R34/Variance_S22.fits')
variance_BR = Table.read(inpath+'BR/e2v/R14/Variance_S22.fits')
variance_TL = Table.read(inpath+'TL/e2v/R30/Variance_S22.fits')
variance_BL_e2v = Table.read(inpath+'BL/e2v/R12/Variance_S22.fits')
variance_BL_itl = Table.read(inpath+'BL/itl/R02/Variance_S22.fits')
variance_BL = Table.read(inpath+'BL/e2v/R12/Variance_S22.fits')

#check points
#for i in range(len(variance_TR['bin_centers'])):
    #print(str(i) + ' | ' + str(variance_TR['bin_centers'][i]) + ' | ' + str(variance_TR['bin_means'][i]))
    #print(str(i) + ' | ' + str(variance_BR['bin_centers'][i]) + ' | ' + str(variance_BR['bin_means'][i]))
    #print(str(i) + ' | ' + str(variance_TL['bin_centers'][i]) + ' | ' + str(variance_TL['bin_means'][i]))
    #print(str(variance_BL_e2v['bin_centers'][i]) + ' | ' + str(variance_BL_e2v['bin_means'][i]))
    #print(str(i) + ' | ' + str(variance_BL_itl['bin_centers'][i]) + ' | ' + str(variance_BL_itl['bin_means'][i]))

#HACK: correct variance_TL for missing values
for i in range(395,405):
    variance_TL['bin_means'][i]=variance_TL['bin_means'][394]

#continuity between e2v and ITL for BL
norm_BL_e2v = np.mean(variance_BL_e2v['bin_means'][490:509])
norm_BL_itl = np.mean(variance_BL_itl['bin_means'][490:509])
scaling_factor = norm_BL_e2v/norm_BL_itl
print(norm_BL_e2v)
print(norm_BL_itl)
print(scaling_factor)
for i in range(len(variance_BL['bin_centers'])):
    if(i>500):
        variance_BL['bin_means'][i]=variance_BL_itl['bin_means'][i]*scaling_factor

#compare corrections        
plt.figure()
plt.scatter(variance_TR['bin_centers'][0:729], variance_TR['bin_means'][0:729], marker='.', color = 'black', s=40, label='TR')
plt.scatter(variance_BR['bin_centers'][0:729], variance_BR['bin_means'][0:729], marker='.', color = 'blue', s=40, label='BR')
plt.scatter(variance_TL['bin_centers'][0:729], variance_TL['bin_means'][0:729], marker='.', color = 'red', s=40, label='TL')
plt.scatter(variance_BL['bin_centers'][0:729], variance_BL['bin_means'][0:729], marker='.', color = 'brown', s=40, label='BL')
#plt.scatter(variance_BL_e2v['bin_centers'][:], variance_BL_e2v['bin_means'][:], marker='.', color = 'brown', s=40, label='BL_e2v')
#plt.scatter(variance_BL_itl['bin_centers'][:], variance_BL_itl['bin_means'][:], marker='.', color = 'purple', s=40, label='BL_ITL')
plt.xlim([0,450])
#plt.xlim([250,350])
plt.ylim([0.,1.2])
plt.xlabel('Radius [mm]')
plt.ylabel('Correction factor')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.title('r band : gradient + vignetting correction')
plt.savefig(outpath+'vignetting_model.png', bbox_inches='tight')

#splines
from scipy.interpolate import make_interp_spline
from scipy.interpolate import splrep, splev
from scipy import interpolate

#first test using all bins
#bspl_TR = make_interp_spline(variance_TR['bin_centers'][0:729], variance_TR['bin_means'][0:729], k=3)
#bspl_BR = make_interp_spline(variance_BR['bin_centers'][0:729], variance_BR['bin_means'][0:729], k=3)
#bspl_TL = make_interp_spline(variance_TL['bin_centers'][0:729], variance_TL['bin_means'][0:729], k=3)
#bspl_BL = make_interp_spline(variance_BL['bin_centers'][0:729], variance_BL['bin_means'][0:729], k=3)

xx = np.linspace(0, 364.75, 1000)
#define nodes by hand
nodes=[10,50,100,150,200,250,300,350,400,450,500,520,540,560,580,600,620,640,660,680,700,720,729]
n_nodes =len(nodes)
#TR
x_TR = np.zeros(n_nodes)
y_TR = np.zeros(n_nodes)
for i_node in range(n_nodes):
    x_TR[i_node]=variance_TR['bin_centers'][nodes[i_node]]
    y_TR[i_node]=variance_TR['bin_means'][nodes[i_node]]
#BR
x_BR = np.zeros(n_nodes)
y_BR = np.zeros(n_nodes)
for i_node in range(n_nodes):
    x_BR[i_node]=variance_BR['bin_centers'][nodes[i_node]]
    y_BR[i_node]=variance_BR['bin_means'][nodes[i_node]]
#BL
x_BL = np.zeros(n_nodes)
y_BL = np.zeros(n_nodes)
for i_node in range(n_nodes):
    x_BL[i_node]=variance_BL['bin_centers'][nodes[i_node]]
    y_BL[i_node]=variance_BL['bin_means'][nodes[i_node]]
#TL
x_TL = np.zeros(n_nodes)
y_TL = np.zeros(n_nodes)
for i_node in range(n_nodes):
    x_TL[i_node]=variance_TL['bin_centers'][nodes[i_node]]
    y_TL[i_node]=variance_TL['bin_means'][nodes[i_node]]

#create splines  
pclk_TR = interpolate.splrep(x_TR, y_TR)
bspl_TR=interpolate.splev(xx, pclk_TR)
pclk_BR = interpolate.splrep(x_BR, y_BR)
bspl_BR=interpolate.splev(xx, pclk_BR)
pclk_BL = interpolate.splrep(x_BL, y_BL)
bspl_BL=interpolate.splev(xx, pclk_BL)
pclk_TL = interpolate.splrep(x_BR, y_TL)
bspl_TL=interpolate.splev(xx, pclk_TL)
    
def bspl_TR(x):
    return interpolate.splev(x, pclk_TR)    
def bspl_BR(x):
    return interpolate.splev(x, pclk_BR)
def bspl_BL(x):
    return interpolate.splev(x, pclk_BL)
def bspl_TL(x):
    return interpolate.splev(x, pclk_TL)

plt.figure()
#plt.scatter(variance_S22['bin_centers'][:], variance_S22['bin_means'][:], marker='.', color = 'black', s=40, label='Data from 4 rafts')
plt.xlim([0,450])
#plt.xlim([250,350])
#plt.xlim([300,370])
#plt.plot(xx, bspl_TR,color = 'black', label='TR model')
#plt.plot(xx, bspl_BR,color = 'blue', label='BR model')
#plt.plot(xx, bspl_BL,color = 'red', label='BL model')
plt.plot(xx, bspl_TR(xx),color = 'black', label='TL model')
plt.plot(xx, bspl_BR(xx),color = 'blue', label='BR model')
plt.plot(xx, bspl_TL(xx),color = 'red', label='TL model')
plt.plot(xx, bspl_BL(xx),color = 'brown', label='BL model')
plt.ylim([0.,1.2])
plt.xlabel('Radius [mm]')
plt.ylabel('Correction factor')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.title('r band : gradient + vignetting correction')
plt.savefig(outpath+'vignetting_model_combined.png', bbox_inches='tight')

#save splines
file = open(outpath+'vignetting_splines.pkl', 'wb')
pickle.dump(pclk_TR, file)
pickle.dump(pclk_BR, file)
pickle.dump(pclk_TL, file)
pickle.dump(pclk_BL, file)
file.close()

sys.exit()
