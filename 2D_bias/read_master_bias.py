#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.io.fits as pyfits
from astropy.io import ascii
import sys
import os
import shutil
import pickle
import gc

#outpath = "/sps/lsst/users/tguillem/web/by_hand/plots/master_bias/2D/13162/"
outpath = "/sps/lsst/users/tguillem/web/debug/mess/"
if os.path.exists(outpath):
    shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

###t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line','var_column', 'mean_column', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D'), meta={'name': 'Variances'})

inpath_base = '/sps/lsst/users/tguillem/web/by_hand/master_bias/fix1/2D/13162/'

rafts_itl = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43']
rafts_itl = ['R01']
rafts_e2v = ['R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
rafts_e2v = ['R11']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
ccds=['S00']

#science image definitions
first_line_itl=1
first_p_over_itl=2000
first_col_itl=3
first_s_over_itl=4000
first_line_e2v=1
first_p_over_e2v=2000
first_col_e2v=10
first_s_over_e2v=4000

#read master bias images
#itl
master_bias_itl=[]
master_bias_corr_2D_itl=[]
for j in range(len(rafts_itl)):
    for k in range(len(ccds)):
        inpath = inpath_base + '/' + rafts_itl[j] + '/' + ccds[k] + '/masterbias.FITS'
        fits=pyfits.open(inpath)
        for i in range(16) :
            imarr = fits[i].data
            master_bias_itl.append(imarr)
        fits.close()
        inpath = inpath_base + '/' + rafts_itl[j] + '/' + ccds[k] + '/masterbias_rms.FITS'
        fits=pyfits.open(inpath)
        for i in range(16) :
            imarr = fits[i].data
            master_bias_corr_2D_itl.append(imarr)
        fits.close()
#e2v
master_bias_e2v=[]
master_bias_corr_2D_e2v=[]
for j in range(len(rafts_e2v)):
    for k in range(len(ccds)):
        inpath = inpath_base + '/' + rafts_e2v[j] + '/' + ccds[k] + '/masterbias.FITS'
        fits=pyfits.open(inpath)
        for i in range(16) :
            imarr = fits[i].data
            master_bias_e2v.append(imarr)
        fits.close()
        inpath = inpath_base + '/' + rafts_e2v[j] + '/' + ccds[k] + '/masterbias_rms.FITS'
        fits=pyfits.open(inpath)
        for i in range(16) :
            imarr = fits[i].data
            master_bias_corr_2D_e2v.append(imarr)
        fits.close()

#itl
means_mb_itl = []
variances_mb_itl = []
for j in range(len(master_bias_itl)):
    imarr = master_bias_itl[j]
    imarr_science = imarr[first_line_itl:first_p_over_itl,first_col_itl:first_s_over_itl]
    mean_mb_itl = np.mean(imarr_science)
    variance_mb_itl = np.sqrt(np.var(imarr_science))
    means_mb_itl.append(mean_mb_itl)
    variances_mb_itl.append(variance_mb_itl)
    del imarr
    del imarr_science
gc.collect()
means_mb_corr_2D_itl = []
variances_mb_corr_2D_itl = []
for j in range(len(master_bias_corr_2D_itl)):
    imarr = master_bias_corr_2D_itl[j]
    imarr_science = imarr[first_line_itl:first_p_over_itl,first_col_itl:first_s_over_itl]
    mean_mb_itl = np.mean(imarr_science)
    variance_mb_itl = np.sqrt(np.var(imarr_science))
    means_mb_corr_2D_itl.append(mean_mb_itl)
    variances_mb_corr_2D_itl.append(variance_mb_itl)
    del imarr
    del imarr_science
gc.collect()

#e2v
means_mb_e2v = []
variances_mb_e2v = []
for j in range(len(master_bias_e2v)):
    imarr = master_bias_e2v[j]
    imarr_science = imarr[first_line_e2v:first_p_over_e2v,first_col_e2v:first_s_over_e2v]
    mean_mb_e2v = np.mean(imarr_science)
    variance_mb_e2v = np.sqrt(np.var(imarr_science))
    means_mb_e2v.append(mean_mb_e2v)
    variances_mb_e2v.append(variance_mb_e2v)
    del imarr
    del imarr_science
gc.collect()
means_mb_corr_2D_e2v = []
variances_mb_corr_2D_e2v = []
for j in range(len(master_bias_corr_2D_e2v)):
    imarr = master_bias_corr_2D_e2v[j]
    imarr_science = imarr[first_line_e2v:first_p_over_e2v,first_col_e2v:first_s_over_e2v]
    mean_mb_e2v = np.mean(imarr_science)
    variance_mb_e2v = np.sqrt(np.var(imarr_science))
    means_mb_corr_2D_e2v.append(mean_mb_e2v)
    variances_mb_corr_2D_e2v.append(variance_mb_e2v)
    del imarr
    del imarr_science
gc.collect()

#mean_mb plot
bin_range = [-0.5,0.5]
nbins = 100
plt.figure()
plt.hist(means_mb_corr_2D_e2v, range=bin_range, bins=nbins, label='e2v', histtype='step', color = 'blue')
plt.hist(means_mb_corr_2D_itl, range=bin_range, bins=nbins, label='itl', histtype='step', color = 'red')
#plt.ylim([0.8, 1.1])
plt.xlabel('mean_mb')
plt.ylabel('n_amp')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.savefig(outpath+'mean_mb.png', bbox_inches='tight')

#var_mb plot
bin_range = [0,3]
nbins = 100
plt.figure()
plt.hist(variances_mb_corr_2D_e2v, range=bin_range, bins=nbins, label='e2v', histtype='step', color = 'blue')
plt.hist(variances_mb_corr_2D_itl, range=bin_range, bins=nbins, label='itl', histtype='step', color = 'red')
#plt.ylim([0.8, 1.1])
plt.xlabel('var_mb')
plt.ylabel('n_amp')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.savefig(outpath+'var_mb.png', bbox_inches='tight')
