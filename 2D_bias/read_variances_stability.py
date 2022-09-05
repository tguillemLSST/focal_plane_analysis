#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
import sys
import os
import shutil
import pickle

outpath = "/sps/lsst/users/tguillem/web/debug/test1/2D/"
#outpath = "/sps/lsst/users/tguillem/web/variances/2D_corr/itl_e2v/"


if os.path.exists(outpath):
    shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

###t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line','var_column', 'mean_column', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D'), meta={'name': 'Variances'})

#input file selection
#2D row-col correction
#inpath_base = '/sps/lsst/users/tguillem/web/batch/variances_2D_corr/13162/bias_bias_'
#inpath_base = '/sps/lsst/users/tguillem/web/batch/master_bias/v1/13162/bias_bias_'
inpath_base = '/sps/lsst/users/tguillem/web/batch/master_bias/v1/13162/after_master_bias/bias_bias_'

exposures = ['000', '001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019']
#exposures = ['015','016']
#rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
#rafts=['R14']
rafts_itl = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43']
rafts_e2v = ['R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
#rafts_itl = []
#rafts_e2v = ['R14']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
#ccds=['S22']

#define numbers
n_e2v = len(rafts_e2v)*len(ccds)
n_amp_e2v = n_e2v*16
n_itl = len(rafts_itl)*len(ccds)
n_amp_itl = n_itl*16
n_exposures = len(exposures)

variances_itl_exposures = []
for i in range(len(exposures)):
    variances_itl = []
    for j in range(len(rafts_itl)):
        for k in range(len(ccds)):
            inpath = inpath_base + exposures[i] + '/' + rafts_itl[j] + '/' + ccds[k] + '/Variances_2D_corr.fits'
            variance_itl = Table.read(inpath)
            variances_itl.append(variance_itl)
    variances_itl_exposures.append(variances_itl)

variances_e2v_exposures = []
for i in range(len(exposures)):
    variances_e2v = []
    for j in range(len(rafts_e2v)):
        for k in range(len(ccds)):
            inpath = inpath_base + exposures[i] + '/' + rafts_e2v[j] + '/' + ccds[k] + '/Variances_2D_corr.fits'
            variance_e2v = Table.read(inpath)
            variances_e2v.append(variance_e2v)
    variances_e2v_exposures.append(variances_e2v)
    
#variances_..._exposures[exposure][radt_cdd_id]['variable'][amp_id]
#print(variances_itl_exposures[1][3]['mean_total'][5])

#compute variance over exposures
#itl
#e2v
n_itl = len(rafts_itl)*len(ccds)
n_amp_itl = n_itl*16
var_exposures_itl = []
var_exposures_corr_2D_itl = []
ratio_var_exposures_corr_2D_itl = []

for i in range(n_itl):
    for j in range(16) :
        mean_exposures_tmp = np.zeros(n_exposures)
        for k in range(n_exposures):
            mean_exposures_tmp[k] = variances_itl_exposures[k][i]['mean_total_corr_2D'][j]
        var_exposures_tmp = np.sqrt(np.var(mean_exposures_tmp[:]))
        var_exposures_itl.append(var_exposures_tmp)
#print(var_exposures_itl)

#e2v
n_e2v = len(rafts_e2v)*len(ccds)
n_amp_e2v = n_e2v*16
var_exposures_e2v = []
var_exposures_corr_2D_e2v = []
ratio_var_exposures_corr_2D_e2v = []

for i in range(n_e2v):
    for j in range(16) :
        mean_exposures_tmp = np.zeros(n_exposures)
        for k in range(n_exposures):
            mean_exposures_tmp[k] = variances_e2v_exposures[k][i]['mean_total_corr_2D'][j]
        var_exposures_tmp = np.sqrt(np.var(mean_exposures_tmp[:]))
        var_exposures_e2v.append(var_exposures_tmp)
#print(var_exposures_e2v)

#var_exposures plot
bin_range = [0,0.2]
nbins = 80
plt.figure()
plt.hist(var_exposures_e2v, range=bin_range, bins=nbins, label='e2v', histtype='step', color = 'blue')
plt.hist(var_exposures_itl, range=bin_range, bins=nbins, label='itl', histtype='step', color = 'red')
#plt.ylim([0.8, 1.1])
plt.xlabel('var_exposures')
plt.ylabel('n_amp')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.savefig(outpath+'var_exposures.png', bbox_inches='tight')

sys.exit()
