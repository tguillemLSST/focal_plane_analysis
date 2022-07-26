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

outpath = "/sps/lsst/users/tguillem/web/debug/"
#outpath = "/sps/lsst/users/tguillem/web/variances/2D_corr/itl_e2v/"


if os.path.exists(outpath):
    shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

###t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line','var_column', 'mean_column', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D'), meta={'name': 'Variances'})

inpath_base = '/sps/lsst/users/tguillem/web/batch/variances_2D_corr/13162/bias_bias_'
#inpath_base = '/sps/lsst/users/tguillem/web/batch/variances_fix3/13162/bias_bias_'
#exposures = ['005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019']
exposures = ['015']
#rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
#rafts=['R14']
rafts_itl = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43']
rafts_e2v = ['R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
#rafts_itl = []
#rafts_e2v = ['R14']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
#ccds=['S22']

variances_itl = []
for i in range(len(exposures)):
    for j in range(len(rafts_itl)):
        for k in range(len(ccds)):
            inpath = inpath_base + exposures[i] + '/' + rafts_itl[j] + '/' + ccds[k] + '/Variances_2D_corr.fits'
            variance_itl = Table.read(inpath)
            variances_itl.append(variance_itl)
                    
variances_e2v = []
for i in range(len(exposures)):
    for j in range(len(rafts_e2v)):
        for k in range(len(ccds)):
            inpath = inpath_base + exposures[i] + '/' + rafts_e2v[j] + '/' + ccds[k] + '/Variances_2D_corr.fits'
            variance_e2v = Table.read(inpath)
            variances_e2v.append(variance_e2v)            

#1-image analysis
#e2v
var_one_exposure_e2v = []
var_one_exposure_corr_2D_e2v = []
ratio_var_one_exposure_corr_2D_e2v = []
var_one_exposure_row_e2v = []
var_one_exposure_row_corr_2D_e2v = []
ratio_var_one_exposure_row_corr_2D_e2v = []
var_one_exposure_col_e2v = []
var_one_exposure_col_corr_2D_e2v = []
ratio_var_one_exposure_col_corr_2D_e2v = []
for j in range(len(variances_e2v)):
    var_one_exposure = np.zeros(16)
    var_one_exposure_corr_2D = np.zeros(16)
    ratio_var_one_exposure_corr_2D = np.zeros(16)
    var_one_exposure_row = np.zeros(16)
    var_one_exposure_row_corr_2D = np.zeros(16)
    ratio_var_one_exposure_row_corr_2D = np.zeros(16)
    var_one_exposure_col = np.zeros(16)
    var_one_exposure_col_corr_2D = np.zeros(16)
    ratio_var_one_exposure_col_corr_2D = np.zeros(16)
    for i in range(16) :
        var_one_exposure[i]=np.sqrt(variances_e2v[j]['var_total'][i])
        var_one_exposure_corr_2D[i]=np.sqrt(variances_e2v[j]['var_total_corr_2D'][i])
        ratio_var_one_exposure_corr_2D[i]=var_one_exposure_corr_2D[i]/var_one_exposure[i]
        #row
        var_one_exposure_row[i]=np.sqrt(np.var(variances_e2v[j]['mean_line'][i]))
        var_one_exposure_row_corr_2D[i]=np.sqrt(np.var(variances_e2v[j]['mean_line_corr_2D'][i]))
        ratio_var_one_exposure_row_corr_2D[i]=var_one_exposure_row_corr_2D[i]/var_one_exposure_row[i]
        #col
        var_one_exposure_col[i]=np.sqrt(np.var(variances_e2v[j]['mean_column'][i]))
        var_one_exposure_col_corr_2D[i]=np.sqrt(np.var(variances_e2v[j]['mean_column_corr_2D'][i]))
        ratio_var_one_exposure_col_corr_2D[i]=var_one_exposure_col_corr_2D[i]/var_one_exposure_col[i]
    var_one_exposure_e2v.append(var_one_exposure)
    var_one_exposure_corr_2D_e2v.append(var_one_exposure_corr_2D)
    ratio_var_one_exposure_corr_2D_e2v.append(ratio_var_one_exposure_corr_2D)
    var_one_exposure_row_e2v.append(var_one_exposure_row)
    var_one_exposure_row_corr_2D_e2v.append(var_one_exposure_row_corr_2D)
    ratio_var_one_exposure_row_corr_2D_e2v.append(ratio_var_one_exposure_row_corr_2D)
    var_one_exposure_col_e2v.append(var_one_exposure_col)
    var_one_exposure_col_corr_2D_e2v.append(var_one_exposure_col_corr_2D)
    ratio_var_one_exposure_col_corr_2D_e2v.append(ratio_var_one_exposure_col_corr_2D)

#duplicate
#itl
var_one_exposure_itl = []
var_one_exposure_corr_2D_itl = []
ratio_var_one_exposure_corr_2D_itl = []
var_one_exposure_row_itl = []
var_one_exposure_row_corr_2D_itl = []
ratio_var_one_exposure_row_corr_2D_itl = []
var_one_exposure_col_itl = []
var_one_exposure_col_corr_2D_itl = []
ratio_var_one_exposure_col_corr_2D_itl = []
for j in range(len(variances_itl)):
    var_one_exposure = np.zeros(16)
    var_one_exposure_corr_2D = np.zeros(16)
    ratio_var_one_exposure_corr_2D = np.zeros(16)
    var_one_exposure_row = np.zeros(16)
    var_one_exposure_row_corr_2D = np.zeros(16)
    ratio_var_one_exposure_row_corr_2D = np.zeros(16)
    var_one_exposure_col = np.zeros(16)
    var_one_exposure_col_corr_2D = np.zeros(16)
    ratio_var_one_exposure_col_corr_2D = np.zeros(16)
    for i in range(16) :
        var_one_exposure[i]=np.sqrt(variances_itl[j]['var_total'][i])
        var_one_exposure_corr_2D[i]=np.sqrt(variances_itl[j]['var_total_corr_2D'][i])
        ratio_var_one_exposure_corr_2D[i]=var_one_exposure_corr_2D[i]/var_one_exposure[i]
        #row
        var_one_exposure_row[i]=np.sqrt(np.var(variances_itl[j]['mean_line'][i]))
        var_one_exposure_row_corr_2D[i]=np.sqrt(np.var(variances_itl[j]['mean_line_corr_2D'][i]))
        ratio_var_one_exposure_row_corr_2D[i]=var_one_exposure_row_corr_2D[i]/var_one_exposure_row[i]
        #col
        var_one_exposure_col[i]=np.sqrt(np.var(variances_itl[j]['mean_column'][i]))
        var_one_exposure_col_corr_2D[i]=np.sqrt(np.var(variances_itl[j]['mean_column_corr_2D'][i]))
        ratio_var_one_exposure_col_corr_2D[i]=var_one_exposure_col_corr_2D[i]/var_one_exposure_col[i]
    var_one_exposure_itl.append(var_one_exposure)
    var_one_exposure_corr_2D_itl.append(var_one_exposure_corr_2D)
    ratio_var_one_exposure_corr_2D_itl.append(ratio_var_one_exposure_corr_2D)
    var_one_exposure_row_itl.append(var_one_exposure_row)
    var_one_exposure_row_corr_2D_itl.append(var_one_exposure_row_corr_2D)
    ratio_var_one_exposure_row_corr_2D_itl.append(ratio_var_one_exposure_row_corr_2D)
    var_one_exposure_col_itl.append(var_one_exposure_col)
    var_one_exposure_col_corr_2D_itl.append(var_one_exposure_col_corr_2D)
    ratio_var_one_exposure_col_corr_2D_itl.append(ratio_var_one_exposure_col_corr_2D)
    
#var_total plot
plt.figure()
plt.scatter(var_one_exposure_e2v,ratio_var_one_exposure_corr_2D_e2v, marker='.',color = 'blue', s=10, alpha=0.3, label='e2v')
plt.scatter(var_one_exposure_itl,ratio_var_one_exposure_corr_2D_itl, marker='.',color = 'red', s=10, alpha=0.3, label='itl')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([3, 7])
plt.ylim([0.8, 1.1])
plt.xlabel('var_total')
plt.ylabel('ratio_corr_2D')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.savefig(outpath+'ratio_var_total_corr_2D.png', bbox_inches='tight')
#plt.savefig('ratio_var_total_corr_2D.png', bbox_inches='tight')

#var_row plot
plt.figure()
plt.scatter(var_one_exposure_row_e2v,ratio_var_one_exposure_row_corr_2D_e2v, marker='.',color = 'blue', s=10, alpha=0.3, label='e2v')
plt.scatter(var_one_exposure_row_itl,ratio_var_one_exposure_row_corr_2D_itl, marker='.',color = 'red', s=10, alpha=0.3, label='itl')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0, 7])
plt.ylim([0, 3])
plt.xlabel('var_mean_row')
plt.ylabel('ratio_corr_2D')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.savefig(outpath+'ratio_var_row_corr_2D.png', bbox_inches='tight')
#plt.savefig('ratio_var_total_corr_2D.png', bbox_inches='tight')

#var_col plot
plt.figure()
plt.scatter(var_one_exposure_col_e2v,ratio_var_one_exposure_col_corr_2D_e2v, marker='.',color = 'blue', s=10, alpha=0.3, label='e2v')
plt.scatter(var_one_exposure_col_itl,ratio_var_one_exposure_col_corr_2D_itl, marker='.',color = 'red', s=10, alpha=0.3, label='itl')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0, 7])
plt.ylim([0, 3])
plt.xlabel('var_mean_col')
plt.ylabel('ratio_corr_2D')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
plt.savefig(outpath+'ratio_var_col_corr_2D.png', bbox_inches='tight')

sys.exit()

#compute variance over exposures
var_exposures = np.zeros(16)
var_exposures_corr_2D = np.zeros(16)
ratio_var_exposures_corr_2D = np.zeros(16)
for i in range(16) :
    var_exposures_tmp =  np.zeros(len(exposures))
    var_exposures_corr_2D_tmp =  np.zeros(len(exposures))
    for j in range(len(exposures)):
        var_exposures_tmp[j] = variances[j]['mean_total'][i]
        var_exposures_corr_2D_tmp[j] = variances[j]['mean_total_corr_2D'][i]
    var_exposures[i]=np.sqrt(np.var(var_exposures_tmp))
    var_exposures_corr_2D[i]=np.sqrt(np.var(var_exposures_corr_2D_tmp))
    ratio_var_exposures_corr_2D[i]=var_exposures_corr_2D[i]/var_exposures[i]
print(var_exposures)
print(var_exposures_corr_2D)
print(ratio_var_exposures_corr_2D)

#variance per line (debugging)
#line
#x=np.sqrt(np.var(variances_e2v[j]['mean_line'][i]))
#y=np.sqrt(np.var(variances_e2v[j]['mean_line_corr_2D'][i]))
#print(x)
#print(y)
