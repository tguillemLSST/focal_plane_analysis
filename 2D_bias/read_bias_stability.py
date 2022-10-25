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

#outpath = "/sps/lsst/users/tguillem/web/debug/test1/2D/"
outpath = "/sps/lsst/users/tguillem/web/bias_stability/v4_1D/"

if os.path.exists(outpath):
    shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

###t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line','var_column', 'mean_column', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D'), meta={'name': 'Variances'})

#input file selection
inpath_base = '/sps/lsst/users/tguillem/web/batch/dark/ccd_assembly_fix_1/v4_1D/13161/after_master_bias/'

#exposures = ['000', '001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019']
exposures = ['000157','000158','000159','000160','000161','000162','000163','000164','000165','000166','000167','000168','000169','000170','000171','000172','000173','000174','000175','000176']
#exposures = ['000157','000158','000159','000160','000161']

rafts = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43', 'R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
#rafts = ['R01']
#ccds=['S00']

colors = ['black', 'red', 'blue', 'purple', 'darkorange', 'darkgreen', 'saddlebrown', 'navy', 'hotpink', 'dimgrey', 'gold', 'cyan', 'lightgreen', 'indigo', 'coral', 'fuchsia', 'tan', 'mediumaquamarine', 'darkkhaki', 'crimson'] 

#loop over all amplifiers
for i in range(len(rafts)):
    for j in range(len(ccds)):
        outpath_final = outpath + str(rafts[i])+'/'+str(ccds[j])+ '/'
        if os.path.exists(outpath_final):
            shutil.rmtree(outpath_final)
        os.makedirs(outpath_final)
        variances_exposures = []
        for k in range(len(exposures)):
            variances = []
            inpath = inpath_base + exposures[k] + '/' + rafts[i] + '/' + ccds[j] + '/Variances_2D_corr.fits'
            variances = Table.read(inpath)
            variances_exposures.append(variances)
        for l in range(16) :
                plt.figure()
                for m in range(len(exposures)):
                    line=variances_exposures[m]['mean_column'][l]
                    n_bins=len(line)
                    bin_x=np.empty(n_bins)
                    xedges=np.empty(n_bins+1)
                    for ix in range(n_bins):
                        bin_x[ix]=ix+0.5
                        xedges[ix]=ix+1
                        xedges[n_bins]=n_bins+1
                    plt.hist(bin_x, weights=line, bins=xedges, histtype='step', color=colors[m], label=exposures[m])
                plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
                plt.ylim([-5, +5])
                plt.xlabel('Column index')
                plt.ylabel('ADU')
                plt.title(str(rafts[i])+' '+str(ccds[j])+' amp '+str(l+1))
                plt.savefig(outpath_final+'mean_column_exposures_amp_'+str(l+1)+'.png')
                                    
sys.exit()
