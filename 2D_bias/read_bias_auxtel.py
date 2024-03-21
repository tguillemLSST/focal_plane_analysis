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

run = '13161'

outpath = "/sps/lsst/users/tguillem/web/run6/reference/run_comparison/"
#outpath = "/sps/lsst/users/tguillem/web/stack/waves/run_13161/stability/"
#outpath = '/sps/lsst/users/tguillem/web/bias_stability_fix/' + run + '/v4_1D_no_mb/'

#if os.path.exists(outpath):
#    shutil.rmtree(outpath)
#os.makedirs(outpath)
#print('outpath = ' + outpath)

###t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line','var_column', 'mean_column', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D'), meta={'name': 'Variances'})

#input file selection
#inpath_base = '/sps/lsst/users/tguillem/web/batch/dark/ccd_assembly_fix_1/v4_1D/13161/'#after_master_bias/'
#inpath_base = '/sps/lsst/users/tguillem/web/batch/dark/ccd_assembly_fix_1/v4_2D/' + run + '/'
#inpath_base = '/sps/lsst/users/tguillem/web/debug/auxtel/160523/'
inpath_base = '/sps/lsst/users/tguillem/web/run6/reference/run_13372/bias/1D/exposures/'

#configuration selection
#exposures = ['biasGen.20221107b','biasGen.20230428a']
exposures = ['exposure_3023061800021','exposure_3023061800021']
#label_exposures = ['DM-36719','DM-38946']
label_exposures = ['3023061800021','3023061800021']
print('Configuration arguments: ', str(sys.argv))
run = str(sys.argv[1])
rafts = [str(sys.argv[2])]
ccds = [str(sys.argv[3])]

colors = ['black', 'red', 'blue', 'purple', 'darkorange', 'darkgreen', 'saddlebrown', 'navy', 'hotpink', 'dimgrey', 'gold', 'cyan', 'lightgreen', 'indigo', 'coral', 'fuchsia', 'tan', 'mediumaquamarine', 'darkkhaki', 'crimson'] 

amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']

#loop over all amplifiers
for i in range(len(rafts)):
    for j in range(len(ccds)):
        outpath_final = outpath + str(rafts[i])+'/'+str(ccds[j])+ '/'
        #if os.path.exists(outpath_final):
        #    shutil.rmtree(outpath_final)
        os.makedirs(outpath_final,exist_ok=True)
        shutil.copy2('index_files/index_stability_4_plots.html',outpath_final+'index.html')
        variances_exposures = []
        for k in range(len(exposures)):
            variances = []
            #by hand
            #inpath = inpath_base + exposures[k] + '/' + rafts[i] + '/' + ccds[j] + '/Variances_2D_corr.fits'
            #DM
            inpath = inpath_base + exposures[k] + '/results/Variances_' + rafts[i] + '_' + ccds[j] + '.fits' 
            variances = Table.read(inpath)
            variances_exposures.append(variances)
        for l in range(16) :
                plt.figure()
                for m in range(len(exposures)):
                    line=variances_exposures[m]['mean_column'][l]
                    #print(line)
                    n_bins=len(line)
                    bin_x=np.empty(n_bins)
                    xedges=np.empty(n_bins+1)
                    for ix in range(n_bins):
                        bin_x[ix]=ix+0.5
                        xedges[ix]=ix+1
                        xedges[n_bins]=n_bins+1
                    plt.hist(bin_x, weights=line, bins=xedges, histtype='step', color=colors[m], label=label_exposures[m])
                plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
                #plt.ylim([23000, 25000])
                #line=variances_exposures[0]['mean_column'][l]
                #mean=np.mean(line)
                plt.ylim([-20,40])
                plt.xlabel('Column index')
                plt.ylabel('ADU')
                plt.title('Master bias ' + str(rafts[i])+' '+str(ccds[j])+ ' ' + str(amps[l]))
                plt.legend(loc='upper left')
                plt.savefig(outpath_final+'mean_column_exposures_amp_'+str(l+1)+'.png')
                plt.close()

sys.exit()
