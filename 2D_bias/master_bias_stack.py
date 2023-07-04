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

outpath = "/sps/lsst/users/tguillem/web/stack/run_6/comparison/test_1/"

if os.path.exists(outpath):
    shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

###t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, mean_corner, var_corner, variance_q95, variance_mean, variance_rms, variance_skewness], \
#names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column', 'mean_corner', 'var_corner', 'variance_q95', 'variance_mean', 'variance_rms', 'variance_skewness'), meta={'name': 'Variances'})

inpath_base = '/sps/lsst/users/tguillem/web/stack/run_6/reference/'
inpath_base_1 = inpath_base + '13391' + '/results/'
inpath_base_2 = inpath_base + '13392' + '/results/'

rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']

variances_1 = []
variances_2 = []
for j in range(len(rafts)):
    for k in range(len(ccds)):
        inpath_1 = inpath_base_1 + '/Variances_'+ rafts[j] + '_' + ccds[k] + '.fits'
        variance_1 = Table.read(inpath_1)
        variances_1.append(variance_1)
        inpath_2 = inpath_base_2 + '/Variances_'+ rafts[j] + '_' + ccds[k] + '.fits'
        variance_2 = Table.read(inpath_2)
        variances_2.append(variance_2)

#mean_col
amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
for i in range(len(rafts)):
    for j in range(len(ccds)):
        #outpath_final = outpath+rafts[i]+'/'+ccds[j]        
        outpath_final = outpath + '/Waves/'
        os.makedirs(outpath_final,exist_ok=True)
        inpath_1 = inpath_base_1 + '/Variances_'+ rafts[j] + '_' + ccds[k] + '.fits'
        variance_1 = Table.read(inpath_1)
        inpath_2 = inpath_base_2 + '/Variances_'+ rafts[j] + '_' + ccds[k] + '.fits'
        variance_2 = Table.read(inpath_2)
            
        for l in range(16) :
            #diff = variance['mean_corner'][l]-variance['mean_total'][l]
            #diff = round(diff,2)
            #if(diff>2 and diff<20):
            #    print(rafts[i]+' '+ccds[j]+' '+amps[l] + ' : ' + str(diff))
            #continue    
            #print(amps[l])
            plt.figure()
            line_1=variance_1['mean_column'][l]
            line_2=variance_2['mean_column'][l]
            n_bins=len(line_1)
            bin_x=np.empty(n_bins)
            xedges=np.empty(n_bins+1)
            for ix in range(n_bins):
                bin_x[ix]=ix+0.5
                xedges[ix]=ix+1
                xedges[n_bins]=n_bins+1
            plt.hist(bin_x, weights=line_1, bins=xedges, histtype='step', color = 'black', label='13391')
            plt.hist(bin_x, weights=line_2, bins=xedges, histtype='step', color = 'red', label='13392')
            plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
            plt.ylim([-10, +10])
            plt.xlabel('Column index')
            plt.ylabel('ADU')
            plt.title(rafts[i]+' '+ccds[j]+' '+amps[l])
            plt.legend()
            #plt.savefig(outpath_final+'/mean_column_amp_'+amps[l]+'.png')
            plt.savefig(outpath_final+'/mean_column_'+rafts[i]+'_'+ccds[j]+'_'+amps[l]+'.png')
            plt.close()
            #print n_total
            #n_total = 0
            #for ix in range(len(line)):
            #    if(abs(line[ix])>1):
            #        n_total +=1
            #if(n_total>0):
            #    print(rafts[i]+' '+ccds[j]+' '+amps[l] + ' : ' + str(n_total))        
            #extra requirement to store all unstable CCDs
            #if(n_total<10):
            #    plt.close()
            #    continue
            #print(rafts[i]+' '+ccds[j]+' '+amps[l] + ' : ' + str(n_total))
            #outpath_final_waves = outpath + 'Waves/'
            #plt.savefig(outpath_final_waves+'mean_column_'+rafts[i]+'_'+ccds[j]+'_'+amps[l])
            #plt.close()

sys.exit()
