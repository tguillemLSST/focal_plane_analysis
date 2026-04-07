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

outpath = '/sdf/home/t/tguillem/public_html/cbp/vignetting/efficiencies/'

#inpath_base = '/sdf/home/a/amouroux/public_html/cbp/lsstcam/telescope_response/20250518/cbp_qe_none/results/'
inpath_base = '/sdf/home/a/amouroux/public_html/cbp/lsstcam/telescope_response/20250516/cbp_none_all_2/results/'
#inpath_base = '/sdf/home/a/amouroux/public_html/cbp/lsstcam/telescope_response/20250518/cbp_qe_none/results/'

###jobs
print('Configuration arguments: ', str(sys.argv))
#run = str(sys.argv[1])
rafts = [str(sys.argv[1])]
ccds = [str(sys.argv[2])]

rafts=['R01','R02','R03','R10','R11','R12','R13','R14','R20','R21','R22','R23','R24','R30','R31','R32','R33','R34','R41','R42','R43']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']

amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
colors = ['black', 'red', 'blue', 'purple', 'darkorange', 'darkgreen', 'saddlebrown', 'navy', 'hotpink', 'dimgrey', 'gold', 'cyan', 'lightgreen', 'indigo', 'coral', 'fuchsia', 'tan', 'mediumaquamarine', 'darkkhaki', 'crimson'] 

detectors_diagonal=['R22_S11','R23_S20','R33_S01','R34_S10','R34_S22']

variances_all=[]
variances_all_ratio=[]
variances_names=[]

#choose a reference
variances_reference = Table.read(inpath_base+'R22_S11.fits')
ratio_reference=[]
for k in range(0,len(variances_reference)):
    ratio_reference.append(variances_reference['telescope_response'][k])

#loop over all amplifiers
for i in range(len(rafts)):
    for j in range(len(ccds)):
        detector_full_name = rafts[i] + '_' +ccds[j]
        if(detector_full_name not in detectors_diagonal):
            continue
        inpath = inpath_base + rafts[i] + '_' + ccds[j] + '.fits' 
        variances = Table.read(inpath)
        variances_all.append(variances)
        variances_names.append(detector_full_name)
        ratios=[]
        for k in range(0,len(variances)):
            print(variances['wavelength'][k])
            ratio=variances['telescope_response'][k]/ratio_reference[k]
            #print(variances_all)
            ratios.append(ratio)
        variances_all_ratio.append(ratios)
        
#plot
plt.figure()
for i in range(0,4):
    #plt.scatter(variances_all[i]['wavelength'][:], variances_all[i]['telescope_response'][:], marker='.', color = colors[i], s=30, label=variances_names[i])
    #plt.scatter(variances_all[i]['wavelength'][:], variances_all[i]['tabulated_response'][:], marker='*', color = colors[i], s=30, label=variances_names[i])
    plt.scatter(variances_all[i]['wavelength'][:], variances_all_ratio[i][:], marker='.', color = colors[i], s=30, label=variances_names[i])
plt.xlim([405,930])
#plt.ylim([0,0.8])
#plt.xlim([490,910])
plt.ylim([0.9,1.05])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Transmission')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.legend()
#plt.title(rafts[i] + '_' + ccds[j])
plt.savefig(outpath+'comparison.png', bbox_inches='tight')

sys.exit()
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
            plt.hist(bin_x, weights=line, bins=xedges, histtype='step', color=colors[m], label=exposures[m])
            plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
            #plt.ylim([23000, 25000])
        line=variances_exposures[1]['mean_column'][l]
        #mean=np.mean(line)
        #plt.ylim([mean-5,mean+5])
        max_line=np.max(line)
        min_line=np.min(line)
        plt.ylim([min_line-5,max_line+5])
        plt.ylim([-10,10])
        plt.xlabel('Column index')
        plt.ylabel('ADU')
        plt.title('Run ' + run + ' ' + str(rafts[i])+' '+str(ccds[j])+' amp '+ amps[l])
        #plt.legend()
        plt.savefig(outpath_final+'mean_column_exposures_amp_'+str(l+1)+'.png')
        plt.close()
        
        #mean stability
        for l in range(16) :
            mean_total_exp=[]
            for m in range(len(exposures)):
                mean_total_exp.append(variances_exposures[m]['mean_total'][l])
                #plot
                
            
        #fix to get the mean
        mean_fix=np.zeros(16)
        for l in range(16) :
            mean_exposures=np.zeros(20)
            for m in range(len(exposures)):
                mean_exposures[m]=variances_exposures[m]['mean_total'][l]
                mean_fix[l]=np.mean(mean_exposures)                
                
                for l in range(16) :
                    for m in range(len(exposures)):
                        mean_exposures_all_amp[l][m] =  variances_exposures[m]['mean_total'][l] - mean_fix[l] 
                        
                        #print(mean_exposures_all_amp)    
                        t_variance = Table([mean_exposures_all_amp],names=['mean_exposures_all_amp'],meta={'name': 'Variances'})
                        t_variance.write(outpath_final + 'bias_jumps.fits', overwrite=True)
                        
                        
sys.exit()
#print(mean_exposures_all_amp)         

#var_total plot
plt.figure()
index_amp = np.full(20,10.5)
plt.xlim([-5, 3030])
plt.ylim([-5,5])
for i in range(len(rafts)):
    for j in range(len(ccds)):
        z = i*9*16+j*16+0.5
        index_amp = np.full(20,z)
        for k in range(len(exposures)):
            plt.scatter(index_amp[k],mean_exposures_all_amp[0][k], marker='.',color = colors[k], s=10, alpha=0.3, label=exposures[m])
plt.xlabel('bias shift')
plt.ylabel('amplifier index')
plt.grid(which='major', axis='y', linestyle='-', linewidth='0.5', color='grey')
#plt.legend()
plt.savefig(outpath+'bias_jumps.png', bbox_inches='tight')
sys.exit()
