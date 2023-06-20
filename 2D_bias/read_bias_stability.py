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

run = '13372'

outpath = "/sps/lsst/users/tguillem/web/run6/reference/run_13372/bias/1D/stability_waves/"

###t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line','var_column', 'mean_column', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D'), meta={'name': 'Variances'})

#input file selection
#inpath_base = '/sps/lsst/users/tguillem/web/stack/waves/run_13161/'
inpath_base = '/sps/lsst/users/tguillem/web/run6/reference/run_13372/bias/1D/exposures/'

#exposure selection
#exposures = ['000', '001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019']
###run 13161
#bias
#exposures = ['000137','000138','000139','000140','000141','000142','000143','000144','000145','000146','000147','000148','000149','000150','000151','000152','000153','000154','000155','000156']
#exposures = ['3021121200145','3021121200137','3021121200141','3021121200150','3021121200151','3021121200156','3021121200152','3021121200148','3021121200139','3021121200138','30211212001\
#54','3021121200149','3021121200153','3021121200147','3021121200140','3021121200146','3021121200143','3021121200155','3021121200142','3021121200144']
#exposures = ['3021121200137','3021121200138','3021121200139','3021121200140','3021121200141','3021121200142','3021121200143','3021121200144','3021121200145','3021121200146','3021121200147','3021121200148','3021121200149','3021121200150','3021121200151','3021121200152','3021121200153','3021121200154','3021121200155','3021121200156']

#dark
#exposures = ['000157','000158','000159','000160','000161','000162','000163','000164','000165','000166','000167','000168','000169','000170','000171','000172','000173','000174','000175','000176']
#exposures = ['000157','000158','000159','000160','000161']
###run 13018
#exposures = ['000269','000268','000267','000266','000265','000264','000263','000262','000261','000260','000259','000258','000257','000256','000255','000254','000253','000252','000251']

###Run 6
exposures = [ '3023061800079','3023061800100','3023061800058','3023061800025','3023061800106','3023061800024','3023061800049','3023061800023','3023061800021','3023061800073','3023061800052','3023061800027','3023061800070','3023061800076','3023061800064','3023061800091','3023061800082','3023061800094','3023061800103','3023061800030']

#rafts = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43', 'R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
#ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
#rafts = ['R20']
#ccds=['S20']
###jobs
print('Configuration arguments: ', str(sys.argv))
run = str(sys.argv[1])
rafts = [str(sys.argv[2])]
ccds = [str(sys.argv[3])]

colors = ['black', 'red', 'blue', 'purple', 'darkorange', 'darkgreen', 'saddlebrown', 'navy', 'hotpink', 'dimgrey', 'gold', 'cyan', 'lightgreen', 'indigo', 'coral', 'fuchsia', 'tan', 'mediumaquamarine', 'darkkhaki', 'crimson'] 

mean_exposures_all_amp = np.zeros((16,20))

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
            inpath = inpath_base + 'exposure_'+exposures[k]+'/results/Variances_'+ rafts[i] + '_' + ccds[j] + '.fits' 
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
                    plt.hist(bin_x, weights=line, bins=xedges, histtype='step', color=colors[m], label=exposures[m])
                plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
                #plt.ylim([23000, 25000])
                line=variances_exposures[10]['mean_column'][l]
                mean=np.mean(line)
                #plt.ylim([mean-10,mean+30])
                plt.ylim([0,5])
                plt.xlabel('Column index')
                plt.ylabel('ADU')
                plt.title('Run ' + run + ' ' + str(rafts[i])+' '+str(ccds[j])+' amp '+str(l+1))
                #plt.legend()
                plt.savefig(outpath_final+'mean_column_exposures_amp_'+str(l+1)+'.png')
                plt.close()

        #extra plot
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
