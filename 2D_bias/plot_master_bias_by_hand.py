#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: plot master biases produced by hand
##########

# system imports
import time
import sys
from sys import exit
import glob
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table
import matplotlib
import os
import shutil
from astropy.visualization import (ImageNormalize,PercentileInterval,HistEqStretch)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc 
import lsst.daf.butler as dafButler
from conversion import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import skew

#science image definitions
first_line_itl=0
first_p_over_itl=2000
first_col_itl=3
first_s_over_itl=509
first_line_e2v=0
first_p_over_e2v=2000
first_col_e2v=10
first_s_over_e2v=512
################
def SingleImageIR(image,first_col,first_cover,first_line,first_lower,is_e2v):
    # Display an IR2 image , with amplifiers set at the right place ...there is a DM version which does this better...
    # but here you are in stand alone
    # the default associated to the image area (pre-overscan excluded) are for e2v IR2 files
    #
    col_size=first_cover-first_col
    line_size=first_lower-first_line
    #
    spf=np.zeros((line_size*2,col_size*8))
    if(is_e2v==True):
        for i in range(16) :
            if i<8 :
                xx=i*col_size-1
                for x in range(first_col,first_cover) :
                    for y in range(first_line,first_lower) :
                        spf[2*line_size-(y-first_line)-1,xx+col_size-(x-first_col)]=image[i][y,x]
            else :
                xx=(15-i)*col_size
                for y in range(first_line,first_lower) :
                    spf[y-first_line,xx:xx+col_size]=image[i][y,first_col:first_cover]
        return spf
    else:
        for i in range(16) :
            if i<8 :
                xx=i*col_size-1
                for x in range(first_col,first_cover) :
                    for y in range(first_line,first_lower) :
                        spf[2*line_size-(y-first_line)-1,xx+col_size-(x-first_col)]=image[i][y,x]
            else :
                xx=(15-i)*col_size
                for x in range(first_col,first_cover) :
                    for y in range(first_line,first_lower) :
                        spf[y-first_line,xx+col_size-(x-first_col)-1]=image[i][y,x]
        return spf
################
        
#runs=('13162')
#rafts=['R01','R02','R03','R10','R11','R12','R13','R14','R20','R21','R22','R23','R24','R30','R31','R32','R33','R34','R41','R42','R43']
rafts=['R02']
#ccds=['S00','S01','S02','S10','S11','S12','S20','S21','S22']
ccds=['S01']
#for a set of batch jobs
str_raft = str(sys.argv[1])
str_all_sensors = str(sys.argv[2])
rafts=[str_raft]
ccds=[str_all_sensors]

#files to access: /sps/lsst/users/tguillem/web/by_hand/master_bias/1D/13162/R41/S10/masterbias_rms.FITS
path_input = '/sps/lsst/users/tguillem/web/by_hand/master_bias/1D/13162/'
#path_input2 = '/sps/lsst/users/tguillem/web/by_hand/master_bias/fix1/2D/13162/'
#path_input = '/sps/lsst/users/tguillem/web/debug/waves/metrics/13162/'
path_output = '/sps/lsst/users/tguillem/web/by_hand/plots/master_bias/waves/skewness/1D/'
#path_output = '/sps/lsst/users/tguillem/web/by_hand/plots/master_bias_details/overscan/1D/13162/'
#path_output = '/sps/lsst/users/tguillem/web/debug/rms/'
#select raw or corr_2D files 
#image_file = 'masterbias.FITS'
#image_rms_file = 'masterbias_rms.FITS'
#corr_2D
image_file = 'masterbias_corr_2D.FITS'
image_rms_file = 'masterbias_corr_2D_rms.FITS'
###read files from Pierre
#path_input = '/sps/lsst/users/astier/slac/13144/'
#path_output = '/sps/lsst/users/tguillem/web/by_hand/plots/master_bias/from_pierre/13144/'
#read a second file in order to do a variance difference image
#path_input2 = '/sps/lsst/users/tguillem/web/by_hand/master_bias/2D/13162/'
path_input2 = '/sps/lsst/users/tguillem/web/by_hand/master_bias/fix1/1D/13162/'
#path_input2 = '/sps/lsst/users/tguillem/web/run6/metrics/test_130623/13162/'
image_file2 = 'masterbias_corr_2D.FITS'
image_rms_file2 = 'masterbias_corr_2D_rms.FITS'
#image_file2 = 'masterbias.FITS'
#image_rms_file2 = 'masterbias_rms.FITS'

os.makedirs(path_output,exist_ok=True)
os.makedirs(path_output+'variance/',exist_ok=True)
os.makedirs(path_output+'results/',exist_ok=True)

#loop over all CCDs
for i_raft in range(len(rafts)):
    #create raft directory
    outpath_final = path_output+rafts[i_raft]+'/'
    outpath_final2 = path_output+'variance/'+rafts[i_raft]+'/'
    outpath_final3 = path_output+'variance_difference/'+rafts[i_raft]+'/'
    outpath_final4 = path_output+'overscan/'+rafts[i_raft]+'/'
    #shutil.copy2('index_files/index_raft_plots.html',outpath_final+'/index.html')

    first_line=first_line_e2v
    first_p_over=first_p_over_e2v
    first_col=first_col_e2v
    first_s_over=first_s_over_e2v
    ccd_e2v=True
    if rafts[i_raft] in ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43']:
        first_line=first_line_itl
        first_p_over=first_p_over_itl
        first_col=first_col_itl
        first_s_over=first_s_over_itl
        ccd_e2v=False
    im_y_size=first_p_over-first_line
    im_x_size=first_s_over-first_col

    #loop over all CCDs
    for i_ccd in range(len(ccds)):
        outpath_final = outpath_final+ccds[i_ccd]+'/'
        outpath_final2 = outpath_final2 +ccds[i_ccd]+'/'
        outpath_final3 = outpath_final3 +ccds[i_ccd]+'/'
        outpath_final4 = outpath_final4 +ccds[i_ccd]+'/'
        os.makedirs(outpath_final,exist_ok=True)
        os.makedirs(outpath_final2,exist_ok=True)
        os.makedirs(outpath_final3,exist_ok=True)
        os.makedirs(outpath_final4,exist_ok=True)
        title = rafts[i_raft]+' '+ccds[i_ccd]
        print('Processing ' + title)
        #mean
        fits=pyfits.open(path_input2+'/'+rafts[i_raft]+'/'+ccds[i_ccd]+'/'+image_file2)
        image_tmp=[]
        for i in range(16) :
           imarr = fits[i].data
           mean = np.mean(imarr[first_line:first_p_over,first_col:first_s_over])
           imarr = imarr - mean
           image_tmp.append(imarr)
        image = SingleImageIR(image_tmp,first_line,first_s_over,first_col,first_p_over,ccd_e2v)    
        plt.figure()
        plt.imshow(image, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_final+'/image.png')
        plt.close()

        #rms
        fits=pyfits.open(path_input2+'/'+rafts[i_raft]+'/'+ccds[i_ccd]+'/'+image_rms_file2)
        image_tmp=[]
        for i in range(16) :
            imarr = fits[i].data
            image_tmp.append(imarr)
            #print('check')
            #print(imarr)
        image = SingleImageIR(image_tmp,first_line,first_s_over,first_col,first_p_over,ccd_e2v)    
        plt.figure()
        #rms
        plt.imshow(image, origin='lower', vmin=0, vmax=10, cmap = 'hot')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_final2+'/image.png')
        plt.close()
        #1D histogram of variance per pixel
        rms_flat = image.flat
        #for i in range(len(rms_flat)):
        #    print(rms_flat[i])
        plt.figure()
        plt.hist(rms_flat, range=[0,10], bins=1000, label='RMS/pixel', histtype='step', color = 'black')
        #plt.ylim([0.8, 1.1])
        plt.xlabel('RMS')
        plt.ylabel('pixels')
        plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
        plt.legend()
        plt.savefig(outpath_final2+'/RMS.png', bbox_inches='tight')
        
        #rms ratio
        fits=pyfits.open(path_input2+'/'+rafts[i_raft]+'/'+ccds[i_ccd]+'/'+image_rms_file2)
        fits2=pyfits.open(path_input2+'/'+rafts[i_raft]+'/'+ccds[i_ccd]+'/'+image_rms_file2)
        image_tmp=[]
        for i in range(16) :
            imarr1 = fits[i].data**2
            imarr2 = fits2[i].data**2
            imarr = imarr1 - imarr2
            #imarr = imarr2 / imarr1
            image_tmp.append(imarr)
        image = SingleImageIR(image_tmp,first_line,first_s_over,first_col,first_p_over,ccd_e2v)
        plt.figure()
        #rms
        #plt.imshow(image, origin='lower', vmin=-1, vmax=1, cmap = 'hot')
        plt.imshow(image, origin='lower', vmin=-2, vmax=10, cmap = 'hot')
        #plt.colorbar()
        clb = plt.colorbar()
        #clb.ax.set_title('$ADU^{2}$')
        plt.title(title)
        plt.savefig(outpath_final3+'/image.png')
        plt.close()

        #metrics per amp
        amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
        fits=pyfits.open(path_input2+'/'+rafts[i_raft]+'/'+ccds[i_ccd]+'/'+image_rms_file2)
        image_tmp=[]
        mean = np.zeros(16)
        rms = np.zeros(16)
        q95 = np.zeros(16)
        skewness = np.zeros(16)
        for m in range(16) :
            imarr1 = fits[m].data
            values = imarr1.flatten()
            #print(values)
            q95[m]=round(np.percentile(values,95),3)
            #print('====q95')
            #print(q95)
            mean[m] = round(np.mean(values),3)
            #print(mean[m])
            rms[m] = round(np.std(values),3)
            #print(rms[m])
            #skewness
            skewness[m] = round(skew(values),3)
            #print(skewness[m])
            plt.figure()
            bin_range = [0,10]
            nbins = 1000
            plt.figure()
            plt.hist(values, range=bin_range, bins=nbins, label=amps[m], histtype='step', color = 'blue')
            plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
            #plt.ylim([23000, 25000])
            #line=variances_exposures[10]['mean_column'][l]
            #mean=np.mean(line)
            #plt.ylim([mean-10,mean+30])
            plt.xlabel('RMS')
            plt.ylabel('Pixels')
            #plt.title('Run ' + run + ' ' + str(rafts[i])+' '+str(ccds[j])+' amp '+str(l+1))
            plt.legend(["mean={m}\nRMS={r}\nq95={q}\nskew={sk}".format(m=mean[m],r=rms[m],q=q95[m],sk=skewness[m])])
            plt.savefig(outpath_final2+'/RMS_amp_'+str(m+1)+'.png')
            plt.close()
            #yellow corner
            
            
        t_variance = Table([mean, rms, q95, skewness], names=('mean', 'rms', 'q95', 'skewness'), meta={'name': 'Variances'})
        print(t_variance)
        t_variance.write(path_output + '/results/Variances_' + rafts[i_raft] + '_'+ ccds[i_ccd] + '.fits', overwrite=True)
        continue

        #overscan studies    
        fits=pyfits.open(path_input+'/'+rafts[i_raft]+'/'+ccds[i_ccd]+'/'+image_rms_file)
        image_tmp=[]
        var_line = np.zeros((16,im_y_size))
        mean_line = np.zeros((16,im_y_size))
        var_column =  np.zeros((16,im_x_size))
        mean_column = np.zeros((16,im_x_size))
        var_line_overscan = np.zeros((16,im_y_size))
        mean_line_overscan = np.zeros((16,im_y_size))
        var_column_overscan =  np.zeros((16,im_x_size))
        mean_column_overscan = np.zeros((16,im_x_size))
        for i in range(16) :
            imarr = fits[i].data
            imarr_science = imarr[first_line:first_p_over,first_col:first_s_over]
            imarr_overscan_serial = imarr[first_line:first_p_over,first_s_over+1:amp_x_size]
            imarr_overscan_parallel = imarr[first_p_over+1:amp_y_size,first_col:first_s_over]
            #imarr = imarr2 / imarr1
            var_line[i,:] = np.var(imarr_science,axis=1)
            mean_line[i,:] = np.mean(imarr_science,axis=1)
            var_column[i,:] = np.var(imarr_science,axis=0)
            mean_column[i,:] = np.mean(imarr,axis=0)
            var_line_overscan[i,:] = np.var(imarr_overscan_serial,axis=1)
            mean_line_overscan[i,:] = np.mean(imarr_overscan_serial,axis=1)
            var_column_overscan[i,:] = np.var(imarr_overscan_parallel,axis=0)
            mean_column_overscan[i,:] = np.mean(imarr_overscan_parallel,axis=0)
        image_tmp.append(imarr)
        image = SingleImageIR(image_tmp,first_line,first_s_over,first_col,first_p_over,ccd_e2v)
        plt.figure()
        #rms 
        plt.imshow(image, origin='lower', vmin=-1, vmax=1, cmap = 'hot')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_final3+'/image_'+ccds[i_ccd]+'.png')
        plt.close()
        #debug
        print('========overscan')
        print(mean_line[0])
        print(var_line_overscan[0])

sys.exit()
