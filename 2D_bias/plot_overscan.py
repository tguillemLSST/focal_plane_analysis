#!/usr/bin/env python
# coding: utf-8

##########
# Author: P. Antilogus
# Adapted by T. Guillemin
# Goal: plot images and overscans per amplifier, per CCD and per raft
##########

from eochar.bot_frame_op import *
from eochar.GetPhotoFlux import *

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
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
#correct line: from lsst.eotest.fitsTools import fitsWriteto
#hack to work in batch
from eotest.fitsTools import fitsWriteto
#imutils copied for now
import imutils as iu

# ==config
run='13202'
raft='R13'
ccd='S11'
directory='/xtalk_-117.0_116.8_*' 
file='MC_C_*_'+raft+'_'+ccd+'.fits'
# Method for Bias correction
#UnBias='1D'
UnBias='2D'

#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run6_per_run/'+run+'/'+directory+'/'+file

#HACK
#select one file by hand
#file_to_read='R14_S22_13159_pca_bias.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/13161/dark_dark_032/MC_C_20211212_000169_R14_S22.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/storage/20211212/MC_C_20211212_000157/MC_C_20211212_000157_R14_S22.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/13161/bias_bias_013/MC_C_20211212_000150_R14_S22.fits'
#list of files
run='13161'
raft='R14'
ccd='S21'
#directory='dark_dark_*'
directory1='dark_dark_xxxxxxxxxxxxx*'
directory2='bias_bias_*'
#directory2='bias_*'
#file='MC_C_*_'+raft+'_'+ccd+'.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/'+run+'/'+directory+'/'+file

#output_data='/sps/lsst/users/tguillem/web/debug/images_dark_run/batch/'

print('Configuration arguments: ', str(sys.argv))
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
repo_path=str(sys.argv[4])
output_data=str(sys.argv[5])
overscan_correction=str(sys.argv[6])

run=str_run_all
raft=str_raft
ccd=str_all_sensors
file='MC_C_*_'+raft+'_'+ccd+'.fits'
file_to_read1='/sps/lsst/groups/FocalPlane/SLAC/run6_per_run/'+run+'/'+directory1+'/'+file
file_to_read2='/sps/lsst/groups/FocalPlane/SLAC/run6_per_run/'+run+'/'+directory2+'/'+file
file_list=glob.glob(file_to_read1)+glob.glob(file_to_read2)

file_list=['/sps/lsst/groups/FocalPlane/SLAC/run6_per_run/13555/dark_dark_0381/MC_C_20231116_000781_R23_S12.fits','/sps/lsst/groups/FocalPlane/SLAC/run6_per_run/13555/dark_dark_0688/MC_C_20231116_001088_R23_S12.fits']
print(file_list)
##########################

# compute unbiased image
#not working: FileRaw=InFile(dirall=file_list,Slow=False,verbose=False,Bias='NoCorr')
FileUnBias=InFile(dirall=file_list,Slow=False,verbose=False,Bias=UnBias) 

# get the image area ( overscan excluded ) 
first_col=FileUnBias.all_file[0].first_col
first_s_over=FileUnBias.all_file[0].first_s_over
first_line=FileUnBias.all_file[0].first_line
first_p_over=FileUnBias.all_file[0].first_p_over
amp_y_size=len(FileUnBias.all_file[0].Image[0][:,0])
amp_x_size=len(FileUnBias.all_file[0].Image[0][0,:])
im_y_size=first_p_over-first_line
im_x_size=first_s_over-first_col
#define the amplifier corner
n_line_amp_corner = first_line + 200
n_col_amp_corner = first_col + 200
print('Number of lines read=',amp_y_size,' Number of colomns read=',amp_x_size)
print('Number of lines in Image area=',im_y_size,' Number of colomns in Image area=',im_x_size)
print('First line in Image area=',first_line,' First column in Image area=',first_col)
print('First line in overscan area=',first_p_over,' First column in overscan area=',first_s_over)

fits=pyfits.open(file_list[0])
# print the image main information
print(fits[0].header.tostring(sep='\n', endcard=True, padding=True))
print(fits)
#sys.exit()
#define CCD type 
#rafts_e2v = ['R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
rafts_itl = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43']
ccd_e2v = True
if raft in rafts_itl:
     ccd_e2v = False
print('ccd_ee2v: ' + str(ccd_e2v))
#select overscan correction method
overscan_1D = False
if(overscan_correction=='1D'):
     overscan_1D = True
print('overscan_1D: ' + str(overscan_1D))

###functions
def SaveFig(fig,rawPlotFile,run_cur='',raft_cur='',ccd_cur='',hdu=0):
    if hdu<1 : 
        root_plt=os.path.join(output_data,run_cur,raft_cur,ccd_cur)
    else : 
        hdu_cur='%d' % (hdu)
        root_plt=os.path.join(output_data,run_cur,raft_cur,ccd_cur,hdu_cur)
    # 
    PlotFile=rawPlotFile.replace('.','_')
    PlotFile=PlotFile.replace(' ','_')
    os.makedirs(root_plt,exist_ok=True)
    plotfile=os.path.join(root_plt,PlotFile)
    print ('PlotFile=',plotfile)
    fig.savefig(plotfile,bbox_inches='tight')
    plt.close(fig) 
    return

def SingleImageIR(image,first_col=first_col,first_cover=first_s_over,first_line=first_line,first_lower=first_p_over,is_e2v=ccd_e2v):
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
        
def SingleImageIR_old(image,first_col=first_col,first_cover=first_s_over,first_line=first_line,first_lower=first_p_over):
        # Display an IR2 image , with amplifiers set at the right place ...there is a DM version which does this better...
        # but here you are in stand alone 
        # the default associated to the image area (pre-overscan excluded) are for e2v IR2 files 
        #
        col_size=first_cover-first_col
        line_size=first_lower-first_line
        #
        spf=np.zeros((line_size*2,col_size*8))
        for i in range(16) :
             if i<8 :
                xx=i*col_size-1
                yy=0
                for x in range(first_col,first_cover) :  
                    spf[yy:yy+line_size,xx+col_size-(x-first_col)]=image[i][first_line:first_lower,x]
             else :
                xx=(15-i)*col_size
                yy=-1
                for y in range(first_line,first_lower) :  
                    spf[yy+2*line_size-(y-first_line),xx:xx+col_size]=image[i][y,first_col:first_cover]
                    
        return spf   
###end functions

images_corr_2D = []
for i in range(16):
     images_corr_2D.append([])

#plots
print(len(file_list))
run_base = run
#ifile=0
for ifile in range(len(file_list)):
    # access directly the fits file image 
    fits=pyfits.open(file_list[ifile])
    # print the image main information 
    print(fits[0].header.tostring(sep='\n', endcard=True, padding=True))
    exp_time = str(fits[0].header['DARKTIME'])
    print(exp_time)
    img_type = str(fits[0].header['IMGTYPE'])
    print(img_type)
    #continue
    
    #add image number to run
    print(file_list[ifile])
    frame = (file_list[ifile].partition('/MC')[0])[-4:]
    #frame = file_list[ifile][-19:-13]
    run = run_base + '/' + frame
    print(run)
    
    # plot the raw images per amplifier 
    fig=plt.figure(figsize=[25,20])
    title='Raw image per amplifier (mean per amp substracted)'#  :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    image_txt='Amp_RawImage'
    plt.suptitle(title, fontsize=20)
    for i in range(16) :
        #norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
        plt.subplot(2,8,i+1,title=i+1)
        imarr = fits[i+1].data
        mean = np.mean(imarr[first_line:first_p_over,first_col:first_s_over])
        imarr = imarr - mean
        #plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        plt.imshow(imarr, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        if not(i%8 ==0) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        plt.colorbar()
    #save also plots separately
    #for i in range(16) :
    #    norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
    #    fig_amp=plt.figure()
    #    plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
    #    plt.colorbar()
    #    #SaveFig(fig_amp,image_txt+'_amp_'+str(i+1),run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    #plt.close(fig)
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)

    # plot all the amplifiers on one image
    nb_amp=16
    #try one image per CCD
    #raw image
    fig=plt.figure(figsize=[25,20])
    title='Raw image (mean per amp substracted)'#\nFile: %s\nExposure time: %s\n%s %s' % (os.path.basename(file_list[ifile]),exp_time,raft,ccd)
    plt.suptitle(title, fontsize=20)
    #create an image (fix the 0/1 index issue of the object fits)
    image_tmp=[]
    for i in range(nb_amp):
        imarr = fits[i+1].data
        mean = np.mean(imarr[first_line:first_p_over,first_col:first_s_over])
        imarr = imarr - mean
        image_tmp.append(imarr)
        #if(i==9):
        #    #print(imarr)
        #    for i in range(amp_x_size):
        #        print('========column ' + str(i))
        #        for j in range(amp_y_size):
        #            print('column ' + str(i) + ' row ' + str(j) + ' :' + str(imarr[j][i]))
        #    sys.exit()
    image = SingleImageIR(image_tmp)
    norm = ImageNormalize(image, interval=PercentileInterval(70.))
    #plt.imshow(image,cmap = 'hot',origin='lower',norm=norm)
    plt.imshow(image, vmin=-5, vmax=5, cmap = 'hot', origin='lower') 
    plt.colorbar()
    image_txt='CCD_RawImage'
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    
    #plot overscans
    # plot the raw image / amplifiers 
    fig=plt.figure(figsize=[25,20])
    title='Overscan serial per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    image_txt='Amp_RawOverscanSerial'
    plt.suptitle(title, fontsize=20)
    for i in range(16) :
        #norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_s_over+1:amp_x_size], interval=PercentileInterval(70.))
        ax = plt.subplot(2,8,i+1,title=i+1)
        imarr=fits[i+1].data[first_line:first_p_over,first_s_over+1:amp_x_size]-np.mean(fits[i+1].data[first_line:first_p_over,first_s_over+1:amp_x_size])
        #plt.imshow(fits[i+1].data[first_line:first_p_over,first_s_over+1:amp_x_size],cmap = 'hot',origin='lower',norm=norm)
        plt.imshow(imarr, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        if not(i%8 ==0) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        plt.colorbar()
        ax.set_aspect('auto')
    #save also plots separately
    for i in range(16) :
        norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
        fig_amp=plt.figure()
        ax = plt.subplot(1,1,1)
        ax.set_aspect('auto')
        plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        plt.colorbar()
        #SaveFig(fig_amp,image_txt+'_amp_'+str(i+1),run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
        plt.close(fig)
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    #print('OVERSCAN-----------')
    #continue

    # plot the raw image / amplifiers 
    fig=plt.figure(figsize=[25,20])
    title='Overscan parallel per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    image_txt='Amp_RawOverscanParallel'
    plt.suptitle(title, fontsize=20)
    for i in range(16):
        #norm = ImageNormalize(fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over], interval=PercentileInterval(70.))
        imarr=fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over]-np.mean(fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over])
        ax = plt.subplot(2,8,i+1,title=i+1)
        #plt.imshow(fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        plt.imshow(imarr, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        if not(i%8 ==0) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        plt.colorbar()
        ax.set_aspect('auto')
    #save also plots separately
    for i in range(16) :
        norm = ImageNormalize(fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over], interval=PercentileInterval(70.))
        fig_amp=plt.figure()
        ax = plt.subplot(1,1,1)
        ax.set_aspect('auto')
        plt.imshow(fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        plt.colorbar()
        #SaveFig(fig_amp,image_txt+'_amp_'+str(i+1),run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
        plt.close(fig_amp)
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)

    ######debugging mean/median
    for i in range(1,2) :
         imarr = fits[i+1].data
         mean_image = np.mean(imarr[first_line:first_p_over,first_col:first_s_over])
         mean_serial = np.mean(imarr[first_line:first_p_over,first_s_over+1:amp_x_size])
         diff = mean_image - mean_serial
         median_image = np.median(imarr[first_line:first_p_over,first_col:first_s_over])
         median_serial = np.median(imarr[first_line:first_p_over,first_s_over+1:amp_x_size])
         diff_median = median_image - median_serial
         print('amp ' + str(i))
         print('mean_image = ' + str(mean_image))
         print('mean_serial = ' + str(mean_serial))
         print('diff = ' + str(diff))
         print('median_image = ' + str(median_image))
         print('median_serial = ' + str(median_serial))
         print('diff_median = ' + str(diff_median))
         #debug lines
         mean_line_overscan = np.zeros(im_y_size)
         median_line_overscan = np.zeros(im_y_size)
         imarr_science = imarr[first_line:first_p_over,first_col:first_s_over]
         imarr_overscan_serial = imarr[first_line:first_p_over,first_s_over+1:amp_x_size]
         mean_line_overscan[:] = np.mean(imarr_overscan_serial,axis=1)
         median_line_overscan[:] = np.median(imarr_overscan_serial,axis=1)
         print(mean_line_overscan)
         print(median_line_overscan)
         print('Comparison')
         print('mean_serial = ' + str(mean_serial))
         average_median = 1/len(median_line_overscan)*np.sum(median_line_overscan)
         print('1/n sum(median) = ' + str(average_median))
         average_mean = 1/len(mean_line_overscan)*np.sum(mean_line_overscan)
         print('1/n sum(mean) = ' + str(average_mean))
         print('Final')
         print('diff = ' + str(diff))
         diff_average_median = mean_image - average_median
         print('diff_average_median = ' + str(diff_average_median))
         #sys.exit()    
    continue
    ###########################################################################
    
    ###test direct 2D correction
    print('++++++++++++++++2D')
    fig=plt.figure(figsize=[25,20])
    #title='Image per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    title='Image per amplifier: direct 2D correction from overscan\nFile: %s\nExposure time: %s\n%s %s' % (os.path.basename(file_list[ifile]),exp_time,raft,ccd)
    image_txt='Amp_Image_corr'
    last_l=len(fits[1].data[:,0])
    last_s=len(fits[1].data[0,:])
    plt.suptitle(title, fontsize=20)
    #create an image (fix the 0/1 index issue of the object fits)
    image_tmp=[]
    #save variances in a table
    var_total =  np.zeros(16)
    mean_total = np.zeros(16)
    var_line = np.zeros((16,im_y_size))
    mean_line = np.zeros((16,im_y_size))
    var_column =  np.zeros((16,im_x_size))
    mean_column = np.zeros((16,im_x_size))
    var_line_overscan = np.zeros((16,im_y_size))
    mean_line_overscan = np.zeros((16,im_y_size))
    var_column_overscan =  np.zeros((16,im_x_size))
    mean_column_overscan = np.zeros((16,im_x_size))
    var_amp_corner = np.zeros(16)
    mean_amp_corner = np.zeros(16)
    var_total_corr_2D = np.zeros(16)
    mean_total_corr_2D = np.zeros(16)
    var_line_corr_2D = np.zeros((16,im_y_size))
    mean_line_corr_2D = np.zeros((16,im_y_size))
    var_column_corr_2D = np.zeros((16,im_x_size))
    mean_column_corr_2D = np.zeros((16,im_x_size))
    var_line_corr_2D_overscan = np.zeros((16,im_y_size))
    mean_line_corr_2D_overscan = np.zeros((16,im_y_size))
    var_column_corr_2D_overscan = np.zeros((16,im_x_size))
    mean_column_corr_2D_overscan = np.zeros((16,im_x_size))
    var_amp_corner_corr_2D = np.zeros(16)
    mean_amp_corner_corr_2D = np.zeros(16)
    for i in range(16) :
        imarr = fits[i+1].data    
        #variances for image before correction
        imarr_science = imarr[first_line:first_p_over,first_col:first_s_over]
        imarr_overscan_serial = imarr[first_line:first_p_over,first_s_over+1:amp_x_size]
        imarr_overscan_parallel = imarr[first_p_over+1:amp_y_size,first_col:first_s_over]
        #imarr_amp_corner = imarr[first_line:n_line_amp_corner,first_col:n_col_amp_corner]
        #compute variances
        var_total[i] = np.var(imarr_science)
        mean_total[i] = np.mean(imarr_science)
        var_line[i,:] = np.var(imarr_science,axis=1)
        mean_line[i,:] = np.mean(imarr_science,axis=1)
        var_column[i,:] = np.var(imarr_science,axis=0)
        mean_column[i,:] = np.mean(imarr_science,axis=0)
        var_line_overscan[i,:] = np.var(imarr_overscan_serial,axis=1)
        mean_line_overscan[i,:] = np.mean(imarr_overscan_serial,axis=1)
        var_column_overscan[i,:] = np.var(imarr_overscan_parallel,axis=0)
        mean_column_overscan[i,:] = np.mean(imarr_overscan_parallel,axis=0)
        var_amp_corner[i] = np.var(imarr_science[first_line:n_line_amp_corner,first_col:n_col_amp_corner])
        mean_amp_corner[i] = np.mean(imarr_science[first_line:n_line_amp_corner,first_col:n_col_amp_corner])
        #print('------Checking overscan')
        #print(mean_line_overscan[i])
        #print(mean_column_overscan[i])
        ###2D correction
        mean_over_per_line=np.mean(imarr[:,first_s_over+2:],axis=1)
        rawl=np.zeros((last_l-first_p_over-2,last_s))
        ###HACK
        #mean_over_per_line=np.mean(imarr[:,first_s_over+2:],axis=1)-np.mean(imarr[:,first_s_over+2:],axis=1)
        for l in range(first_p_over+2,last_l):
            rawl[l-first_p_over-2,:]=imarr[l,:]-mean_over_per_line[l]
            #####################HACK to apply here line correction#####################
            mean_over_per_column=np.mean(rawl[:,:],axis=0)
            #next line is the hack
            if(overscan_1D==True):
                 mean_over_per_column=np.mean(rawl[:,:],axis=0)-np.mean(rawl[:,:],axis=0)
            linef=mean_over_per_line[:,np.newaxis]
        # generate the 2D correction (thank's to numpy)
        over_cor_mean=mean_over_per_column+linef
        imarr = imarr - over_cor_mean
        ###
        #norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
        norm = ImageNormalize(imarr[first_line:first_p_over,first_col:first_s_over],interval=PercentileInterval(70.))
        plt.subplot(2,8,i+1,title=i+1)
        #plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        #plt.imshow(imarr[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        plt.imshow(imarr[first_line:first_p_over,first_col:first_s_over], vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        if not(i%8 ==0) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        plt.colorbar()
        image_tmp.append(imarr)
        imarr_science_corr_2D = imarr[first_line:first_p_over,first_col:first_s_over]
        imarr_amp_corner = imarr[first_line:n_line_amp_corner,first_col:n_col_amp_corner]
        #compute variances
        var_total_corr_2D[i] = np.var(imarr_science_corr_2D)
        mean_total_corr_2D[i] = np.mean(imarr_science_corr_2D)
        var_line_corr_2D[i,:] = np.var(imarr_science_corr_2D,axis=1)
        mean_line_corr_2D[i,:] = np.mean(imarr_science_corr_2D,axis=1)
        var_column_corr_2D[i,:] = np.var(imarr_science_corr_2D,axis=0)
        mean_column_corr_2D[i,:] = np.mean(imarr_science_corr_2D,axis=0)
        var_amp_corner_corr_2D[i] = np.var(imarr_science_corr_2D[first_line:n_line_amp_corner,first_col:n_col_amp_corner])
        mean_amp_corner_corr_2D[i] = np.mean(imarr_science_corr_2D[first_line:n_line_amp_corner,first_col:n_col_amp_corner])

        imarr_np = np.array(imarr, dtype=np.float32)
        images_corr_2D[i].append(afwImage.ImageF(imarr_np))
        
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    plt.close(fig)
    t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_amp_corner, mean_amp_corner, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D, var_amp_corner_corr_2D, mean_amp_corner_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column', 'var_amp_corner', 'mean_amp_corner', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D', 'var_amp_corner_corr_2D', 'mean_amp_corner_corr_2D'), meta={'name': 'Variances'})
    #print(t_variance)
    t_variance.write(output_data + run + '/' + raft + '/' + ccd + '/Variances_2D_corr.fits', overwrite=True)
    #continue

    #try one image per CCD
    #2D-corrected image
    fig=plt.figure(figsize=[25,20])
    title='2D-corrected CCD image: direct 2D correction from overscan\nFile: %s\nExposure time: %s\n%s %s' % (os.path.basename(file_list[ifile]),exp_time,raft,ccd) 
    image_txt='CCD_Image_corr'
    plt.suptitle(title, fontsize=20)
    #create an image (fix the 0/1 index issue of the object fits)
    image = SingleImageIR(image_tmp)
    norm = ImageNormalize(image, interval=PercentileInterval(70.))
    #plt.imshow(image,cmap = 'hot',origin='lower',norm=norm)
    plt.imshow(image, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
    plt.colorbar()
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0) 

    fits.close()
    
    print('free memory')
    print(dir())
    del image
    del image_tmp
    del imarr
    del imarr_np
    del imarr_science
    del imarr_science_corr_2D
    gc.collect()

sys.exit()
