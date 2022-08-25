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
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/'+run+'/'+directory+'/'+file

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
directory1='dark_dark_*'
directory2='bias_bias_*'
#file='MC_C_*_'+raft+'_'+ccd+'.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/'+run+'/'+directory+'/'+file

#output_data='/sps/lsst/users/tguillem/web/debug/images_dark_run/batch/'

print('Configuration arguments: ', str(sys.argv))
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
repo_path=str(sys.argv[4])
output_data=str(sys.argv[5])

run=str_run_all
raft=str_raft
ccd=str_all_sensors
file='MC_C_*_'+raft+'_'+ccd+'.fits'
file_to_read1='/sps/lsst/groups/FocalPlane/SLAC/run5/'+run+'/'+directory1+'/'+file
file_to_read2='/sps/lsst/groups/FocalPlane/SLAC/run5/'+run+'/'+directory2+'/'+file
file_list=glob.glob(file_to_read1)+glob.glob(file_to_read2)
file_list.sort()
print(file_list)

###HACK
#file_list=['/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/test_master_dark_3/u/tguillem/master_dark_0/20220621T090913Z/dark/dark_LSSTCam_R14_S22_u_tguillem_master_dark_0_20220621T090913Z.fits']
#file_list=['/sps/lsst/groups/FocalPlane/SLAC/run5/13161/dark_dark_032/MC_C_20211212_000169_R14_S22.fits']
#file_list=['/sps/lsst/groups/FocalPlane/SLAC/run5/13162/bias_bias_000/MC_C_20211212_000177_R41_S21.fits', '/sps/lsst/groups/FocalPlane/SLAC/run5/13162/bias_bias_001/MC_C_20211212_000178_R41_S21.fits']

print(file_list)

fits=pyfits.open(file_list[0])
# print the image main information
print(fits[0].header.tostring(sep='\n', endcard=True, padding=True))

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

print('Number of lines read=',amp_y_size,' Number of colomns read=',amp_x_size)
print('Number of lines in Image area=',im_y_size,' Number of colomns in Image area=',im_x_size)
print('First line in Image area=',first_line,' First column in Image area=',first_col)

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

def SingleImageIR(image,first_col=first_col,first_cover=first_s_over,first_line=first_line,first_lower=first_p_over):
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

#compute a masterbias (adapted from eotest/python/lsst/eotest/sensor/ccd_bias_pca.py)
outfile = 'masterbias.FITS'
path_output = os.path.join(output_data,run,raft,ccd)
os.makedirs(path_output)

images = []
for i in range(16):
            images.append([])
medianed_images = dict()
statistic=afwMath.MEAN
for ifile in range(len(file_list)):
    fits=pyfits.open(file_list[ifile])
    for i in range(0,16):
        imarr = fits[i+1].data
        imarr_np = np.array(imarr, dtype=np.float32)
        images[i].append(afwImage.ImageF(imarr_np))
    fits.close()
for i in range(0,16):
    medianed_images[i] = afwMath.statisticsStack(images[i], statistic)
with pyfits.open(file_list[0]) as hdus:
    for amp, image in medianed_images.items():
        hdus[amp].data = image.array
        hdus[0].header['FILENAME'] = os.path.basename(outfile)
        fitsWriteto(hdus, path_output+'/'+outfile, overwrite=True)
print('Master bias after correction produced')
#plot
fig=plt.figure(figsize=[25,20])
title='CCD image: master bias'
image_txt='CCD_Image_master_bias'
plt.suptitle(title, fontsize=20)
fits=pyfits.open(path_output+'/'+outfile)
image_tmp=[]
for i in range(16) :
    imarr = fits[i].data
    mean = np.mean(imarr[first_line:first_p_over,first_col:first_s_over])
    imarr = imarr - mean
    image_tmp.append(imarr)
print(image_tmp)
image = SingleImageIR(image_tmp)
norm = ImageNormalize(image, interval=PercentileInterval(70.))
plt.imshow(image, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
plt.colorbar()
SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
fits.close()
#sys.exit()

#prepare corrected masterbias
outfile = 'masterbias_corr_2D.FITS'
images_corr_2D = []
for i in range(16):
        images_corr_2D.append([])
medianed_images_corr_2D = dict() 

run_base = run
#plots
print(len(file_list))
#ifile=0
for ifile in range(len(file_list)): 
    #if(ifile==1):
    #    break
    
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
    frame = (file_list[ifile].partition('/MC')[0])[-13:]
    run = run_base + '/' + frame
    print(run)

    # plot the raw image / amplifiers 
    fig=plt.figure(figsize=[25,20])
    title='Image per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    image_txt='RawImagePerAmp'
    plt.suptitle(title)
    for i in range(16) :
        norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
        plt.subplot(2,8,i+1,title=i+1)
        plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        if not(i%8 ==0) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        plt.colorbar()
    #save also plots separately
    for i in range(16) :
        norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
        fig_amp=plt.figure()
        plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        plt.colorbar()
        #SaveFig(fig_amp,image_txt+'_amp_'+str(i+1),run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    plt.close(fig)
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)

    # plot all the amplifiers on one image
    nb_amp=16

    fig=plt.figure(figsize=[25,20])
    title='Master bias for PCA. Unbiased image per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    plt.suptitle(title)
    first_col=FileUnBias.all_file[ifile].first_col
    first_s_over=FileUnBias.all_file[ifile].first_s_over
    first_line=FileUnBias.all_file[ifile].first_line
    first_p_over=FileUnBias.all_file[ifile].first_p_over
    
    for i in range(nb_amp) :
        norm = ImageNormalize(FileUnBias.all_file[ifile].Image[i][first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
        plt.subplot(2,8,i+1,title=i+1)
        plt.imshow(FileUnBias.all_file[ifile].Image[i][first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
        if not(i%8 ==0) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        plt.colorbar()
    image_txt='UnbiasedImagePerAmp'
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)

    #try one image per CCD
    #raw image
    fig=plt.figure(figsize=[25,20])
    title='Master bias for PCA. Raw image per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    plt.suptitle(title)
    #create an image (fix the 0/1 index issue of the object fits)
    image_tmp=[]
    for i in range(nb_amp):
        imarr = fits[i+1].data
        mean = np.mean(imarr[first_line:first_p_over,first_col:first_s_over])
        imarr = imarr - mean
        image_tmp.append(imarr)
    image = SingleImageIR(image_tmp)
    norm = ImageNormalize(image, interval=PercentileInterval(70.))
    plt.imshow(image,cmap = 'hot',origin='lower',norm=norm)
    plt.colorbar()
    image_txt='CCD_RawImage'
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0) 

    fig=plt.figure(figsize=[25,20])
    title='Master bias for PCA. Unbiased image per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    plt.suptitle(title) 
    image = SingleImageIR(FileUnBias.all_file[ifile].Image)
    norm = ImageNormalize(image, interval=PercentileInterval(70.))
    plt.imshow(image,cmap = 'hot',origin='lower',norm=norm)
    plt.colorbar()
    image_txt='CCD_UnbiasedImage'
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0) 

    #plot overscans
    # plot the raw image / amplifiers 
    fig=plt.figure(figsize=[25,20])
    title='Master bias for PCA. Overscan serial per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    image_txt='RawOverscanSerialPerAmp'
    plt.suptitle(title)
    for i in range(16) :
        norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_s_over+1:amp_x_size], interval=PercentileInterval(70.))
        ax = plt.subplot(2,8,i+1,title=i+1)
        plt.imshow(fits[i+1].data[first_line:first_p_over,first_s_over+1:amp_x_size],cmap = 'hot',origin='lower',norm=norm)
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

    # plot the raw image / amplifiers 
    fig=plt.figure(figsize=[25,20])
    title='Master bias for PCA. Overscan parallel per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    image_txt='RawOverscanParallelPerAmp'
    plt.suptitle(title)
    for i in range(16):
        norm = ImageNormalize(fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over], interval=PercentileInterval(70.))
        ax = plt.subplot(2,8,i+1,title=i+1)
        plt.imshow(fits[i+1].data[first_p_over+1:amp_y_size,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
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

    ###test direct 2D correction
    print('++++++++++++++++2D')
    fig=plt.figure(figsize=[25,20])
    #title='Image per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
    title='Image per amplifier: direct 2D correction from overscan\nFile: %s\nExposure time: %s' % (os.path.basename(file_list[ifile]),exp_time)
    image_txt='RawImagePerAmp_2D_corr'
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
    var_total_corr_2D =  np.zeros(16)
    mean_total_corr_2D = np.zeros(16)
    var_line_corr_2D = np.zeros((16,im_y_size))
    mean_line_corr_2D = np.zeros((16,im_y_size))
    var_column_corr_2D =  np.zeros((16,im_x_size))
    mean_column_corr_2D = np.zeros((16,im_x_size))
    for i in range(16) :
        imarr = fits[i+1].data    
        #variances for image before correction
        imarr_science = imarr[first_line:first_p_over,first_col:first_s_over]
        #compute variances
        var_total[i] = np.var(imarr_science)
        mean_total[i] = np.mean(imarr_science)
        var_line[i,:] = np.var(imarr_science,axis=1)
        mean_line[i,:] = np.mean(imarr_science,axis=1)
        var_column[i,:] = np.var(imarr_science,axis=0)
        mean_column[i,:] = np.mean(imarr_science,axis=0)
        ###2D correction
        mean_over_per_line=np.mean(imarr[:,first_s_over+2:],axis=1)
        rawl=np.zeros((last_l-first_p_over-2,last_s))
        for l in range(first_p_over+2,last_l):
            rawl[l-first_p_over-2,:]=imarr[l,:]-mean_over_per_line[l]
            #####################HACK to apply here line correction#####################
            mean_over_per_column=np.mean(rawl[:,:],axis=0)
            #next line is the hack
            #mean_over_per_column=np.mean(rawl[:,:],axis=0)-np.mean(rawl[:,:],axis=0)
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
        #compute variances
        var_total_corr_2D[i] = np.var(imarr_science_corr_2D)
        mean_total_corr_2D[i] = np.mean(imarr_science_corr_2D)
        var_line_corr_2D[i,:] = np.var(imarr_science_corr_2D,axis=1)
        mean_line_corr_2D[i,:] = np.mean(imarr_science_corr_2D,axis=1)
        var_column_corr_2D[i,:] = np.var(imarr_science_corr_2D,axis=0)
        mean_column_corr_2D[i,:] = np.mean(imarr_science_corr_2D,axis=0)

        imarr_np = np.array(imarr, dtype=np.float32)
        images_corr_2D[i].append(afwImage.ImageF(imarr_np))

    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    plt.close(fig)
    t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, var_total_corr_2D, mean_total_corr_2D, var_line_corr_2D, mean_line_corr_2D, var_column_corr_2D, mean_column_corr_2D], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column', 'var_total_corr_2D', 'mean_total_corr_2D', 'var_line_corr_2D', 'mean_line_corr_2D', 'var_column_corr_2D', 'mean_column_corr_2D'), meta={'name': 'Variances'})
    print(t_variance)
    t_variance.write(output_data + run + '/' + raft + '/' + ccd + '/Variances_2D_corr.fits', overwrite=True)
    #continue

    #try one image per CCD
    #2D-corrected image
    fig=plt.figure(figsize=[25,20])
    title='2D-corrected CCD image: direct 2D correction from overscan\nFile: %s\nFrame: %s\nExposure time: %s' % (os.path.basename(file_list[ifile]),img_type,exp_time)
    image_txt='CCD_Image_2D_corr'
    plt.suptitle(title, fontsize=20)
    #create an image (fix the 0/1 index issue of the object fits)
    image = SingleImageIR(image_tmp)
    norm = ImageNormalize(image, interval=PercentileInterval(70.))
    #plt.imshow(image,cmap = 'hot',origin='lower',norm=norm)
    plt.imshow(image, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
    plt.colorbar()
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0) 

    ###test imutils correction
    print('++++++++++++++++imutils')
    fig=plt.figure(figsize=[25,20])
    title='Image per amplifier: imutils correction\nFile: %s\nExposure time: %s' % (os.path.basename(file_list[ifile]),exp_time)
    image_txt='RawImagePerAmp_imutils_corr'
    plt.suptitle(title, fontsize=20)
    pbias = "bycol"
    sbias = "byrowe2v"
    #create an image (fix the 0/1 index issue of the object fits)
    image_tmp=[]
    for i in range(16) :
        iu.subtract_bias(sbias,pbias,fits[i+1])
        imarr = fits[i+1].data
        image_tmp.append(imarr)
        plt.subplot(2,8,i+1,title=i+1)
        plt.imshow(imarr[first_line:first_p_over,first_col:first_s_over], vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        if not(i%8 ==0) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        plt.colorbar()
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    plt.close(fig)

    #try one image per CCD
    #imutils-corrected image
    fig=plt.figure(figsize=[25,20])
    title='imutils-corrected CCD image\nFile: %s\nFrame: %s\nExposure time: %s' % (os.path.basename(file_list[ifile]),img_type,exp_time)
    image_txt='CCD_Image_imutils_corr'
    plt.suptitle(title, fontsize=20)
    #create an image (fix the 0/1 index issue of the object fits)
    image = SingleImageIR(image_tmp)
    plt.imshow(image, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
    plt.colorbar()
    SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0) 

    fits.close()
    
    print('free memory')
    print(dir())
    del image
    del image_tmp
    del imarr
    del imarr_science
    del imarr_science_corr_2D
    del imarr_np
    gc.collect()

#master bias after correction
for amp in range(0,16):
    medianed_images_corr_2D[amp] = afwMath.statisticsStack(images_corr_2D[amp], statistic)
with pyfits.open(file_list[0]) as hdus:
    for amp, image in medianed_images_corr_2D.items():
        hdus[amp].data = image.array
        hdus[0].header['FILENAME'] = os.path.basename(outfile)
        fitsWriteto(hdus, path_output+'/'+outfile, overwrite=True) 
print('Master bias after correction produced')

#plot
#2D corr
fig=plt.figure(figsize=[25,20])
title='2D-corrected CCD image: master bias'
image_txt='CCD_Image_master_bias_2D_corr'
plt.suptitle(title, fontsize=20)
fits=pyfits.open(path_output+'/'+outfile)
image_tmp=[]
for i in range(16) :
    imarr = fits[i].data
    image_tmp.append(imarr)
#print(image_tmp)
image = SingleImageIR(image_tmp)
norm = ImageNormalize(image, interval=PercentileInterval(70.))
plt.imshow(image, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
plt.colorbar()
SaveFig(fig,image_txt,run_cur=run_base,raft_cur=raft,ccd_cur=ccd,hdu=0)
fits.close()
    
sys.exit()














#try one image per raft ==> NOT YET WORKING
fig,ax = plt.subplots(3,3,figsize=(15,15))
line=2
for row in ax:
    column=0
    for col in row:
        i=column+3*line
        print(i)
        out1=SingleImageIR(FileUnBias.all_file[ifile].Image[i])
        norm = ImageNormalize(out1, interval=PercentileInterval(80.))
        #norm = ImageNormalize(out1,  stretch=HistEqStretch(out1))
        im = col.imshow(out1,cmap = 'hot',origin='lower',norm=norm)
        col.set_title(os.path.basename(file_list[i]))
        divider = make_axes_locatable(col)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        if not(i%3 ==0) :
            #figure=plt.gca()
            y_axis = col.axes.get_yaxis()
            y_axis.set_visible(False)
            
        fig.colorbar(im, cax=cax, orientation='vertical')
        column+=1
        del(out1)
        del(norm)
    line-=1
plt.show()
SaveFig(fig,'CCDUnbiasedImage',run_cur=run,raft_cur=raft)
