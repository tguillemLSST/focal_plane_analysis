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
import matplotlib
import os
from astropy.visualization import (ImageNormalize,PercentileInterval,HistEqStretch)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc 

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
#file_to_read='R14_S22_13159_pca_bias.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/13161/dark_dark_032/MC_C_20211212_000169_R14_S22.fits'
file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/13161/bias_bias_013/MC_C_20211212_000150_R14_S22.fits'
run='13159'
raft='R14'
ccd='S22'

output_data='/sps/lsst/users/tguillem/web/debug/images_dark_run/' 

file_list=glob.glob(file_to_read)
file_list.sort()

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

#plots
ifile=0
# acces directly the fits file image 
fits=pyfits.open(file_list[ifile])
# print the image main information 
print(fits[0].header.tostring(sep='\n', endcard=True, padding=True))
#print(fits)

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
    SaveFig(fig_amp,image_txt+'_amp_'+str(i+1),run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    plt.close()
SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)

# plot all the amplifiers on one image
ifile=0
nb_amp=16
# 
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
    image_tmp.append(fits[i+1].data)
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
    SaveFig(fig_amp,image_txt+'_amp_'+str(i+1),run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    plt.close()
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
    SaveFig(fig_amp,image_txt+'_amp_'+str(i+1),run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
    plt.close()
SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)

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
