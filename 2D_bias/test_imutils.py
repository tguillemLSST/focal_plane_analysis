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
from astropy.io import fits
import matplotlib
import os
from astropy.visualization import (ImageNormalize,PercentileInterval,HistEqStretch)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc 
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
#file_to_read='R14_S22_13159_pca_bias.fits'
file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/13161/dark_dark_032/MC_C_20211212_000169_R14_S22.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/storage/20211212/MC_C_20211212_000157/MC_C_20211212_000157_R14_S22.fits'
#file_to_read='/sps/lsst/groups/FocalPlane/SLAC/run5/13161/bias_bias_013/MC_C_20211212_000150_R14_S22.fits'
run='13159'
raft='R14'
ccd='S22'

output_data='/sps/lsst/users/tguillem/web/debug/images_dark_run/169/' 

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
#fits=pyfits.open(file_list[ifile])
# print the image main information 
#print(fits[0].header.tostring(sep='\n', endcard=True, padding=True))
#print(fits)

###test imutils
#hdulist = fits.open(ffile)
print(file_list[ifile])
fits = fits.open(file_list[ifile])
print(fits)
fits.info()
#bad_segs = iu.get_union_of_bad_column_segs(hdulist)
#hdu_s = hdulist[5]
pbias = "bycol"
sbias = "byrowe2v"
#iu.subtract_bias(sbias, pbias, hdu_s)
fig=plt.figure(figsize=[25,20])
title='Image per amplifier :  (70%s percentile) \n%s' % ('%',os.path.basename(file_list[ifile]))
image_txt='ImagePerAmp_imutils'
plt.suptitle(title)
for i in range(16) :
    iu.subtract_bias(sbias, pbias,fits[i+1])
    norm = ImageNormalize(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], interval=PercentileInterval(70.))
    plt.subplot(2,8,i+1,title=i+1)
    #plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over],cmap = 'hot',origin='lower',norm=norm)
    plt.imshow(fits[i+1].data[first_line:first_p_over,first_col:first_s_over], vmin=-5, vmax=5, cmap = 'hot', origin='lower')
    if not(i%8 ==0) :
        figure=plt.gca()
        y_axis = figure.axes.get_yaxis()
        y_axis.set_visible(False)
    plt.colorbar()
SaveFig(fig,image_txt,run_cur=run,raft_cur=raft,ccd_cur=ccd,hdu=0)
plt.close()
print('=======DONE')
#sys.exit()
###end test imutils
