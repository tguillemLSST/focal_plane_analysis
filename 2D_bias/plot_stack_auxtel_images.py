#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: plot images directly from the butler
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
from scipy.stats import skew
import matplotlib
import os
import shutil
from astropy.visualization import (ImageNormalize,PercentileInterval,HistEqStretch)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc 
import lsst.daf.butler as dafButler
from conversion import *

print('Configuration arguments: ', str(sys.argv))
str_raft = str(sys.argv[1])
str_all_sensors = str(sys.argv[2])
collection_1 = str(sys.argv[3])
collection_2 = str(sys.argv[4])
output_data=str(sys.argv[5])

os.makedirs(output_data,exist_ok=True)

#butler access
#repo='/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/'
#repo='/sps/lsst/groups/auxtel/softs/shared/auxteldm_gen3/data/'
repo='/sps/lsst/groups/FocalPlane/SLAC/run6/butler/test_auxtel/main_240311/'
print('=========butler: ' + repo)

###calibrations
plot_calibrations = False
calibType='bias'
#calibType='dark'
###user collections
plot_images = True
datasetType='cpBiasProc'
#datasetType='raw'

print('=========collections: ' + collection_1 + ' and ' + collection_2)
butler_1 = dafButler.Butler(repo,collections=collection_1)
registry_1 = butler_1.registry
butler_2 = dafButler.Butler(repo,collections=collection_2)
registry_2 = butler_2.registry

#for a set of batch jobs
rafts=[str_raft]
ccds=[str_all_sensors]

ref_dir='/sps/lsst/users/tguillem/web/stack/run_6/reference/13391/focal_plane/'
shutil.copytree(ref_dir,output_data+'/focal_plane', dirs_exist_ok = True)
shutil.copytree(ref_dir,output_data+'/variance/focal_plane', dirs_exist_ok = True)
shutil.copytree(ref_dir,output_data+'/difference/focal_plane', dirs_exist_ok = True)

############################calibrations
os.makedirs(output_data+'/results',exist_ok=True)
if plot_calibrations:
    #loop over all CCDs
    for i_raft in range(len(rafts)):
        #create raft directory
        outpath_final = output_data+rafts[i_raft]
        os.makedirs(outpath_final,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final+'/index.html')
        #for variance
        outpath_final_variance = output_data+'variance/'+rafts[i_raft]
        os.makedirs(outpath_final_variance,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final_variance+'/index.html')
        #for difference
        outpath_final_difference = output_data+'difference/'+rafts[i_raft]
        os.makedirs(outpath_final_difference,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final_difference+'/index.html')
        #variance results
        os.makedirs(outpath_final+'/results',exist_ok=True)
        for i_ccd in range(len(ccds)):
            detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
            #to read official bias/dark calibrations
            #dataId = {'exposure': 2021121200150, 'detector': detector_id, "instrument": 'LATISS'}
            #calib_1 = butler_1.get(calibType,dataId=dataId)
            calib_1 = butler_1.get(calibType, instrument='LATISS', detector=detector_id)
            title = rafts[i_raft]+' '+ccds[i_ccd]
            # get full array
            arr_1 = calib_1.getImage().getArray()
            plt.figure()
            #plt.gcf().subplots_adjust(left = 0.05, bottom = 0.05,right = 0.95, top = 0.95)
            #plt.imshow(arr_1, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
            plt.imshow(arr_1, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final+'/image_'+ccds[i_ccd]+'.png')
            plt.close()
            # get variance array
            arr_var_1 = calib_1.getVariance().getArray()#calib.variance.array
            plt.figure()
            plt.imshow(arr_var_1, origin='lower', vmin=0.5, vmax=2, cmap = 'hot')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_variance+'/image_'+ccds[i_ccd]+'.png')
            plt.close()
            # variance difference collection_1 - collection_2
            calib_2 = butler_2.get(calibType, instrument='LATISS', detector=detector_id)
            arr_var_2 = calib_2.getVariance().getArray()
            arr = arr_var_1 - arr_var_2
            plt.figure()
            plt.imshow(arr, origin='lower', vmin=-0.5, vmax=0.5, cmap = 'hot')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_difference+'/image_'+ccds[i_ccd]+'.png')
            plt.close()
            #get per amplifier results
            detector = calib_1.getDetector()
            #get size
            amplifier = detector['C00']
            sub_im1 = calib_1.getMaskedImage()[amplifier.getBBox()]
            arr_amp1 = sub_im1.getImage().getArray()
            im_x_size = arr_amp1.shape[1]
            im_y_size = arr_amp1.shape[0]
            var_total =  np.zeros(16)
            mean_total = np.zeros(16)
            mean_corner = np.zeros(16)
            var_corner = np.zeros(16)
            var_line = np.zeros((16,im_y_size))
            mean_line = np.zeros((16,im_y_size))
            var_column =  np.zeros((16,im_x_size))
            mean_column = np.zeros((16,im_x_size))
            variance_mean = np.zeros(16)
            variance_rms = np.zeros(16)
            variance_q95 = np.zeros(16)
            variance_skewness = np.zeros(16)
            amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
            detector = calib_1.getDetector()
            up_line = im_y_size - 100
            up_column = im_x_size - 100
            #print(up_line)
            #print(up_column)
            for i_amp in range(len(amps)):
                amplifier = detector[amps[i_amp]]
                sub_im1 = calib_1.getMaskedImage()[amplifier.getBBox()]
                arr_amp_1 = sub_im1.getImage().getArray()
                plt.figure()
                #print(arr_amp_1.shape)
                arr_amp_1_small = arr_amp_1[0:100,0:100]
                arr_amp_1_corner = arr_amp_1[0:20,0:40]
                if(i_amp>7):
                    arr_amp_1_small = arr_amp_1[up_line:im_y_size,up_column:im_x_size]
                    arr_amp_1_corner = arr_amp_1[up_line+80:im_y_size,up_column+60:im_x_size]
                #title = rafts[i_raft]+' '+ccds[i_ccd]+' '+amps[i_amp]+': corner'    
                #plt.imshow(arr_amp_1_small, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
                #plt.colorbar()
                #plt.title(title)
                #plt.savefig(outpath_final+amps[i_amp]+'.png')
                var_total[i_amp] = np.var(arr_amp_1)
                mean_total[i_amp] = np.mean(arr_amp_1)
                var_line[i_amp,:] = np.var(arr_amp_1,axis=1)
                mean_line[i_amp,:] = np.mean(arr_amp_1,axis=1)
                var_column[i_amp,:] = np.var(arr_amp_1,axis=0)
                mean_column[i_amp,:] = np.mean(arr_amp_1,axis=0)
                #check yellow corner
                mean_corner[i_amp] = np.mean(arr_amp_1_corner)
                var_corner[i_amp] = np.var(arr_amp_1_corner)
                #print(amps[i_amp])
                diff = mean_corner[i_amp]-mean_total[i_amp]
                #print(diff)
                #variance_metrics
                arr_var_amp_1 = sub_im1.getVariance().getArray()
                sub_im2 = calib_2.getMaskedImage()[amplifier.getBBox()]
                arr_var_amp_2 = sub_im2.getVariance().getArray()
                values_1 = arr_var_amp_1.flatten()
                values_2 = arr_var_amp_2.flatten()
                variance_q95[i_amp]=np.percentile(values_1,95)
                variance_mean[i_amp] = np.nanmean(values_1)
                variance_rms[i_amp] = np.nanstd(values_1)
                variance_skewness[i_amp] = skew(values_1,nan_policy='omit')
                print('------')
                print(mean_total[i_amp])
                print(variance_mean[i_amp])
                print(values_1)
            #write results    
            t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, mean_corner, var_corner, variance_q95, variance_mean, variance_rms, variance_skewness], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column', 'mean_corner', 'var_corner', 'variance_q95', 'variance_mean', 'variance_rms', 'variance_skewness'), meta={'name': 'Variances'})
            #print(t_variance)
            t_variance.write(output_data + '/results/Variances_' + rafts[i_raft] + '_'+ ccds[i_ccd] + '.fits', overwrite=True)

            #print('free memory')
            #print(dir())
            #del arr
            #del arr_amp0
            #gc.collect()

############################images            
if plot_images:
    #get the list of exposures
    datasetRefs=list(registry_1.queryDatasets(datasetType=datasetType, instrument='LATISS', collections=collection_1))
    #bug with: where="detector.full_name=R14_S22"
    #exposures = []
    #for i in range(len(datasetRefs)):
    #    exposure = datasetRefs[i].dataId['exposure']
    #    if(exposure not in exposures):
    #        exposures.append(exposure)
    #print(exposures)
    #sys.exit()
    #by hand
    #exposures = ['3023062100512','3023062100513','3023062100514']
    #exposures = ['3023062100285','3023062100286','3023062100287']
    exposures = ['2024030800199']
    for i_exp in range(len(exposures)):
        #create exposure directory
        str_exp = str(exposures[i_exp])
        outpath_exp = output_data+'exposures/exposure_'+ str_exp +'/'
        os.makedirs(outpath_exp+'/results',exist_ok=True)
        #print(outpath_exp)
        #loop over all CCDs
        for i_raft in range(len(rafts)):
            #create raft directory
            outpath_final = outpath_exp+rafts[i_raft]+''
            os.makedirs(outpath_final,exist_ok=True)
            #shutil.copy2('index_files/index_raft_plots.html',outpath_final+'/index.html')
        
            for i_ccd in range(len(ccds)):
                #detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=3021121200153'
                detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=' + str_exp
                #without seq_num
                #detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\''
                print(detector)
                datasetRefs=list(registry_1.queryDatasets(datasetType=datasetType, instrument='LATISS', collections=collection_1, where=f"detector.full_name={detector}"))
                # get image
                exposure = datasetRefs[0].dataId['exposure']
                raw = butler_1.get(datasetRefs[0])
                #title = 'Raft '+rafts[i_raft]+' CCD '+ccds[i_ccd]+' exposure '+str(exposure)
                title = rafts[i_raft]+' '+ccds[i_ccd]+' exp='+str(exposure)
                plt.figure(figsize=(20, 18))
                # get full array
                arr = raw.getImage().getArray()
                print(arr.shape)
                print(arr)
                #plt.imshow(arr, origin='lower', vmin=-1, vmax=1, cmap = 'hot')
                #plt.imshow(arr, origin='lower', vmin=-2, vmax=2, cmap = 'hot')
                plt.imshow(arr, origin='lower', vmin=73500, vmax=74300, cmap='gist_gray')
                #display the image directly
                #afw_display.scale('linear', 'zscale')
                plt.colorbar()
                plt.title(title)
                plt.savefig(outpath_final+'/image_'+ccds[i_ccd]+'.png')
                plt.close()
                
                #get per amplifier results
                detector = raw.getDetector()
                #get size
                amplifier = detector['C00']
                sub_im0 = raw.getMaskedImage()[amplifier.getBBox()]
                arr_amp0 = sub_im0.getImage().getArray()
                im_x_size = arr_amp0.shape[1]
                im_y_size = arr_amp0.shape[0]
                var_total =  np.zeros(16)
                mean_total = np.zeros(16)
                mean_corner = np.zeros(16)
                var_corner = np.zeros(16)
                var_line = np.zeros((16,im_y_size))
                mean_line = np.zeros((16,im_y_size))
                var_column =  np.zeros((16,im_x_size))
                mean_column = np.zeros((16,im_x_size))
                #20*20 pixels
                mean_local_20 = np.zeros((16,100,25))
                #50*50 pixels
                mean_local_50 = np.zeros((16,40,10))
                amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
                detector = raw.getDetector()
                for i_amp in range(len(amps)):
                    amplifier = detector[amps[i_amp]]
                    sub_im0 = raw.getMaskedImage()[amplifier.getBBox()]
                    arr_amp = sub_im0.getImage().getArray()
                    #print(arr_amp0.shape)
                    plt.figure(figsize=(5,20))
                    plt.imshow(arr_amp, origin='lower', vmin=73500, vmax=74300, cmap = 'gist_gray')
                    plt.colorbar()
                    plt.savefig(outpath_final+'/raw_amp_'+str(i_amp)+'.png')
                    var_total[i_amp] = np.var(arr_amp)
                    mean_total[i_amp] = np.mean(arr_amp)
                    var_line[i_amp,:] = np.var(arr_amp,axis=1)
                    mean_line[i_amp,:] = np.mean(arr_amp,axis=1)
                    var_column[i_amp,:] = np.var(arr_amp,axis=0)
                    mean_column[i_amp,:] = np.mean(arr_amp,axis=0)
                    #check yellow corner
                    mean_corner[i_amp] = np.mean(arr_amp[0:19,0:19])
                    var_corner[i_amp] = np.var(arr_amp[0:19,0:19])
                    #print(mean_corner[i_amp])
                    #print(mean_total[i_amp])
                    #if(i_amp==15):
                    #    print(mean_column[i_amp])
                    #20*20 pixels
                    for j in range(100):
                        xmin=20*j
                        xmax=20*(j+1)
                        for k in range(25):
                            ymin=20*k
                            ymax=20*(k+1)
                            mean_local_20[i_amp,j,k]=np.mean(arr_amp[xmin:xmax,ymin:ymax])
                    #50*50 pixels
                    for j in range(40):
                        xmin=50*j
                        xmax=50*(j+1)
                        for k in range(10):
                            ymin=50*k
                            ymax=50*(k+1)
                            mean_local_50[i_amp,j,k]=np.mean(arr_amp[xmin:xmax,ymin:ymax])
                            
                    del arr_amp
                #np.set_printoptions(threshold=np.inf)
                #print(mean_local_20)
                #print('===50')
                #print(mean_local_50)
                
                #write results
                t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, mean_local_20, mean_local_50], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column', 'mean_local_20', 'mean_local_50'), meta={'name': 'Variances'})
                #print(t_variance)
                t_variance.write(outpath_exp + '/results/Variances_' + rafts[i_raft] + '_'+ ccds[i_ccd] + '.fits', overwrite=True)

                #print('free memory')
                #print(dir())
                del arr
                del arr_amp0
                gc.collect()
                
print('DONE')
sys.exit()
