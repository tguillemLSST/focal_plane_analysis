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
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
collection_1 = str(sys.argv[4])
collection_2 = str(sys.argv[5])
output_data=str(sys.argv[6])

os.makedirs(output_data,exist_ok=True)

#butler access
#repo='/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/'
repo = '/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/'
print('=========butler: ' + repo)

###calibrations
plot_calibrations = True
calibType='bias'
#calibType='dark'
###user collections
plot_images = False
datasetType='cpBiasProc'
print('=========collections: ' + collection_1 + ' and ' + collection_2)
butler_1 = dafButler.Butler(repo,collections=collection_1)
registry_1 = butler_1.registry
butler_2 = dafButler.Butler(repo,collections=collection_2)
registry_2 = butler_2.registry

#for a set of batch jobs
rafts=[str_raft]
ccds=[str_all_sensors]

ref_dir='/sps/lsst/users/tguillem/web/stack/yellow_corner_fix1/MB/focal_plane'
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
            #dataId = {'exposure': 2021121200150, 'detector': detector_id, "instrument": 'LSSTCam'}
            #calib_1 = butler_1.get(calibType,dataId=dataId)
            calib_1 = butler_1.get(calibType, instrument='LSSTCam', detector=detector_id)
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
            arr_1_var = calib_1.getVariance().getArray()#calib.variance.array
            print('variance--')
            print(arr_1_var)
            plt.figure()
            plt.imshow(arr_1_var, origin='lower', vmin=0.2, vmax=1.0, cmap = 'hot')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_variance+'/image_'+ccds[i_ccd]+'.png')
            plt.close()
            # difference collection_2 - collection_1
            calib_2 = butler_2.get('bias', instrument='LSSTCam', detector=detector_id)
            arr_2 = calib_2.getImage().getArray()
            arr = arr_2 - arr_1
            plt.figure()
            plt.imshow(arr, origin='lower', vmin=-2, vmax=2, cmap = 'hot')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_difference+'/image_'+ccds[i_ccd]+'.png')
            plt.close()
            #get per amplifier results
            detector = calib_1.getDetector()
            #get size
            amplifier = detector['C00']
            sub_im0 = calib_1.getMaskedImage()[amplifier.getBBox()]
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
                sub_im0 = calib_1.getMaskedImage()[amplifier.getBBox()]
                arr_amp = sub_im0.getImage().getArray()
                plt.figure()
                print(arr_amp.shape)
                arr_amp_small = arr_amp[0:100,0:100]
                arr_amp_corner = arr_amp[0:20,0:40]
                if(i_amp>7):
                    arr_amp_small = arr_amp[up_line:im_y_size,up_column:im_x_size]
                    arr_amp_corner = arr_amp[up_line+80:im_y_size,up_column+60:im_x_size]
                title = rafts[i_raft]+' '+ccds[i_ccd]+' '+amps[i_amp]+': corner'    
                plt.imshow(arr_amp_small, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
                plt.colorbar()
                plt.title(title)
                plt.savefig(outpath_final+amps[i_amp]+'.png')
                var_total[i_amp] = np.var(arr_amp)
                mean_total[i_amp] = np.mean(arr_amp)
                var_line[i_amp,:] = np.var(arr_amp,axis=1)
                mean_line[i_amp,:] = np.mean(arr_amp,axis=1)
                var_column[i_amp,:] = np.var(arr_amp,axis=0)
                mean_column[i_amp,:] = np.mean(arr_amp,axis=0)
                #check yellow corner
                mean_corner[i_amp] = np.mean(arr_amp_corner)
                var_corner[i_amp] = np.var(arr_amp_corner)
                print(amps[i_amp])
                diff = mean_corner[i_amp]-mean_total[i_amp]
                print(diff)
                #variance_metrics
                arr_var_amp = sub_im0.getVariance().getArray()
                values = arr_var_amp.flatten()
                variance_q95[i_amp]=np.percentile(values,95)
                variance_mean[i_amp] = np.mean(values)
                variance_rms[i_amp] = np.std(values)
                variance_skewness[i_amp] = skew(values)
                
            #write results    
            t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, mean_corner, var_corner, variance_q95, variance_mean, variance_rms, variance_skewness], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column', 'mean_corner', 'var_corner', 'variance_q95', 'variance_mean', 'variance_rms', 'variance_skewness'), meta={'name': 'Variances'})
            print(t_variance)
            t_variance.write(output_data + '/results/Variances_' + rafts[i_raft] + '_'+ ccds[i_ccd] + '.fits', overwrite=True)

            #print('free memory')
            #print(dir())
            #del arr
            #del arr_amp0
            #gc.collect()

############################images            
if plot_images:
    #loop over exposures
    #exposures = ['3021121200137','3021121200138','3021121200139','3021121200140','3021121200141','3021121200142','3021121200143','3021121200144','3021121200145','3021121200146','3021121200147','3021121200148','3021121200149','3021121200150','3021121200151','3021121200152','3021121200153','3021121200154','3021121200155','3021121200156','3021121200157','3021121200158','3021121200159','3021121200160','3021121200161','3021121200162','3021121200163','3021121200164','3021121200165','3021121200166','3021121200167','3021121200168','3021121200169','3021121200170','3021121200171','3021121200172','3021121200173','3021121200174','3021121200175','3021121200176']
    #exposures = ['3021121200146','3021121200158','3021121200174']
    #exposures = ['3021121200145','3021121200137','3021121200141','3021121200150','3021121200151','3021121200156','3021121200152','3021121200148','3021121200139','3021121200138','3021121200154','3021121200149','3021121200153','3021121200147','3021121200140','3021121200146','3021121200143','3021121200155','3021121200142','3021121200144']
    ###Run 6
    #exposures = ['3023061800079']
    exposures = ['3023061800147']
    #Run 13372
    #exposures = [ '3023061800079','3023061800100','3023061800058','3023061800025','3023061800106','3023061800024','3023061800049','3023061800023','3023061800021','3023061800073','3023061800052','3023061800027','3023061800070','3023061800076','3023061800064','3023061800091','3023061800082','3023061800094','3023061800103','3023061800030']
    for i_exp in range(len(exposures)):
        #create exposure directory
        outpath_exp = output_data+'exposure_'+exposures[i_exp]+'/'
        os.makedirs(outpath_exp+'/results',exist_ok=True)
        #print(outpath_exp)
        #loop over all CCDs
        for i_raft in range(len(rafts)):
            #create raft directory
            outpath_final = outpath_exp+rafts[i_raft]+''
            os.makedirs(outpath_final,exist_ok=True)
            shutil.copy2('index_files/index_raft_plots.html',outpath_final+'/index.html')
        
            for i_ccd in range(len(ccds)):
                #HACK: pass the seq_num here
                #Run 13161 bias: exposure.id=3021121200146
                #Run 13161 dark: exposure.id=3021121200160
                #Run 13161: bias_bias_000 <=> 000137
                #detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=3021121200153'
                detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=' + exposures[i_exp]
                #without seq_num
                #detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\''
                print(detector)
                datasetRefs=list(registry_1.queryDatasets(datasetType=datasetType, instrument='LSSTCam', collections=collection_1, where=f"detector.full_name={detector}"))
                print(datasetRefs)
            
                # get image
                exposure = datasetRefs[0].dataId['exposure']
                raw = butler_1.get(datasetRefs[0])
                #title = 'Raft '+rafts[i_raft]+' CCD '+ccds[i_ccd]+' exposure '+str(exposure)
                title = rafts[i_raft]+' '+ccds[i_ccd]+' exp='+str(exposure)
                plt.figure()
                # get full array
                arr = raw.getImage().getArray()
                #plt.imshow(arr, origin='lower', vmin=-1, vmax=1, cmap = 'hot')
                plt.imshow(arr, origin='lower', vmin=25, vmax=100, cmap = 'hot')
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
                amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
                detector = raw.getDetector()
                for i_amp in range(len(amps)):
                    amplifier = detector[amps[i_amp]]
                    sub_im0 = raw.getMaskedImage()[amplifier.getBBox()]
                    arr_amp = sub_im0.getImage().getArray()
                    #plt.figure()
                    #plt.imshow(arr_amp, origin='lower', vmin=-5, vmax=5, cmap = 'hot')
                    #plt.colorbar()
                    #plt.savefig(outpath_final+'/raw_amp.png')
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
                    if(i_amp==15):
                        print(mean_column[i_amp])
                    del arr_amp
                #write results    
                t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column'), meta={'name': 'Variances'})
                #print(t_variance)
                t_variance.write(outpath_exp + '/results/Variances_' + rafts[i_raft] + '_'+ ccds[i_ccd] + '.fits', overwrite=True)

                #print('free memory')
                #print(dir())
                del arr
                del arr_amp0
                gc.collect()
                
print('DONE')
sys.exit()
