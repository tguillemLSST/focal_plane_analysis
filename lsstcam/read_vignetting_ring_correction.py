#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: apply vignetting ring correction arrays
##########

# system imports
import time
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import lsst.daf.butler as dafButler
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from conversion import *
from scipy import stats
from scipy.interpolate import make_interp_spline
from scipy.interpolate import splrep, splev

print('Configuration arguments: ', str(sys.argv))
str_raft = str(sys.argv[1])
str_all_sensors = str(sys.argv[2])

#butler access
repo='dp2_prep'
print('=========butler: ' + repo)
butler = dafButler.Butler(repo)
registry = butler.registry

final_collection='LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage2'
final_dataset='post_isr_image'
exposures = ['2025120200381']

#for a set of batch jobs
rafts=[str_raft]
ccds=[str_all_sensors]

edge_detectors=['R01_S00','R01_S01','R01_S02','R01_S10','R02_S00','R02_S01','R02_S02','R03_S00','R03_S01','R03_S02','R03_S12','R14_S01','R14_S02','R14_S12','R14_S22','R24_S02','R24_S12','R24_S22','R34_S02','R34_S12','R34_S22','R34_S21','R34_S12','R34_S22','R43_S12','R43_S22','R43_S21','R43_S20','R42_S22','R42_S21','R42_S20','R41_S22','R42_S21','R41_S20','R41_S10','R30_S21','R30_S20','R30_S10','R30_S00','R20_S20','R20_S21','R20_S00','R10_S20','R10_S10','R10_S00','R10_S01']

#loop over exposures and detectors
for i_exp in range(len(exposures)):
    str_exp = str(exposures[i_exp])
    for i_raft in range(len(rafts)):
        for i_ccd in range(len(ccds)):
            detector_id=detector(rafts[i_raft],ccds[i_ccd])[1]
            detector_full_name = rafts[i_raft] + '_' +ccds[i_ccd]
            if(detector_full_name=='R30_S12'):
                continue
            if(detector_full_name not in edge_detectors):
                continue

            detector_condition = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=' + str_exp

            datasetRefs=list(registry.queryDatasets(datasetType=final_dataset, instrument='LSSTCam', collections=final_collection, where=f"detector.full_name={detector_condition}"))
            post_isr_image = butler.get(datasetRefs[0])
            arr_post_isr_image=post_isr_image.getMaskedImage().getImage().getArray()
            print(arr_post_isr_image[1300,2800])
            
            ################METHOD 1 : derive a static correction from flattened skyflats
            ##read VRCA
            #with open('/sdf/data/rubin/user/tguillem/shared_data/vignetting_ring_correction_arrays/VRCA_'+detector_full_name+'.npy', 'rb') as f:
            #    arr_final_correction = np.load(f)
            #    
            ##just divide the arr_post_isr_image by the correction array
            #arr_post_isr_image_corrected = np.copy(arr_post_isr_image)
            #arr_post_isr_image_corrected /= arr_final_correction
            #print(arr_post_isr_image_corrected[1300,2800])
            #
            ####bonus part : images before/after
            #plots='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/dp2/tests/vignetting/ring_geometry/validation/'+str_exp+'/'
            #os.makedirs(plots,exist_ok=True)
            ##before correction
            #median_scale= np.median(arr_post_isr_image)
            #vmin_scale=median_scale-5
            #vmax_scale=median_scale+5
            #plt.figure(figsize=(25, 20))
            #plt.imshow(arr_post_isr_image,origin='lower', vmin=vmin_scale, vmax=vmax_scale, cmap = 'grey')
            #plt.colorbar()
            #plt.savefig(plots+'image_'+detector_full_name+'.png', bbox_inches='tight')
            #plt.close()
            ##after correction
            #median_scale= np.median(arr_post_isr_image_corrected)
            #vmin_scale=median_scale-5
            #vmax_scale=median_scale+5
            #plt.figure(figsize=(25, 20))
            #plt.imshow(arr_post_isr_image_corrected,origin='lower', vmin=vmin_scale, vmax=vmax_scale, cmap = 'grey')
            #plt.colorbar()
            #plt.savefig(plots+'image_corrected_'+detector_full_name+'.png', bbox_inches='tight')
            #plt.close()
            
            ################METHOD 2 : derive a dynamic radial correction for each exposure
            #read geometry arrays
            with open('/sdf/data/rubin/user/tguillem/shared_data/geometry/radius_'+detector_full_name+'.npy', 'rb') as f:
                arr_radius2 = np.load(f)
                
            #extra vignetting profile versus radius
            vignetting_radius=[]
            vignetting_value=[]
            n_bins_x=len(arr_post_isr_image[0,:])
            n_bins_y=len(arr_post_isr_image[:,0])
            arr_2nd = np.copy(arr_post_isr_image)
            arr_2nd[:,:]=0
            for i_x in range(n_bins_x):
                for i_y in range(n_bins_y):
                    radius=arr_radius2[i_y,i_x]
                    value=arr_post_isr_image[i_y,i_x]
                    vignetting_radius.append(radius)
                    vignetting_value.append(value)

            #compute median versus radius
            bin_means, bin_edges, binnumber = stats.binned_statistic(vignetting_radius, vignetting_value, statistic='median', range=[(0, 370)], bins=740)
            bin_width = (bin_edges[1] - bin_edges[0])
            bin_centers = bin_edges[1:] - bin_width/2

            #renormalize bin_means (to start at 1 at r=0)
            average_vignetting = []
            for i in range(len(bin_means)):
                if(np.isnan(bin_means[i])):
                    bin_means[i]=0

            for i in range(len(bin_means)):
                if(bin_means[i]>0.8):
                    average_vignetting.append(bin_means[i])
                if(len(average_vignetting)>10):
                    break

            average_vignetting_v = np.mean(average_vignetting[1:10])
            arr_vignetting_value = np.array(vignetting_value)
            arr_vignetting_value = 1/average_vignetting_v*arr_vignetting_value
            bin_means = 1/average_vignetting_v*bin_means
            
            #find first and last points automatically
            index_min=0
            index_max=0
            for i in range(len(bin_means)):
                if(bin_means[i]>0.5):
                    index_min=i
                    break
            for i in range(len(bin_means)):
                j=len(bin_means)-1-i
                if(bin_means[j]>0.5):
                    index_max=j
                    break

            #spline
            bspl_1 = make_interp_spline(bin_centers[index_min:index_max], bin_means[index_min:index_max],k=3)
            xx_min=bin_centers[index_min]
            xx_max=bin_centers[index_max]
            #write array built from spline
            arr_final_correction = np.copy(arr_post_isr_image)
            arr_final_correction[:,:]=0
            for i_x in range(n_bins_x):
                for i_y in range(n_bins_y):
                    radius=arr_radius2[i_y,i_x]
                    value=bspl_1(radius)
                    #be careful with the first/last bins
                    if(radius<xx_min):
                        value=bspl_1(xx_min)
                    if(radius>xx_max):
                        value=bspl_1(xx_max)
                    arr_final_correction[i_y,i_x]=value

            #just divide the arr_post_isr_image by the correction array
            arr_post_isr_image_corrected = np.copy(arr_post_isr_image)
            arr_post_isr_image_corrected /= arr_final_correction
            print(arr_post_isr_image_corrected[1300,2800])
            
            ###bonus part : images before/after
            plots='/sdf/home/t/tguillem/public_html/camera_tma/data_validation/dp2/tests/vignetting/ring_geometry/validation_tmp/'+str_exp+'/'
            os.makedirs(plots,exist_ok=True)
            #before correction
            median_scale= np.median(arr_post_isr_image)
            vmin_scale=median_scale-5
            vmax_scale=median_scale+5
            plt.figure(figsize=(25, 20))
            plt.imshow(arr_post_isr_image,origin='lower', vmin=vmin_scale, vmax=vmax_scale, cmap = 'grey')
            plt.colorbar()
            plt.savefig(plots+'image_'+detector_full_name+'.png', bbox_inches='tight')
            plt.close()
            #after correction
            median_scale= np.median(arr_post_isr_image_corrected)
            vmin_scale=median_scale-5
            vmax_scale=median_scale+5
            plt.figure(figsize=(25, 20))
            plt.imshow(arr_post_isr_image_corrected,origin='lower', vmin=vmin_scale, vmax=vmax_scale, cmap = 'grey')
            plt.colorbar()
            plt.savefig(plots+'image_corrected_'+detector_full_name+'.png', bbox_inches='tight')
            plt.close()
                        
print('DONE')
sys.exit()
