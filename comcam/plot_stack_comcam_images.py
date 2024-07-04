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
#from conversion import *
from lsst.ip.isr import AssembleCcdTask
assembleTask = AssembleCcdTask()

output_data='/sdf/home/t/tguillem/public_html/Operations_Rehearsal/OR4/20240625/'
os.makedirs(output_data,exist_ok=True)

#butler access
#repo_1='/sdf/data/rubin/repo/ops-rehearsal-3-prep/butler.yaml'
#print('=========butler 1: ' + repo_1)
#repo='embargo_or4.yaml'
#print('=========butler : ' + repo)

#collections
#collection_1='LSSTComCamSim/raw/all'#runs/nightlyvalidation/20240403/d_2024_03_29/DM-43612'
collection_1='LSSTComCamSim/raw/all'
datasetType_1_1='raw'
#OR3
#collection_2='LSSTComCamSim/runs/nightlyvalidation/20240403/d_2024_03_29/DM-43612'
#OR4
collection_2='LSSTComCamSim/nightlyValidation'
#collection_2='LSSTComCamSim/prompt/output-2024-06-25'

datasetType_2_1='postISRCCD'
datasetType_2_2='calexp'
datasetType_2_3='calexpBackground'

#print('=========collections: ' + collection_1 + ' and ' + collection_2)
#butler_1 = dafButler.Butler(repo_1)
#registry_1 = butler_1.registry
repo = 'embargo_or4'
butler = dafButler.Butler(repo,writeable=False)
registry = butler.registry
#old
#butler = dafButler.Butler(repo)
#registry = butler.registry

#CCD
rafts=['R22']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
ccds=['S01', 'S02']

#get the list of exposures
#exposures=[]
#datasetRefs=list(registry_1.queryDatasets(datasetType=datasetType, instrument='LSSTComCamSim', collections=collection_1))
#raw_exps = sorted(registry_1.queryDatasets('raw', instrument='LSSTComCamSim',where="exposure.day_obs=20240117",collections=collection_1))
#print(raw_exps)

#by hand
#exposures = ['7024062500131']
exposures = ['7024062500500']
arr_raw_ccd=[]
arr_postISRCCD_ccd=[]
arr_calexp_ccd=[]
arr_calexpBackground_ccd=[]
for i_exp in range(len(exposures)):
    #create exposure directory
    str_exp = str(exposures[i_exp])
    #loop over all CCDs
    for i_raft in range(len(rafts)):
        outpath_exp_raw = output_data+'exposures/raw/exposure_'+ str_exp +'/'
        os.makedirs(outpath_exp_raw+'/results',exist_ok=True)
        outpath_final_raw = outpath_exp_raw+rafts[i_raft]+''
        os.makedirs(outpath_final_raw,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final_raw+'/index.html')

        outpath_exp_postISRCCD = output_data+'exposures/postISRCCD/exposure_'+ str_exp +'/'
        os.makedirs(outpath_exp_postISRCCD+'/results',exist_ok=True)
        outpath_final_postISRCCD = outpath_exp_postISRCCD+rafts[i_raft]+''
        os.makedirs(outpath_final_postISRCCD,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final_postISRCCD+'/index.html')
        
        outpath_exp_calexp = output_data+'exposures/calexp/exposure_'+ str_exp +'/'
        os.makedirs(outpath_exp_calexp+'/results',exist_ok=True)
        outpath_final_calexp = outpath_exp_calexp+rafts[i_raft]+''
        os.makedirs(outpath_final_calexp,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final_calexp+'/index.html')

        outpath_exp_calexpBackground = output_data+'exposures/calexpBackground/exposure_'+ str_exp +'/'
        os.makedirs(outpath_exp_calexpBackground+'/results',exist_ok=True)
        outpath_final_calexpBackground = outpath_exp_calexpBackground+rafts[i_raft]+''
        os.makedirs(outpath_final_calexpBackground,exist_ok=True)
        shutil.copy2('index_files/index_raft_plots.html',outpath_final_calexpBackground+'/index.html')

        os.makedirs(outpath_exp_raw+'/raft_level/',exist_ok=True)
        os.makedirs(outpath_exp_postISRCCD+'/raft_level/',exist_ok=True)
        os.makedirs(outpath_exp_calexp+'/raft_level/',exist_ok=True)
        os.makedirs(outpath_exp_calexpBackground+'/raft_level/',exist_ok=True)
        
        #CCD-level analysis
        for i_ccd in range(len(ccds)):
            detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND exposure.id=' + str_exp
            #print(detector)

            #raw
            datasetRefs=list(registry.queryDatasets(datasetType=datasetType_1_1, instrument='LSSTComCamSim', collections=collection_1, where=f"detector.full_name={detector}"))
            # get image
            exposure = datasetRefs[0].dataId['exposure']
            raw = butler.get(datasetRefs[0])
            title = rafts[i_raft]+' '+ccds[i_ccd]+' exp='+str(exposure)
            # get full array
            arr = raw.getImage().getArray()
            #print(arr.shape)
            #print(arr[300,300])
            raw_science = assembleTask.assembleCcd(raw)
            arr_raw = raw_science.getImage().getArray()
            arr_raw_ccd.append(arr_raw)
            mean_raw=np.mean(arr_raw)
            max_raw=mean_raw+100
            min_raw=mean_raw-100
            plt.figure(figsize=(20, 18))
            plt.imshow(arr_raw, origin='lower', vmin=min_raw, vmax=max_raw, cmap='gist_gray')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_raw+'/image_'+ccds[i_ccd]+'.png')
            plt.close()
            #sys.exit()

            #postISRCCD
            datasetRefs=list(registry.queryDatasets(datasetType=datasetType_2_1, instrument='LSSTComCamSim', collections=collection_2, where=f"detector.full_name={detector}"))
            postISRCCD = butler.get(datasetRefs[0])
            #print(datasetRefs)
            # get full array
            arr_postISRCCD = postISRCCD.getImage().getArray()
            arr_postISRCCD_ccd.append(arr_postISRCCD)
            mean_postISRCCD=np.mean(arr_postISRCCD)
            max_postISRCCD=mean_postISRCCD+50
            min_postISRCCD=mean_postISRCCD-50
            plt.figure(figsize=(20, 18))
            plt.imshow(arr_postISRCCD, origin='lower', vmin=min_postISRCCD, vmax=max_postISRCCD, cmap='gist_gray')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_postISRCCD+'/image_'+ccds[i_ccd]+'.png')
            plt.close()

            detector = '\'' + rafts[i_raft] + '_' +ccds[i_ccd] + '\' AND visit.id=' + str_exp
            #calexp
            datasetRefs=list(registry.queryDatasets(datasetType=datasetType_2_2, instrument='LSSTComCamSim', collections=collection_2, where=f"detector.full_name={detector}"))
            calexp = butler.get(datasetRefs[0])
            #print(datasetRefs)
            # get full array
            arr_calexp = calexp.getImage().getArray()
            arr_calexp_ccd.append(arr_calexp)
            mean_calexp=np.mean(arr_calexp)
            max_calexp=mean_calexp+50
            min_calexp=mean_calexp-50
            plt.figure(figsize=(20, 18))
            plt.imshow(arr_calexp, origin='lower', vmin=min_calexp, vmax=max_calexp, cmap='gist_gray')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_calexp+'/image_'+ccds[i_ccd]+'.png')
            plt.close()

            #calexpBackground
            datasetRefs=list(registry.queryDatasets(datasetType=datasetType_2_3, instrument='LSSTComCamSim', collections=collection_2, where=f"detector.full_name={detector}"))
            calexpBackground = butler.get(datasetRefs[0])
            #print(datasetRefs)
            # get full array
            arr_calexpBackground = calexpBackground.getImage().getArray()
            arr_calexpBackground_ccd.append(arr_calexpBackground)
            mean_calexpBackground=np.mean(arr_calexpBackground)
            max_calexpBackground=mean_calexpBackground+10
            min_calexpBackground=mean_calexpBackground-10
            #max_calexpBackground=1220
            #min_calexpBackground=1180
            plt.figure(figsize=(20, 18))
            plt.imshow(arr_calexpBackground, origin='lower', vmin=min_calexpBackground, vmax=max_calexpBackground, cmap='gist_gray')
            plt.colorbar()
            plt.title(title)
            plt.savefig(outpath_final_calexpBackground+'/image_'+ccds[i_ccd]+'.png')
            plt.close()
                        
            #continue

            ##get per amplifier results
            #detector = raw.getDetector()
            ##get size
            #amplifier = detector['C00']
            #sub_im0 = raw.getMaskedImage()[amplifier.getBBox()]
            #arr_amp0 = sub_im0.getImage().getArray()
            #im_x_size = arr_amp0.shape[1]
            #im_y_size = arr_amp0.shape[0]
            #var_total =  np.zeros(16)
            #mean_total = np.zeros(16)
            #mean_corner = np.zeros(16)
            #var_corner = np.zeros(16)
            #var_line = np.zeros((16,im_y_size))
            #mean_line = np.zeros((16,im_y_size))
            #var_column =  np.zeros((16,im_x_size))
            #mean_column = np.zeros((16,im_x_size))
            ##20*20 pixels
            #mean_local_20 = np.zeros((16,100,25))
            ##50*50 pixels
            #mean_local_50 = np.zeros((16,40,10))
            #amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']
            #detector = raw.getDetector()
            #for i_amp in range(len(amps)):
            #    amplifier = detector[amps[i_amp]]
            #    sub_im0 = raw.getMaskedImage()[amplifier.getBBox()]
            #    arr_amp = sub_im0.getImage().getArray()
            #    #print(arr_amp0.shape)
            #    var_total[i_amp] = np.var(arr_amp)
            #    mean_total[i_amp] = np.mean(arr_amp)
            #    var_line[i_amp,:] = np.var(arr_amp,axis=1)
            #    mean_line[i_amp,:] = np.mean(arr_amp,axis=1)
            #    var_column[i_amp,:] = np.var(arr_amp,axis=0)
            #    mean_column[i_amp,:] = np.mean(arr_amp,axis=0)
            #    #check yellow corner
            #    mean_corner[i_amp] = np.mean(arr_amp[0:19,0:19])
            #    var_corner[i_amp] = np.var(arr_amp[0:19,0:19])
            #    #print(mean_corner[i_amp])
            #    #print(mean_total[i_amp])
            #    #if(i_amp==15):
            #    #    print(mean_column[i_amp])
            #    #20*20 pixels
            #    for j in range(100):
            #        xmin=20*j
            #        xmax=20*(j+1)
            #        for k in range(25):
            #            ymin=20*k
            #            ymax=20*(k+1)
            #            mean_local_20[i_amp,j,k]=np.mean(arr_amp[xmin:xmax,ymin:ymax])
            #    #50*50 pixels
            #    for j in range(40):
            #        xmin=50*j
            #        xmax=50*(j+1)
            #        for k in range(10):
            #            ymin=50*k
            #            ymax=50*(k+1)
            #            mean_local_50[i_amp,j,k]=np.mean(arr_amp[xmin:xmax,ymin:ymax])
            #                
            #    del arr_amp
            #    #np.set_printoptions(threshold=np.inf)
            #    #print(mean_local_20)
            #    #print('===50')
            #    #print(mean_local_50)
                
            ##write results
            #t_variance = Table([var_total, mean_total, var_line, mean_line, var_column, mean_column, mean_local_20, mean_local_50], names=('var_total', 'mean_total', 'var_line', 'mean_line', 'var_column', 'mean_column', 'mean_local_20', 'mean_local_50'), meta={'name': 'Variances'})
            ##print(t_variance)
            #t_variance.write(outpath_exp + '/results/Variances_' + rafts[i_raft] + '_'+ ccds[i_ccd] + '.fits', overwrite=True)

            #print('free memory')
            #print(dir())
            del arr
            #del arr_amp0
            gc.collect()

        #sys.exit()    
        #Raft-level analysis            
        title = rafts[i_raft]+' exp='+str(exposure)
        #raw
        arr_raw_raft_S0=np.concatenate((arr_raw_ccd[0],arr_raw_ccd[1],arr_raw_ccd[2]),axis=1)
        arr_raw_raft_S1=np.concatenate((arr_raw_ccd[3],arr_raw_ccd[4],arr_raw_ccd[5]),axis=1)
        arr_raw_raft_S2=np.concatenate((arr_raw_ccd[6],arr_raw_ccd[7],arr_raw_ccd[8]),axis=1)
        arr_raw_raft=np.concatenate((arr_raw_raft_S0,arr_raw_raft_S1,arr_raw_raft_S2),axis=0)
        #print(arr_raw_raft)
        plt.figure(figsize=(20, 18))
        mean_raw=np.mean(arr_raw_raft)
        max_raw=np.max(arr_raw_raft)
        min_raw=np.min(arr_raw_raft)
        max_raw=mean_raw+100
        min_raw=min_raw-100
        plt.imshow(arr_raw_raft, origin='lower', vmin=min_raw, vmax=max_raw, cmap='gist_gray')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_exp_raw+'/raft_level/'+'image_raft.png')
        plt.close()

        #1D histogram
        histo_arr_raw_raft = arr_raw_raft.flat
        plt.figure()
        plt.hist(histo_arr_raw_raft, range=[min_raw,max_raw], bins=1000, label='ADU', histtype='step', color = 'black')
        #plt.ylim([0.8, 1.1])
        plt.xlabel('ADU')
        plt.ylabel('pixels')
        plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
        plt.legend()
        plt.savefig(outpath_exp_raw+'/raft_level/'+'histo_raft.png')

        #postISRCCD
        arr_postISRCCD_raft_S0=np.concatenate((arr_postISRCCD_ccd[0],arr_postISRCCD_ccd[1],arr_postISRCCD_ccd[2]),axis=1)
        arr_postISRCCD_raft_S1=np.concatenate((arr_postISRCCD_ccd[3],arr_postISRCCD_ccd[4],arr_postISRCCD_ccd[5]),axis=1)
        arr_postISRCCD_raft_S2=np.concatenate((arr_postISRCCD_ccd[6],arr_postISRCCD_ccd[7],arr_postISRCCD_ccd[8]),axis=1)
        arr_postISRCCD_raft=np.concatenate((arr_postISRCCD_raft_S0,arr_postISRCCD_raft_S1,arr_postISRCCD_raft_S2),axis=0)
        #print(arr_postISRCCD_raft)
        plt.figure(figsize=(20, 18))
        mean_postISRCCD=np.mean(arr_postISRCCD_raft)
        max_postISRCCD=np.max(arr_postISRCCD_raft)
        min_postISRCCD=np.min(arr_postISRCCD_raft)
        max_postISRCCD=mean_postISRCCD+100
        min_postISRCCD=mean_postISRCCD-100
        plt.imshow(arr_postISRCCD_raft, origin='lower', vmin=min_postISRCCD, vmax=max_postISRCCD, cmap='gist_gray')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_exp_postISRCCD+'/raft_level/'+'image_raft.png')
        plt.close()

        #1D histogram
        histo_arr_postISRCCD_raft = arr_postISRCCD_raft.flat
        plt.figure()
        plt.hist(histo_arr_postISRCCD_raft, range=[min_postISRCCD,max_postISRCCD], bins=1000, label='ADU', histtype='step', color = 'black')
        #plt.ylim([0.8, 1.1])
        plt.xlabel('ADU')
        plt.ylabel('pixels')
        plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
        plt.legend()
        plt.savefig(outpath_exp_postISRCCD+'/raft_level/'+'histo_raft.png')

        #calexp
        arr_calexp_raft_S0=np.concatenate((arr_calexp_ccd[0],arr_calexp_ccd[1],arr_calexp_ccd[2]),axis=1)
        arr_calexp_raft_S1=np.concatenate((arr_calexp_ccd[3],arr_calexp_ccd[4],arr_calexp_ccd[5]),axis=1)
        arr_calexp_raft_S2=np.concatenate((arr_calexp_ccd[6],arr_calexp_ccd[7],arr_calexp_ccd[8]),axis=1)
        arr_calexp_raft=np.concatenate((arr_calexp_raft_S0,arr_calexp_raft_S1,arr_calexp_raft_S2),axis=0)
        #print(arr_calexp_raft)
        plt.figure(figsize=(20, 18))
        mean_calexp=np.mean(arr_calexp_raft)
        max_calexp=np.max(arr_calexp_raft)
        min_calexp=np.min(arr_calexp_raft)
        max_calexp=mean_calexp+100
        min_calexp=mean_calexp-100
        plt.imshow(arr_calexp_raft, origin='lower', vmin=min_calexp, vmax=max_calexp, cmap='gist_gray')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_exp_calexp+'/raft_level/'+'image_raft.png')
        plt.close()

        #1D histogram
        histo_arr_calexp_raft = arr_calexp_raft.flat
        plt.figure()
        plt.hist(histo_arr_calexp_raft, range=[min_calexp,max_calexp], bins=1000, label='ADU', histtype='step', color = 'black')
        #plt.ylim([0.8, 1.1])
        plt.xlabel('ADU')
        plt.ylabel('pixels')
        plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
        plt.legend()
        plt.savefig(outpath_exp_calexp+'/raft_level/'+'histo_raft.png')

        #calexpBackground
        arr_calexpBackground_raft_S0=np.concatenate((arr_calexpBackground_ccd[0],arr_calexpBackground_ccd[1],arr_calexpBackground_ccd[2]),axis=1)
        arr_calexpBackground_raft_S1=np.concatenate((arr_calexpBackground_ccd[3],arr_calexpBackground_ccd[4],arr_calexpBackground_ccd[5]),axis=1)
        arr_calexpBackground_raft_S2=np.concatenate((arr_calexpBackground_ccd[6],arr_calexpBackground_ccd[7],arr_calexpBackground_ccd[8]),axis=1)
        arr_calexpBackground_raft=np.concatenate((arr_calexpBackground_raft_S0,arr_calexpBackground_raft_S1,arr_calexpBackground_raft_S2),axis=0)
        #print(arr_calexpBackground_raft)
        plt.figure(figsize=(20, 18))
        mean_calexpBackground=np.mean(arr_calexpBackground_raft)
        max_calexpBackground=np.max(arr_calexpBackground_raft)
        min_calexpBackground=np.min(arr_calexpBackground_raft)
        plt.imshow(arr_calexpBackground_raft, origin='lower', vmin=min_calexpBackground, vmax=max_calexpBackground, cmap='gist_gray')
        plt.colorbar()
        plt.title(title)
        plt.savefig(outpath_exp_calexpBackground+'/raft_level/'+'image_raft.png')
        plt.close()

        #1D histogram
        histo_arr_calexpBackground_raft = arr_calexpBackground_raft.flat
        plt.figure()
        plt.hist(histo_arr_calexpBackground_raft, range=[min_calexpBackground,max_calexpBackground], bins=1000, label='ADU', histtype='step', color = 'black')
        #plt.ylim([0.8, 1.1])
        plt.xlabel('ADU')
        plt.ylabel('pixels')
        plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
        plt.legend()
        plt.savefig(outpath_exp_calexpBackground+'/raft_level/'+'histo_raft.png')

        
print('DONE')
sys.exit()
