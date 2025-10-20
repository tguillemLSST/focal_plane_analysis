###import eo_pipe and other stack utils
import pandas as pd
import lsst.daf.butler as daf_butler
from lsst.obs.lsst import LsstCam, LsstTS8
import lsst.eo.pipe as eo_pipe

###import various packages
from collections import defaultdict
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

#prepare per amp data table
amp_data = defaultdict(lambda: defaultdict(dict))
amp_data_relative = defaultdict(lambda: defaultdict(dict))
rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
#e2v
#rafts=['R11' ,'R12' ,'R13' ,'R14' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']

#QE from CBP campaigns
inpath_base = '/sdf/home/a/amouroux/public_html/cbp/lsstcam/telescope_response/20250518/cbp_qe_none/results/'
#inpath_base = '/sdf/home/a/amouroux/public_html/cbp/lsstcam/telescope_response/20250516/cbp_none_all_1/results/'
#inpath_base = '/sdf/home/a/amouroux/public_html/cbp/lsstcam/telescope_response/20250516/cbp_none_all_2/results/'

outpath = '/sdf/home/t/tguillem/public_html/cbp/mosaics/'
os.makedirs(outpath,exist_ok=True)
os.makedirs(outpath+'QE/',exist_ok=True)
os.makedirs(outpath+'QE_relative/',exist_ok=True)

#get a reference
variance_reference=Table.read(inpath_base+'R22_S11.fits')
#for k in range(0,len(variance_reference)):
#            print(str(k) + ' ' + str(variance_reference['wavelength'][k]))
            
#index_lambda=9
#loop over wavelegths
for index_lambda in range(0,len(variance_reference)):
            for i in range(len(rafts)):
	                for j in range(len(ccds)):
                                    detector_full_name = rafts[i] + '_' +ccds[j]
                                    if(detector_full_name=='R30_S10' or detector_full_name=='R30_S11' or detector_full_name=='R30_S12'):
                                                continue
                                    inpath = inpath_base + rafts[i] + '_' + ccds[j] + '.fits'
                                    variance = Table.read(inpath)
                                    for m in range(16) :	
                                                amp_data['mean'][rafts[i] + '_' + ccds[j]][amps[m]]=variance['telescope_response'][index_lambda]
                                                amp_data_relative['mean'][rafts[i] + '_' + ccds[j]][amps[m]]=variance['telescope_response'][index_lambda]/variance_reference['telescope_response'][index_lambda]
                                                
            #mosaic plot
            #print(amp_data['mean'])
            my_z_range = [0.4,0.72]
            my_title = 'QE ' + str(variance_reference['wavelength'][index_lambda]) + ' nm'
            plt.figure(figsize=(9,9))
            ax = plt.gca()
            eo_pipe.plotting.plot_focal_plane(ax, amp_data['mean'], title = my_title, z_range=my_z_range)
            plt.savefig(outpath+'QE/mosaic_'+str(variance_reference['wavelength'][index_lambda])+'.png', bbox_inches='tight')
            plt.close()
            
            #relative QE
            my_z_range = [0.8,1.05]
            my_title = 'Relative QE ' + str(variance_reference['wavelength'][index_lambda]) + ' nm'
            plt.figure(figsize=(9,9))
            ax = plt.gca()
            eo_pipe.plotting.plot_focal_plane(ax, amp_data_relative['mean'], title = my_title, z_range=my_z_range)
            plt.savefig(outpath+'QE_relative/mosaic_'+str(variance_reference['wavelength'][index_lambda])+'.png', bbox_inches='tight')
            plt.close()
            
print('DONE')
sys.exit()
