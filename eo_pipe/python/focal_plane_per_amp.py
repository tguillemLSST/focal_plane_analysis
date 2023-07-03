###import eo_pipe and other stack utils
import pandas as pd
import lsst.daf.butler as daf_butler
from lsst.obs.lsst import LsstCam, LsstTS8
import lsst.eo.pipe as eo_pipe

###import
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

#var_total plot
#prepare per amp data table
amp_data = defaultdict(lambda: defaultdict(dict))
rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
amps=['C00','C01','C02','C03','C04','C05','C06','C07','C10','C11','C12','C13','C14','C15','C16','C17']

inpath_base = '/sps/lsst/users/tguillem/web/debug/stack/run_6/13392/results/'
outpath = inpath_base+'mosaics/'
os.makedirs(outpath,exist_ok=True)
print('outpath = ' + outpath)

for i in range(len(rafts)):
	for j in range(len(ccds)):
		inpath = inpath_base + '/Variances_'+ rafts[i] + '_' + ccds[j] + '.fits'
		variance = Table.read(inpath)
		#print(variance)
		for m in range(16) :	
			amp_data['mean'][rafts[i] + '_' + ccds[j]][amps[m]]=variance['variance_mean'][m]
			amp_data['rms'][rafts[i] + '_' + ccds[j]][amps[m]]=variance['variance_rms'][m]
			#amp_data['q95'][rafts[i] + '_' + ccds[j]][amps[m]]=variance['q95'][m]
			#amp_data['skewness'][rafts[i] + '_' + ccds[j]][amps[m]]=variance['skewness'][m]

#print(amp_data['mean'])

#mosaic plot
my_z_range = [-2,10]
my_title = 'mean over 20 biases'
plt.figure(figsize=(9,9))
ax = plt.gca()
eo_pipe.plotting.plot_focal_plane(ax, amp_data['mean'], title = my_title, z_range=my_z_range)
plt.savefig(outpath+'mean.png', bbox_inches='tight')

my_z_range = [0.2,1.0]
my_title = 'variance over 20 biases'
plt.figure(figsize=(9,9))
ax = plt.gca()
eo_pipe.plotting.plot_focal_plane(ax, amp_data['rms'], title = my_title, z_range=my_z_range)
plt.savefig(outpath+'variance.png', bbox_inches='tight')

sys.exit()

###function used
#def plot_focal_plane(ax, amp_data, camera=None, cm=plt.cm.hot,
#		     x_range=None, y_range=None,
#		     z_range=None, use_log10=False, scale_factor='1',
#		     title='', nsigma=4):
