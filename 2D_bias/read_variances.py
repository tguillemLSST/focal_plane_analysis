#!/usr/bin/env python
# coding: utf-8

###import
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

outpath = "/sps/lsst/users/tguillem/web/debug/variances/summary/"

if os.path.exists(outpath):
    shutil.rmtree(outpath)
    os.makedirs(outpath)
    print('outpath = ' + outpath)

inpath_base = '/sps/lsst/users/tguillem/web/debug/variances/13162/bias_bias_'
exposures = ['005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019']
#rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
rafts=['R14']
#ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']
ccds=['S22']

variances = []
for i in range(len(exposures)):
    for j in range(len(rafts)):
        for k in range(len(ccds)):
            inpath = inpath_base + exposures[i] + '/' + rafts[j] + '/' + ccds[k] + '/Variances_2D_corr.fits'
            variance = Table.read(inpath)
            #print(variances)
            variances.append(variance)

#print(variances)
print(variances[2]['mean_total_corr_2D'])
