#!/usr/bin/env python
# coding: utf-8

##########
# Adaptation of a DM 2019 bootcamp notebook
# Author: T. Guillemin
# Date: May 2022
# Goal: access CCD images of a bias run using the butler
##########

# system imports
from matplotlib import pylab as plt
import numpy as np
import os
import sys

# LSST stack imports
from lsst.daf.persistence import Butler
import lsst.afw.display as afwDisplay
from lsst.ip.isr import IsrTask
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage

# Access images from the Butler
repo_path='/sps/lsst/groups/FocalPlane/SLAC/run5/butler/all_runs/13144/'
butler = Butler(repo_path)

#check the butler content
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "PRAGMA table_info(raw);"')
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "select * from raw limit 10;"')
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "select distinct testType from raw;"')
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "select distinct imageType from raw;"')

#select visits for BIAS images
visits_bias=butler.queryMetadata('raw', ['visit'], dataId={'testType': 'BIAS','imageType': 'BIAS'})
#to select a sequence: 'seqNum': 812
print(visits_bias)

## Specifiy a visit/sensor/amplifier
visit = visits_bias[0]
sensor = 97
amp = 'C16'

## We will do a bias and offset correction on a bias image
#dId = {'visit': visit, 'detectorName': 'S22'}
dId = {'visit': visit, 'detector': sensor}
print(dId)
raw1 = butler.get('raw', **dId)
#bias1 = butler.get('bias', **dId)

#define output folder (useful to be on the web)
output_data='/sps/lsst/users/tguillem/web/tutorial/butler/'

## What do these images look like?
# raw image
plt.figure()
# get full array
arr = raw1.getImage().getArray()
plt.imshow(arr, origin='lower', vmin=26168, vmax=26172)
plt.colorbar()
plt.savefig(output_data+'raw.png') 
plt.close()
print('Full CCD image (includes overscan pixels): (columns,rows) = ' + str(arr.shape))

# per amp
detector = raw1.getDetector()
amplifier = detector[amp]
sub_im0 = raw1.getMaskedImage()[amplifier.getBBox()]
arr_amp = sub_im0.getImage().getArray()
#np.set_printoptions(threshold=sys.maxsize)
#print(arr_amp)
plt.figure()
plt.imshow(arr_amp, origin='lower', vmin=26168, vmax=26170)
plt.colorbar() 
plt.savefig(output_data+'raw_amp.png')
print('Science image of one amplifier: (columns,rows) = ' + str(arr_amp.shape))
