#!/usr/bin/env python
# coding: utf-8

##########
# Started from a DM 2019 bootcamp notebook
# Author: T. Guillemin

# Goal: do some bias plots and access overscan/image pixels
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

# ### Access images from the Butler
repo_path='/sps/lsst/groups/FocalPlane/SLAC/run5/butler/all_runs/13144/'
butler = Butler(repo_path)

#check the butler content
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "PRAGMA table_info(raw);"')
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "select * from raw limit 10;"')
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "select distinct testType from raw;"')
print('---')
os.system('sqlite3 '+repo_path+'/registry.sqlite3 "select distinct imageType from raw;"')

#select visits for BIAS images
visits_bias=butler.queryMetadata('raw', ['visit'], dataId={'testType': 'FLAT','imageType': 'BIAS'})
print(visits_bias)

## Specifiy a visit/sensor/amplifier
visit = visits_bias[3]
sensor = 97
amp = 'C16'

## We will do a bias and offset correction on a bias image
#dId = {'visit': visit, 'detectorName': 'S22'}
dId = {'visit': visit, 'detector': sensor}
print(dId)
raw1 = butler.get('raw', **dId)
#bias1 = butler.get('bias', **dId)

## What do these images look like?
# raw image
plt.figure()
# get full array
arr = raw1.getImage().getArray()
plt.imshow(arr, origin='lower', vmin=24000, vmax=26000)
plt.colorbar()
plt.savefig('raw.png') 
plt.close()
print('Full image (includes overscan pixels): (columns,rows) = ' + str(arr.shape))

# per amp
detector = raw1.getDetector()
amplifier = detector[amp]
sub_im0 = raw1.getMaskedImage()[amplifier.getBBox()]
arr_amp = sub_im0.getImage().getArray()
#np.set_printoptions(threshold=sys.maxsize)
#print(arr_amp)
plt.figure()
plt.imshow(arr_amp, origin='lower', vmin=26165, vmax=26175)
plt.colorbar() 
plt.savefig('raw_amp.png')
print('1 amp: (columns,rows) = ' + str(arr_amp.shape))
sys.exit()
###+++++++++++++++++++STOP HERE FOR NOW+++++++++++++++++++

#bias image
#plt.figure()
#plt.imshow(bias1.getImage().getArray(), origin='lower', vmin=-15, vmax=15)
#plt.colorbar()
#plt.savefig('bias.png')
#plt.close()


# ### Stack ISR bias/offset correction
# What does the superbias look like?

detector = bias1.getDetector()
amplifier = detector[amp]

## Imaging section only 
sub_im0 = bias1.getMaskedImage()[amplifier.getBBox()]

plt.imshow(sub_im0.getImage().getArray(), origin='lower', vmin=-15, vmax=15)
plt.colorbar()
plt.show()


# What does the super bias look like along rows and columns?
## Looks strange...
project(bias1.getImage())


# ##### Run ISR bias and offset correction
result1 = isr.run(raw1.clone(), bias=bias1.clone())


# #### Stack bias and offset-corrected image
# 
# Since our raw image is a bias image, we expect ~zero counts in our final image.

# In[13]:


## Get subimage for a specific amplifier
detector = result1.exposure.getDetector()
amplifier = detector[amp]

sub_im1 = result1.exposure.getMaskedImage()[amplifier.getBBox()]


# In[14]:


plt.imshow(sub_im1.getImage().getArray(), origin='lower', vmin=-15, vmax=15)
plt.colorbar()
plt.show()


# ### EOTest bias/offset correction

# In[15]:


superbias = get_superbias(visits, sensor, amp)


# In[16]:


## Why is mean bias level around -3 along the rows?
project(superbias)


# In[17]:


plt.imshow(superbias.getArray(), origin='lower', vmin=-15, vmax=15)
plt.colorbar()
plt.show()


# In[18]:


unbiased = eotest_unbias(raw1.clone(), amp, superbias)


# In[19]:


plt.imshow(unbiased.getArray(), origin='lower', vmin=-15, vmax=15)
plt.colorbar()
plt.show()


# ### Subtract the two overscan and bias corrected images

# In[20]:


diff = unbiased.clone().getArray() - sub_im1.getImage().clone().getArray()


# In[21]:


plt.imshow(diff, origin='lower', vmin=-5, vmax=5)
plt.colorbar()
plt.show()


# In[ ]:




