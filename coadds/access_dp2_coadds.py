#!/usr/bin/env python
# coding: utf-8

##########
# Author: T. Guillemin
# Goal: access DP2 coadds
##########

import lsst.daf.butler as dafButler
import lsst.geom as geom

#butler
butler=dafButler.Butler('dp2_prep',collections='LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage3',skymap='lsst_cells_v2')
registry = butler.registry

#get tract and patch from on-sky geometrical info
ra, dec = 53.8, -28.1
radec = geom.SpherePoint(ra, dec, geom.degrees)
skymap = butler.get("skyMap")
tractInfo = skymap.findTract(radec)
tract = tractInfo.getId()
patchInfo = tractInfo.findPatch(radec)
patch = tractInfo.getSequentialPatchIndex(patchInfo)
print('tract = ' + str(tract) + ' and patch = ' + str(patch))

bands = ['u','g','r','i','z','y']
for i_band in range(len(bands)):
        coaddId = {'tract': tract, 'patch': patch, 'band': bands[i_band]}
        coadd = butler.get('deep_coadd', dataId=coaddId)
        print(bands[i_band] + ' : ' + str(len(coadd.getInfo().getCoaddInputs().visits['visit'])))
        image = coadd.getImage().getArray()
        print(image)
