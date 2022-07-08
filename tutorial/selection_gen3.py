import pylab as plt
import numpy as np
import pandas as pd
import sys
import os
import time
from astropy.io import fits
import lsst.afw.display as afwDisplay
# LSST Science Pipelines (Stack) packages
import lsst.daf.butler as dafButler
import lsst.afw.display as afwDisplay

#options
afwDisplay.setDefaultBackend('matplotlib')
#pd.set_option("display.max_rows", None, "display.max_columns", None)

start = time.time()

#butler path
repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/all_runs/13032/butler.yaml'
#LSSTCam
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/shared/butler.yaml'

butler = dafButler.Butler(repo)
registry = butler.registry

for c in sorted(registry.queryCollections()):
	print(c)
for x in sorted(registry.queryDatasetTypes()):
	print(x)

collection = "LSSTCam/raw/all"
butler = dafButler.Butler(repo, collections=collection)
registry = butler.registry

# build a DataFrame from exposures
df_exposure = pd.DataFrame(columns=['instrument','id','physical_filter','obs_id','exposure_time','dark_time','observation_type','observation_reason','day_obs','seq_num','group_name','group_id','target_name','science_program','tracking_ra','tracking_dec','sky_angle','zenith_angle','timespan'])

for i, ref in enumerate(registry.queryDimensionRecords('exposure')):
	df_exposure.loc[i] = [ref.instrument,ref.id,ref.physical_filter,ref.obs_id,ref.exposure_time,ref.dark_time,ref.observation_type,ref.observation_reason,ref.day_obs,ref.seq_num,ref.group_name,ref.group_id,ref.target_name,ref.science_program,ref.tracking_ra,ref.tracking_dec,ref.sky_angle,ref.zenith_angle,ref.timespan]

df_bias = df_exposure[df_exposure.observation_type == 'bias']
run_all=['13151','13159','13032']
visits_all=[]
for run_cur in run_all :
	df_selected = df_bias[df_bias.science_program == run_cur]
	visits = df_selected['obs_id'].tolist()
	if len(visits) < 1 :
		print('No image in run ',run_cur)
	else :
		visits.sort()
		visits_all.append(visits)
		print('For run ',run_cur,' we identified ',len(visits_all[-1]),' FLAT images')

print('Selected flats:')
print(visits_all)
