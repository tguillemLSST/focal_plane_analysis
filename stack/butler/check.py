import sys
import pylab as plt
import pandas as pd

# Set a standard figure size to use
plt.rcParams['figure.figsize'] = (15.0, 15.0)

# LSST Science Pipelines (Stack) packages
import lsst.daf.butler as dafButler
import lsst.afw.display as afwDisplay

afwDisplay.setDefaultBackend('matplotlib')

#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/test_master_dark_3/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/all_runs/13162/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/shared/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/copy_from_usdf/sdf/group/rubin/repo/main/butler.yaml'
repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/butler.yaml'
butler = dafButler.Butler(repo)
registry = butler.registry

print('=======Collections')
for c in sorted(registry.queryCollections()):
	print(c)

#print('=======Datasets')
#for x in sorted(registry.queryDatasetTypes()):
#	print(x)
#sys.exit()
#print(registry.dimensions.names)
#print(butler.registry.dimensions["exposure"].RecordClass.fields)

#print('=================Specific collection')
collection = "LSSTCam/raw/all"
#collection="LSSTCam/calib/DM-37455"
#collection = "u/tguillem/DM-37455"
butler = dafButler.Butler(repo, collections=collection)
registry = butler.registry

#datasetRefs = registry.queryDatasets(datasetType='raw', collections=collection)
datasetRefs = registry.queryDatasets(datasetType='*')

#for i, ref in enumerate(datasetRefs):
#	print(ref.dataId.full)

#sys.exit()

print('=================List of runs')
#print the exposure fields
print(registry.dimensions["exposure"].RecordClass.fields)
df_exposure = pd.DataFrame(columns=['science_program'])
list_runs = []
for i, ref in enumerate(registry.queryDimensionRecords('exposure')):
	run = ref.science_program
	if run not in list_runs:
		list_runs.append(run)
	else:
		continue
print(list_runs)

#using datasetRefs
#print the exposure fields
#print(registry.dimensions["exposure"].RecordClass.fields)
datasetRefs = list(registry.queryDatasets(datasetType='raw', instrument='LSSTCam', where="exposure.science_program='13162' AND exposure.observation_type='bias' AND detector.full_name='R14_S12' AND exposure.seq_num>176 AND exposure.seq_num<197"))
for i in range(len(datasetRefs)):
	dataId = datasetRefs[i].dataId
	print(dataId['exposure'])
sys.exit()	
print(datasetRefs)
dataId = datasetRefs[0].dataId
print(dataId)
print(dataId.full)
print(dataId['exposure'])
dataset = butler.get(datasetRefOrType='raw',dataId=dataId,collections=collection)
#dataset = butler.get(datasetRefs[0])
print(dataset)
