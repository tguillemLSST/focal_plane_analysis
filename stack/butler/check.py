import sys
import pylab as plt
import pandas as pd

# Set a standard figure size to use
plt.rcParams['figure.figsize'] = (15.0, 15.0)

# LSST Science Pipelines (Stack) packages
import lsst.daf.butler as dafButler
import lsst.afw.display as afwDisplay
from uuid import UUID

afwDisplay.setDefaultBackend('matplotlib')

def print_calibrations(my_dataset,my_collection):
	#for dataset_type in butler.registry.queryDatasetTypes(...):
	#	if not dataset_type.isCalibration():
	#	continue
	for assoc in butler.registry.queryDatasetAssociations(my_dataset, collections=my_collection):
		print(
			my_dataset,
			assoc.ref.run,
			assoc.collection,
			assoc.ref.dataId,
			assoc.timespan,
		)
		#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/test_master_dark_3/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/all_runs/13162/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/shared/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/copy_from_usdf/sdf/group/rubin/repo/main/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/butler.yaml'
#repo = '/sps/lsst/groups/auxtel/softs/shared/auxteldm_gen3/data/'
repo = '/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run6/butler/test_comcam/main_231023/butler.yaml'
print('Butler: ' + repo)
butler = dafButler.Butler(repo,collections="LSSTCam/raw/all")
registry = butler.registry
#print(registry.getDataset(UUID("05c2f841-9b64-49f7-a177-964e98c49578")))
#sys.exit()
#print_calibrations('bias',['LSSTComCam/calib/DM-33657'])
#print('YES')
#sys.exit()
#dataId = {'exposure': 2022031600563, 'detector': 0, "instrument": 'LATISS'}
#calib = butler.get('bias', collections=['LATISS/calib/DM-38946/noRGseq/biasGen.20230428a/20230428T210637Z'], dataId=dataId)
#dataId = {'exposure': 2021121200150, 'detector': 53, "instrument": 'LSSTCam'}
#calib_1 = butler.get('dark',dataId=dataId,collections=['LSSTCam/calib/DM-36442/dark.20221026a'])
	    
#arr_1 = calib.getImage().getArray()
#print(arr_1[12][18])
#calib = butler.get('bias', collections=['LATISS/calib/DM-36719/biasGen.20221107b/20221107T213306Z'], dataId=dataId)
#arr_1 = calib.getImage().getArray()
#print(arr_1[12][18])

#print('=======Collections')
#for c in sorted(registry.queryCollections()):
#	print(c)
#print('=======Datasets')	
#for x in sorted(registry.queryDatasetTypes()):
#	print(x)
#sys.exit()
#print('Full list')
#collection='LATISS/calib/DM-36719/flatGen-SDSSg.20221107a/20221108T002737Z'
#collection='LATISS/calib/DM-38946/noRGseq/flatGen-i.20230501a/20230501T211541Z'
#'LATISS/calib/DM-38946/noRGseq/flatGen-r.20230501a/20230501T205909Z'
#collection = 'LSSTCam/raw/all'
#collection = 'u/tguillem/run_6_validation/run_13392_bias_exp_0D_20230912a'
#butler = dafButler.Butler(repo, collections=collection)
#registry = butler.registry

#datasetRefs = registry.queryDatasets(datasetType='raw', collections=collection)

#for i, ref in enumerate(datasetRefs):
#	print(ref.dataId.full)

#sys.exit()
#exp = butler.get("raw", instrument="LSSTCam", day_obs=20231115, seq_num=30, full_name="R13_S12", collections="LSSTCam/raw/all")
#expMetadata = butler.get('raw.metadata',instrument='LSSTCam', full_name="R13_S12", science_program="13550", exposure=2023111600008, collections="LSSTCam/raw/all")
#print(expMetadata)
#for i in range(2023111600400,2023111602223):
#for i in range(2023111600008,2023111600273):
#	expMetadata = butler.get('raw.metadata',instrument='LSSTCam', full_name="R13_S12", exposure = i, collections="LSSTCam/raw/all")
#	print(str(i) + ' | ' + expMetadata['TESTTYPE'] + ' | ' + expMetadata['IMGTYPE'] + ' | ' + expMetadata['DATE-BEG'] + ' | ' + expMetadata['DATE-END'] + ' | ' +  str(expMetadata['MJD-BEG']))
	#print([(k,v) for k,v in expMetadata.items() if "DATE-BEG" in k])
#sys.exit()
#check
#my_collections = ['u/tguillem/DM-37455/run_13161_master_dark30_2D_20230315c_v3']
#butler = dafButler.Butler(repo,collections=my_collections)
#registry = butler.registry
#calibType = 'bias'
#detectorId = 120
#calib = butler.get(calibType, instrument='LSSTCam', detector=detectorId)
#print('OK')
#test
#my_collections = ['LATISS/raw/all', 'LATISS/calib']
#butler = dafButler.Butler(repo)#,collections=my_collections)
#registry = butler.registry

#debug
#datasetRefs = registry.queryDatasetAssociations(datasetType='bias',collections='LATISS/calib/DM-38946')#noRGseq/biasGen.20230428a')
#print(datasetRefs[0].ref())
#dataset = butler.get(datasetRefs)
#print(dataset)
#sys.exit()

#calibType = 'bias'
#detectorId = 0
#dataId = {'exposure': 2022031600563, 'detector': 0, "instrument": 'LATISS'}
#calib = butler.get(calibType,  collections=['LATISS/raw/all', 'LATISS/calib/DM-38946/noRGseq/biasGen.20230428a'], dataId=dataId)
#print(calib)
#calib = butler.get(calibType, instrument='LATISS', detector=detectorId)
#sys.exit()
#print('=======Datasets')
#for x in sorted(registry.queryDatasetTypes()):
#	print(x)
#sys.exit()
#print(registry.dimensions.names)
#print(butler.registry.dimensions["exposure"].RecordClass.fields)

#print('=================Specific collection')
collection = "LSSTCam/raw/all"
#collection="LSSTCam/calib/DM-37455"
#collection = "u/tguillem/DM-37455/master_bias_1D_20230111a"
#collection="LSSTCam/calib/DM-36442/dark.20221026b"
#collection="LSSTCam/calib"
#butler = dafButler.Butler(repo, collections=collection)
#butler = dafButler.Butler(repo)
#registry = butler.registry

#exp = butler.get("raw", instrument="LSSTCam", day_obs=20231115, seq_num=30, full_name="R13_S12", collections="LSSTCam/raw/all")
#print(exp)
#sys.exit()
#datasetRefs = registry.queryDatasets(datasetType='bias', collections=collection)
#datasetRefs = registry.queryDatasets(datasetType='*')

#for x in sorted(registry.queryDatasetTypes()):
#	       print(x)
#print('===1')
#for i, ref in enumerate(datasetRefs):
#	print(ref.dataId.full)
#sys.exit()

#try to get a dataset
#calib_1 = butler.get('d', instrument='LSSTCam', detector=48)
#sys.exit()

#to build the list of runs and dump the fields of exposures
#print('=================List of runs')
#print the exposure fields
#print(registry.dimensions["exposure"].RecordClass.fields)
#df_exposure = pd.DataFrame(columns=['science_program'])
#list_runs = []
#for i, ref in enumerate(registry.queryDimensionRecords('exposure')):
#	run = str(ref.science_program)
	#print(ref.id)
	#if(ref.id==3023062100485):
	#	print(ref)
	#'302306210026
#	if(ref.science_program=='13513'): #and ref.observation_type=='bias' and ref.observation_reason=='bias'):
#		print(ref)
#		print('+++')
#		print(str(ref.id) + ' ' + str(ref.timespan))
	#	#print('++++++ '+str(ref))
	#if(ref.science_program=='13391' and ref.observation_type=='dark'):
	#	print('exposure ' + str(ref.id) + ' exposure time = ' + str(ref.exposure_time))
	#continue
	
	#print('reason = ' + ref.observation_reason + 'type = ' + ref.observation_type)

#	if(run.startswith('13')):
#		if run not in list_runs:
#			list_runs.append(run)
#		else:
#			continue
#print(sorted(list_runs))
#sys.exit()
#using datasetRefs
#print the exposure fields
#print(registry.dimensions["exposure"].RecordClass.fields)
#sys.exit()
#datasetRefs = list(registry.queryDatasets(datasetType='raw', instrument='LSSTCam', where="exposure.science_program='13162' AND exposure.observation_type='bias' AND detector.full_name='R14_S12' AND exposure.seq_num>176 AND exposure.seq_num<197"))
#bias
#datasetRefs = list(registry.queryDatasets(datasetType='raw', instrument='LSSTCam', where="exposure.science_program='13391' AND exposure.observation_type='bias' AND exposure.observation_reason='bias' AND detector.full_name='R14_S02'"))
#dark
#datasetRefs = list(registry.queryDatasets(datasetType='raw', instrument='LSSTCam', where="exposure.science_program='13391' AND exposure.observation_type='dark' AND detector.full_name='R14_S02'"))

#dark list
#df_exposure = pd.DataFrame(columns=['science_program'])
#exposures = []
#for i, ref in enumerate(registry.queryDimensionRecords('exposure')):
#	if(ref.science_program=='13550'): #and ref.observation_type=='bias' and ref.observation_reason=='bias'):
		#print(ref.full)
		#print('+++')
#		print(str(ref.id) + ' | observation_type = ' + str(ref.observation_type) + ' | seq_start = ' + str(ref.seq_start) + ' | seq_end = ' + str(ref.seq_end) + ' | exposure_time = ' + str(ref.exposure_time))
#		exposures.append(str(ref.id))
#exposures_sorted = sorted(exposures)
#2nd loop to get results with sorted exposures
#for j in exposures_sorted:
#	for i, ref in enumerate(registry.queryDimensionRecords('exposure')):
#		if(ref.id==exposures_sorted[j]):
#			print(str(ref.id) + ' | exposure time = ' + str(ref.exposure_time))
#sys.exit()

#datasets
datasetRefs = list(registry.queryDatasets(datasetType='raw', instrument='LSSTCam', where="exposure.science_program='13592' AND exposure.observation_type='dark' AND detector.full_name='R14_S02'"))
print(len(datasetRefs))
with open('exposures_13592.txt', 'w') as f:
	for i in range(len(datasetRefs)):
		dataId = datasetRefs[i].dataId
#		#print(str(dataId['exposure'].exposure_time))
		#print("\'" + str(dataId['exposure']) + "\'")
		if(i<len(datasetRefs)-1):
			f.write(str(dataId['exposure']) + "\n")
		else:
			f.write(str(dataId['exposure']))
#		#print(dataId.full)
sys.exit()

print(datasetRefs)
dataId = datasetRefs[0].dataId
print(dataId)
#print(dataId.full)
print(dataId['exposure'])
dataset = butler.get(datasetRefOrType='raw',dataId=dataId,collections=collection)
#dataset = butler.get(datasetRefs[0])
print(dataset)
