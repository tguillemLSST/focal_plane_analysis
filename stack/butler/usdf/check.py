import warnings
warnings.simplefilter("ignore", category=UserWarning)
import pylab as plt
import sys
import pandas as pd

# Set a standard figure size to use
plt.rcParams['figure.figsize'] = (15.0, 15.0)

# LSST Science Pipelines (Stack) packages
import lsst.daf.butler as dafButler
import lsst.afw.display as afwDisplay

afwDisplay.setDefaultBackend('matplotlib')

#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/test_master_dark_3/butler.yaml'
#repo = '/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/all_runs/13161/butler.yaml'
#repo = '/sdf/group/rubin/repo/main/butler.yaml'
#repo = '/sdf/group/rubin/repo/embargo/butler.yaml'
#repo = '/sdf/data/rubin/repo/ir2/butler.yaml'
#repo = '/sdf/data/rubin/repo/ops-rehearsal-3-prep/butler.yaml'
repo = 'embargo_or4'

#collections=["LATISS/raw/all"]#, "LATISS/calib"]
#collections=['LATISS/runs/AUXTEL_DRP_IMAGING_2023-09A-08ABC-07AB-05AB/d_2023_09_25/PREOPS-3780']
#butler = dafButler.Butler(repo,collections=collections)
butler = dafButler.Butler(repo,writeable=False)
registry = butler.registry
#method 1
#defects = butler.get('defects', instrument='LATISS', detector=0, collections='u/czw/calibX.20220608')
#registry.findDataset('defects', instrument='LATISS', detector=0, collections='LATISS/calib', timespan=raw.dataId.timespan)
#method 2
#raws = registry.queryDatasets("raw", where="exposure.day_obs=20240308", instrument='LATISS', detector=0, collections='LATISS/raw/all')
#for exp in raws:
#       print(exp.dataId.full)


#defects = registry.findDataset('defects', instrument='LATISS', detector=0, collections='LATISS/calib', timespan=raw.dataId.timespan)
#print(defects.dataId.full)
#AND exposure.seq_num=26

#calibType = 'bias'
#detectorId = 0
#dataId = {'exposure': 2023042500038, 'detector': 0, "instrument": 'LATISS'}
#calib = butler.get(calibType,  dataId=dataId)
#print(calib)

#check IOV
#results = registry.findDataset("bias", collections="LATISS/calib")
#for result in results:
#    print(result.timespan)

#sys.exit()
#for c in sorted(registry.queryCollections('CALIBRATION')):
#	print(c)

#for x in sorted(registry.queryDatasetTypes()):
#	print(x)

#results = butler.registry.queryDatasetAssociations("bias", collections="LATISS/calib")
#for result in results:
#    print(result.timespan, result.ref.run,)
#sys.exit()
#print(registry.dimensions.names)
#print(butler.registry.dimensions["exposure"].RecordClass.fields)
#raw_exps = registry.queryDatasets('postISRCCD', instrument='LATISS',where="exposure.science_program='AUXTEL_PHOTO_IMAGING' AND exposure.day_obs>20230731 AND exposure.day_obs<20230901",collections=['LATISS/runs/AUXTEL_DRP_IMAGING_20230509_20231207/w_2023_49/PREOPS-4648'])
#with butler.export(filename="export.yaml", transfer="copy", directory="auxtel_exposures_3") as export:
#       export.saveDatasets(raw_exps.expanded(), elements=["exposure"])
#print(raw_exps)
#DRP
#raw_exps = registry.queryDatasets(datasetType='postISRCCD',instrument='LATISS', where="exposure.science_program='AUXTEL_PHOTO_IMAGING'", collections=['LATISS/runs/AUXTEL_DRP_IMAGING_2023-05A/w_2023_19/PREOPS-3444'])
#print(raw_exps.expanded())
#with butler.export(filename="export.yaml", transfer="copy", directory="auxtel_drp_exposures") as export:
#       export.saveDatasets(raw_exps.expanded(), elements=["exposure"])

###to get metadata from an exposure
expMetadata = butler.get('postISRCCD.metadata', instrument='LSSTComCamSim',detector=0,exposure=7024062500429, collections=['LSSTComCamSim/nightlyValidation'])
print(expMetadata)
###end
sys.exit()

#collections=["LSSTComCam/raw/all", "LSSTComCam/calib"]
#butler = dafButler.Butler(repo,collections=collections)
#butler = dafButler.Butler(repo)
#registry = butler.registry
#for c in sorted(registry.queryCollections()):
#       print(c)
#OR3
#raw_exps = sorted(registry.queryDatasets('raw', instrument='LSSTComCamSim',where="day_obs=20240403 AND detector.id=8",collections=['LSSTComCamSim/raw/all']))
#OR4
#raw_exps = sorted(registry.queryDatasets('raw', instrument='LSSTComCamSim',where="day_obs=20240626 and detector.id=8",collections=['LSSTComCamSim/raw/all']))
#raw_exps = sorted(registry.queryDatasets('raw', instrument='LSSTCam',where="exposure.day_obs > 20231001 AND detector.id=69",collections=['LSSTCam/raw/all']))
#raw_exps = registry.queryDatasets('raw', instrument='LSSTComCam',where="exposure.day_obs=20220608 AND exposure.observation_type='bias'",collections=['LSSTComCam/raw/all'])
#raw_exps = registry.queryDatasets('raw', instrument='LSSTComCam',where="exposure.day_obs>20230630 AND exposure.day_obs<20230801",collections=['LSSTComCam/raw/all'])
#raw_exps = registry.queryDatasets('raw', instrument='LSSTComCamSim',where="day_obs=20240529",collections=['LSSTComCamSim/raw/all'])
#postISRCCD
#raw_exps = sorted(registry.queryDatasets('calexp', instrument='LSSTComCamSim',where="day_obs=20240627 and detector.id=8",collections=['LSSTComCamSim/nightlyValidation']))
#calibration exposures
#raw_exps = sorted(registry.queryDatasets('raw', instrument='LSSTComCamSim',where="detector.id=8 AND exposure.observation_type='bias'",collections=['LSSTComCamSim/raw/all']))
#with butler.export(filename="export.yaml", transfer="copy", directory="comcamsim_exposures") as export:
#       export.saveDatasets(raw_exps.expanded(), elements=["exposure"])
#for exp in raw_exps:
#        print(exp)
#sys.exit()
#collection = "LSSTCam/raw/all"
#collection = "LATISS/raw/all"
#butler = dafButler.Butler(repo, collections=collection)
#registry = butler.registry

###TEST
# load the butler
#butler_repo = "/sdf/data/rubin/repo/ir2/"
#butler_repo = "/sdf/group/rubin/repo/embargo/"
#butler_repo = "/sdf/group/rubin/repo/main/"
#butler = dafButler.Butler(butler_repo)
#day_obs = 20230524
#LSSTCam or LSST-TS8
#datasetRefs = butler.registry.queryDatasets(
#    datasetType="raw",
#    collections="LSSTCam/raw/all",
#    where=f"instrument='LSSTCam' AND exposure.day_obs > {day_obs}"
#).expanded()
#for i, ref in enumerate(datasetRefs):
#       print(ref.dataId.full)

#for c in sorted(butler.registry.queryCollections()):
#       print(c)
#sys.exit()
###

#collection='LATISS/calib/DM-43022'
collection='LATISS/calib/DM-43022/refactorCalibs/verifyBias.20240227b/20240228T173415Z'
########
#collection_ref='LATISS/calib/DM-43022'
#with butler.registry.caching_context():
    # Prime the cache with one bulk fetch of collection-summary information.
#    butler.registry.queryDatasets(..., collections=collection_ref).any(execute=False, exact=False)
#    for collection in butler.registry.queryCollections(collection_ref):
#for dataset_type_results in butler.registry.queryDatasets(..., collections=collection_ref).byParentDatasetType():
#    if dataset_type_results.any(execute=False, exact=False):
#        print(collection, dataset_type_results.parentDatasetType)
########

for datasetType in registry.queryDatasetTypes():
    if butler.registry.queryDatasets(datasetType, collections=collection).any(execute=False, exact=False):
        # Limit search results to the data products
           #if ('_config' not in datasetType.name) and ('_log' not in datasetType.name) and ('_metadata' not in datasetType.name) and ('_resource_usage' not in datasetType.name):
           print(datasetType)
sys.exit()
print('Exporting=========')
with butler.export(filename="export.yaml", transfer="copy", directory="./test" ) as export:
    for datasetTypeName in ("bias"):
        print('Exporting dataset: ',datasetTypeName)
        export.saveDatasets(butler.registry.queryDatasets(datasetTypeName, collections=collection),
                            elements=["exposure"])
    #print('bug')
    #export.saveDatasets("verifyBiasProc")
    #for coll in butler.registry.queryCollections(collection, includeChains=True,
    #                                             flattenChains=True):
    #if coll.split('/')[1] != 'calib':
        
    #    print('Exporting collection: ',coll)
    #    export.saveCollection(coll)
sys.exit()

datasetRefs = registry.queryDatasets(datasetType='*', collections=collection)
#datasetRefs = butler.registry.queryDatasets(datasetType='defects', collections=collection)

for i, ref in enumerate(datasetRefs):
	print(ref.dataId.full)
sys.exit()
#method 1
#visits=[]
#for i,dataId in enumerate(registry.queryDimensionRecords('exposure', where='exposure.science_program = \'13162\' AND exposure.observation_type = \'bias\' AND exposure.seq_num>177 AND expos#ure.seq_num<197')):
#for i,dataId in enumerate(butler.registry.queryDimensionRecords('exposure', where='exposure.detector=69 AND exposure.day_obs > 20230523')):
#        print(dataId.full)
#sys.exit()        
#for dataId in enumerate(registry.queryDimensionRecords('exposure', where="exposure.observation_type='bias'")):
        #print(dataId)
        #visits.append(dataId.)
        
#print(visits)

#to build the list of runs and dump the fields of exposures
print('=================List of runs')
#print the exposure fields
print(registry.dimensions["exposure"].RecordClass.fields)
df_exposure = pd.DataFrame(columns=['science_program'])
list_runs = []
for i, ref in enumerate(registry.queryDimensionRecords('exposure')):
        run = str(ref.science_program)
        #print(ref.id)
        #if(ref.id==3023062100485):
        #       print(ref)
        #'302306210026
        #if(ref.science_program=='13391' and ref.observation_type=='bias' and ref.observation_reason=='bias'):
                #print(ref)
                #print('+++')
        #       print(str(ref.id) + ' ' + str(ref.timespan))
        #       #print('++++++ '+str(ref))
        #if(ref.science_program=='13391' and ref.observation_type=='dark'):
        #       print('exposure ' + str(ref.id) + ' exposure time = ' + str(ref.exposure_time))
        #continue

        #print('reason = ' + ref.observation_reason + 'type = ' + ref.observation_type)

        if(run.startswith('12') or run.startswith('13')):
               if run not in list_runs:
                      list_runs.append(run)
               else:
                      continue
print(sorted(list_runs))
#print(list_runs)
sys.exit()

#metadata
expId = 2023051000565
#bias
#EXPOSURES=[2023042500023,2023042500028,2023042500033,2023042500038,2023042500043,2023042500048,2023042500053,2023042500058,2023042500063,2023042500068,2023042500323,2023042500328,2023042500333,2023042500338,2023042500343,2023042500348,2023042500353,2023042500358,2023042500363,2023042500368,2023042700002,2023042700004,2023042700006,2023042700008,2023042700010,2023042700012,2023042700014,2023042700016,2023042700018,2023042700020,2023042700105,2023042700107,2023042700109,2023042700111,2023042700113,2023042700115,2023042700117,2023042700119,2023042700121,2023042700123]
#dark
EXPOSURES=[2023042500121,2023042500125,2023042500129,2023042500133,2023042500137,2023042500141,2023042500145,2023042500149,2023042500153,2023042500157,2023042500161,2023042500165,2023042500169,2023042500173,2023042500177,2023042500181,2023042500185,2023042500189,2023042500193,2023042500197,2023042500201,2023042500205,2023042500209,2023042500213,2023042500217,2023042500221,2023042500225,2023042500229,2023042700030,2023042700033, 2023042700036, 2023042700039, 2023042700042]
for i in range(len(EXPOSURES)):
        mData = butler.get('raw.metadata', detector=0, exposure=EXPOSURES[i])
        #print(mData)
        #break
        print(str(EXPOSURES[i]) + ' : ' + str(mData['DARKTIME']))




#for dim in ['exposure', 'detector']:
#print(list(registry.queryDimensionRecords('exposure', where='exposure.science_program=1316 AND detector=166'))[0])

#data = registry.queryDatasets(datasetType='raw', instrument='LSSTCam', collections=collection,
#                              where="exposure.science_program='13162' AND exposure.observation_type='bias' AND detector=166")
#print(data)

###to get paths of files
#data = registry.queryDatasets(datasetType='raw', instrument='LSSTCam', collections=collection,
#                              where="exposure.science_program='13162' AND exposure.observation_type='bias' AND detector=166 AND exposure.id=3021121200177")
#
#with butler.export(filename="test.yaml") as export:
#            export.saveDatasets(data, elements=())
###
