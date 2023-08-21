# LSST Science Pipelines (Stack) packages
import lsst.daf.butler as dafButler
import lsst.afw.display as afwDisplay

repo = '/sdf/group/rubin/repo/ir2/butler.yaml'
collection='u/lsstccs/defects_13391_w_2023_24'
butler = dafButler.Butler(repo)
registry = butler.registry
#for x in sorted(registry.queryDatasets(datasetType='*',collections=collection)):
#    print(x)
##################
    
with butler.export(filename="export.yaml", transfer="copy", directory="/sdf/data/rubin/user/tguillem/stack/w_2022_48/copy_defects_13391_w_2023_24/") as export:
    #export.saveCollection('u/tguillem/DM-37455/masterbias_1D_.20230103.v1')
    #export.saveCollection('u/jchiang/ptc_13162_w_2022_39')
                                                 
    #for datasetTypeName in ("bias", "dark", "flat", "sky", "cpBiasProc"):
    #    print('Exporting dataset: ',datasetTypeName)
    #    export.saveDatasets(butler.registry.queryDatasets(datasetTypeName, collections=collection),
    #                        elements=["exposure"])
    #collections=
    #for coll in butler.registry.queryCollections(collection, includeChains=True,
    #                                                   flattenChains=True):
    #    if coll.split('/')[1] != 'calib':
    #        print('Exporting collection: ',coll)
    #        export.saveCollection(coll)
    #export.saveDatasets('u/jchiang/ptc_13162_w_2022_39/20221117T200148Z')
    export.saveDatasets(butler.registry.queryDatasets('*', collections=collection))
    #                  elements=["exposure"])
    
