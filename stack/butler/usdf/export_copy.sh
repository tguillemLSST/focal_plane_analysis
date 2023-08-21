#butler export-calibs /sdf/group/rubin/repo/main_20210215 ./calibs_DM-36484 LATISS/calib/DM-36484/bias.20221005a
#scp -r calibs_DM-36719 tguillem@cca.in2p3.fr:/sps/lsst/groups/FocalPlane/SLAC/run5/butler/export_060103/.
#scp -r export_user tguillem@cca.in2p3.fr:/sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/.

#defects
#butler export-calibs /sdf/group/rubin/repo/main ./defects_20230817/1 u/jchiang/defects_13162_w_2023_19
#butler export-calibs /sdf/group/rubin/repo/main ./defects_20230817/2 u/jchiang/eo_bright_defects_13162_w_2023_19
#butler export-calibs /sdf/group/rubin/repo/main ./defects_20230817/3 u/jchiang/eo_dark_defects_13162_w_2023_19

#scp -r defects_20230817 tguillem@cca.in2p3.fr:/sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/.
scp -r copy_defects_13391_w_2023_24 tguillem@cca.in2p3.fr:/sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/.


