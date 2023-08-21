export REPO=/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/
#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/
#export REPO=/sps/lsst/groups/auxtel/softs/shared/auxteldm_gen3/data/
#butler register-instrument $REPO 'lsst.obs.lsst.Latiss'

#calibs
#butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/butler/export_060103/calibs --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run5/butler/export_060103/calibs/export.yaml -s instrument -s detector -s physical_filter
#butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/butler/export_060103/calibs --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run5/butler/export_060103/calibs/export.yaml -s instrument -s detector -s physical_filter
#butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/butler/auxtel_calib_validation/calibs_DM-36719 --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run5/butler/auxtel_calib_validation/calibs_DM-36719/export.yaml -s instrument #-s detector -s physical_filter
#butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/defects_20230817/1 --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/defects_20230817/1/export.yaml -s instrument
#butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/defects_20230817/2 --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/defects_20230817/2/export.yaml -s instrument
#butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/defects_20230817/3 --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/defects_20230817/3/export.yaml -s instrument

#user
#butler import $REPO /sps/lsst/users/boutigny/export --transfer copy --export-file /sps/lsst/users/boutigny/export/export.yaml
butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/copy_defects_13391_w_2023_24/ --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/copy_defects_13391_w_2023_24/export.yaml -s instrument
butler import $REPO /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/copy_defects_13392_w_2023_24/ --transfer copy --export-file /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/copy_defects_13392_w_2023_24/export.yaml -s instrument

#remove calibration collections
#butler remove-collections $REPO LSSTCam/calib/DM-28636
#butler remove-collections $REPO LSSTCam/calib/DM-28636
#butler remove-collections $REPO LSSTCam/calib/DM-28636/unbounded
#butler remove-collections $REPO LSSTCam/calib/DM-36442/bias.20221026a
#butler remove-collections $REPO LSSTCam/calib/DM-36442/bias.20221026b
#butler remove-collections $REPO LSSTCam/calib/DM-36442/dark.20221026a
#butler remove-collections $REPO LSSTCam/calib/DM-36442/dark.20221026b
#butler remove-collections $REPO LSSTCam/calib/DM-36442/defects.20221026a
#butler remove-collections $REPO LSSTCam/calib/DM-36442/defects.20221026b
#butler remove-collections $REPO LSSTCam/calib/DM-36442/flat.20221026a
#butler remove-collections $REPO LSSTCam/calib/DM-36442/flat.20221026b
#butler remove-collections $REPO LSSTCam/calib/DM-36442/linearity.20221026a
#butler remove-collections $REPO LSSTCam/calib/DM-36442/linearity.20221026b
#butler remove-collections $REPO LSSTCam/calib/unbounded
#butler remove-collections $REPO u/jchiang/defects_13162_w_2023_19
#butler remove-collections $REPO u/jchiang/eo_bright_defects_13162_w_2023_19
#butler remove-collections $REPO u/jchiang/eo_dark_defects_13162_w_2023_19
