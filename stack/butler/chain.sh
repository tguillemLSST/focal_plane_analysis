export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/

#full command
#butler collection-chain $REPO LSSTCam/calib --mode=prepend LSSTCam/calib/DM-36442/bias.20221026a LSSTCam/calib/DM-36442/dark.20221026a LSSTCam/calib/DM-36442/flat.20221026a LSSTCam/calib/DM-36442/defects.20221026a LSSTCam/calib/DM-36442/linearity.20221026a LSSTCam/calib/DM-36442/bias.20221026b LSSTCam/calib/DM-36442/dark.20221026b LSSTCam/calib/DM-36442/flat.20221026b LSSTCam/calib/DM-36442/defects.20221026b LSSTCam/calib/DM-36442/linearity.20221026b

#only the current valid calibrations
butler collection-chain $REPO LSSTCam/calib --mode=prepend LSSTCam/calib/DM-37455/bias.20230119a LSSTCam/calib/DM-36442/dark.20221026b LSSTCam/calib/DM-36442/flat.20221026b LSSTCa\
m/calib/DM-36442/defects.20221026b LSSTCam/calib/DM-36442/linearity.20221026b
