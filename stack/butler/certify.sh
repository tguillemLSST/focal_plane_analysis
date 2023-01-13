export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/

#butler certify-calibrations $REPO u/tguillem/DM-37455/master_bias_1D_20230111a/20230111T094034Z LSSTCam/calib/DM-37455/bias.20230111a --begin-date 2023-01-01 --end-date 2050-01-01 bias
butler certify-calibrations $REPO u/tguillem/DM-37455/master_bias_1D_20230113c LSSTCam/calib/DM-37455/bias.20230113c --begin-date 2023-01-01 --end-date 2050-01-01 bias

#check certification
#butler query-collections $REPO LSSTCam/calib
