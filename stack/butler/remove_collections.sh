export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/

#get list of CHAINS
#butler query-collections /sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/ | grep tguillem | grep chained

#remove chains first
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/biasCorr.20230117a

#remove runs after
butler remove-runs --no-confirm $REPO u/tguillem/DM-37455/*

#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/biasCorr.20230117a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/biasCorr.20230117a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/biasCorr.20230117a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/biasCorr.20230117c
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/biasCorr_1D.20230118a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/biasCorr_2D.20230118a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/master_bias_1D_20230111a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/master_bias_1D_20230113b
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/master_bias_1D_20230113c
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/master_bias_1D_20230116a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13161_1D_20230118b
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13161_2D_20230118b
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_1D.20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_1D_MB.20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_1D_MB.20230119b
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D.20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D.20230119b
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D_MB.20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D_MB.20230119b
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D_MB_fix.20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D_fix.20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D_my_MB.20230119b
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_2D_my_MB.20230119c
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_master_bias_2D_20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_median.20230119a
#butler remove-collections --no-confirm $REPO u/tguillem/DM-37455/run_13162_median_MB.20230119a
