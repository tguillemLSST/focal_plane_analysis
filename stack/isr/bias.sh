###BOT  data
#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/butler.yaml
export REPO=/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/butler.yaml
#export REPO=/sps/lsst/groups/auxtel/softs/shared/auxteldm_gen3/data/butler.yaml
#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run6/butler/test_comcam/main_231023/butler.yaml

###from Dominique
#3021121200146,3021121200158,3021121200174
#exposure_id=3023062100508
#EXPOSURES='3023062100497, 3023062100498, 3023062100499, 3023062100500, 3023062100501, 3023062100502, 3023062100503, 3023062100504, 3023062100505, 3023062100506, 3023062100507, 3023062100508, 3023062100509, 3023062100510, 3023062100511'
#EXPOSURES='3023062100512, 3023062100513, 3023062100514'
#EXPOSURES='3023062100512'
EXPOSURES='3023102600140'
pipetask --log-level DEBUG --long-log run -b $REPO -p cp_pipe/pipelines/cpBias_corr.yaml \
    -j 8 \
    -i LSSTCam/raw/all,LSSTCam/calib,u/lsstccs/defects_13392_w_2023_24 -o u/tguillem/run_6_validation/debug_run_6b_0D_20231027a -c isr:doDefect=False \
    -d "instrument='LSSTCam' AND exposure IN ($EXPOSURES) AND detector.full_name='R13_S11'" \
    --register-dataset-types

#AuxTel
#EXPOSURES='2023091400065'
#pipetask --log-level DEBUG --long-log run -b $REPO -p cp_pipe/pipelines/cpBias_corr.yaml \
#    -j 8 \
#    -i LATISS/raw/all,LATISS/calib -o u/tguillem/auxtel_isr_validation/raw_exp_2D_20231023a -c isr:doDefect=False \
#    -d "instrument='LATISS' AND exposure IN ($EXPOSURES)" \
#    --register-dataset-types

#ComCam
#EXPOSURES='2023071900023'
#pipetask --log-level DEBUG --long-log run -b $REPO -p cp_pipe/pipelines/cpBias_corr.yaml \
#    -j 8 \
#    -i LSSTComCam/raw/all -o u/tguillem/comcam_isr_validation/raw_exp_0D_20231025a -c isr:doDefect=False \
#    -d "instrument='LSSTComCam' AND exposure IN ($EXPOSURES) AND detector=0" \
#    --register-dataset-types

# --instrument lsst.obs.lsst.LsstComCam 
#-d "instrument='LSSTComCam'" \
#--log-level DEBUG
#20211207/MC_C_20211207_000484
#AND exposure.id=3021120700063" \
#=3021121200146 20211207_000634
#-i u/tguillem/DM-37455/run_13162_master_bias_2D_20230119a
#--log-level DEBUG
#-d "instrument='LSSTCam' AND exposure.observation_type='dark' AND exposure.id=3021121200161" \
#AND exposure.seq_num > 0 AND exposure.seq_num < 20
#-d "instrument='LSSTCam' AND detector.full_name='R14_S22' AND exposure.observation_type='dark' AND exposure.dark_time > 119.0 AND exposure.dark_time < 121.0" \    
#detector.full_name='R14_S22'
#detector.raft='R14'
#exposure.seq_num=161
#AND exposure.id=$exposure_id
##pipetask --long-log run -b $REPO -p $CP_PIPE_DIR/pipelines/LSSTCam/cpBias.yaml \
##     -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/darkGen_0 \
##   -c isr:doDefect=False
##     -d "instrument='LSSTCam' AND detector=0 AND exposure IN ($EXPOSURES) \
##     -c isr:doDefect=False -c isr:doEmpiricalReadNoise=True >& ./dark.log

###AUXTEL data
#export REPO=/sps/lsst/groups/auxtel/softs/shared/auxteldm_gen3/data/
#
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-30000/biasGen.full.20220705a -c isr:doDefect=False \
