###BOT  data 
export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/

pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpDark.yaml \
    -i LSSTCam/raw/all,LSSTCam/calib \
    -o u/tguillem/DM-37455/darkGen.20230516a \
    -c isr:doDefect=False \
    -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='dark' AND detector.full_name='R14_S22'" \
    --register-dataset-types

#works for test_master_dark_2
#-i LSSTCam/raw/all,LSSTCam/calib \
#-i LSSTCam/raw/all,LSSTCam/calib,u/tguillem/DM-30001/biasGen.full.20220919a \

#does not work
#-i LSSTCam/raw/all,u/tguillem/master_bias_0/20220620T084256Z,LSSTCam/calib/DM-30000 \

##############
#pipetask --long-log run -b $REPO -p $CP_PIPE_DIR/pipelines/cpBias.yaml \
#     -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/darkGen_0 \
#     -c isr:doDefect=False

#     -d "instrument='LSSTCam' AND detector=0 AND exposure IN ($EXPOSURES) \
#     -c isr:doDefect=False -c isr:doEmpiricalReadNoise=True >& ./dark.log

#RERUN=20210707a
#pipetask --long-log run -b $BUTLER_REPO -p $CP_PIPE_DIR/pipelines/Latiss/cpDark.yaml \
#    -i LATISS/raw/all,u/czw/DM-28920/defectGen.20210706h,u/czw/DM-28920/biasGen.20210702a,LATISS/calib \
#    -o u/czw/DM-28920/darkGen
#    -d "instrument='LATISS' AND detector=0 AND exposure IN ($EXPOSURES) \
#    >& dark.$RERUN.log

#-d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='dark' AND detector.full_name=$detector AND exposure.seq_num>156 AND exposure.seq_num<162"
#-d "instrument='LSSTCam' AND detector.full_name='R14_S22' AND exposure.observation_type='dark' AND exposure.dark_time > 29.0 AND exposure.dark_time < 31.0" \
