###BOT  data
#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/test_master_dark_3
#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/all_runs/13162/
export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/butler.yaml

#working
export detector_1='R14_S22'
export detector=\'$detector_1\'
echo $detector

###from Dominique
pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
    -j 8 \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-37455/run_13161_master_bias_1D_20230314bfix -c isr:doDefect=False -c isr:overscan.doParallelOverscan=False\
    -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='dark' AND detector.full_name=$detector AND exposure.seq_num>156 AND exposure.seq_num<162" \
    --register-dataset-types

#-d "instrument='LSSTCam' AND exposure.observation_type='bias' AND detector.full_name='R14_S22' AND exposure.seq_num>176 AND exposure.seq_num<197" \
#AND exposure.seq_num > 0 AND exposure.seq_num < 20
#-d "instrument='LSSTCam' AND detector.full_name='R14_S22' AND exposure.observation_type='dark' AND exposure.dark_time > 119.0 AND exposure.dark_time < 121.0" \    
#detector.full_name='R14_S22'
#dark 30 s: exposure.seq_num>156 AND exposure.seq_num<162
##pipetask --long-log run -b $REPO -p $CP_PIPE_DIR/pipelines/LSSTCam/cpBias.yaml \
##     -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/darkGen_0 \
##   -c isr:doDefect=False
##     -d "instrument='LSSTCam' AND detector=0 AND exposure IN ($EXPOSURES) \
##     -c isr:doDefect=False -c isr:doEmpiricalReadNoise=True >& ./dark.log

#pipetask --log-level DEBUG

###AUXTEL data
#export REPO=/sps/lsst/groups/auxtel/softs/shared/auxteldm_gen3/data/
#
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-30000/biasGen.full.20220705a -c isr:doDefect=False \
