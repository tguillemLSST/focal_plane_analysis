#!/usr/bin/bash
 
#WARNING: LSST must not be already set up! 
echo "LSST load" 
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2023_05/loadLSST.bash  
echo "LSST setup" 
setup lsst_distrib  
 
echo "Job configuration:"  
echo ${run}  
echo ${raft}  
echo ${ccd}  
echo ${exposure}

export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/
export det=\'${raft}_${ccd}\'
export exp=${exposure}

###master bias
#cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-37455/run_13161_master_dark30_2D_20230314b -c isr:doDefect=False \
#   -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='dark' AND detector.full_name=$det AND exposure.seq_num>156 AND exposure.seq_num<162" \
#    --register-dataset-types

###bias
cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias_corr.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,u/tguillem/DM-37455/master_bias_2D_20230206a -o u/tguillem/waves/run_13161_raw_20230412d -c isr:doDefect=False \
    -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.id=$exp" \
    --register-dataset-types

###copy of useful lines
#-d "instrument='LSSTCam' AND exposure.science_program='13162' AND exposure.observation_type='bias' AND exposure.seq_num=192" \

###produce several master biases
#cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-37455/run_13161_master_dark240_1D_20230315d -c isr:doDefect=False -c isr:overscan.doParallelOverscan=False \
#    -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='dark' AND exposure.seq_num>171 AND exposure.seq_num<177" \
#    --register-dataset-types

#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-37455/run_13161_master_dark240_2D_20230315d -c isr:doDefect=False -c isr:overscan.doParallelOverscan=True \
#    -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='dark' AND exposure.seq_num>171 AND exposure.seq_num<177" \
#    --register-dataset-types

#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-37455/run_13161_master_bias_1D_20230315d -c isr:doDefect=False -c isr:overscan.doParallelOverscan=False \
#    -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='bias' AND exposure.seq_num>137 AND exposure.seq_num<143" \
#    --register-dataset-types

#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-37455/run_13161_master_bias_2D_20230315d -c isr:doDefect=False -c isr:overscan.doParallelOverscan=True \
#    -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='bias' AND exposure.seq_num>137 AND exposure.seq_num<143" \
#    --register-dataset-types
