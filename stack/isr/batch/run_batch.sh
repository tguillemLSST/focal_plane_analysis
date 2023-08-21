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

#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111/
export REPO=/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/butler.yaml
export det=\'${raft}_${ccd}\'
export run_number=\'${run}\'

###master bias
#cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_13372_master_bias_2D_20230619a -c isr:doDefect=False \
#    -d "instrument='LSSTCam' AND exposure.science_program='13372' AND exposure.observation_type='bias'" \
#    --register-dataset-types

###bias
#cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias_corr.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_13378_1D_20230621a -c isr:doDefect=False \
#    -d "instrument='LSSTCam' AND exposure.science_program='13378' AND exposure.observation_type='bias' AND exposure.id=3023061900157" \
 #   --register-dataset-types

###copy of useful lines
#-d "instrument='LSSTCam' AND exposure.science_program='13162' AND exposure.observation_type='bias' AND exposure.seq_num=192" \

###produce several master biases
cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_master_bias_0D_20230817a -c isr:doDefect=True \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=False -c isr:overscan.fitType='MEDIAN' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.observation_type='bias' AND exposure.observation_reason='bias'" \
    --register-dataset-types

cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_master_bias_1D_20230817a -c isr:doDefect=True \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=False -c isr:overscan.fitType='MEDIAN_PER_ROW' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.observation_type='bias' AND exposure.observation_reason='bias'" \
    --register-dataset-types

cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_master_bias_2D_20230817a -c isr:doDefect=True \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=True -c isr:overscan.fitType='MEDIAN_PER_ROW' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.observation_type='bias' AND exposure.observation_reason='bias'" \
    --register-dataset-types

#doOverscan: true
#     overscan.doParallelOverscan: true
#     overscan.fitType: 'MEDIAN_PER_ROW'

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
