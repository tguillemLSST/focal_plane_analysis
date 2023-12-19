#!/usr/bin/bash
 
#WARNING: LSST must not be already set up! 
echo "LSST load" 
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2023_34/loadLSST.bash  
echo "LSST setup" 
setup lsst_distrib  
 
echo "Job configuration:"
echo ${run}
echo ${raft}
echo ${ccd}
#echo ${exposure}

export REPO=/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/butler.yaml
export det=\'${raft}_${ccd}\'
export run_number=\'${run}\'
#export exposure_id=${exposure}
#export EXPOSURES=\'${exposure_1}\'\,\'${exposure_2}\'\,\'${exposure_3}\'\,\'${exposure_4}\'\,\'${exposure_5}\'
export exposure_id_1=${exposure_1}
export exposure_id_2=${exposure_2}
export exposure_id_3=${exposure_3}
export exposure_id_4=${exposure_4}
export exposure_id_5=${exposure_5}
#echo $EXPOSURES

###master bias
#cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_05/
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/DM-37455/run_13161_master_dark30_2D_20230314b -c isr:doDefect=False \
#   -d "instrument='LSSTCam' AND exposure.science_program='13161' AND exposure.observation_type='dark' AND detector.full_name=$det AND exposure.seq_num>156 AND exposure.seq_num<162" \
#    --register-dataset-types

###bias
#cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_34/
#pipetask --long-log run -b $REPO -p cp_pipe/pipelines/cpBias_corr.yaml \
#    -j $SLURM_CPUS_PER_TASK \
#    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/debug_run_6b_0D_20231027d/exp_$exp -c isr:doDefect=False \
#    -d "instrument='LSSTCam' AND exposure.id=$exp" \
#    --register-dataset-types

#test many exposures
cd /sps/lsst/users/tguillem/Rubin/stack/w_2023_34/
pipetask --long-log run -b $REPO -p cp_pipe/pipelines/LsstCam/cpBias_corr.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_0D_exp_${exposure_id_1}_20231127d -c isr:doDefect=False \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=False -c isr:overscan.fitType='MEDIAN' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.id=$exposure_id_1" \
    --register-dataset-types

pipetask --long-log run -b $REPO -p cp_pipe/pipelines/LsstCam/cpBias_corr.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_0D_exp_${exposure_id_2}_20231127d -c isr:doDefect=False \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=False -c isr:overscan.fitType='MEDIAN' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.id=$exposure_id_2" \
    --register-dataset-types

pipetask --long-log run -b $REPO -p cp_pipe/pipelines/LsstCam/cpBias_corr.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_0D_exp_${exposure_id_3}_20231127d -c isr:doDefect=False \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=False -c isr:overscan.fitType='MEDIAN' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.id=$exposure_id_3" \
    --register-dataset-types

pipetask --long-log run -b $REPO -p cp_pipe/pipelines/LsstCam/cpBias_corr.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_0D_exp_${exposure_id_4}_20231127d -c isr:doDefect=False \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=False -c isr:overscan.fitType='MEDIAN' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.id=$exposure_id_4" \
    --register-dataset-types

pipetask --long-log run -b $REPO -p cp_pipe/pipelines/LsstCam/cpBias_corr.yaml \
    -j $SLURM_CPUS_PER_TASK \
    -i LSSTCam/raw/all,LSSTCam/calib -o u/tguillem/run_6_validation/run_${run}_0D_exp_${exposure_id_5}_20231127d -c isr:doDefect=False \
    -c isr:doOverscan=True -c isr:overscan.doParallelOverscan=False -c isr:overscan.fitType='MEDIAN' \
    -d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.id=$exposure_id_5" \
    --register-dataset-types

#-d "instrument='LSSTCam' AND exposure.science_program=$run_number AND (exposure.id=$exposure_id_1 OR exposure.id=$exposure_id_2 OR exposure.id=$exposure_id_3 OR exposure.id=$exposure_id_4 OR exposure.id=$exposure_id_5)"
#-d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure IN ($EXPOSURES)" \
#-d "instrument='LSSTCam' AND exposure.science_program=$run_number AND exposure.id=$exposure_id" \
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
