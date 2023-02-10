#!/usr/bin/zsh 
 
#WARNING: LSST must not be already set up! 
echo "LSST load" 
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2022_01/loadLSST.zsh  
echo "LSST setup" 
setup lsst_distrib  
 
echo "setup eotest" 
cd /sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/eotest 
source setup.sh 
 
cd /sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/focal_plane_analysis/2D_bias 
echo "Starting to run eotest" 
 
echo "Job configuration:"  
echo ${run}  
echo ${raft}  
echo ${ccd}  
echo ${input_path}  
echo ${output_path} 
#python test_PCA.py 
#python test_PCA.py ${run} ${raft} ${ccd} ${input_path} ${output_path} 
python plot_images.py ${run} ${raft} ${ccd} ${input_path} ${output_path} 1D
echo "DONE" 
