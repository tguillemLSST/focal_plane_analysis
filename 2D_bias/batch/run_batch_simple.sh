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
echo "Starting to run job"
#python read_bias_stability.py
source debug.sh
echo "DONE"
