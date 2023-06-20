#!/usr/bin/zsh

#cd /sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3

#WARNING: LSST must not be already set up!
echo "LSST load"
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2022_52/loadLSST.bash  

echo "LSST setup"
setup lsst_distrib 

echo "Starting butler ingestion"
#echo ${run}
#rm -rf all_runs/${run}
#mkdir -p all_runs/${run}
#cp /sps/lsst/groups/FocalPlane/SLAC/run5/butler/all_runs/${run}_list.txt all_runs/${run}/.
#cd all_runs/${run}

#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/all_runs/${run}
#export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/all_runs/13162
#butler create $REPO 
#butler register-instrument $REPO lsst.obs.lsst.LsstCam   

#mix 
#butler ingest-raws --transfer symlink -j 8 $REPO ${run}_list.txt
#butler ingest-raws --transfer symlink -j 8 $REPO 13151_list.txt
#butler ingest-raws --transfer symlink -j 8 $REPO `cat /sps/lsst/groups/FocalPlane/SLAC/run5/butler/all_runs/13151/13151_list.txt` 
#butler ingest-raws --transfer symlink -j 8 $REPO `cat 13151_list.txt`
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/13151
#butler ingest-raws --transfer direct -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/13151
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/storage/20211209/

#NEW METHOD
export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111
find -L /sps/lsst/groups/FocalPlane/SLAC/run5/${run} -name "*.fits" | xargs butler ingest-raws --transfer symlink -j 8 $REPO

butler write-curated-calibrations $REPO lsst.obs.lsst.LsstCam

echo "DONE"
