#!/usr/bin/bash

#to debug
#export input_path_batch=/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/data_PTC
#general case
export input_path_batch=/sps/lsst/groups/FocalPlane/SLAC/run5/
export output_path_batch=/sps/lsst/users/tguillem/web/batch/dark/ccd_assembly_fix_2/v4_1D/
export log_path=/sps/lsst/users/tguillem/batch_logs

#full focal plane
#runs=('13159' '13144' '13162' '13161' '13018')
runs=('13161')
rafts=('R01' 'R02' 'R03' 'R10' 'R11' 'R12' 'R13' 'R14' 'R20' 'R21' 'R22' 'R23' 'R24' 'R30' 'R31' 'R32' 'R33' 'R34' 'R41' 'R42' 'R43')
ccds=('S00' 'S01' 'S02' 'S10' 'S11' 'S12' 'S20' 'S21' 'S22')
#exposures=('3021120700078' '3021120700333')
#DOES NOT WORK: exposures=\($(cat /sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/file.txt)\)
exposures=('3021121200145' '3021121200137' '3021121200141' '3021121200150' '3021121200151' '3021121200156' '3021121200152' '3021121200148' '3021121200139' '3021121200138' '3021121200154' '3021121200149' '3021121200153' '3021121200147' '3021121200140' '3021121200146' '3021121200143' '3021121200155' '3021121200142' '3021121200144')
#just one job
rafts=('R01')
ccds=('S00')
for run in "${runs[@]}" ;do
    for raft in "${rafts[@]}" ;do
	for ccd in "${ccds[@]}" ;do
	    for exposure in "${exposures[@]}" ;do
		echo "---Run: $run"
		echo "Raft: $raft"
		echo "CCD: $ccd"
		echo "Exposure: $exposure"
	        #qsub -P P_lsst -o ${log_path} -e ${log_path} -v run=${run} -v raft=${raft} -v ccd=${ccd} -v input_path=${input_path_batch} -v output_path=${output_path_batch} run_batch.sh
	        #--mem=64G
	        #sbatch --mem=64G --export=run=${run},raft=${raft},ccd=${ccd},input_path=${input_path_batch},output_path=${output_path_batch} run_batch.sh 
		sbatch -t 12:00:00 --cpus-per-task=8 --mem=16G --export=run=${run},raft=${raft},ccd=${ccd},exposure=${exposure} run_batch_exposure.sh
	        #sleep 2s
	    done
	done	
    done
done
#qsub -P P_lsst -v run=${run_batch} -v raft='R14' -v ccd='S22' -v input_path=${input_path_batch} -v output_path=${output_path_batch} run_batch.sh
