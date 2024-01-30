export stamp=w_2023_49_20240129a
#master bias and master dark
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_bias_0D.yaml w_2023_49
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_bias_1D.yaml w_2023_49
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_bias_2D.yaml w_2023_49
export EXPOSURES='2023111600660, 2023111600661, 2023111600662, 2023111600663, 2023111600664, 2023111600665, 2023111600666, 2023111600667, 2023111600668, 2023111600669, 2023111600670, 2023111600671, 2023111600672, 2023111600673, 2023111600674, 2023111600675, 2023111600676, 2023111600677, 2023111600678, 2023111600679'
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_dark15_0D.yaml w_2023_49
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_dark15_1D.yaml w_2023_49
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_dark15_2D.yaml w_2023_49
export EXPOSURES='2023111600760, 2023111600761, 2023111600762, 2023111600763, 2023111600764, 2023111600765, 2023111600766, 2023111600767, 2023111600768, 2023111600769, 2023111600770, 2023111600771, 2023111600772, 2023111600773, 2023111600774, 2023111600775, 2023111600776, 2023111600777, 2023111600778, 2023111600779'
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_dark30_0D.yaml w_2023_49
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_dark30_1D.yaml w_2023_49
sbatch /pbs/throng/lsst/software/parsl/tools/bps_submit.sh bps_master_dark30_2D.yaml w_2023_49
