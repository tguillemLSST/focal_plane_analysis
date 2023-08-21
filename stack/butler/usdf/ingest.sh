export REPO=/sps/lsst/groups/FocalPlane/SLAC/run6/butler/main/

#bias
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/13162/bias_bias*/*R14*S22*.fits
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/storage/run6/LSSTCam/20230601/MC_C_*/*.fits

#debug
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/storage/20211209/MC_C_20211209_001704/MC_C_20211209_001704_R21_S22.fits
#butler ingest-raws --transfer symlink -j 8 $REPO all_runs/13161_short_list.txt

#dark
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/13161/dark_dark_*/*R14*S22*.fits

#bias
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/13161/bias*/*R14*S22*.fits 

#mix
#butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/13161/*/*R14*S22*.fits

#butler write-curated-calibrations $REPO lsst.obs.lsst.LsstCam

#copy of user collections from usdf
butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run6/butler/import/250523/copy_defects_1339*_w_2023_24/u/lsstccs/defects_1339*_w_2023_24/*/defects/defects_LSSTCam_R*.f
