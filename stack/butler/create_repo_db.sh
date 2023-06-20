export REPO=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/gen3/main_20230111
#mkdir -p $REPO

#butler db configuration
#cat $REPO/butler-seed.yaml
#registry:
#    db: "postgresql://ccpglsstdev.in2p3.fr:6553/messier"
#    namespace: "main_20230111"

#butler create --seed-config $REPO/butler-seed.yaml --override $REPO

#instrument registration
butler register-instrument $REPO 'lsst.obs.lsst.LsstCam'

#test ingestion
butler ingest-raws --transfer symlink -j 8 $REPO /sps/lsst/groups/FocalPlane/SLAC/run5/13162/bias\
*/*R14*S22*.fits
