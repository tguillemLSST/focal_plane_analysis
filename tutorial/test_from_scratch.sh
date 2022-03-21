mkdir /tmp/focal_plane_210324
cd /tmp/focal_plane_210324
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2022_01/loadLSST.zsh 
setup lsst_distrib
git clone https://github.com/tguillemLSST/focal_plane_analysis.git
cd focal_plane_analysis/tutorial
python play_with_images.py
