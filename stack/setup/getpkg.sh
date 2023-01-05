export version=w.2022.52

#cp_pipe
git clone https://github.com/lsst/cp_pipe.git
cd cp_pipe
#git tag
git checkout version
setup -k -r .
scons -j 4
eups declare -r . -t tguillem
cd ..
setup ip_isr -t tguillem

#ip_isr
git clone https://github.com/lsst/ip_isr.git
cd ip_isr
#git tag
git checkout version
setup -k -r .
scons -j 4
eups declare -r . -t tguillem
cd ..
setup ip_isr -t tguillem
