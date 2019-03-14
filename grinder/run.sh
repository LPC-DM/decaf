#!/usr/bin/env bash
export USER=matteoc
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
xrdcp -s root://cmseos.fnal.gov//store/user/$USER/decaf.tgz .
echo "Decaf correctly copied"
xrdcp -s root://cmseos.fnal.gov//store/user/$USER/pylocal.tgz .
echo "Python correctly copied" 
tar -zxvf decaf.tgz
tar -zxvf pylocal.tgz
rm decaf.tgz
rm pylocal.tgz
cd decaf
source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc62-opt/setup.sh
#wget https://github.com/xrootd/xrootd/archive/v4.8.4.tar.gz
#tar zxf v4.8.4.tar.gz && rm v4.8.4.tar.gz
#cp xrd_setup.py xrootd-4.8.4/bindings/python/
#cd xrootd-4.8.4/bindings/python/
##python xrd_setup.py install --install-option="--prefix=${_CONDOR_SCRATCH_DIR}/uscms/home/matteoc/.local/lib/python3.6/site-packages"
#python xrd_setup.py install --user
#cd ${_CONDOR_SCRATCH_DIR}/decaf
#rm -rf xrootd-4.8.4
export PATH=${_CONDOR_SCRATCH_DIR}/uscms/home/matteoc/.local/bin:$PATH
export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/uscms/home/matteoc/.local/lib/python3.6/site-packages:$PYTHONPATH
export PYTHONPATH=$(find ${_CONDOR_SCRATCH_DIR}/uscms/home/matteoc/.local/lib/python3.6/site-packages/ -name *.egg |tr '\n' ':')$PYTHONPATH
echo "Updated python path: " $PYTHONPATH
#source setup_lcg.sh  ## if a bash script, use .sh instead of .csh
cd grinder
echo "year: " $year
echo "lumi: " $lumi
echo "selection: " $selection
echo "dataset: " $sample
echo "python run.py --year $year --lumi $lumi --selection $selection --dataset $sample"
python run.py --year $year --lumi $lumi --selection $selection --dataset $sample
ls ../pods/2017
#cd ${_CONDOR_SCRATCH_DIR}
#rm -rf decaf
