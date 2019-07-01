#!/usr/bin/env bash
export USER=${3}
echo "User is: ${3}"
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
export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/site-packages:$PYTHONPATH
export PYTHONPATH=$(find ${_CONDOR_SCRATCH_DIR}/site-packages/ -name *.egg |tr '\n' ':')$PYTHONPATH
echo "Updated python path: " $PYTHONPATH
cd grinder
echo "python run.py --year ${1} --dataset ${2}"
python run.py --year ${1} --dataset ${2}
ls pods/${1}/${2}.pkl.gz
cp pods/${1}/${2}.pkl.gz ${_CONDOR_SCRATCH_DIR}/${1}_${2}.pkl.gz
