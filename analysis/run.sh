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
if uname -r | grep -q el6; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-slc6-gcc8-opt/setup.sh
else
  source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
fi
export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/site-packages:$PYTHONPATH
export PYTHONPATH=$(find ${_CONDOR_SCRATCH_DIR}/site-packages/ -name *.egg |tr '\n' ':')$PYTHONPATH
export PYTHONWARNINGS="ignore"
echo "Updated python path: " $PYTHONPATH
cd analysis
echo "python run.py --year ${1} --dataset ${2} --analysis ${4}"
python run.py --year ${1} --dataset ${2} --analysis ${4}
ls hists/${4}${1}/${2}.futures
cp hists/${4}${1}/${2}.futures ${_CONDOR_SCRATCH_DIR}/${4}${1}_${2}.futures
