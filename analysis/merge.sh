#!/usr/bin/env bash

export USER=${4}
echo "User is: ${4}"
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
echo $(hostname)

if [ "${3}" == "kisti" ]; then
    env
    /usr/bin/voms-proxy-info -exists
    if [ $? -eq 0 ]; then
	echo "No need to copy"
	ls -l /tmp/x509up_u$(id -u)
	/usr/bin/voms-proxy-info -all
    else
	cp ./x509up_u* /tmp
	ls -l /tmp/x509up_u$(id -u)
	/usr/bin/voms-proxy-info -all
    fi
    xrdcp -s root://cms-xrdr.private.lo:2094//xrd/store/user/$USER/decaf.tgz .
    echo "Decaf correctly copied"
    xrdcp -s root://cms-xrdr.private.lo:2094//xrd/store/user/$USER/pylocal.tgz .
    echo "Python correctly copied" 
else
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/decaf.tgz .
    echo "Decaf correctly copied"
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/pylocal.tgz .
    echo "Python correctly copied"
fi
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
ls
echo "python merge.py --folder ${1} --variable ${2}"
python merge.py --folder ${1} --variable ${2}
ls ${1}/${2}.merged
cp ${1}/${2}.merged ${_CONDOR_SCRATCH_DIR}/${2}.merged
