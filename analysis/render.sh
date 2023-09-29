#!/usr/bin/env bash
export USER=${3}
echo "User is: ${3}"
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
echo $(hostname)
source /cvmfs/cms.cern.ch/cmsset_default.sh

if [ "${2}" == "kisti" ]; then
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
    xrdcp -s root://cms-xrdr.private.lo:2094//xrd/store/user/$USER/cmssw.tgz .
    echo "cmssw correctly copied"
else
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/cmssw.tgz .
    echo "cmssw correctly copied"
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/py2local.tgz .
    echo "py2local correctly copied"
    tar -zxvf py2local.tgz
    rm py2local.tgz
    export PYTHONPATH=${_CONDOR_SCRATCH_DIR}/site-packages:$PYTHONPATH
    export PYTHONPATH=$(find ${_CONDOR_SCRATCH_DIR}/site-packages/ -name *.egg |tr '\n' ':')$PYTHONPATH
    export PYTHONWARNINGS="ignore"
    echo "Updated python path: " $PYTHONPATH
fi
tar -zxvf cmssw.tgz
rm cmssw.tgz
setenv SCRAM_ARCH slc7_amd64_gcc700
cd CMSSW_10_2_13/src
scramv1 b ProjectRename
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
cd decaf/analysis
echo "python render.py --model ${1}"
python render.py --model ${1}
ls datacards/${1}/*
tar -czvf ${_CONDOR_SCRATCH_DIR}/${1}.tgz datacards/${1}/*

