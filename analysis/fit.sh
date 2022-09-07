#!/usr/bin/env bash
export USER=${6}
echo "User is: ${6}"
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
echo $(hostname)
source /cvmfs/cms.cern.ch/cmsset_default.sh

if [ "${5}" == "kisti" ]; then
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
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/cmssw.tgz .
    echo "cmssw correctly copied"
fi
echo "untar cmssw"
tar -zxvf cmssw.tgz
echo "cmssw untarred"
rm cmssw.tgz
export SCRAM_ARCH=slc7_amd64_gcc700
cd CMSSW_10_2_13/src
scramv1 b ProjectRename
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
cd decaf/analysis
if [ "${3}" == "None"  ]; then
    echo "python fit.py -M ${2} -w ${1}"
    python fit.py -M ${2} -w ${1}
else
    export spaces=$( echo ${3} | tr '+' ' ' )
    export arguments=$( echo $spaces | tr 'X' '"')
    echo "python fit.py -M ${2} -w ${1} -a '$arguments'"
    python fit.py -M ${2} -w ${1} -a '$arguments'
fi
ls ${2}
tar -czvf ${_CONDOR_SCRATCH_DIR}/${4}.tgz results/${4}/*${2}*
