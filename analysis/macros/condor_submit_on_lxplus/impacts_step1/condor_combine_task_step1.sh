#!/bin/sh
ulimit -s unlimited
set -e
cd /afs/cern.ch/work/j/${USER}/CMSSW_10_2_13/src
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/j/${USER}/CMSSW_10_2_13/src

echo "cp ${3} ./"
cp ${3} ./

echo "python CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ${1}.root -m 125 --doInitialFit --robustFit 1 -t -1 --setParameters ${2}_r=1 -n R1_${2}"
python CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ${1}.root -m 125 --doInitialFit --robustFit 1 -t -1 --setParameters ${2}_r=1 -n R1_${2}

ls -ltrh
