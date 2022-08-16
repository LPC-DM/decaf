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
    xrdcp -s root://cmseos.fnal.gov//store/user/$USER/cmssw.tgz .
    echo "cmssw correctly copied"
fi
tar -zxvf cmssw.tgz
rm cmssw.tgz
export SCRAM_ARCH=slc7_amd64_gcc700
cd CMSSW_10_2_13/src
scramv1 b ProjectRename
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
cd decaf/datacards/${1}
#echo "text2workspace.py ${1}.txt --channel-masks"
echo "text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO map=.*/Mz1000_mhs50_Mdm1000:Mz1000_mhs50_Mdm1000_r[1,0,2] --PO map=.*/Mz1000_mhs50_Mdm150:Mz1000_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz1000_mhs50_Mdm500:Mz1000_mhs50_Mdm500_r[1,0,2] --PO map=.*/Mz1000_mhs70_Mdm1000:Mz1000_mhs70_Mdm1000_r[1,0,2] --PO map=.*/Mz1000_mhs70_Mdm150:Mz1000_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz1000_mhs70_Mdm500:Mz1000_mhs70_Mdm500_r[1,0,2] --PO map=.*/Mz1000_mhs90_Mdm1000:Mz1000_mhs90_Mdm1000_r[1,0,2] --PO map=.*/Mz1000_mhs90_Mdm150:Mz1000_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz1000_mhs90_Mdm500:Mz1000_mhs90_Mdm500_r[1,0,2] --PO map=.*/Mz2000_mhs50_Mdm1000:Mz2000_mhs50_Mdm1000_r[1,0,2] --PO map=.*/Mz2000_mhs50_Mdm1500:Mz2000_mhs50_Mdm1500_r[1,0,2] --PO map=.*/Mz2000_mhs50_Mdm500:Mz2000_mhs50_Mdm500_r[1,0,2] --PO map=.*/Mz2000_mhs70_Mdm1000:Mz2000_mhs70_Mdm1000_r[1,0,2] --PO map=.*/Mz2000_mhs70_Mdm1500:Mz2000_mhs70_Mdm1500_r[1,0,2] --PO map=.*/Mz2000_mhs70_Mdm500:Mz2000_mhs70_Mdm500_r[1,0,2] --PO map=.*/Mz2000_mhs90_Mdm1000:Mz2000_mhs90_Mdm1000_r[1,0,2] --PO map=.*/Mz2000_mhs90_Mdm1500:Mz2000_mhs90_Mdm1500_r[1,0,2] --PO map=.*/Mz2000_mhs90_Mdm500:Mz2000_mhs90_Mdm500_r[1,0,2] --PO map=.*/Mz200_mhs50_Mdm100:Mz200_mhs50_Mdm100_r[1,0,2] --PO map=.*/Mz200_mhs50_Mdm150:Mz200_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz200_mhs70_Mdm100:Mz200_mhs70_Mdm100_r[1,0,2] --PO map=.*/Mz200_mhs70_Mdm150:Mz200_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz200_mhs90_Mdm100:Mz200_mhs90_Mdm100_r[1,0,2] --PO map=.*/Mz200_mhs90_Mdm150:Mz200_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz2500_mhs50_Mdm1250:Mz2500_mhs50_Mdm1250_r[1,0,2] --PO map=.*/Mz2500_mhs50_Mdm750:Mz2500_mhs50_Mdm750_r[1,0,2] --PO map=.*/Mz2500_mhs70_Mdm1250:Mz2500_mhs70_Mdm1250_r[1,0,2] --PO map=.*/Mz2500_mhs70_Mdm750:Mz2500_mhs70_Mdm750_r[1,0,2] --PO map=.*/Mz2500_mhs90_Mdm1250:Mz2500_mhs90_Mdm1250_r[1,0,2] --PO map=.*/Mz2500_mhs90_Mdm750:Mz2500_mhs90_Mdm750_r[1,0,2] --PO map=.*/Mz3000_mhs50_Mdm1000:Mz3000_mhs50_Mdm1000_r[1,0,2] --PO map=.*/Mz3000_mhs50_Mdm1500:Mz3000_mhs50_Mdm1500_r[1,0,2] --PO map=.*/Mz3000_mhs70_Mdm1000:Mz3000_mhs70_Mdm1000_r[1,0,2] --PO map=.*/Mz3000_mhs70_Mdm1500:Mz3000_mhs70_Mdm1500_r[1,0,2] --PO map=.*/Mz3000_mhs90_Mdm1000:Mz3000_mhs90_Mdm1000_r[1,0,2] --PO map=.*/Mz3000_mhs90_Mdm1500:Mz3000_mhs90_Mdm1500_r[1,0,2] --PO map=.*/Mz300_mhs50_Mdm100:Mz300_mhs50_Mdm100_r[1,0,2] --PO map=.*/Mz300_mhs50_Mdm150:Mz300_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz300_mhs70_Mdm100:Mz300_mhs70_Mdm100_r[1,0,2] --PO map=.*/Mz300_mhs70_Mdm150:Mz300_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz300_mhs90_Mdm100:Mz300_mhs90_Mdm100_r[1,0,2] --PO map=.*/Mz300_mhs90_Mdm150:Mz300_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs50_Mdm150:Mz500_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs50_Mdm250:Mz500_mhs50_Mdm250_r[1,0,2] --PO map=.*/Mz500_mhs50_Mdm500:Mz500_mhs50_Mdm500_r[1,0,2] --PO map=.*/Mz500_mhs70_Mdm150:Mz500_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs70_Mdm250:Mz500_mhs70_Mdm250_r[1,0,2] --PO map=.*/Mz500_mhs70_Mdm500:Mz500_mhs70_Mdm500_r[1,0,2] --PO map=.*/Mz500_mhs90_Mdm150:Mz500_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs90_Mdm250:Mz500_mhs90_Mdm250_r[1,0,2] --PO map=.*/Mz500_mhs90_Mdm500:Mz500_mhs90_Mdm500_r[1,0,2] ${1}.txt -o ${1}.root"

text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO map=.*/Mz1000_mhs50_Mdm1000:Mz1000_mhs50_Mdm1000_r[1,0,2] --PO map=.*/Mz1000_mhs50_Mdm150:Mz1000_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz1000_mhs50_Mdm500:Mz1000_mhs50_Mdm500_r[1,0,2] --PO map=.*/Mz1000_mhs70_Mdm1000:Mz1000_mhs70_Mdm1000_r[1,0,2] --PO map=.*/Mz1000_mhs70_Mdm150:Mz1000_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz1000_mhs70_Mdm500:Mz1000_mhs70_Mdm500_r[1,0,2] --PO map=.*/Mz1000_mhs90_Mdm1000:Mz1000_mhs90_Mdm1000_r[1,0,2] --PO map=.*/Mz1000_mhs90_Mdm150:Mz1000_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz1000_mhs90_Mdm500:Mz1000_mhs90_Mdm500_r[1,0,2] --PO map=.*/Mz2000_mhs50_Mdm1000:Mz2000_mhs50_Mdm1000_r[1,0,2] --PO map=.*/Mz2000_mhs50_Mdm1500:Mz2000_mhs50_Mdm1500_r[1,0,2] --PO map=.*/Mz2000_mhs50_Mdm500:Mz2000_mhs50_Mdm500_r[1,0,2] --PO map=.*/Mz2000_mhs70_Mdm1000:Mz2000_mhs70_Mdm1000_r[1,0,2] --PO map=.*/Mz2000_mhs70_Mdm1500:Mz2000_mhs70_Mdm1500_r[1,0,2] --PO map=.*/Mz2000_mhs70_Mdm500:Mz2000_mhs70_Mdm500_r[1,0,2] --PO map=.*/Mz2000_mhs90_Mdm1000:Mz2000_mhs90_Mdm1000_r[1,0,2] --PO map=.*/Mz2000_mhs90_Mdm1500:Mz2000_mhs90_Mdm1500_r[1,0,2] --PO map=.*/Mz2000_mhs90_Mdm500:Mz2000_mhs90_Mdm500_r[1,0,2] --PO map=.*/Mz200_mhs50_Mdm100:Mz200_mhs50_Mdm100_r[1,0,2] --PO map=.*/Mz200_mhs50_Mdm150:Mz200_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz200_mhs70_Mdm100:Mz200_mhs70_Mdm100_r[1,0,2] --PO map=.*/Mz200_mhs70_Mdm150:Mz200_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz200_mhs90_Mdm100:Mz200_mhs90_Mdm100_r[1,0,2] --PO map=.*/Mz200_mhs90_Mdm150:Mz200_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz2500_mhs50_Mdm1250:Mz2500_mhs50_Mdm1250_r[1,0,2] --PO map=.*/Mz2500_mhs50_Mdm750:Mz2500_mhs50_Mdm750_r[1,0,2] --PO map=.*/Mz2500_mhs70_Mdm1250:Mz2500_mhs70_Mdm1250_r[1,0,2] --PO map=.*/Mz2500_mhs70_Mdm750:Mz2500_mhs70_Mdm750_r[1,0,2] --PO map=.*/Mz2500_mhs90_Mdm1250:Mz2500_mhs90_Mdm1250_r[1,0,2] --PO map=.*/Mz2500_mhs90_Mdm750:Mz2500_mhs90_Mdm750_r[1,0,2] --PO map=.*/Mz3000_mhs50_Mdm1000:Mz3000_mhs50_Mdm1000_r[1,0,2] --PO map=.*/Mz3000_mhs50_Mdm1500:Mz3000_mhs50_Mdm1500_r[1,0,2] --PO map=.*/Mz3000_mhs70_Mdm1000:Mz3000_mhs70_Mdm1000_r[1,0,2] --PO map=.*/Mz3000_mhs70_Mdm1500:Mz3000_mhs70_Mdm1500_r[1,0,2] --PO map=.*/Mz3000_mhs90_Mdm1000:Mz3000_mhs90_Mdm1000_r[1,0,2] --PO map=.*/Mz3000_mhs90_Mdm1500:Mz3000_mhs90_Mdm1500_r[1,0,2] --PO map=.*/Mz300_mhs50_Mdm100:Mz300_mhs50_Mdm100_r[1,0,2] --PO map=.*/Mz300_mhs50_Mdm150:Mz300_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz300_mhs70_Mdm100:Mz300_mhs70_Mdm100_r[1,0,2] --PO map=.*/Mz300_mhs70_Mdm150:Mz300_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz300_mhs90_Mdm100:Mz300_mhs90_Mdm100_r[1,0,2] --PO map=.*/Mz300_mhs90_Mdm150:Mz300_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs50_Mdm150:Mz500_mhs50_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs50_Mdm250:Mz500_mhs50_Mdm250_r[1,0,2] --PO map=.*/Mz500_mhs50_Mdm500:Mz500_mhs50_Mdm500_r[1,0,2] --PO map=.*/Mz500_mhs70_Mdm150:Mz500_mhs70_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs70_Mdm250:Mz500_mhs70_Mdm250_r[1,0,2] --PO map=.*/Mz500_mhs70_Mdm500:Mz500_mhs70_Mdm500_r[1,0,2] --PO map=.*/Mz500_mhs90_Mdm150:Mz500_mhs90_Mdm150_r[1,0,2] --PO map=.*/Mz500_mhs90_Mdm250:Mz500_mhs90_Mdm250_r[1,0,2] --PO map=.*/Mz500_mhs90_Mdm500:Mz500_mhs90_Mdm500_r[1,0,2] ${1}.txt -o ${1}.root

ls ${1}.root
cp ${1}.root ${_CONDOR_SCRATCH_DIR}/${1}.root
