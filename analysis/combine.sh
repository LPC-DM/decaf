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
setenv SCRAM_ARCH slc7_amd64_gcc700
cd CMSSW_10_2_13/src
scramv1 b ProjectRename
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
cd decaf/analysis/datacards/${1}
rm fitDiagnostics.root

############################ FitDiagnostics #########################################
echo ""
echo "First step: Generate toy sample"
echo "command:"
echo "combine -M GenerateOnly -t -1 --saveToys --toysFrequentist --bypassFrequentistFit --expectSignal 0 -n ${1} -d ${1}.root"
combine -M GenerateOnly -t -1 --saveToys --toysFrequentist --bypassFrequentistFit --expectSignal 0 -n ${1} -d ${1}.root

if [ -f "higgsCombine${1}.GenerateOnly.mH120.123456.root" ]; then
    ls -l higgsCombine${1}.GenerateOnly.mH120.123456.root
    echo "The toy sample will be copied to ${_CONDOR_SCRATCH_DIR}"
    cp higgsCombine${1}.GenerateOnly.mH120.123456.root ${_CONDOR_SCRATCH_DIR}/GenerateOnly.${1}.root
fi

echo ""
echo "Second step: Run FitDiagnostics"
echo "command:"
echo "combine -M FitDiagnostics -t -1 --rMin -1.5 --rMax +1.6 --saveWorkspace --toysFile=higgsCombine${1}.GenerateOnly.mH120.123456.root -d ${1}.root"
combine -M FitDiagnostics -t -1 --rMin -1.5 --rMax +1.6 --saveWorkspace --toysFile=higgsCombine${1}.GenerateOnly.mH120.123456.root -d ${1}.root

if [ -f "./fitDiagnosticsTest.root" ]; then
    ls -l fitDiagnosticsTest.root
    echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
    cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.root
fi
#####################################################################################

############################### AsymptoticLimits ####################################
#echo ""
#echo "Run AsymptoticLimits"
#echo "command:"
#echo "combine -M AsymptoticLimits ${1}.txt"
#combine -M AsymptoticLimits ${1}.txt
#
#if [ -f "higgsCombineTest.AsymptoticLimits.mH120.root" ]; then
#    ls -l higgsCombineTest.AsymptoticLimits.mH120.root 
#    echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
#    cp higgsCombineTest.AsymptoticLimits.mH120.root ${_CONDOR_SCRATCH_DIR}/AsymptoticLimits.${1}.root
#fi
#####################################################################################

############################ FitDiagnostics #########################################
#echo "combine -M FitDiagnostics -d ${1}.root --expectSignal 0 --forceRecreateNLL --cminDefaultMinimizerType Minuit --ignoreCovWarning --saveWithUncertainties -t -1"
#combine -M FitDiagnostics -d ${1}.root --expectSignal 0 --forceRecreateNLL --cminDefaultMinimizerType Minuit --cminDefaultMinimizerStrategy 0 --rMin 0 --rMax 2 --ignoreCovWarning --saveWithUncertainties -t -1
#
#### Check the combine fit succeeded
#if [ $? -eq 0 ]; then
#    echo ""
#    echo "==================================="
#    echo "The combine fit has been succeeded"
#    echo "==================================="
#
#    if [ -f "./fitDiagnosticsTest.root" ]; then
#        ls -l fitDiagnosticsTest.root 
#        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
#        cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.root
#    else
#        echo ""
#        echo "==================================="
#        echo "However, the fitDiagonsticsTest.root is not exist!"
#        echo "==================================="
#    fi
#else
#    echo ""
#    echo "==================================="
#    echo "The combine fit has been failed"
#    echo "==================================="
#    echo ""
#fi
#####################################################################################
