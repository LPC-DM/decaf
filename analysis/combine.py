#!/usr/bin/env python
import os
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', '--signal', help='signal', dest='signal')
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    (options, args) = parser.parse_args()
    
    rootfiles:[]
    for filename in os.listdir('datacards/'+options.folder):
        if '.root' not in filename: continue
        if folder not in filename: continue
        if options.signal:
            if not any(_signal in filename for _signal in options.signal.split(',')): continue
        else:
            try:
                if filename.split('_',1)[1]
                    continue
            except:
                pass
        
        rootfiles.append(filename)
        
    print(filename)
          

'''
if [ "${4}" == "limit" ]; then

    echo ""
    echo "Run AsymptoticLimits"
    #echo "command: combine -M AsymptoticLimits ${1}.root --cminDefaultMinimizerStrategy 0 --rMax +1.6 -v 3"
    #combine -M AsymptoticLimits ${1}.root --cminDefaultMinimizerStrategy 0 --rMax +1.6 -v 3
    echo "combine -M AsymptoticLimits -d darkhiggs.root --cminDefaultMinimizerStrategy 0 --setParameters 'rgx{^(?!${5}.)Mz.*$}=0' --setParameterRanges 'rgx{^(?!${5}.)Mz.*$}'=0,0 --freezeParameters var{'^(?!${5}.)Mz.*$'} -v 3"
    combine -M AsymptoticLimits -d darkhiggs.root --cminDefaultMinimizerStrategy 0 --setParameters "rgx{^(?!${5}.)Mz.*$}=0" --setParameterRanges "rgx{^(?!${5}.)Mz.*$}"=0,0 --freezeParameters var{"^(?!${5}.)Mz.*$"} -v 3

    if [ -f "higgsCombineTest.AsymptoticLimits.mH120.root" ]; then
        ls -l higgsCombineTest.AsymptoticLimits.mH120.root 
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp higgsCombineTest.AsymptoticLimits.mH120.root ${_CONDOR_SCRATCH_DIR}/AsymptoticLimits.result.${5}.root
    fi

elif [ "${4}" == "shape" ]; then

    echo ""
    echo "Run PostFitShapesFromWorkspace"
    echo "command: PostFitShapesFromWorkspace -w ${1}.root -d ${1}.txt -f fitDiagnostics.cminresult.${1}.root:fit_s --postfit --sampling --samples 300 --skip-proc-errs -o postfitshapes.result.${1}.root"
     PostFitShapesFromWorkspace -w ${1}.root -d ${1}.txt -f fitDiagnostics.cminresult.${1}.root:fit_s --postfit --sampling --samples 300 --skip-proc-errs -o postfitshapes.result.${1}.root

    if [ -f "postfitshapes.result.${1}.root" ]; then
        ls -l postfitshapes.result.${1}.root 
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp postfitshapes.result.${1}.root ${_CONDOR_SCRATCH_DIR}/postfitshapes.result.${1}.root
    fi

elif [ "${4}" == "fit" ]; then
    
    echo ""
    echo "Run FitDiagnostics"
    echo "command: combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d ${1}.root -v 3"
    combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d ${1}.root -v 3 
    #combine -M FitDiagnostics -d ${1}.root 

    if [ -f "fitDiagnosticsTest.root" ]; then
        ls -l fitDiagnosticsTest.root
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.result.${1}.root
    fi

elif [ "${4}" == "cminfit" ]; then
    
    echo ""
    echo "Run FitDiagnostics including --cminDefaultMinimizerStrategy 0"
    #echo "command: combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d ${1}.root --cminDefaultMinimizerStrategy 0 -v 3"
    #combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d ${1}.root --cminDefaultMinimizerStrategy 0 -v 3
    echo "combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d darkhiggs.root --cminDefaultMinimizerStrategy 0 --setParameters 'rgx{^(?!${5}.)Mz.*$}=0' --setParameterRanges 'rgx{^(?!${5}.)Mz.*$}'=0,0 --freezeParameters var{'^(?!${5}.)Mz.*$'} -v 3"
    combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d darkhiggs.root --cminDefaultMinimizerStrategy 0 --setParameters "rgx{^(?!${5}.)Mz.*$}=0" --setParameterRanges "rgx{^(?!${5}.)Mz.*$}"=0,0 --freezeParameters var{"^(?!${5}.)Mz.*$"} -v 3

    if [ -f "fitDiagnosticsTest.root" ]; then
        ls -l fitDiagnosticsTest.root
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.cminresult.${5}.root
    fi

else

    echo "combine -M FitDiagnostics -d ${1}.root --expectSignal 0 --forceRecreateNLL --cminDefaultMinimizerType Minuit --ignoreCovWarning --saveWithUncertainties -t -1"
    combine -M FitDiagnostics -d ${1}.root --expectSignal 0 --forceRecreateNLL --cminDefaultMinimizerType Minuit --rMin 0 --rMax 2 --ignoreCovWarning --saveWithUncertainties -t -1
    
    if [ -f "fitDiagnosticsTest.root" ]; then
        ls -l fitDiagnosticsTest.root
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.result2.${1}.root
    fi
