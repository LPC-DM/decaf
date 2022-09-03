#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', '--signal', help='signal', dest='signal')
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-o', '--option', help='option', dest='option', default='')
    parser.add_option('--multiSignal', action='store_true', dest='multiSignal')
    parser.add_option('--AsymptoticLimits', action='store_true', dest='AsymptoticLimits')
    parser.add_option('--FitDiagnostics', action='store_true', dest='FitDiagnostics')
    parser.add_option('--PostFitShapesFromWorkspace', action='store_true', dest='PostFitShapesFromWorkspace')
    (options, args) = parser.parse_args()
    
    signals=[]
    for k,v in processes.items():
        process = k
        if not isinstance(k, str):
            process = k[0]
        #if options.signal.split(':')[0] not in process: continue
        if not any(_signal in process for _signal in options.signal.split(',')): continue
        if process not in signals:
            signals.append(process)
    print(signals)

    rootfiles=[]
    for filename in os.listdir('datacards/'+options.folder):
        if options.multiSignal and options.folder+'.root' not in filename: continue
        if not any(signal in filename for signal in signals): continue
        rootfiles.append(filename)
    print(rootfiles)

    commands=[]
    if (options.AsymptoticLimits): 
        method='combine -M AsymptoticLimits -d '
        standard_options='--cminDefaultMinimizerStrategy 0 -v 3 '


    for signal in signals:
        special_options=options.option.replace('SIGNAL',signal)
        if options.multiSignal:
            command=method+'datacards/'+options.folder+'/'+options.folder+'.root '+standard_options+special_options
            if command not in commands: commands.append(command)
        else:
            for datacard in rootfiles:
                if signal not in datacard: continue
                command=method+'datacards/'+options.folder+'/'+datacard+' '+standard_options+special_options
                if command not in commands: commands.append(command)



    print(len(commands))
    for command in commands:
        print(command)
'''
if [ "${4}" == "limit" ]; then
    combine -M AsymptoticLimits -d darkhiggs.root --cminDefaultMinimizerStrategy 0 --setParameters "rgx{^(?!${5}.)Mz.*$}=0" --setParameterRanges "rgx{^(?!${5}.)Mz.*$}"=0,0 --freezeParameters var{"^(?!${5}.)Mz.*$"} -v 3

    if [ -f "higgsCombineTest.AsymptoticLimits.mH120.root" ]; then
        ls -l higgsCombineTest.AsymptoticLimits.mH120.root 
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp higgsCombineTest.AsymptoticLimits.mH120.root ${_CONDOR_SCRATCH_DIR}/AsymptoticLimits.result.${5}.root
    fi

elif [ "${4}" == "shape" ]; then

     PostFitShapesFromWorkspace -w ${1}.root -d ${1}.txt -f fitDiagnostics.cminresult.${1}.root:fit_s --postfit --sampling --samples 300 --skip-proc-errs -o postfitshapes.result.${1}.root

    if [ -f "postfitshapes.result.${1}.root" ]; then
        ls -l postfitshapes.result.${1}.root 
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp postfitshapes.result.${1}.root ${_CONDOR_SCRATCH_DIR}/postfitshapes.result.${1}.root
    fi

elif [ "${4}" == "fit" ]; then

    combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d ${1}.root -v 3 

    if [ -f "fitDiagnosticsTest.root" ]; then
        ls -l fitDiagnosticsTest.root
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.result.${1}.root
    fi

elif [ "${4}" == "cminfit" ]; then
    
    combine -M FitDiagnostics --rMax +1.6 --saveWorkspace -d darkhiggs.root --cminDefaultMinimizerStrategy 0 --setParameters "rgx{^(?!${5}.)Mz.*$}=0" --setParameterRanges "rgx{^(?!${5}.)Mz.*$}"=0,0 --freezeParameters var{"^(?!${5}.)Mz.*$"} -v 3

    if [ -f "fitDiagnosticsTest.root" ]; then
        ls -l fitDiagnosticsTest.root
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.cminresult.${5}.root
    fi

else

    combine -M FitDiagnostics -d ${1}.root --expectSignal 0 --forceRecreateNLL --cminDefaultMinimizerType Minuit --rMin 0 --rMax 2 --ignoreCovWarning --saveWithUncertainties -t -1
    
    if [ -f "fitDiagnosticsTest.root" ]; then
        ls -l fitDiagnosticsTest.root
        echo "The output will be copied to ${_CONDOR_SCRATCH_DIR}"
        cp fitDiagnosticsTest.root ${_CONDOR_SCRATCH_DIR}/fitDiagnostics.result2.${1}.root
    fi
'''
