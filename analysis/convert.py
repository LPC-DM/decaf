#!/usr/bin/env python
import os
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d', '--datacard', help='datacard', dest='datacard')
    parser.add_option('-o', '--outfile', help='outfile', dest='outfile')
    parser.add_option('-m', '--maps', help='maps', dest='maps')
    parser.add_option('-a', '--arguments', help='arguments', dest='arguments', default='')
    (options, args) = parser.parse_args()
    
    command = 'text2workspace.py '+options.datacard
    
    datacard=open(options.datacard,'r')
    process_lines=[]
    for line in datacard.readlines():
        if not line.startswith('process'): continue
        process_lines.append(line.split())

    signal_indices = [i for i in range(1, len(process_lines[1])) if int(process_lines[1][i]) <= 0]      
    signals = set([process_lines[0][i] for i in signal_indices if process_lines[0][i]])

    def add_maps(command, options):
        #command += ' -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose'
        command += ' -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel'
        maps= ''
        for option in options.split('--PO '):
            if not option: continue
            if 'SIGNAL' in option:
                for signal in signals:
                    if not all(piece in signal for piece in option.split(':')[0].split('/')[1].split('*')): continue
                    #if option.split(':')[0].split('/')[1].replace('*','') not in signal: continue
                    maps += ' --PO '+option.replace(option.split(':')[0].split('/')[1],signal).replace('SIGNAL',signal)
            else:
                maps += ' --PO '+option
        command += maps
        return command
    
    commands=[]
    if options.maps: 
        if 'SIGNAL:' in options.maps:
            for signal in signals:
                commands.append(add_maps(command, options.maps.replace('SIGNAL',signal))+' -o '+options.outfile.replace('SIGNAL',signal)+' '+options.arguments)
        else:
            commands.append(add_maps(command, options.maps)+' -o '+options.outfile+' '+options.arguments)
    else:
        commands.append(command+' -o '+options.outfile+' '+options.arguments)

    for command in commands:
        #print(command)
        os.system(command)
