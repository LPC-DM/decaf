#!/usr/bin/env python
import os
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-w', '--workspace', help='workspace', dest='workspace')
    parser.add_option('-M', '--method', help='method', dest='method')
    parser.add_option('-a', '--arguments', help='arguments', dest='arguments')
    (options, args) = parser.parse_args()
    
    command='combine -M '+options.method

    rootfile = options.workspace.split('/')[-1]
    folder = options.workspace.replace(rootfile, '')
    
    datacard=''
    for filename in os.listdir(folder):
        if '.txt' in filename: datacard=folder+'/'+filename
          
    process_lines=[]
    for line in open(datacard,'r').readlines():
        if not line.startswith('process'): continue
        process_lines.append(line.split())

    signal_indices = [i for i in range(1, len(process_lines[1])) if int(process_lines[1][i]) <= 0]      
    signals = set([process_lines[0][i] for i in signal_indices if process_lines[0][i]])
    
    workspaces=[]
    for workspace in os.listdir(folder):
        if '.root' not in workspace: continue
        if not all(piece in workspace for piece in rootfile.split('*')): continue
        workspaces.append(workspace)

    commands=[]
    for workspace in workspaces:
        if options.arguments:
            if 'SIGNAL' in options.arguments:
                for signal in signals:
                    if signal not in workspace: continue
                    commands.append(command+' -d '+folder+'/'+workspace+' '+
                                    '-n .'+options.method+'Results ' +
                                    options.arguments.replace('SIGNAL',signal).replace('\\"','\''))
            else:
                commands.append(command+' -d '+folder+'/'+workspace+' '+
                                '-n .'+options.method+'Results ' +
                                options.arguments.replace('\\"','\''))
        else:
            commands.append(command+' -d '+folder+'/'+workspace+' '+
                            '-n .'+options.method+'Results ')
                
    for command in commands:
        os.system(command)
        os.system('mkdir -p results/'+options.method+'Results_'+command.split('-d ')[1].split('.root')[0].split('/')[-1])
        os.system('mv *.'+options.method+'Results.* results/'+options.method+'Results_'+command.split('-d ')[1].split('.root')[0].split('/')[-1])
    
