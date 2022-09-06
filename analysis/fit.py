#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-w', '--workspace', help='workspace', dest='workspace')
    parser.add_option('-M', '--method', help='method', dest='method')
    parser.add_option('-o', '--options', help='options', dest='options')
    
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
    for workspace in in os.listdir(folder)
        if '.root' not in workspace: continue
        if rootfile.replace('.root','').replace('*','') not in workspace: continue
        workspaces.append(workspace)
    
    commands=[]
    for workspace in workspaces:
        if options.options:
            if 'SIGNAL' in options.options
                for signal in signals:
                    if signal not in option: continue
                    commands.append(command+' -d '+folder+'/'+workspace+' '+options.options.replace('SIGNAL',signal))
            else:
                commands.append(command+' -d '+folder+'/'+workspace+' '+options.options)
                
    for command in commands:
        print(command)
        #os.system(command)
    
