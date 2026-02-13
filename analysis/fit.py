#!/usr/bin/env python
import os
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-w', '--workspace', help='workspace', dest='workspace')
    parser.add_option('-M', '--method', help='method', dest='method')
    parser.add_option('-n', '--name', help='name', dest='name', default='')
    parser.add_option('-a', '--arguments', help='arguments', dest='arguments')
    (options, args) = parser.parse_args()
    
    command='combine -M '+options.method

    rootfile = options.workspace.split('/')[-1]
    folder = options.workspace.replace(rootfile, '')
    
    datacards=[]
    for filename in os.listdir(folder):
        if '_' in filename and '.root' in filename: 
            #print(filename.split('_',1)[1].replace('.root',''))
            datacards.append(filename.split('_',1)[1].replace('.root',''))

    signals = set(datacards)
    
    workspaces=[]
    for workspace in os.listdir(folder):
        if '.root' not in workspace: continue
        if not all(piece in workspace for piece in rootfile.split('*')): continue
        workspaces.append(workspace)

    commands=[]
    tag=options.name+options.method
    for workspace in workspaces:
        if options.arguments:
            if 'SIGNAL' in options.arguments:
                for signal in signals:
                    if len(workspaces)>1 and signal not in workspace: continue
                    commands.append(command+' -d '+folder+workspace+' ' +
                                    '-n _'+workspace.replace('.root','')+'_'+tag+' ' +
                                    options.arguments.replace('SIGNAL',signal).replace('\\"','\''))
            else:
                commands.append(command+' -d '+folder+workspace+' ' +
                                '-n _'+workspace.replace('.root','')+'_'+tag+' ' +
                                options.arguments.replace('\\"','\''))
        else:
            commands.append(command+' -d '+folder+workspace+' ' +
                            '-n _'+workspace.replace('.root','')+'_'+tag)
                
    for command in commands:
        os.system(command)
        folder='results/'+command.split('-n ')[1].split(' ')[0].split('_',1)[1]
        os.system('mkdir -p '+folder)
        os.system('mv *'+tag+'* '+folder)
