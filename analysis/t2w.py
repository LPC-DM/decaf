#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-m', '--map', help='maps', dest='maps')
    parser.add_option('-d', '--datacard', help='datacard', dest='datacard')
    parser.add_option('-o', '--output', help='output', dest='output')
    parser.add_option('-p', '--process', help='process', dest='process')
    (options, args) = parser.parse_args()

    datacard = ' datacards/'+options.folder+'/'+options.folder+'.txt'
    outfile = ' -o datacards/'+options.folder+'/'+options.output+'.root'

    process_list=[]
    for k,v in processes.items():
        process = k
        if not isinstance(k, str):
            process = k[0]
        if process not in process_list:
            process_list.append(process)

    def write_command():
        command = 'text2workspace.py'+datacard+outfile
        option = ''
        if options.maps:
            command += ' -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose'
            for _map in options.maps.split(';'):
                if _map.split(':')[0].split('/')[1] in _map.split(':')[1]:
                    for process in process_list:
                        if _map.split(':')[0].split('/')[1].replace('*','') not in process: continue
                        option += ' --PO map='+_map.replace(_map.split(':')[0].split('/')[1], process)
                else:
                    option += ' --PO map='+_map
    
        command += option
        return command
    
    commands=[]
    for process in process_list:
        if process not in options.process: continue
        commands.append(write_command().replace('PROCESS',process))
    print(write_command())
    #os.system(command)
