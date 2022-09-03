#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', '--signal', help='signal', dest='signal')
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-o', '--output', help='output', dest='output')
    (options, args) = parser.parse_args()

    process_list=[]
    for k,v in processes.items():
        process = k
        if not isinstance(k, str):
            process = k[0]
        if process not in process_list:
            process_list.append(process)

    option = ''
    for signal in options.signal.split(';'):
        if 'SIGNAL' in signal.split(':')[1]:
            for process in process_list:
                if signal.split(':')[0].replace('*','') not in process: continue
                option += ' --PO map=.*/'+process+':'+signal.split(':')[1].replace('SIGNAL',process)
        else:
            option += ' --PO map=.*/'+signal.split(':')[0]+':'+signal.split(':')[1]
    
    command = 'text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose'+option+' datacards/'+options.folder+'/'+options.folder+'.txt -o datacards/'+options.folder+'/'+options.output+'.root'
    print(command)
    #os.system(command)
