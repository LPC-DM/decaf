#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', '--signal', help='signal', dest='signal')
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    (options, args) = parser.parse_args()

    '''
    signals=[]
    for k,v in processes.items():
        process = k
        if not isinstance(k, str):
            process = k[0]
        if options.signal.split(':')[0] not in process: continue
        if process not in signals:
            signals.append(process)
    '''
    option = ''
    for signal in options.signal.split(','):
        option += ' --PO map=.*/'+signal.split(:)[0]+':'+signal.split(:)[1]
    
    command = 'text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose'+option+' datacards/'+options.folder+'/'+options.folder+'.txt -o datacards/'+options.folder+'/'+options.folder+name+'.root'
    return command
    
    print(command)
    #os.system(command)
