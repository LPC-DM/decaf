#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', '--signal', help='signal', dest='signal')
    parser.add_option('-d', '--datacard', help='datacard', dest='datacard')
    (options, args) = parser.parse_args()
    
    signals=[]
    for k,v in processes.items():
        process = k
        if not isinstance(k, str):
            process = k[0]
        if options.signal.split(':')[0] not in process: continue
        signals.append(process)
    
    tmp = []
    for signal in signals:
        try:
            if options.signal.split(':')[1] in signal:
                label = ' --PO map=.*/'+signal+':'+signal+'_r[1,0,10]'
            else:
                label = ' --PO map=.*/'+signal+':0'
        except:
            label = ' --PO map=.*/'+signal+':'+signal+'_r[1,0,10]'
        tmp.append(label)
    option = ''.join(tmp)
    
    command = "text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose"+option+' '+options.datacard+' -o '+options.datacard.split('.')[0]+'.root'
    #print(command)
    os.system(command)
