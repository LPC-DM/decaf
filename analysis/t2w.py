#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', '--signal', help='signal', dest='signal')
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('--multiSignal', action='store_true', dest='multiSignal')
    (options, args) = parser.parse_args()

    def multiSignalModel():
        
        signals=[]
        for k,v in processes.items():
            process = k
            if not isinstance(k, str):
                process = k[0]
            if options.signal.split(':')[0] not in process: continue
            if process not in signals:
                signals.append(process)
    
        tmp = []
        name = ''
        for signal in signals:
            try:
                if options.signal.split(':')[1] in signal:
                    name='_'+signal
                    label = ' --PO map=.*/'+signal+':'+signal+'_r[1,0,10]'
                else:
                    label = ' --PO map=.*/'+signal+':0'
            except:
                label = ' --PO map=.*/'+signal+':'+signal+'_r[1,0,10]'
            tmp.append(label)
        option = ''.join(tmp)
    
        command = 'text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose'+option+' datacards/'+options.folder+'/'+options.folder+'.txt -o datacards/'+options.folder+'/'+options.folder+name+'.root'
        return command
    
    
    command = 'text2workspace.py datacards/'+options.folder+'/'+options.folder+'.txt -o datacards/'+options.folder+'/'+options.folder+'.root'
    if options.multiSignal: 
        command = multiSignalModel()
    print(command)
    #os.system(command)
