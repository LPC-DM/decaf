#!/usr/bin/env python
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='darkhiggs')
(options, args) = parser.parse_args()

folder_list = []
for tgz in os.listdir('datacards/'):
    if '.tgz' not in tgz: continue
    folder_list.append(tgz.split('.')[0])

txt_list = []
root_list = []
for folder in folder_list:
    for txt in os.listdir('datacards/'+folder+'/'):
        if '.txt' not in txt: continue
        txt_list.append('datacards/'+folder+'/'+txt)
    for root in os.listdir('datacards/'+folder+'/'):
            if '.root' not in root: continue
            root_list.append('datacards/'+folder+'/'+root)
os.system('mkdir -p datacards/'+options.analysis)
os.system('rm datacards/'+options.analysis+'/*')

command='combineCards.py '
for card in txt_list:
    #if 'sr' not in card: continue
    #if 'monojet' not in card: continue
    if options.analysis not in card: continue
    filename=card.strip()
    print(filename)
    os.system('cp '+filename+' .')
    binname=filename.split(".")[0].split('/')[len(filename.split(".")[0].split('/'))-1].replace('-','')
    command=command+binname+'='+filename.split('/')[len(filename.split('/'))-1]+' '
command=command+' > datacards/'+options.analysis+'/'+options.analysis+'.txt'  
os.system(command)
for rootfile in root_list:
    if options.analysis not in rootfile: continue
    filename=rootfile.strip()
    os.system('cp '+filename+' datacards/'+options.analysis)
