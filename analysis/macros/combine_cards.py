#!/usr/bin/env python
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='darkhiggs')
(options, args) = parser.parse_args()
os.system('mkdir -p datacards/'+options.analysis)
os.system('rm datacards/'+options.analysis+'/*')
command='combineCards.py '
cards = 'datacards/*'+options.analysis+'* -name \'*.txt\''
rootfiles = 'datacards/*'+options.analysis+'* -name \'*.root\''
os.system('find '+cards+' > cards.txt')
for card in open('cards.txt'):
    #if 'sr' not in card: continue
    #if 'monojet' not in card: continue
    filename=card.strip()
    os.system('cp '+filename+' .')
    binname=filename.split(".")[0].split('/')[len(filename.split(".")[0].split('/'))-1].replace('-','')
    command=command+binname+'='+filename.split('/')[len(filename.split('/'))-1]+' '
command=command+' > datacards/'+options.analysis+'/'+options.analysis+'.txt'  
os.system(command)
os.system('rm *.txt')
#os.system('mv *.txt datacards/'+options.analysis)
os.system('find '+rootfiles+' > rootfiles.txt')
for rootfile in open('rootfiles.txt'):
    filename=rootfile.strip()
    os.system('cp '+filename+' datacards/'+options.analysis)
os.system('rm rootfiles.txt')
#print(command)
