#!/usr/bin/env python
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='darkhiggs')
(options, args) = parser.parse_args()
os.system('rm datacards/'+options.analysis+'.txt')
command='combineCards.py '
filenames = 'datacards/*'+options.analysis+'* -name \'*.txt\''
os.system('find '+filenames+' > cards.txt')
for card in open('cards.txt'):
    #if 'sr' not in card: continue
    #if 'monojet' not in card: continue
    filename=card.strip()
    binname=filename.split(".")[0].split('/')[len(filename.split(".")[0].split('/'))-1].replace('-','')
    command=command+binname+'='+filename+' '
command=command+' > datacards/'+options.analysis+'.txt'  
os.system(command)
os.system('rm cards.txt')
#print(command)
