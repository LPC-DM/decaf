#!/usr/bin/env python
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='darkhiggs')
(options, args) = parser.parse_args()

for filename in os.listdir('results'):
    if '.tgz' not in filename: continue
    if options.analysis not in filename: continue
    print('Removing')
    os.system('ls results/'+filename.split('.')[0]+'/*')
    os.system('rm results/'+filename.split('.')[0]+'/*')
    print('Untarring',filename)
    os.system('tar -zxvf results/'+filename)
