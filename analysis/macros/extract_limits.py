#!/usr/bin/env python
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-l', '--logs', help='logs', dest='logs', default='logs')
(options, args) = parser.parse_args()

logs=options.logs.split('/')[-1]
folder=options.logs.replace(logs,'')

for log in os.listdir(folder):
  if logs.replace('*','') not in log: continue
  for line in open(folder+'/'+log,'r').readlines():
    if 'Observed Limit:' not in line: continue
    print(line.strip().split('<')[0].split(':')[1],line.strip().split('<')[1])
