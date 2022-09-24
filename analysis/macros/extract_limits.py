#!/usr/bin/env python
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-l', '--logs', help='logs', dest='logs', default='logs')
(options, args) = parser.parse_args()

command = 'tail '+ options.logs
results=os.popen(command).read()
for line in results.split():
  if 'Observed Limit:' not in line: continue
  print(line.strip().split('<')[0].split(:))[1], line.strip().split('<')[1])
