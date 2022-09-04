#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d', '--datacard', help='datacard', dest='datacard')
    (options, args) = parser.parse_args()
    
    datacard=open(options.datacard,'r')
    processes = []
    for line in datacard:
      if 'process' not in line: continue
      processes.append(line)
      
    print(len(processes))
    print(processes)
    
