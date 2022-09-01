#!/usr/bin/env python
import os
from optparse import OptionParser
from data.process import *

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-s', '--signal', help='signal', dest='signal')
    parser.add_option('-d', '--datacard', help='datacard', dest='datacard')
    (options, args) = parser.parse_args()
    
    for k,v in processes.items():
      print(k)
