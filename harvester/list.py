#!/usr/bin/env python
import uproot
import json
import pyxrootd.client
import fnmatch
import numpy as np
import numexpr
import subprocess
import concurrent.futures
import warnings
import os
import difflib
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-y', '--year', help='year', dest='year')
(options, args) = parser.parse_args()
filenames = []

slist = open("samples.txt")
for sample in slist:
    s = sample.strip()
    os.system("dasgoclient --query=\"dataset dataset=/"+s+"*/RunIIFall17MiniAODv2*/MINIAODSIM\" >"+s+".txt")
    filenames.append(s+".txt")

with open('miniaod'+options.year+'.txt', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
        os.system('rm '+fname)
