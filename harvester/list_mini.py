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
campaigns ={}
campaigns['2016'] = 'RunIISummer16MiniAODv3'
campaigns['2017'] = 'RunIIFall17MiniAODv2'
campaigns['2018'] = 'RunIIAutumn18MiniAOD'

slist = open("samples.txt")
for sample in slist:
    s = sample.strip()
    os.system("dasgoclient --query=\"dataset dataset=/"+s+"*/"+campaigns[options.year]+"*/MINIAODSIM\" >"+s+".txt")
    filenames.append(s+".txt")

with open('miniaod'+options.year+'.txt', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
        os.system('rm '+fname)
