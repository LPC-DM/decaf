#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import numexpr
import os
from optparse import OptionParser
import uproot, uproot_methods
import numpy as np
from fnal_column_analysis_tools import hist
#from saiyan import Builder
from analysis.darkhiggs import analysis,samples

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-s', '--selection', help='selection', dest='selection')
parser.add_option('-l', '--lumi', help='lumi', dest='lumi')
(options, args) = parser.parse_args()

with open("../beans/"+options.year+".json") as fin:
    datadef = json.load(fin)

lumi = 1000.
if options.lumi: lumi=lumi*float(options.lumi)
os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz ../../decaf')
os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')
os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
for selection,v in samples.items():
    if options.selection and options.selection not in selection: continue    
    os.system("mkdir -p pods/"+options.year+"/"+selection+"/condor/out pods/"+options.year+"/"+selection+"/condor/err pods/"+options.year+"/"+selection+"/condor/log")
    os.system("cp run run.sh pods/"+options.year+"/"+selection+"/")
    print(os.environ['PWD'])
    for i in range (0,len(v)):
        for dataset, info in datadef.items():
            if options.dataset and options.dataset not in dataset: continue
            if v[i] not in dataset: continue
            os.system("./submit.sh "+options.year+" "+options.lumi+" "+selection+" "+dataset)
