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
parser.add_option('-t', '--tar', action="store_true", dest="tar")
(options, args) = parser.parse_args()

with open("../harvester/beans/"+options.year+".json") as fin:
    datadef = json.load(fin)

lumi = 1000.
if options.lumi: lumi=lumi*float(options.lumi)
if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')
    os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
    os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')

for selection,v in samples.items():
    if options.selection and options.selection not in selection: continue    
    os.system("mkdir -p pods/"+options.year+"/"+selection+"/condor/out pods/"+options.year+"/"+selection+"/condor/err pods/"+options.year+"/"+selection+"/condor/log")
    os.system("rm -rf pods/"+options.year+"/"+selection+"/condor/out/* pods/"+options.year+"/"+selection+"/condor/err/* pods/"+options.year+"/"+selection+"/condor/log/*")
    for i in range (0,len(v)):
        for dataset, info in datadef.items():
            if options.dataset and options.dataset not in dataset: continue
            if v[i] not in dataset: continue
            os.environ['SAMPLE'] = dataset
            os.environ['YEAR']   = options.year
            os.environ['LUMI']   = options.lumi
            os.environ['SELECTION'] = selection
            os.system("condor_submit run")
