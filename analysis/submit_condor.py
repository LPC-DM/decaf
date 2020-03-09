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
from coffea import hist

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-e', '--exclude', help='exclude', dest='exclude')
parser.add_option('-p', '--analysis', help='analysis', dest='analysis')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-t', '--tar', action="store_true", dest="tar")
(options, args) = parser.parse_args()

year=''
if options.year: year=options.year

os.system("mkdir -p hists/"+options.analysis+year+"/condor/out hists/"+options.analysis+year+"/condor/err hists/"+options.analysis+year+"/condor/log")
os.system("rm -rf hists/"+options.analysis+year+"/condor/out/* hists/"+options.analysis+year+"/condor/err/* hists/"+options.analysis+year+"/condor/log/*")

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')
    os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
    os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')

with open("metadata/"+year+".json") as fin:
    datadef = json.load(fin)

for dataset, info in datadef.items():
    if options.dataset and options.dataset not in dataset: continue
    if options.exclude and options.exclude in dataset: continue
    os.environ['SAMPLE'] = dataset
    os.environ['ANALYSIS']   = options.analysis
    os.environ['YEAR']   = options.year
    os.system("condor_submit run.jdl")
