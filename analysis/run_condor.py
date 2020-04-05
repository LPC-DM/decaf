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
parser.add_option('-d', '--dataset', help='dataset', dest='dataset', default='')
parser.add_option('-e', '--exclude', help='exclude', dest='exclude', default='')
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='')
parser.add_option('-y', '--year', help='year', dest='year', default='')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-k', '--kisti', action='store_true', dest='kisti')
(options, args) = parser.parse_args()

os.system('mkdir -p hists/'+options.analysis+options.year+'/condor/out hists/'+options.analysis+options.year+'/condor/err hists/'+options.analysis+options.year+'/condor/log')
os.system('rm -rf hists/'+options.analysis+options.year+'/condor/err/'+options.dataset+'*')
os.system('rm -rf hists/'+options.analysis+options.year+'/condor/log/'+options.dataset+'*')

jdl = 'run.jdl'
if options.kisti: jdl = 'run_kisti.jdl'
os.system('cp jdls/'+jdl+' .')
print('Using',jdl)

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*____*\' --exclude=\'analysis/hists/*/condor/*/*\' ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')
    if options.kisti: 
        os.system('xrdcp -f ../../decaf.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal.tgz')
    else:
        os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')

with open('metadata/'+options.year+'.json') as fin:
    datadef = json.load(fin)

for dataset, info in datadef.items():
    if options.dataset and options.dataset not in dataset: continue
    if options.exclude and options.exclude in dataset: continue
    os.environ['SAMPLE'] = dataset
    os.environ['ANALYSIS']   = options.analysis
    os.environ['YEAR']   = options.year
    os.system('condor_submit '+jdl)
os.system('rm '+jdl)
