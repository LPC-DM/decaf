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

parser = OptionParser()
parser.add_option('-m', '--model', help='model', dest='model', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw.tgz '
              '--exclude=\'src/decaf/analysis/logs\' '
              '--exclude=\'src/decaf/analysis/plots\' '
              '--exclude=\'src/decaf/analysis/hists\' '
              '--exclude=\'src/decaf/analysis/results\' '
              '--exclude=\'src/decaf/analysis/datacards/*\' '
              '--exclude=\'src/decaf.tgz\' '
              '--exclude=\'src/pylocal.tgz\' '
              '../../../../CMSSW_10_2_13')

if options.cluster == 'kisti':
    if options.copy:
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
        print('cmssw removed')
        os.system('xrdcp -f ../../../../cmssw.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = render.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = render.sh, /tmp/x509up_u556950957, /cvmfs/cms.cern.ch/cmsset_default.sh
Output = logs/condor/render/out/$ENV(MODEL)_$(Cluster)_$(Process).stdout
Error = logs/condor/render/err/$ENV(MODEL)_$(Cluster)_$(Process).stderr
Log = logs/condor/render/log/$ENV(MODEL)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(MODEL).tgz=$ENV(PWD)/datacards/$ENV(MODEL).tgz"
Arguments = $ENV(MODEL) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(MODEL)
request_memory = 8000
request_cpus = 16
Queue 1"""

if options.cluster == 'lpc':
    if options.tar:
        os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../py2local.tgz -C ~/.local/lib/python2.7/ site-packages')
    if options.copy:
        os.system('xrdcp -f ../../../../cmssw.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw.tgz')
        os.system('xrdcp -f ../../py2local.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/py2local.tgz')
    jdl = """universe = vanilla
Executable = render.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = render.sh
Output = logs/condor/render/out/$ENV(MODEL)_$(Cluster)_$(Process).stdout
Error = logs/condor/render/err/$ENV(MODEL)_$(Cluster)_$(Process).stderr
Log = logs/condor/render/log/$ENV(MODEL)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(MODEL).tgz=$ENV(PWD)/datacards/$ENV(MODEL).tgz"
Arguments = $ENV(MODEL) $ENV(CLUSTER) $ENV(USER)
request_memory = 8000
request_cpus = 16
Queue 1"""

jdl_file = open("render.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

for filename in os.listdir('data/models'):
    if '.model' not in filename: continue
    if options.model:
        if not any(model in filename for model in options.model.split(',')): continue
    print('Preparing job for model', filename.split('.')[0])
    os.system('mkdir -p datacards/'+filename.split('.')[0])
    os.system('mkdir -p logs/condor/render/err/')
    os.system('rm -rf logs/condor/render/err/*'+filename.split('.')[0]+'*')
    os.system('mkdir -p logs/condor/render/log/')
    os.system('rm -rf logs/condor/render/run/*'+filename.split('.')[0]+'*')
    os.system('mkdir -p logs/condor/render/out/')
    os.system('rm -rf logs/condor/render/out/*'+filename.split('.')[0]+'*')
    os.environ['MODEL']   = filename.split('.')[0]
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit render.submit')
os.system('rm render.submit')
