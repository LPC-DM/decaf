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
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
(options, args) = parser.parse_args()

os.system('mkdir -p datacards/condor/out datacards/condor/err datacards/condor/log')

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw.tgz --exclude=\'src/decaf/analysis/hists/*\' ../../../../CMSSW_10_2_13')

if options.cluster == 'kisti':
    os.system('xrdcp -f ../../../../cmssw.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = render.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = render.sh, /tmp/x509up_u556950957
Output = datacards/condor/out/$ENV(MODEL)_$(Cluster)_$(Process).stdout
Error = datacards/condor/err/$ENV(MODEL)_$(Cluster)_$(Process).stderr
Log = datacards/condor/log/$ENV(MODEL)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(MODEL).tgz=$ENV(PWD)/datacards/$ENV(MODEL).tgz"
Arguments = $ENV(MODEL) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
request_memory = 8000
request_cpus = 16
Queue 1"""

if options.cluster == 'lpc':
    os.system('xrdcp -f ../../../../cmssw.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = render.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = render.sh
Output = datacards/condor/out/$ENV(MODEL)_$(Cluster)_$(Process).stdout
Error = datacards/condor/err/$ENV(MODEL)_$(Cluster)_$(Process).stderr
Log = datacards/condor/log/$ENV(MODEL)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(MODEL).tgz=$ENV(PWD)/datacards/$ENV(MODEL).tgz"
Arguments = $ENV(MODEL) $ENV(CLUSTER) $ENV(USER)
request_memory = 8000
request_cpus = 16
Queue 1"""

jdl_file = open("render.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

if options.analysis:
    for filename in os.listdir('data'):
        if '.model' not in filename: continue
        if options.analysis not in filename: continue
        print('Preparing job for model', filename.split('.')[0])
        os.system('mkdir -p datacards/'+filename.split('.')[0])
        os.system('rm -rf datacards/condor/err/'+filename.split('.')[0]+'*')
        os.system('rm -rf datacards/condor/log/'+filename.split('.')[0]+'*')
        os.system('rm -rf datacards/condor/out/'+filename.split('.')[0]+'*')
        os.environ['MODEL']   = filename.split('.')[0]
        os.environ['CLUSTER'] = options.cluster
        os.system('condor_submit render.submit')
elif options.model:
    print('Preparing job for model', options.model)
    os.system('mkdir -p datacards/'+options.model)
    os.system('rm -rf datacards/condor/err/'+options.model+'*')
    os.system('rm -rf datacards/condor/log/'+options.model+'*')
    os.system('rm -rf datacards/condor/out/'+options.model+'*')
    os.environ['MODEL']   = options.model
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit render.submit')
os.system('rm render.submit')
