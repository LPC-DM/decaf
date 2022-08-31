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
parser.add_option('-p', '--processor', help='processor', dest='processor', default='')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

os.system('mkdir -p hists/'+options.processor+'/run_condor/out hists/'+options.processor+'/run_condor/err hists/'+options.processor+'/run_condor/log')

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*condor/*/*\' --exclude=\'analysis/hists/*/*.futures\' --exclude=\'analysis/hists/*/*.merged\' --exclude=\'analysis/hists/*/*.reduced\' --exclude=\'analysis/plots\' ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')

if options.cluster == 'kisti':
    if options.copy:
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
        print('decaf removed')
        os.system('xrdcp -f ../../decaf.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/pylocal.tgz') 
        print('pylocal removed')
        os.system('xrdcp -f ../../pylocal.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh, /tmp/x509up_u556950957
Output = hists/$ENV(PROCESSOR)/run_condor/out/$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = hists/$ENV(PROCESSOR)/run_condor/err/$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = hists/$ENV(PROCESSOR)/run_condor/log/$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(PROCESSOR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(PROCESSOR)/$ENV(SAMPLE).futures"
Arguments = $ENV(METADATA) $ENV(SAMPLE) $ENV(PROCESSOR) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(BTCN)
request_cpus = 8
request_memory = 7000
Queue 1"""

if options.cluster == 'lpc':
    if options.copy:
        os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh
Output = hists/$ENV(PROCESSOR)/run_condor/out/$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = hists/$ENV(PROCESSOR)/run_condor/err/$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = hists/$ENV(PROCESSOR)/run_condor/log/$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(PROCESSOR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(PROCESSOR)/$ENV(SAMPLE).futures"
Arguments = $ENV(METADATA) $ENV(SAMPLE) $ENV(PROCESSOR) $ENV(CLUSTER) $ENV(USER) 
request_cpus = 8
request_memory = 5700
Queue 1"""

jdl_file = open("run.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

with open('metadata/'+options.metadata+'.json') as fin:
    datadef = json.load(fin)

for dataset, info in datadef.items():
    if options.dataset:
        if not any(_dataset in dataset for _dataset in options.dataset.split(',')): continue
    if options.exclude:
        if any(_dataset in dataset for _dataset in options.exclude.split(',')): continue
    os.system('rm -rf hists/'+options.processor+'/run_condor/err/'+dataset+'*')
    os.system('rm -rf hists/'+options.processor+'/run_condor/log/'+dataset+'*')
    os.system('rm -rf hists/'+options.processor+'/run_condor/out/'+dataset+'*')
    os.environ['SAMPLE'] = dataset
    os.environ['BTCN'] = dataset.split('____')[0]
    os.environ['PROCESSOR']   = options.processor
    os.environ['METADATA']   = options.metadata
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit run.submit')
os.system('rm run.submit')
