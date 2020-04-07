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
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
(options, args) = parser.parse_args()

os.system('mkdir -p hists/'+options.analysis+options.year+'/run_condor/out hists/'+options.analysis+options.year+'/run_condor/err hists/'+options.analysis+options.year+'/run_condor/log')
os.system('rm -rf hists/'+options.analysis+options.year+'/run_condor/err/'+options.dataset+'*')
os.system('rm -rf hists/'+options.analysis+options.year+'/run_condor/log/'+options.dataset+'*')
os.system('rm -rf hists/'+options.analysis+options.year+'/run_condor/out/'+options.dataset+'*')

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*____*\' --exclude=\'analysis/hists/*/*condor/*/*\' ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')

if options.cluster == 'kisti':
    os.system('xrdcp -f ../../decaf.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
    os.system('xrdcp -f ../../pylocal.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh, /tmp/x509up_u556950957
Output = hists/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/out/$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = hists/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/err/$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = hists/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/log/$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(ANALYSIS)$ENV(YEAR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(ANALYSIS)$ENV(YEAR)/$ENV(SAMPLE).futures"
Arguments = $ENV(YEAR) $ENV(SAMPLE) $ENV(ANALYSIS) $ENV(CLUSTER) $ENV(USER)
JobBatchName = $ENV(SAMPLE)
accounting_group=group_cms
request_cpus = 8
request_memory = 6000
Queue 1"""

if options.cluster == 'lpc':
    os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
    os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = run.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = run.sh
Output = hists/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/out/$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = hists/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/err/$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = hists/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/log/$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(ANALYSIS)$ENV(YEAR)_$ENV(SAMPLE).futures=$ENV(PWD)/hists/$ENV(ANALYSIS)$ENV(YEAR)/$ENV(SAMPLE).futures"
Arguments = $ENV(YEAR) $ENV(SAMPLE) $ENV(ANALYSIS) $ENV(CLUSTER) $ENV(USER) 
request_cpus = 8
request_memory = 5700
Queue 1"""

jdl_file = open("run.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

with open('metadata/'+options.year+'.json') as fin:
    datadef = json.load(fin)

for dataset, info in datadef.items():
    if options.dataset and options.dataset not in dataset: continue
    if options.exclude and options.exclude in dataset: continue
    os.environ['SAMPLE'] = dataset
    os.environ['ANALYSIS']   = options.analysis
    os.environ['YEAR']   = options.year
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit run.submit')
os.system('rm run.submit')
