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
parser.add_option('-m', '--mass', help='mass', dest='mass', default='')
parser.add_option('-g', '--category', help='category', dest='category', default='')
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='')
parser.add_option('-y', '--year', help='year', dest='year', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*____*\' --exclude=\'analysis/hists/*/*condor/*/*\' ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')

if options.cluster == 'kisti':
    os.system('xrdcp -f ../../decaf.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
    os.system('xrdcp -f ../../pylocal.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = limit.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = limit.sh, /tmp/x509up_u556950957
Output = datacards/$ENV(ANALYSIS)$ENV(YEAR)/condor/out/$ENV(MASS)$ENV(CATEGORY)_$(Cluster)_$(Process).stdout
Error = datacards/$ENV(ANALYSIS)$ENV(YEAR)/condor/err/$ENV(MASS)$ENV(CATEGORY)_$(Cluster)_$(Process).stderr
Log = datacards/$ENV(ANALYSIS)$ENV(YEAR)/condor/log/$ENV(MASS)$ENV(CATEGORY)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(ANALYSIS)$ENV(YEAR)_$ENV(MASS).tgz=$ENV(PWD)/datacards/$ENV(ANALYSIS)$ENV(YEAR)/$ENV(MASS).tgz"
Arguments = $ENV(YEAR) $ENV(SAMPLE) $ENV(ANALYSIS) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
request_cpus = 8
Queue 1"""

if options.cluster == 'lpc':
    os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
    os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = limit.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = limit.sh
Output = datacards/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/out/$ENV(SAMPLE)_$(Cluster)_$(Process).stdout
Error = datacards/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/err/$ENV(SAMPLE)_$(Cluster)_$(Process).stderr
Log = datacards/$ENV(ANALYSIS)$ENV(YEAR)/run_condor/log/$ENV(SAMPLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(ANALYSIS)$ENV(YEAR)_$ENV(MASS).tgz=$ENV(PWD)/datacards/$ENV(ANALYSIS)$ENV(YEAR)/$ENV(MASS).tgz"
Arguments = $ENV(YEAR) $ENV(SAMPLE) $ENV(ANALYSIS) $ENV(CLUSTER) $ENV(USER) 
request_cpus = 8
Queue 1"""

jdl_file = open("limit.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

with open('metadata/'+options.year+'.json') as fin:
    datadef = json.load(fin)

for mass in ['mass0','mass1','mass2','mass3','mass4']:
    for category in ['monojet','monohs']:
        if options.mass and options.mass not in mass: continue
        if options.category and options.category in category: continue
        if options.year and options.year not in year: continue
        os.system('mkdir -p datacards/'+options.analysis+options.year+'/condor/out')
        os.system('mkdir -p datacards/'+options.analysis+options.year+'/condor/err')
        os.system('mkdir -p datacards/'+options.analysis+options.year+'/condor/log')
        os.system('rm -rf datacards/'+options.analysis+options.year+'/condor/err/'options.mass+options.category+'*')
        os.system('rm -rf datacards/'+options.analysis+options.year+'/condor/log/'options.mass+options.category+'*')
        os.system('rm -rf datacards/'+options.analysis+options.year+'/condor/out/'options.mass+options.category+'*')
        os.environ['ANALYSIS']   = options.analysis
        os.environ['YEAR']   = options.year
        os.environ['MASS']   = options.mass
        os.environ['CATEGORY']   = options.category
        os.environ['CLUSTER'] = options.cluster
        os.system('condor_submit limit.submit')
os.system('rm limit.submit')
