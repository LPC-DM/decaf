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
parser.add_option('-f', '--folder', help='folder', dest='folder')
parser.add_option('-v', '--variable', help='variable', dest='variable', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

os.system('mkdir -p '+options.folder+'/merge_condor/out '+options.folder+'/merge_condor/err '+options.folder+'/merge_condor/log')
os.system('rm -rf '+options.folder+'/merge_condor/err/'+options.variable+'*')
os.system('rm -rf '+options.folder+'/merge_condor/log/'+options.variable+'*')
os.system('rm -rf '+options.folder+'/merge_condor/out/'+options.variable+'*')

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*condor/*/*\' --exclude=\'analysis/hists/*/*.futures\' --exclude=\'analysis/hists/*/*.merged\' --exclude=\'analysis/plots\' ../../decaf')
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
Executable = merge.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = merge.sh, /tmp/x509up_u556950957
Output = logs/condor/merge/out/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/merge/err/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/merge/log/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE).merged=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE).merged"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(CLUSTER) $ENV(USER)
JobBatchName = $ENV(VARIABLE)
accounting_group=group_cms
request_cpus = 16
Queue 1"""

if options.cluster == 'lpc':
    if options.copy:
        os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = merge.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = merge.sh
Output = logs/condor/merge/out/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/merge/err/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/merge/log/$ENV(TAG)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE).merged=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE).merged"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(CLUSTER) $ENV(USER)
request_cpus = 16
Queue 1"""

jdl_file = open("merge.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

tag=options.folder.split('/')[-1]

variables = []
for filename in os.listdir(options.folder):
          if '.reduced' not in filename: continue
          if filename.split('--')[0] not in variables: variables.append(filename.split('--')[0])

for variable in variables:
    if options.variable and options.variable not in variable: continue
    if options.variable:
        if not any(_variable==variable for _variable in options.variable.split(',')): continue
    os.system('mkdir -p logs/condor/merge/err/')
    os.system('rm -rf logs/condor/merge/err/*'+tag+'*'+variable+'*')
    os.system('mkdir -p logs/condor/merge/log/')
    os.system('rm -rf logs/condor/merge/run/*'+tag+'*'+variable+'*')
    os.system('mkdir -p logs/condor/merge/out/')
    os.system('rm -rf logs/condor/merge/out/*'+tag+'*'+variable+'*')
    os.environ['TAG'] = tag
    os.environ['FOLDER'] = options.folder
    os.environ['VARIABLE'] = variable
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit merge.submit')
os.system('rm merge.submit')
