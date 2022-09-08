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
from coffea.util import load

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset', default='')
parser.add_option('-e', '--exclude', help='exclude', dest='exclude', default='')
parser.add_option('-f', '--folder', help='folder', dest='folder')
parser.add_option('-v', '--variable', help='variable', dest='variable')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*condor/*/*\' --exclude=\'analysis/hists/*/*.reduced\' --exclude=\'analysis/hists/*/*.merged\' --exclude=\'analysis/plots\' ../../decaf')
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
Executable = reduce.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = reduce.sh, /tmp/x509up_u556950957
Output = logs/condor/reduce/out/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/reduce/err/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/reduce/log/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE)_$ENV(SAMPLE).reduced=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE)--$ENV(SAMPLE).reduced"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(SAMPLE) $ENV(CLUSTER) $ENV(USER)
JobBatchName = $ENV(VARIABLE)
accounting_group=group_cms
request_cpus = 16
request_disk = 10G
Queue 1"""

if options.cluster == 'lpc':
    if options.copy:
        os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = reduce.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = reduce.sh
Output = logs/condor/reduce/out/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = logs/condor/reduce/err/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = logs/condor/reduce/log/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE)_$ENV(SAMPLE).reduced=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE)--$ENV(SAMPLE).reduced"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(SAMPLE) $ENV(CLUSTER) $ENV(USER)
request_cpus = 16
Queue 1"""

jdl_file = open("reduce.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

pd = []
futurefile=''
for filename in os.listdir(options.folder):
    if '.futures' not in filename: continue
    futurefile=filename
    if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])

variables=load(options.folder+'/'+futurefile).keys()
for pdi in pd:
    if options.dataset:
        if not any(_dataset in pdi for _dataset in options.dataset.split(',')): continue
    if options.exclude:
        if any(_dataset in pdi for _dataset in options.exclude.split(',')): continue
    for variable in variables:
        if options.variable and options.variable not in variable: continue
        os.environ['FOLDER'] = options.folder
        os.environ['SAMPLE'] = pdi
        os.environ['VARIABLE'] = variable
        os.environ['CLUSTER'] = options.cluster
        os.system('condor_submit reduce.submit')
os.system('rm reduce.submit')
