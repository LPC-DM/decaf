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
parser.add_option('-f', '--folder', help='folder', dest='folder')
parser.add_option('-v', '--variable', help='variable', dest='variable')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
(options, args) = parser.parse_args()

os.system('mkdir -p '+options.folder+'/reduce_condor/out '+options.folder+'/reduce_condor/err '+options.folder+'/reduce_condor/log')
os.system('rm -rf '+options.folder+'/reduce_condor/err/'+options.dataset+'*')
os.system('rm -rf '+options.folder+'/reduce_condor/log/'+options.dataset+'*')
os.system('rm -rf '+options.folder+'/reduce_condor/out/'+options.dataset+'*')

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*condor/*/*\' ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')

if options.cluster == 'kisti':
    if options.tar:
        os.system('xrdcp -f ../../decaf.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = reduce.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = reduce.sh, /tmp/x509up_u556950957
Output = $ENV(FOLDER)/reduce_condor/out/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = $ENV(FOLDER)/reduce_condor/err/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = $ENV(FOLDER)/reduce_condor/log/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE)_$ENV(SAMPLE).reduced=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE)--$ENV(SAMPLE).reduced"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(SAMPLE) $ENV(CLUSTER) $ENV(USER)
JobBatchName = $ENV(VARIABLE)_$ENV(SAMPLE)
accounting_group=group_cms
request_cpus = 16
Queue 1"""

if options.cluster == 'lpc':
    if options.tar:
        os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = reduce.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = reduce.sh
Output = $ENV(FOLDER)/reduce_condor/out/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = $ENV(FOLDER)/reduce_condor/err/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = $ENV(FOLDER)/reduce_condor/log/$ENV(SAMPLE)_$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE)_$ENV(SAMPLE).reduced=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE)--$ENV(SAMPLE).reduced"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(SAMPLE) $ENV(CLUSTER) $ENV(USER)
request_cpus = 16
Queue 1"""

jdl_file = open("reduce.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

pd = []
for filename in os.listdir(options.folder):
    if '.futures' not in filename: continue
    if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])

variables = [
    'sumw',
    'CaloMinusPfOverRecoil',
    'template',
    'recoil',
    'met',
    'metphi',
    'mindphi',
    'j1pt',
    'j1eta',
    'j1phi',
    'fj1pt',
    'fj1eta',
    'fj1phi',
    'njets',
    'ndcsvL',
    'ndflvL',
    'nfjclean',
    'fjmass',
    'e1pt',
    'e1eta',
    'e1phi',
    'dielemass',
    'dielept',
    'mu1pt',
    'mu1eta',
    'mu1phi',
    'dimumass',
    'dimupt',
    'ZHbbvsQCD'
]

for pdi in pd:
    if options.dataset and options.dataset not in pdi: continue
    if options.exclude and options.exclude in pdi: continue
    for variable in variables:
        if options.variable and options.variable not in variable: continue
        os.environ['FOLDER'] = options.folder
        os.environ['SAMPLE'] = pdi
        os.environ['VARIABLE'] = variable
        os.environ['CLUSTER'] = options.cluster
        os.system('condor_submit reduce.submit')
os.system('rm reduce.submit')
