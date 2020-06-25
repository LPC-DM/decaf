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
(options, args) = parser.parse_args()

os.system('mkdir -p '+options.folder+'/merge_condor/out '+options.folder+'/merge_condor/err '+options.folder+'/merge_condor/log')
os.system('rm -rf '+options.folder+'/merge_condor/err/'+options.variable+'*')
os.system('rm -rf '+options.folder+'/merge_condor/log/'+options.variable+'*')
os.system('rm -rf '+options.folder+'/merge_condor/out/'+options.variable+'*')

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*____*\' --exclude=\'analysis/hists/*/*condor/*/*\' ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')

if options.cluster == 'kisti':
    if options.tar:
        os.system('xrdcp -f ../../decaf.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = merge.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = merge.sh, /tmp/x509up_u556950957
Output = $ENV(FOLDER)/merge_condor/out/$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = $ENV(FOLDER)/merge_condor/err/$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = $ENV(FOLDER)/merge_condor/log/$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE).merged=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE).merged"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(CLUSTER) $ENV(USER)
JobBatchName = $ENV(VARIABLE)
accounting_group=group_cms
request_cpus = 16
Queue 1"""

if options.cluster == 'lpc':
    if options.tar:
        os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')
    jdl = """universe = vanilla
Executable = merge.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = merge.sh
Output = $ENV(FOLDER)/merge_condor/out/$ENV(VARIABLE)_$(Cluster)_$(Process).stdout
Error = $ENV(FOLDER)/merge_condor/err/$ENV(VARIABLE)_$(Cluster)_$(Process).stderr
Log = $ENV(FOLDER)/merge_condor/log/$ENV(VARIABLE)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(VARIABLE).merged=$ENV(PWD)/$ENV(FOLDER)/$ENV(VARIABLE).merged"
Arguments = $ENV(FOLDER) $ENV(VARIABLE) $ENV(CLUSTER) $ENV(USER)
request_cpus = 16
Queue 1"""

jdl_file = open("merge.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

variables = [
    'sumw',
    'CaloMinusPfOverRecoil',
    'recoil',
    'template',
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

for variable in variables:
    if options.variable and options.variable not in variable: continue
    os.environ['FOLDER'] = options.folder
    os.environ['VARIABLE'] = variable
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit merge.submit')
os.system('rm merge.submit')
