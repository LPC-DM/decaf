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
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw.tgz --exclude=\'src/decaf/analysis/hists/*\' --exclude=\'src/decaf/analysis/plots/*\' --exclude=\'src/decaf/analysis/datacards/*-*\' ../../../../CMSSW_10_2_13')

if options.cluster == 'kisti':
    if options.copy:
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
        print('cmssw removed')
        os.system('xrdcp -f ../../../../cmssw.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = t2w.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = t2w.sh, /tmp/x509up_u556950957
Output = datacards/$ENV(FOLDER)/condor/t2w/out/$(Cluster)_$(Process).stdout
Error = datacards/$ENV(FOLDER)/condor/t2w/err/$ENV(FOLDER)_$(Cluster)_$(Process).stderr
Log = datacards/$ENV(FOLDER)/condor/t2w/log/$ENV(FOLDER)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(FOLDER).root=$ENV(PWD)/datacards/$ENV(FOLDER)/$ENV(FOLDER).root"
Arguments = $ENV(FOLDER) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(FOLDER)
request_memory = 8000
Queue 1"""

if options.cluster == 'lpc':
    os.system('xrdcp -f ../../../../cmssw.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = t2w.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = t2w.sh, /tmp/x509up_u556950957
Output = datacards/$ENV(FOLDER)/condor/t2w/out/$(Cluster)_$(Process).stdout
Error = datacards/$ENV(FOLDER)/condor/t2w/err/$ENV(FOLDER)_$(Cluster)_$(Process).stderr
Log = datacards/$ENV(FOLDER)/condor/t2w/log/$ENV(FOLDER)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(FOLDER).root=$ENV(PWD)/datacards/$ENV(FOLDER)/$ENV(FOLDER).root"
Arguments = $ENV(FOLDER) $ENV(CLUSTER) $ENV(USER)
request_memory = 8000
Queue 1"""

jdl_file = open("t2w.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

for folder in os.listdir('datacards/'):
    if options.analysis not in folder: continue
    if '-' in folder: continue
    print('Preparing job for folder', folder)
    os.system('mkdir -p datacards/'+folder+'/condor/t2w/err/')
    os.system('rm -rf datacards/'+folder+'/condor/t2w/err/*')
    os.system('mkdir -p datacards/'+folder+'/condor/t2w/log/')
    os.system('rm -rf datacards/'+folder+'/condor/t2w/log/*')
    os.system('mkdir -p datacards/'+folder+'/condor/t2w/out/')
    os.system('rm -rf datacards/'+folder+'/condor/t2w/out/*')
    os.environ['FOLDER']   = folder
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit t2w.submit')
os.system('rm t2w.submit')
