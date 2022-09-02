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
from data.process import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-s', '--signal', help='signal', dest='signal')
parser.add_option('-f', '--folder', help='folder', dest='folder')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw.tgz --exclude=\'src/decaf/analysis/hists/*\' --exclude=\'src/decaf/analysis/plots/*\' --exclude=\'src/decaf/analysis/datacards/*-*\' --exclude=\'src/decaf.tgz\' --exclude=\'src/pylocal.tgz\' ../../../../CMSSW_10_2_13')

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
Output = datacards/$ENV(FOLDER)/condor/t2w/out/$ENV(FOLDER)_$ENV(SIGNAL)_$(Cluster)_$(Process).stdout
Error = datacards/$ENV(FOLDER)/condor/t2w/err/$ENV(FOLDER)_$ENV(SIGNAL)_$(Cluster)_$(Process).stderr
Log = datacards/$ENV(FOLDER)/condor/t2w/log/$ENV(FOLDER)_$ENV(SIGNAL)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(FOLDER)_$ENV(SIGNAL).root=$ENV(PWD)/datacards/$ENV(FOLDER)/$ENV(FOLDER)_$ENV(SIGNAL).root"
Arguments = $ENV(SIGNALS) $ENV(SIGNAL) $ENV(FOLDER) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(FOLDER)_$ENV(SIGNALS)_$ENV(SIGNAL)
request_memory = 8000
Queue 1"""

if options.cluster == 'lpc':
    os.system('xrdcp -f ../../../../cmssw.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = t2w.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = t2w.sh, /tmp/x509up_u556950957
Output = datacards/$ENV(FOLDER)/condor/t2w/out/$ENV(FOLDER)_$ENV(SIGNAL)_$(Cluster)_$(Process).stdout
Error = datacards/$ENV(FOLDER)/condor/t2w/err/$ENV(FOLDER)_$ENV(SIGNAL)_$(Cluster)_$(Process).stderr
Log = datacards/$ENV(FOLDER)/condor/t2w/log/$ENV(FOLDER)_$ENV(SIGNAL)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(FOLDER)_$ENV(SIGNAL).root=$ENV(PWD)/datacards/$ENV(FOLDER)/$ENV(FOLDER)_$ENV(SIGNAL).root"
Arguments = $ENV(SIGNAL) $ENV(FOLDER) $ENV(CLUSTER) $ENV(USER)
request_memory = 8000
Queue 1"""

jdl_file = open("t2w.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

signals=[]
for k,v in processes.items():
    process = k
    if not isinstance(k, str):
        process = k[0]
    if options.signal.split(':')[0] not in process: continue
    if process not in signals: signals.append(process)

os.system('mkdir -p datacards/'+options.folder+'/condor/t2w/err/')
os.system('rm -rf datacards/'+options.folder+'/condor/t2w/err/*')
os.system('mkdir -p datacards/'+options.folder+'/condor/t2w/log/')
os.system('rm -rf datacards/'+options.folder+'/condor/t2w/log/*')
os.system('mkdir -p datacards/'+options.folder+'/condor/t2w/out/')
os.system('rm -rf datacards/'+options.folder+'/condor/t2w/out/*')
    
for signal in signals:
    try:
        if not any(_signal in signal for _signal in options.signal.split(':')[1].split(',')): continue
    except: 
        pass
    os.environ['FOLDER']   = options.folder
    os.environ['SIGNALS']  = options.signal.split(':')[0]
    os.environ['SIGNAL']  = signal
    os.environ['CLUSTER'] = options.cluster
    os.system('condor_submit t2w.submit')
os.system('rm t2w.submit')
