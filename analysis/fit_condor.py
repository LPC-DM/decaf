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
parser.add_option('-w', '--workspace', help='workspace', dest='workspace')
parser.add_option('-M', '--method', help='method', dest='method')
parser.add_option('-n', '--name', help='name', dest='name', default='')
parser.add_option('-a', '--arguments', help='arguments', dest='arguments')
parser.add_option('-l', '--list', help='list', dest='list')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw.tgz '
              '--exclude=\'src/decaf/analysis/data\' '
              '--exclude=\'src/decaf/analysis/metadata\' '
              '--exclude=\'src/decaf/analysis/logs\' '
              '--exclude=\'src/decaf/analysis/plots\' '
              '--exclude=\'src/decaf/analysis/hists\' '
              '--exclude=\'src/decaf/analysis/notebooks\' '
              '--exclude=\'src/decaf/analysis/results/*\' '
              '--exclude=\'src/decaf/analysis/datacards/*-*\' '
              '--exclude=\'src/decaf/analysis/datacards/*.tgz\' '
              '--exclude=\'src/decaf/analysis/datacards/*/*.txt\' '
              '--exclude=\'src/decaf.tgz\' '
              '--exclude=\'src/pylocal.tgz\' '
              '../../../../CMSSW_10_2_13')
   
if options.cluster == 'kisti':
    if options.copy:
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
        print('cmssw removed')
        os.system('xrdcp -f ../../../../cmssw.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = fit.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = fit.sh, /tmp/x509up_u556950957
Output = logs/condor/fit/out/$ENV(OUTFOLDER)_$(Cluster)_$(Process).stdout
Error = logs/condor/fit/err/$ENV(OUTFOLDER)_$(Cluster)_$(Process).stderr
Log = logs/condor/fit/log/$ENV(OUTFOLDER)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(OUTFOLDER).tgz=$ENV(PWD)/results/$ENV(OUTFOLDER).tgz"
Arguments = $ENV(WORKSPACE) $ENV(METHOD) $ENV(ARGUMENTS) $ENV(OUTFOLDER) $ENV(NAME) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(OUTFOLDER)
request_memory = 8000
Queue 1"""

if options.cluster == 'lpc':
    os.system('xrdcp -f ../../../../cmssw.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = fit.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = fit.sh, /tmp/x509up_u556950957
Output = logs/condor/fit/out/$ENV(OUTFOLDER)_$(Cluster)_$(Process).stdout
Error = logs/condor/fit/err/$ENV(OUTFOLDER)_$(Cluster)_$(Process).stderr
Log = logs/condor/fit/log/$ENV(OUTFOLDER)_$(Cluster)_$(Process).log
TransferOutputRemaps = "$ENV(OUTFOLDER).tgz=$ENV(PWD)/results/$ENV(OUTFOLDER).tgz"
Arguments = $ENV(WORKSPACE) $ENV(METHOD) $ENV(ARGUMENTS) $ENV(OUTFOLDER) $ENV(NAME) $ENV(CLUSTER) $ENV(USER)
request_memory = 8000
Queue 1"""

jdl_file = open("fit.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

rootfile = options.workspace.split('/')[-1]
folder = options.workspace.replace(rootfile, '')

datacard=''
for filename in os.listdir(folder):
    if '.txt' in filename: datacard=folder+filename

process_lines=[]
for line in open(datacard,'r').readlines():
    if not line.startswith('process'): continue
    process_lines.append(line.split())

signal_indices = [i for i in range(1, len(process_lines[1])) if int(process_lines[1][i]) <= 0]      
signals = set([process_lines[0][i] for i in signal_indices if process_lines[0][i]])

workspaces=[]
for workspace in os.listdir(folder):
    if '.root' not in workspace: continue
    if not all(piece in workspace for piece in rootfile.split('*')): continue
    workspaces.append(workspace)

tag=options.name+options.method
for workspace in workspaces:
    if options.arguments:
        if 'SIGNAL' in options.arguments:
            for signal in signals:
                if options.list:
                    if not any(_item in signal for _item in options.list.split(',')): continue
                if len(workspaces)>1 and signal not in workspace: continue
                outfolder = workspace.replace('.root','')+'_'+tag
                os.system('mkdir -p logs/condor/fit/err/')
                os.system('rm -rf logs/condor/fit/err/*'+outfolder+'*')
                os.system('mkdir -p logs/condor/fit/log/')
                os.system('rm -rf logs/condor/fit/log/*'+outfolder+'*')
                os.system('mkdir -p logs/condor/fit/out/')
                os.system('rm -rf logs/condor/fit/out/*'+outfolder+'*')
                os.environ['CLUSTER'] = options.cluster
                os.environ['WORKSPACE'] = folder+workspace
                os.environ['METHOD'] = options.method
                os.environ['NAME'] = options.name
                os.environ['OUTFOLDER']  = outfolder
                os.environ['ARGUMENTS']     = options.arguments.replace('SIGNAL',signal).replace(' ','+').replace('"','X')
                os.system('condor_submit fit.submit')
        else:
            #outfolder = options.workspace.split('/')[-2]+'_'+tag
            outfolder = workspace.replace('.root','')+'_'+tag
            os.system('mkdir -p logs/condor/fit/err/')
            os.system('rm -rf logs/condor/fit/err/*'+outfolder+'*')
            os.system('mkdir -p logs/condor/fit/log/')
            os.system('rm -rf logs/condor/fit/log/*'+outfolder+'*')
            os.system('mkdir -p logs/condor/fit/out/')
            os.system('rm -rf logs/condor/fit/out/*'+outfolder+'*')
            os.environ['CLUSTER'] = options.cluster
            os.environ['WORKSPACE'] = folder+workspace
            os.environ['METHOD'] = options.method
            os.environ['NAME'] = options.name
            os.environ['OUTFOLDER']  = outfolder
            os.environ['ARGUMENTS']     = options.arguments.replace(' ','+').replace('"','X')
            os.system('condor_submit fit.submit')
    else:
        #outfolder = options.workspace.split('/')[-2]+'_'+tag
        outfolder = workspace.replace('.root','')+'_'+tag
        os.system('mkdir -p logs/condor/fit/err/')
        os.system('rm -rf logs/condor/fit/err/*'+outfolder+'*')
        os.system('mkdir -p logs/condor/fit/log/')
        os.system('rm -rf logs/condor/fit/log/*'+outfolder+'*')
        os.system('mkdir -p logs/condor/fit/out/')
        os.system('rm -rf logs/condor/fit/out/*'+outfolder+'*')
        os.environ['CLUSTER'] = options.cluster
        os.environ['WORKSPACE'] = folder+workspace
        os.environ['METHOD'] = options.method
        os.environ['NAME'] = options.name
        os.environ['OUTFOLDER']  = outfolder
        os.environ['ARGUMENTS']  = 'None'
        os.system('condor_submit fit.submit')
os.system('rm fit.submit')
