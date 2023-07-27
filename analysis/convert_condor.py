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
parser.add_option('-d', '--datacard', help='datacard', dest='datacard')
parser.add_option('-o', '--outfile', help='outfile', dest='outfile')
parser.add_option('-m', '--maps', help='maps', dest='maps')
parser.add_option('-l', '--list', help='list', dest='list')
parser.add_option('-c', '--cluster', help='cluster', dest='cluster', default='lpc')
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-x', '--copy', action='store_true', dest='copy')
(options, args) = parser.parse_args()

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../../../cmssw.tgz '
              '--exclude=\'src/decaf/analysis/data\' '
              '--exclude=\'src/decaf/analysis/logs\' '
              '--exclude=\'src/decaf/analysis/plots\' '
              '--exclude=\'src/decaf/analysis/hists\' '
              '--exclude=\'src/decaf/analysis/results\' '
              '--exclude=\'src/decaf/analysis/datacards/*-*\' '
              '--exclude=\'src/decaf/analysis/datacards/*.tgz\' '
              '--exclude=\'src/decaf.tgz\' '
              '--exclude=\'src/pylocal.tgz\' '
              '../../../../CMSSW_10_2_13')

if options.cluster == 'kisti':
    if options.copy:
        os.system('xrdfs root://cms-xrdr.private.lo:2094/ rm /xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
        print('cmssw removed')
        os.system('xrdcp -f ../../../../cmssw.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/cmssw.tgz')
    jdl = """universe = vanilla
Executable = convert.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = convert.sh, /tmp/x509up_u556950957
Output = logs/condor/convert/out/$ENV(TAG)_$(Cluster)_$(Process).stdout
Error = logs/condor/convert/err/$ENV(TAG)_$(Cluster)_$(Process).stderr
Log = logs/condor/convert/log/$ENV(TAG)_$(Cluster)_$(Process).log
TransferOutputRemaps = "outfile.root=$ENV(PWD)/$ENV(OUTFILE)"
Arguments = $ENV(DATACARD) $ENV(OUTFILE) $ENV(MAPS) $ENV(CLUSTER) $ENV(USER)
accounting_group=group_cms
JobBatchName = $ENV(TAG)
request_memory = 8000
request_disk = 100G
Queue 1"""

if options.cluster == 'lpc':
    if options.tar:
        os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../py2local.tgz -C ~/.local/lib/python2.7/ site-packages')
    if options.copy:
        os.system('xrdcp -f ../../../../cmssw.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/cmssw.tgz')
        os.system('xrdcp -f ../../py2local.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/py2local.tgz')
    jdl = """universe = vanilla
Executable = convert.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = convert.sh, /tmp/x509up_u556950957
Output = logs/condor/convert/out/$ENV(TAG)_$(Cluster)_$(Process).stdout
Error = logs/condor/convert/err/$ENV(TAG)_$(Cluster)_$(Process).stderr
Log = logs/condor/convert/log/$ENV(TAG)_$(Cluster)_$(Process).log
TransferOutputRemaps = "outfile.root=$ENV(PWD)/$ENV(OUTFILE)"
Arguments = $ENV(DATACARD) $ENV(OUTFILE) $ENV(MAPS) $ENV(CLUSTER) $ENV(USER)
request_memory = 8000
Queue 1"""

jdl_file = open("convert.submit", "w") 
jdl_file.write(jdl) 
jdl_file.close() 

datacard=open(options.datacard,'r')
process_lines=[]
for line in datacard.readlines():
    if not line.startswith('process'): continue
    process_lines.append(line.split())
signal_indices = [i for i in range(1, len(process_lines[1])) if int(process_lines[1][i]) <= 0]      
signals = set([process_lines[0][i] for i in signal_indices if process_lines[0][i]])

if options.maps: 
        if 'SIGNAL:' in options.maps:
            for signal in signals:
                if options.list:
                    if not any(_item in signal for _item in options.list.split(',')): continue
                tag = options.outfile.split('/')[-1].replace('.root','').replace('SIGNAL',signal)
                print(tag)
                os.system('mkdir -p logs/condor/convert/err/')
                os.system('rm -rf logs/condor/convert/err/*'+tag+'*')
                os.system('mkdir -p logs/condor/convert/log/')
                os.system('rm -rf logs/condor/convert/log/*'+tag+'*')
                os.system('mkdir -p logs/condor/convert/out/')
                os.system('rm -rf logs/condor/convert/out/*'+tag+'*')
                os.environ['TAG'] = tag
                os.environ['CLUSTER'] = options.cluster
                os.environ['DATACARD'] = options.datacard
                os.environ['OUTFILE']  = options.outfile.replace('SIGNAL',signal)
                os.environ['MAPS']     = options.maps.replace('SIGNAL',signal).replace(' ','+')
                os.system('condor_submit convert.submit')
        else:
            tag = options.outfile.split('/')[-1].replace('.root','')
            os.system('mkdir -p logs/condor/convert/err/')
            os.system('rm -rf logs/condor/convert/err/*'+tag+'*')
            os.system('mkdir -p logs/condor/convert/log/')
            os.system('rm -rf logs/condor/convert/log/*'+tag+'*')
            os.system('mkdir -p logs/condor/convert/out/')
            os.system('rm -rf logs/condor/convert/out/*'+tag+'*')
            os.environ['TAG'] = tag
            os.environ['CLUSTER'] = options.cluster
            os.environ['DATACARD'] = options.datacard
            os.environ['OUTFILE']  = options.outfile
            os.environ['MAPS']     = options.maps.replace(' ','+')
            os.system('condor_submit convert.submit')
else:
    tag = options.outfile.split('/')[-1].replace('.root','')
    os.system('mkdir -p logs/condor/convert/err/')
    os.system('rm -rf logs/condor/convert/err/*'+tag+'*')
    os.system('mkdir -p logs/condor/convert/log/')
    os.system('rm -rf logs/condor/convert/log/*'+tag+'*')
    os.system('mkdir -p logs/condor/convert/out/')
    os.system('rm -rf logs/condor/convert/out/*'+tag+'*')
    os.environ['TAG'] = tag
    os.environ['CLUSTER'] = options.cluster
    os.environ['DATACARD'] = options.datacard
    os.environ['OUTFILE']  = options.outfile
    os.environ['MAPS']     = 'None'
    os.system('condor_submit convert.submit')
os.system('rm convert.submit')
