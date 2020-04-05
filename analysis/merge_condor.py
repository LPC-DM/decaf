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
parser.add_option('-t', '--tar', action='store_true', dest='tar')
parser.add_option('-k', '--kisti', action='store_true', dest='kisti')
parser.add_option('-r', '--reduce', action='store_true', dest='reduce')
(options, args) = parser.parse_args()

os.system('mkdir -p '+options.folder+'/merge_condor/out '+options.folder+'/merge_condor/err '+options.folder+'/merge_condor/log')
os.system('rm -rf '+options.folder+'/merge_condor/err/'+options.dataset+'*')
os.system('rm -rf '+options.folder+'/merge_condor/log/'+options.dataset+'*')
os.system('rm -rf '+options.folder+'/merge_condor/out/'+options.dataset+'*')

jdl = 'merge'
if options.kisti:  jdl = 'merge_kisti'
if options.reduce: 
    jdl = 'reduce'
    if options.kisti: jdl = 'reduce_kisti'
jdl=jdl+'.jdl'
os.system('cp jdls/'+jdl+' .')
print('Using',jdl)

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz --exclude=\'analysis/hists/*/*condor/*/*\' ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')
    if options.kisti: 
        os.system('xrdcp -f ../../decaf.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cms-xrdr.private.lo:2094//xrd/store/user/'+os.environ['USER']+'/pylocal.tgz')
    else:
        os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
        os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')

pd = []
for filename in os.listdir(options.folder):
    if '.futures' not in filename: continue
    if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])

variables = [
    'sumw',
    'CaloMinusPfOverRecoil',
    'recoil',
    'met',
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

if options.reduce:
    for variable in variables:
        if options.variable and options.variable not in variable: continue
        os.environ['FOLDER'] = options.folder
        os.environ['VARIABLE'] = variable
        os.system('condor_submit '+jdl)
else:
    for pdi in pd:
        if options.dataset and options.dataset not in pdi: continue
        if options.exclude and options.exclude in pdi: continue
        for variable in variables:
            if options.variable and options.variable not in variable: continue
            #print(variable)
            os.environ['FOLDER'] = options.folder
            os.environ['SAMPLE'] = pdi
            os.environ['VARIABLE'] = variable
            os.system('condor_submit '+jdl)
os.system('rm '+jdl)
