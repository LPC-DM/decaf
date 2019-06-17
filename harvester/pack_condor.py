#!/usr/bin/env python
import json
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-e', '--exclude', help='exclude', dest='exclude')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-t', '--tar', action="store_true", dest="tar")
(options, args) = parser.parse_args()

os.system("mkdir -p condor/"+options.year+"/out condor/"+options.year+"/err condor/"+options.year+"/log")
os.system("rm -rf condor/"+options.year+"/out/* condor/"+options.year+"/err/* condor/"+options.year+"/log/*")

if options.tar:
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../decaf.tgz ../../decaf')
    os.system('tar --exclude-caches-all --exclude-vcs -czvf ../../pylocal.tgz -C ~/.local/lib/python3.6/ site-packages')
    os.system('xrdcp -f ../../decaf.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/decaf.tgz')
    os.system('xrdcp -f ../../pylocal.tgz root://cmseos.fnal.gov//store/user/'+os.environ['USER']+'/pylocal.tgz')

with open("beans/"+options.year+".json") as fin:
    datadef = json.load(fin)

for dataset, info in datadef.items():
    if options.dataset and options.dataset not in dataset: continue
    if options.exclude and options.exclude in dataset: continue
    os.environ['SAMPLE'] = dataset
    os.environ['YEAR']   = options.year
    os.system("condor_submit pack")
