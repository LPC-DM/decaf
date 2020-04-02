#!/usr/bin/env python
import lz4.frame as lz4f
import pickle
import json
import time
import cloudpickle
import gzip
import os
from optparse import OptionParser

import uproot
import numpy as np
from coffea import hist, processor
from coffea.util import load, save

parser = OptionParser()
parser.add_option('-a', '--analysis', help='analysis', dest='analysis')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
(options, args) = parser.parse_args()

year=''
if options.year: year=options.year

processor_file = ''
for filename in os.listdir('data'):
    if '.processor' not in filename: continue
    if options.analysis+year in filename: processor_file = filename

processor_instance=load('data/'+processor_file)

fileslice = slice(None)
with open("metadata/"+options.year+".json") as fin:
    samplefiles = json.load(fin)

filelist = {}
for dataset, info in samplefiles.items():
    if options.dataset and options.dataset not in dataset: continue
    files = []
    for file in info['files'][fileslice]:
        files.append(file)
    filelist[dataset] = files

from distributed import Client
client = Client('coffea-dask.fnal.gov:8786')
tstart = time.time()
output = processor.run_uproot_job(filelist,
                                  treename='Events',
                                  processor_instance=processor_instance,
                                  executor=processor.dask_executor,
                                  executor_args={'client': client},
                                  chunksize=50000,
                                  )

# Pickle is not very fast or memory efficient, will be replaced by something better soon
#    with lz4f.open("pods/"+options.year+"/"+dataset+".pkl.gz", mode="xb", compression_level=5) as fout:
os.system("mkdir -p hists/"+options.analysis+year)
save(output,'hists/'+options.analysis+year+'/'+dataset+'.dask')
dt = time.time() - tstart
print(dt)
