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
from coffea.nanoaod.methods import collection_methods, LorentzVector

collection_methods['AK15Puppi'] = LorentzVector
collection_methods['AK15PuppiSubJet'] = LorentzVector

parser = OptionParser()
parser.add_option('-p', '--processor', help='processor', dest='processor')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-w', '--workers', help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=8)
(options, args) = parser.parse_args()

processor_instance=load('data/'+options.processor+'.processor')

fileslice = slice(None)
with open("metadata/"+options.metadata+".json") as fin:
    samplefiles = json.load(fin)

for dataset, info in samplefiles.items():
    filelist = {}
    if options.dataset and options.dataset not in dataset: continue
    print('Processing:',dataset)
    files = []
    for file in info['files'][fileslice]:
        files.append(file)
    filelist[dataset] = files

    tstart = time.time()
    output = processor.run_uproot_job(filelist,
                                      treename='Events',
                                      processor_instance=processor_instance,
                                      executor=processor.futures_executor,
                                      executor_args={'nano': True, 'workers': options.workers},
                                      )
    
    #nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in output.values() if isinstance(h, hist.Hist))
    #nfilled = sum(sum(np.sum(arr > 0) for arr in h._sumw.values()) for h in output.values() if isinstance(h, hist.Hist))
    #print("Filled %.1fM bins" % (nbins/1e6, ))
    #print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

    os.system("mkdir -p hists/"+options.processor)
    save(output,'hists/'+options.processor+'/'+dataset+'.futures')        
    dt = time.time() - tstart
    nworkers = options.workers
    print("%.2f us*cpu overall" % (1e6*dt*nworkers, ))
