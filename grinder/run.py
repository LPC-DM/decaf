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
from fnal_column_analysis_tools import hist
#from saiyan import Builder
from analysis.darkhiggs import analysis,samples

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-s', '--selection', help='selection', dest='selection')
parser.add_option('-l', '--lumi', help='lumi', dest='lumi')
(options, args) = parser.parse_args()

with open("../beans/"+options.year+".json") as fin:
    datadef = json.load(fin)

#for h in hists.values(): h.clear()
hists = {}
dataset_xs = {k: v['xs'] for k,v in datadef.items()}
lumi = 1000.
if options.lumi: lumi=lumi*float(options.lumi)

tstart = time.time()
def clean(val, default):
    val[np.isnan(val)|(val==-999.)] = default
    return val

nworkers = 8
fileslice = slice(None)
with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = set()
    for dataset, info in datadef.items():
        nevents = 0
        sumw = 0
        if options.dataset and options.dataset not in dataset: continue
        for selection,v in samples.items():
            if options.selection and options.selection not in selection: continue    
            for i in range (0,len(v)):
                if v[i] not in dataset: continue
                print(dataset)
                futures.update(executor.submit(analysis, selection, options.year, dataset_xs[dataset], dataset, file) for file in info['files'][fileslice])
                if(len(futures)==0): continue
                try:
                    total = len(futures)
                    processed = 0
                    while len(futures) > 0:
                        finished = set(job for job in futures if job.done())
                        for job in finished:
                            dataset, sumws, nentries, hout = job.result()
                            nevents += nentries
                            sumw += sumws
                            for k in hout.keys():
                                if k in hists: hists[k] += hout[k]
                                else: hists[k]= hout[k]
                            print("Processing: done with % 4d / % 4d files" % (processed, total))
                            processed += 1
                        futures -= finished
                    del finished
                except KeyboardInterrupt:
                    print("Ok quitter")
                    for job in futures: job.cancel()
                except:
                    for job in futures: job.cancel()
                    raise
                print(dataset,"nevents:",nevents,"sumw:",sumw)
                scale = 1
                if dataset_xs[dataset]!= -1: scale = lumi*dataset_xs[dataset] / sumw
                print("xsec:",dataset_xs[dataset],"xsec weight:",scale)
                for h in hists.values(): h.scale(scale)
                dt = time.time() - tstart

                print("%.2f us*cpu/event" % (1e6*dt*nworkers/nevents, ))
                nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in hists.values())
                nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in hists.values())
                print("Processed %.1fM events" % (nevents/1e6, ))
                print("Filled %.1f bins" % nbins)
                print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))
                os.system("mkdir -p ../pods/"+options.year+"/"+selection)
                # Pickle is not very fast or memory efficient, will be replaced by something better soon
                with gzip.open("../pods/"+options.year+"/"+selection+"/"+dataset+".pkl.gz", "wb") as fout:
                    pickle.dump(hists, fout)

