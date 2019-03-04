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
from optparse import OptionParser
import uproot, uproot_methods
import numpy as np
from fnal_column_analysis_tools import hist
#from saiyan import Builder
from darkhiggs import analysis

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
(options, args) = parser.parse_args()

with open("data/coffeabeans2016.json") as fin:
    datadef = json.load(fin)

dataset_xs = {k: v['xs'] for k,v in datadef.items()}
lumi = 1000.  # [1/pb]
monojet_recoil_binning = [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]

tstart = time.time()
nevents = defaultdict(lambda: 0.)

def clean(val, default):
    val[np.isnan(val)|(val==-999.)] = default
    return val

nworkers = 15
fileslice = slice(None)
with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = set()
    for dataset, info in datadef.items():
        hists = {}
        hists['recoil'] = hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("recoil","Hadronic Recoil",monojet_recoil_binning))
        for h in hists.values(): h.clear()
        if options.dataset:
            if options.dataset in dataset:
                futures.update(executor.submit(analysis, dataset, hists, file) for file in info['files'][fileslice])
            else:
                continue
        else:
            futures.update(executor.submit(analysis, dataset, hists, file) for file in info['files'][fileslice])
        try:
            total = len(futures)
            processed = 0
            while len(futures) > 0:
                finished = set(job for job in futures if job.done())
                for job in finished:
                    dataset, sumw, nentries, hout = job.result()
                    nevents[dataset] += nentries
                    for k in hout.keys():
                        hists[k] += hout[k]
                    processed += 1
                    print("Processing: done with % 4d / % 4d files" % (processed, total))
                futures -= finished
            del finished
        except KeyboardInterrupt:
            print("Ok quitter")
            for job in futures: job.cancel()
        except:
            for job in futures: job.cancel()
            raise


        scale = {}
        print(dataset,"nevents:",nevents[dataset],"sumw:",sumw)
        scale[dataset] = lumi*dataset_xs[dataset] / sumw

        for h in hists.values(): h.scale(scale, axis="dataset")

        dt = time.time() - tstart

        print("%.2f us*cpu/event" % (1e6*dt*nworkers/sum(nevents.values()), ))
        nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in hists.values())
        nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in hists.values())
        print("Processed %.1fM events" % (sum(nevents.values())/1e6, ))
        print("Filled %.1f bins" % nbins)
        print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

        # Pickle is not very fast or memory efficient, will be replaced by something better soon
        with gzip.open("hists/"+dataset+".pkl.gz", "wb") as fout:
            pickle.dump(hists, fout)

