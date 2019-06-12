#!/usr/bin/env python
import lz4.frame as lz4f
import pickle
import json
import time
import cloudpickle
import argparse

import uproot
import numpy as np
from fnal_column_analysis_tools import hist, processor



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run analysis on baconbits files using processor cloudpickle files')
    parser.add_argument('--processor', default='boostedHbbProcessor.cpkl.lz4', help='The name of the compiled processor file')
    parser.add_argument('--output', default='hists.cpkl.lz4', help='Output histogram filename')
    parser.add_argument('--samplejson', default='metadata/samplefiles.json', help='JSON file containing dataset and file locations')
    parser.add_argument('--sample', default='test_skim', help='The sample to use in the sample JSON')
    parser.add_argument('--limit', type=int, default=None, metavar='N', help='Limit to the first N files of each dataset in sample JSON')
    parser.add_argument('--executor', choices=['iterative', 'futures'], default='iterative', help='The type of executor to use')
    parser.add_argument('-j', '--workers', type=int, default=12, help='Number of workers to use for multi-worker executors (e.g. futures or condor)')
    parser.add_argument('--profile-out', dest='profilehtml', default=None, help='Filename for the pyinstrument HTML profile output')
    args = parser.parse_args()

    with open(args.samplejson) as fin:
        samplefiles = json.load(fin)
    sample = samplefiles[args.sample]
    filelist = []
    for dataset, files in sample.items():
        for file in files[:args.limit]:
            filelist.append((dataset, file))

    tstart = time.time()
    output = processor.run_uproot_job(filelist,
                                  treename='Events',
                                  processor_instance=BoostedHbbProcessor(corrections=corrections, columns=allcolumns, debug=args.debug, year=args.year),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 8, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )
    processor_instance.postprocess(output)

    nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in final_accumulator.values() if isinstance(h, hist.Hist))
    nfilled = sum(sum(np.sum(arr > 0) for arr in h._sumw.values()) for h in final_accumulator.values() if isinstance(h, hist.Hist))
    print("Filled %.1fM bins" % (nbins/1e6, ))
    print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

    # Pickle is not very fast or memory efficient, will be replaced by something better soon
    with lz4f.open(args.output, mode="wb", compression_level=5) as fout:
        cloudpickle.dump(final_accumulator, fout)

    dt = time.time() - tstart
    nworkers = 1 if args.executor == 'iterative' else args.workers
    print("%.2f us*cpu overall" % (1e6*dt*nworkers, ))
