#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import pprint
import numpy as np
from coffea import hist, processor
from optparse import OptionParser
from analysis.darkhiggs import AnalysisProcessor

parser = OptionParser()
parser.add_option('-y', '--year', help='year', dest='year')
(options, args) = parser.parse_args()

with open("../harvester/beans/"+options.year+".json") as fin:
    samplefiles = json.load(fin)
    
xsec = {k: v['xs'] for k,v in samplefiles.items()}
filelist = {}
selections = []
fileslice = slice(None)

for dataset, info in samplefiles.items():
    if options.dataset and options.dataset not in dataset: continue
    files = []
    for file in info['files'][fileslice]:
        files.append(file)
    filelist[dataset] = files
    for selection,v in samples.items():
        for i in range (0,len(v)):
            if v[i] not in dataset: continue
            selections.append(selection)

    processor_instance=AnalysisProcessor(selected_regions=selections, year=options.year, xsec=xsec, lumi=lumi)
    tstart = time.time()
    output = processor.run_uproot_job(filelist,
                                      treename='Events',
                                      processor_instance=processor_instance,
                                      executor=processor.futures_executor,
                                      executor_args={'workers': options.workers, 'pre_workers': 1},
                                      chunksize=500000,
                                      )
    processor_instance.postprocess(output)
    
    nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in output.values() if isinstance(h, hist.Hist))
    nfilled = sum(sum(np.sum(arr > 0) for arr in h._sumw.values()) for h in output.values() if isinstance(h, hist.Hist))
    print("Filled %.1fM bins" % (nbins/1e6, ))
    print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

    # Pickle is not very fast or memory efficient, will be replaced by something better soon
    #    with lz4f.open("pods/"+options.year+"/"+dataset+".pkl.gz", mode="xb", compression_level=5) as fout:
    with gzip.open("pods/"+options.year+"/"+dataset+".pkl.gz", "wb") as fout:
        cloudpickle.dump(output, fout)
        
    dt = time.time() - tstart
    nworkers = options.workers
    print("%.2f us*cpu overall" % (1e6*dt*nworkers, ))
