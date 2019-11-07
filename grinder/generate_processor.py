#!/usr/bin/env python
import json
import gzip
import uproot
import numexpr
import numpy as np
from coffea.util import load, save
from coffea import hist, processor
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-p', '--processor', help='processor', dest='processor')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-l', '--lumi', help='lumi', dest='lumi', type=float)
(options, args) = parser.parse_args()

if not options.processor:
        print('Need to specify the processor file to generate')
        exit()

lumis = {} #Values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable                                                    
lumis['2016']=35.92
lumis['2017']=41.53
lumis['2018']=59.97
lumi = 1000.*float(lumis[options.year])
if options.lumi: lumi=1000.*options.lumi

with open("../harvester/beans/"+options.year+".json") as fin:
    samplefiles = json.load(fin)
xsec = {k: v['xs'] for k,v in samplefiles.items()}

if 'darkhiggs' in options.processor: 

        from analysis.darkhiggs import AnalysisProcessor
        corrections = load('corrections.coffea')
        triggers    = load('triggers.coffea')
        ids         = load('ids.coffea')
        metfilters  = load('metfilters.coffea')
        
        processor_instance=AnalysisProcessor(year=options.year, 
                                             xsec=xsec, 
                                             lumi=lumi, 
                                             corrections=corrections, 
                                             triggers=triggers, 
                                             ids=ids, 
                                             metfilters=metfilters)

save(processor_instance, options.processor+'.coffea')
