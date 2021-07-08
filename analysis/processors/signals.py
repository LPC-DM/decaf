#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
import pprint
import numpy as np
import math 
import awkward
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea.arrays import Initialize
from coffea import hist, processor
from coffea.util import load, save
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty, JetTransformer, JetResolution, JetResolutionScaleFactor
from optparse import OptionParser
from uproot_methods import TVector2Array, TLorentzVectorArray

class AnalysisProcessor(processor.ProcessorABC):
    
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            'yields': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('yields', 'Yield', [0, 1])
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator
            
    def process(self, events):
        dataset = events.metadata['dataset']
        hout = self.accumulator.identity()
        isCentral = 'GenModel' in events.columns
        mask=np.zeros(events.size, dtype=np.bool)
        dataset_name=dataset
        if isCentral:
            #mask=events.GenModel.DarkHiggs_MonoHs_LO_MZprime_2000_Mh_50_Mchi_1500_TuneCP5_13TeV-madgraph-pythia8
            column="DarkHiggs_MonoHs_LO_MZprime-500_Mhs-50_Mchi-150_TuneCP5_13TeV-madgraph-pythia8"
            if column in events.GenModel.columns:
                mask=events.GenModel[column]
            dataset_name=dataset
        else:
            mask=np.ones(events.size, dtype=np.bool)
        hout['yields'].fill(dataset=dataset_name, yields=np.zeros(events.size), weight=mask.astype(np.int))
        return hout

    def postprocess(self, accumulator):

        return accumulator

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    (options, args) = parser.parse_args()

    processor_instance=AnalysisProcessor()
    save(processor_instance, 'data/signals'+options.year+'.processor')
