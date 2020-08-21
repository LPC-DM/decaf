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

    lumis = { #Values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable                                                      
        '2016': 35.92,
        '2017': 41.53,
        '2018': 59.74
    }
    
    def __init__(self, year, xsec, corrections):
        self._year = year
        self._lumi = 1000.*float(AnalysisProcessor.lumis[year])
        self._xsec = xsec
        self._corrections = corrections
        self._accumulator = processor.dict_accumulator({
            'sumw': hist.Hist(
                'sumw', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Bin('sumw', 'Weight value', [0.])
            ),
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
        get_nnlo_nlo_weight     = self._corrections['get_nnlo_nlo_weight'][self._year]

        gen = events.GenPart
        gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
        gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
        gen['isA'] = (abs(gen.pdgId)==22)&gen.hasFlags(['isPrompt', 'fromHardProcess', 'isLastCopy'])&(gen.status==1)
        
        ###
        # Calculating gen photon dynamic isolation as in https://arxiv.org/pdf/1705.04664.pdf
        ###

        epsilon_0_dyn = 0.1
        n_dyn = 1
        gen['R_dyn'] = (91.1876/(gen.pt * np.sqrt(epsilon_0_dyn)))*(gen.isA).astype(np.int) + (-999)*(~gen.isA).astype(np.int)
        gen['R_0_dyn'] = gen.R_dyn*(gen.R_dyn<1.0).astype(np.int) + (gen.R_dyn>=1.0).astype(np.int)
        
        def isolation(R):
            hadrons = gen[ #Stable hadrons not in NanoAOD, using quarks/glouns instead
                ((abs(gen.pdgId)<=5)|(abs(gen.pdgId)==21)) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy'])
            ]
            genhadrons = gen.cross(hadrons, nested=True)
            hadronic_et = genhadrons.i1[(genhadrons.i0.delta_r(genhadrons.i1) <= R)].pt.sum()
            return (hadronic_et<=(epsilon_0_dyn * gen.pt * np.power((1 - np.cos(R)) / (1 - np.cos(gen.R_0_dyn)), n_dyn))) | (hadrons.counts==0)

        isIsoA=gen.isA
        iterations = 5.
        for i in range(1, int(iterations) + 1):
            isIsoA=isIsoA&isolation(gen.R_0_dyn*i/iterations)
        gen['isIsoA']=isIsoA

        genWs = gen[gen.isW&(gen.pt>100)]
        genZs = gen[gen.isZ&(gen.pt>100)]
        genIsoAs = gen[gen.isIsoA&(gen.pt>150)] #Based on photon weight distribution

        nnlo_nlo = np.ones(events.size)
        lo = np.ones(events.size)
        if('GJets' in dataset): 
            nnlo_nlo=get_nnlo_nlo_weight['a']['cen'](genIsoAs.pt.max())*(genIsoAs.counts>0)
            lo=(genIsoAs.counts>0).astype(np.int)
        elif('WJets' in dataset): 
            nnlo_nlo=get_nnlo_nlo_weight['w']['cen'](genWs.pt.max())*(genWs.counts>0)
            lo=(genWs.counts>0).astype(np.int)
        elif('DY' in dataset):
            nnlo_nlo=get_nnlo_nlo_weight['dy']['cen'](genZs.pt.max())*(genZs.counts>0)
            lo=(genZs.counts>0).astype(np.int)
        elif('ZJets' in dataset): 
            nnlo_nlo=get_nnlo_nlo_weight['z']['cen'](genZs.pt.max())*(genZs.counts>0)
            lo=(genZs.counts>0).astype(np.int)
        hout['sumw'].fill(dataset='lo--'+dataset, sumw=1, weight=events.genWeight.sum())
        hout['sumw'].fill(dataset='nnlo--'+dataset, sumw=1, weight=events.genWeight.sum())
        hout['yields'].fill(dataset='lo--'+dataset, yields=np.zeros(events.size), weight=events.genWeight*lo)
        hout['yields'].fill(dataset='nnlo--'+dataset, yields=np.zeros(events.size), weight=events.genWeight*nnlo_nlo)
        return hout

    def postprocess(self, accumulator):
        scale = {}
        for d in accumulator['sumw'].identifiers('dataset'):
            print('Scaling:',d.name)
            dataset = d.name
            if '--' in dataset: dataset = dataset.split('--')[1]
            print('Cross section:',self._xsec[dataset])
            if self._xsec[dataset]!= -1: scale[d.name] = self._lumi*self._xsec[dataset]
            else: scale[d.name] = 1

        for histname, h in accumulator.items():
            if histname == 'sumw': continue
            #if histname == 'yields': continue
            if isinstance(h, hist.Hist):
                print(histname)
                h.scale(scale, axis='dataset')

        return accumulator

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    (options, args) = parser.parse_args()


    with open('metadata/'+options.year+'.json') as fin:
        samplefiles = json.load(fin)
        xsec = {k: v['xs'] for k,v in samplefiles.items()}

    corrections = load('data/corrections.coffea')
    processor_instance=AnalysisProcessor(year=options.year,
                                         xsec=xsec,
                                         corrections=corrections)
    save(processor_instance, 'data/kfactors'+options.year+'.processor')
