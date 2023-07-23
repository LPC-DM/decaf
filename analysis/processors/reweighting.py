#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
import pprint
import copy
import numpy as np
import math
import awkward
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea.arrays import Initialize
from coffea import hist, processor
from coffea.util import load, save
from optparse import OptionParser
from uproot_methods import TVector2Array, TLorentzVectorArray

class AnalysisProcessor(processor.ProcessorABC):

    lumis = { #Values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable                                                      
        '2016': 36.31,
        '2017': 41.53,
        '2018': 59.74
    }

    met_filter_flags = {
        '2016': ['goodVertices',
                 'globalSuperTightHalo2016Filter',
                 'HBHENoiseFilter',
                 'HBHENoiseIsoFilter',
                 'EcalDeadCellTriggerPrimitiveFilter',
                 'BadPFMuonFilter'
             ],
        '2017': ['goodVertices',
                 'globalSuperTightHalo2016Filter',
                 'HBHENoiseFilter',
                 'HBHENoiseIsoFilter',
                 'EcalDeadCellTriggerPrimitiveFilter',
                 'BadPFMuonFilter',
                 'ecalBadCalibFilterV2'
             ],
        '2018': ['goodVertices',
                 'globalSuperTightHalo2016Filter',
                 'HBHENoiseFilter',
                 'HBHENoiseIsoFilter',
                 'EcalDeadCellTriggerPrimitiveFilter',
                 'BadPFMuonFilter',
                 'ecalBadCalibFilterV2'
             ]
    }

    def __init__(self, year, xsec, corrections, ids, common):

        self._year = year
        self._lumi = 1000.*float(AnalysisProcessor.lumis[year])
        self._xsec = xsec

        self._gentype_map = {
            'bb':       0,
            'b':        1,
            'cc' :      2,
            'c':        3,
            'other':    4,
            'hs':       5,
        }

        self._btagmu_triggers = {
            '2016': [
                'BTagMu_AK4Jet300_Mu5',
                'BTagMu_AK8Jet300_Mu5',
                #'BTagMu_AK4DiJet170_Mu5'
                ],
            '2017': [
                'BTagMu_AK4Jet300_Mu5',
                'BTagMu_AK8Jet300_Mu5',
                #'BTagMu_AK4DiJet170_Mu5'
                ],
            '2018': [
                'BTagMu_AK4Jet300_Mu5',
                'BTagMu_AK8Jet300_Mu5',
                'BTagMu_AK8Jet300_Mu5_noalgo'
                'HLT_BTagMu_AK4Jet300_Mu5_noalgo'
                #'BTagMu_AK4DiJet170_Mu5'
                ]
        }

        self._corrections = corrections
        self._ids = ids
        self._common = common

        self._accumulator = processor.dict_accumulator({
            'sumw': hist.Hist(
                'sumw',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('sumw', 'Weight value', [0.])
            ),
            'reweighting': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('tau21','tau21', [0.0, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 1.0]),
                hist.Bin('fj1pt','Leading AK15 Jet SoftDrop Pt',[350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950, 1000.0, 2500.0]),
                hist.Bin('fj1eta','Leading AK15 Jet SoftDrop Eta',[-5.0, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 5.0]),
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    def process(self, events):

        dataset = events.metadata['dataset']

        isData = 'genWeight' not in events.columns
        selection = processor.PackedSelection()
        hout = self.accumulator.identity()

        ###
        #Getting ids from .coffea files
        ###

        get_msd_weight  = self._corrections['get_msd_weight']
        get_pu_weight   = self._corrections['get_pu_weight'][self._year]  
        isSoftMuon      = self._ids['isSoftMuon']
        isGoodFatJet    = self._ids['isGoodFatJet']
        isHEMJet        = self._ids['isHEMJet']  

        match = self._common['match']

        ###
        #Initialize physics objects
        ###
        
        mu = events.Muon
        mu['issoft'] = isSoftMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.looseId,self._year)
        mu_soft=mu[mu.issoft.astype(np.bool)]

        j = events.Jet
        j['isHEM'] = isHEMJet(j.pt, j.eta, j.phi)
        j_HEM = j[j.isHEM.astype(np.bool)]
        j_nHEM = j_HEM.counts

        fj = events.AK15Puppi
        fj['sd'] = fj.subjets.sum()
        fj['isgood'] = isGoodFatJet(fj.sd.pt, fj.sd.eta, fj.jetId)
        fj['tau21'] = fj.tau2/fj.tau1
        jetmu = fj.sd.cross(mu_soft, nested=True)
        fj['withmu'] = (mu_soft.counts>0) & ((jetmu.i0.delta_r(jetmu.i1) < 1.5).astype(np.int).sum()>0)
        fj_good = fj[fj.isgood.astype(np.bool)]
        fj_withmu = fj_good[fj_good.withmu.astype(np.bool)]
        fj_nwithmu = fj_withmu.counts
        
        ###
        # Calculating weights
        ###
        
        if not isData:
            pu = get_pu_weight(events.Pileup.nTrueInt)

        #### ak15 jet selection ####
        leading_fj = fj[fj.sd.pt.argmax()]
        leading_fj = leading_fj[leading_fj.isgood.astype(np.bool)]
        
        ###
        # Selections
        ###

        #### trigger selection ####
        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._btagmu_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('btagmu_triggers', triggers)

        #### MET filters ####
        met_filters =  np.ones(events.size, dtype=np.bool)
        if isData:
            met_filters = met_filters & events.Flag['eeBadScFilter'] #this filter is recommended for data only
        for flag in AnalysisProcessor.met_filter_flags[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters',met_filters)

        noHEMj = np.ones(events.size, dtype=np.bool)
        if self._year=='2018': noHEMj = (j_nHEM==0)

        selection.add('noHEMj', noHEMj)
        selection.add('fj_pt', (leading_fj.sd.pt.max() > 350) )
        selection.add('fj_mass', (leading_fj.msd_corr.sum() > 40) ) ## optionally also <130
        selection.add('nwithmu', (fj_nwithmu>0))
        selection.add('fj_withmu', (leading_fj.withmu.sum().astype(np.bool)))

        isFilled = False
        if isData:
            if not isFilled:
                hout['sumw'].fill(dataset=dataset, sumw=1, weight=1)
                isFilled=True

            cut = selection.all(*selection.names)
            ##### template for bb SF #####
            hout['reweighting'].fill(dataset=dataset,
                                    tau21=leading_fj.tau21.sum(),
                                    fj1pt=leading_fj.sd.pt.sum(),
                                    fj1eta=leading_fj.sd.eta.sum(),
                                    weight=np.ones(events.size)*cut)
        else:
            weights = processor.Weights(len(events))
            if 'L1PreFiringWeight' in events.columns: weights.add('prefiring',events.L1PreFiringWeight.Nom)
            weights.add('genw',events.genWeight)
            weights.add('pileup',pu)
            cut = selection.all(*selection.names)

            if 'QCD' in dataset:
            if not isFilled:
                hout['sumw'].fill(dataset=dataset, sumw=1, weight=events.genWeight.sum())
                isFilled=True
            hout['reweighting'].fill(dataset=dataset,
                                    tau21=leading_fj.tau21.sum(),
                                    fj1pt=leading_fj.sd.pt.sum(),
                                    fj1eta=leading_fj.sd.eta.sum(),
                                    weight=weights.weight()*cut)
            
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
            if isinstance(h, hist.Hist):
                h.scale(scale, axis='dataset')

        return accumulator

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
    parser.add_option('-n', '--name', help='name', dest='name')
    (options, args) = parser.parse_args()


    with open('metadata/'+options.metadata+'.json') as fin:
        samplefiles = json.load(fin)
        xsec = {k: v['xs'] for k,v in samplefiles.items()}

    corrections = load('data/corrections.coffea')
    ids         = load('data/ids.coffea')
    common      = load('data/common.coffea')

    processor_instance=AnalysisProcessor(year=options.year,
                                         xsec=xsec,
                                         corrections=corrections,
                                         ids=ids,
                                         common=common)

    save(processor_instance, 'data/reweighting'+options.name+'.processor')
