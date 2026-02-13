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

    def __init__(self, year):
        self._year = year
        self._accumulator = processor.dict_accumulator({
            'yields': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('yields', 'Yield', [0, 1])
            ),
            'met': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('met','MET',30,0,600)
            ),
            'metphi': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('metphi','MET phi',35,-3.5,3.5)
            ),
            'mindphimet': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('mindphimet','Min dPhi(MET,AK4s)',30,0,3.5)
            ),
            'j1pt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('j1pt','AK4 Leading Jet Pt',[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])
            ),
            'j1eta': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('j1eta','AK4 Leading Jet Eta',35,-3.5,3.5)
            ),
            'j1phi': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('j1phi','AK4 Leading Jet Phi',35,-3.5,3.5)
            ),
            'njets': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('njets','AK4 Number of Jets',6,-0.5,5.5)
            ),
            'medmass': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('medmass','Zprime mass', 50, 0, 3000)
            ),
            'medpt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('medpt','Zprime pT', 30, 0, 1500)
            ),
            'medeta': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('medeta','Zprime Eta',35,-3.5,3.5)
            ),
            'dmmass': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('dmmass','Dark matter mass', 25, 0, 1500)
            ),
            'dmpt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('dmpt','Dark matter pT', 30, 0, 1500)
            ),
            'dmeta': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('dmeta','Dark matter Eta',35,-3.5,3.5)
            ),
            'hsmass': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('hsmass','Dark Higgs mass', 10, 0, 100)
            ),
            'hspt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('hspt','Dark Higgs pT', 30, 0, 1500)
            ),
            'hseta': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('hseta','Dark Higgs Eta',35,-3.5,3.5)
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

        ###
        #Initialize global quantities (MET ecc.)
        ###

        met = events.MET
        if self._year == '2017': events.METFixEE2017#Recommended for 2017
        met['T']  = TVector2Array.from_polar(met.pt, met.phi)

        ###
        #Initialize physics objects
        ###

        j = events.Jet
        j['T'] = TVector2Array.from_polar(j.pt, j.phi)
        j_ntot=j.counts
        leading_j = j[j.pt.argmax()]

        ### Gen information
        # pdgId = 52: dark matter, 54: dark Higgs, 55: Zprime
        ###
        gen = events.GenPart
        zp = gen[gen.pdgId==55]
        dm = gen[gen.pdgId==52]
        hs = gen[gen.pdgId==54]

        ### Store variables as set 
        variables = {
            'met':                    met.pt,
            'metphi':                 met.phi,
            'mindphimet':             abs(met.T.delta_phi(j.T)).min(),
            'j1pt':                   leading_j.pt,
            'j1eta':                  leading_j.eta,
            'j1phi':                  leading_j.phi,
            'njets':                  j_ntot,
            'medmass':                zp.mass,
            'medpt':                  zp.pt,
            'medeta':                 zp.eta,
            'dmmass':                 dm.mass,
            'dmpt':                   dm.pt,
            'dmeta':                  dm.eta,
            'hsmass':                 hs.mass,
            'hspt':                   hs.pt,
            'hseta':                  hs.eta,
        }
        print('Variables:',variables.keys())
        flat_variables = {k: v[mask].flatten() for k, v in variables.items()}
        flat_weight = {k: ~np.isnan(v[mask]).flatten() for k, v in variables.items()}

        for histname, h in hout.items():
            if not isinstance(h, hist.Hist):
                continue
            if histname not in variables:
                continue
            elif histname == 'yields':
                hout['yields'].fill(dataset=dataset_name, yields=np.zeros(events.size), weight=mask.astype(np.int))
            else:
                flat_variable = {histname: flat_variables[histname]}
                h.fill(dataset=dataset, **flat_variable, weight=flat_weight[histname])

        return hout

    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    (options, args) = parser.parse_args()

    processor_instance=AnalysisProcessor(year=options.year)
    save(processor_instance, 'data/signals'+options.year+'.processor')
