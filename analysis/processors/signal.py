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
from optparse import OptionParser
from uproot_methods import TVector2Array, TLorentzVectorArray

class AnalysisProcessor(processor.ProcessorABC):

    def __init__(self, year, ids, common):

        self._columns = """
        AK15PuppiSubJet_eta
        AK15PuppiSubJet_mass
        AK15PuppiSubJet_phi
        AK15PuppiSubJet_pt
        AK15Puppi_eta
        AK15Puppi_jetId
        AK15Puppi_msoftdrop
        AK15Puppi_phi
        AK15Puppi_probHbb
        AK15Puppi_probQCDb
        AK15Puppi_probQCDbb
        AK15Puppi_probQCDc
        AK15Puppi_probQCDcc
        AK15Puppi_probQCDothers
        AK15Puppi_probZbb
        AK15Puppi_pt
        AK15Puppi_subJetIdx1
        AK15Puppi_subJetIdx2
        CaloMET_pt
        CaloMET_phi
        Electron_charge
        Electron_cutBased
        Electron_dxy
        Electron_dz
        Electron_eta
        Electron_mass
        Electron_phi
        Electron_pt
        Flag_BadPFMuonFilter
        Flag_EcalDeadCellTriggerPrimitiveFilter
        Flag_HBHENoiseFilter
        Flag_HBHENoiseIsoFilter
        Flag_globalSuperTightHalo2016Filter
        Flag_goodVertices
        GenPart_eta
        GenPart_genPartIdxMother
        GenPart_pdgId
        GenPart_phi
        GenPart_pt
        GenPart_statusFlags
        HLT_Ele115_CaloIdVT_GsfTrkIdT
        HLT_Ele32_WPTight_Gsf
        HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
        HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60
        HLT_Photon200
        Jet_btagDeepB
        Jet_btagDeepFlavB
        Jet_chEmEF
        Jet_chHEF
        Jet_eta
        Jet_hadronFlavour
        Jet_jetId
        Jet_mass
        Jet_neEmEF
        Jet_neHEF
        Jet_phi
        Jet_pt
        Jet_rawFactor
        MET_phi
        MET_pt
        Muon_charge
        Muon_eta
        Muon_looseId
        Muon_mass
        Muon_pfRelIso04_all
        Muon_phi
        Muon_pt
        Muon_tightId
        PV_npvs
        Photon_eta
        Photon_phi
        Photon_pt
        Tau_eta
        Tau_idDecayMode
        Tau_idMVAoldDM2017v2
        Tau_phi
        Tau_pt
        fixedGridRhoFastjetAll
        genWeight
        nAK15Puppi
        nAK15PuppiSubJet
        nElectron
        nGenPart
        nJet
        nMuon
        nPhoton
        nTau
        """.split()

        self._year = year

        self._ids = ids
        self._common = common

        self._accumulator = processor.dict_accumulator({
            'sumw': hist.Hist(
                'sumw',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('sumw', 'Weight value', [0.])
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

    @property
    def columns(self):
        return self._columns

    def process(self, events):

        dataset = events.metadata['dataset']

        hout = self.accumulator.identity()
        match = self._common['match']

        ###
        #Initialize global quantities (MET ecc.)
        ###

        met = events.MET
        if self._year == '2017': events.METFixEE2017#Recommended for 2017
        met['T']  = TVector2Array.from_polar(met.pt, met.phi)
        calomet = events.CaloMET
        puppimet = events.PuppiMET

        ###
        #Initialize physics objects
        ###

        mu = events.Muon
        mu['T'] = TVector2Array.from_polar(mu.pt, mu.phi)
        mu_ntot = mu.counts
        leading_mu = mu[mu.pt.argmax()]

        e = events.Electron
        e['isclean'] = ~match(e,mu_loose,0.3)
        e['T'] = TVector2Array.from_polar(e.pt, e.phi)
        e_clean = e[e.isclean.astype(np.bool)]
        e_ntot = e.counts
        leading_e = e[e.pt.argmax()]
        leading_e = leading_e[leading_e.isclean.astype(np.bool)]

        pho = events.Photon
        pho['isclean']=~match(pho,mu_loose,0.5)&~match(pho,e_loose,0.5)
        _id = 'cutBasedBitmap'
        if self._year=='2016':
            _id = 'cutBased'
        pho['T'] = TVector2Array.from_polar(pho.pt, pho.phi)
        pho_clean=pho[pho.isclean.astype(np.bool)]
        pho_ntot=pho.counts
        leading_pho = leading_pho[leading_pho.isclean.astype(np.bool)]

        j = events.Jet
        j['isclean'] = ~match(j,e_loose,0.4)&~match(j,mu_loose,0.4)&~match(j,pho_loose,0.4)
        j['T'] = TVector2Array.from_polar(j.pt, j.phi)
        j_clean = j[j.isclean.astype(np.bool)]
        j_ntot=j.counts
        j_nclean=j_clean.counts
        leading_j = j[j.pt.argmax()]
        leading_j = leading_j[leading_j.isclean.astype(np.bool)]

        ### Gen information
        # pdgId = 52: dark matter, 54: dark Higgs, 55: Zprime
        ###
        gen = events.GenPart
        zp = gen[gen.pdgId==55]
        dm = gen[gen.pdgId==52]
        hs = gen[gen.pdgId==54]

        ### Official Randomnized parameter sample ####
        # Not required for prtivate sample
        # Need to add GenModel flag
        # Example: GenModel_DarkHiggs_MonoHs_LO_MZprime_300_Mhs_50_Mchi_100_TuneCP5_13TeV_madgraph_pythia8
        ###

        variables = {
            'met':                    met.pt,
            'metphi':                 met.phi,
            'mindphimet':             abs(met.T.delta_phi(j_clean.T)).min(),
            'j1pt':                   leading_j.pt,
            'j1eta':                  leading_j.eta,
            'j1phi':                  leading_j.phi,
            'njets':                  j_nclean,
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

        flat_variables = {k: v.flatten() for k, v in variables.items()}

        for histname, h in hout.items():
            if not isinstance(h, hist.Hist):
                continue
            if histname not in variables:
                continue
            elif histname == 'sumw':
                hout['sumw'].fill(dataset=dataset, sumw=1, weight=1)
                #hout['sumw'].fill(dataset=dataset, sumw=1, weight=events.genWeight.sum())
            else:
                flat_variable = {histname: flat_variables[histname]}
                h.fill(dataset=dataset, **flat_variable)

        return hout

    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    (options, args) = parser.parse_args()

    with open('metadata/'+options.year+'.json') as fin:
        samplefiles = json.load(fin)

    ids         = load('data/ids.coffea')
    common      = load('data/common.coffea')

    processor_instance=AnalysisProcessor(year=options.year,
                                         ids=ids,
                                         common=common)

    save(processor_instance, 'data/signal'+options.year+'.processor')
