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
        GenPart_pdgIdGenPart_phi
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

        self._lumi = 1000.*float(AnalysisProcessor.lumis[year])

        self._xsec = xsec

        self._samples = {
            'sr':('ZJets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET','Mhs_50','Mhs_70','Mhs_90','MonoJet','MonoW','MonoZ','mhs'),
            'wmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'tmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'wecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'tecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'zmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'zecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'gcr':('GJets','QCD','SinglePhoton','EGamma')
        }

        self._gentype_map = {
            'xbb':      1,
            'tbcq':     2,
            'tbqq':     3,
            'zcc':      4,
            'wcq':      5,
            'vqq':      6,
            'bb':       7,
            'bc':       8,
            'b':        9,
            'cc' :     10,
            'c':       11,
            'other':   12
            #'garbage': 13
        }
        
        self._ZHbbvsQCDwp = {
            '2016': 0.53,
            '2017': 0.61,
            '2018': 0.65
        }

        self._met_triggers = {
            '2016': [
                'PFMETNoMu90_PFMHTNoMu90_IDTight',
                'PFMETNoMu100_PFMHTNoMu100_IDTight',
                'PFMETNoMu110_PFMHTNoMu110_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2017': [
                'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2018': [
                'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ]
        }

        self._singlephoton_triggers = {
            '2016': [
                'Photon175',
                'Photon165_HE10'
            ],
            '2017': [
                'Photon200'
            ],
            '2018': [
                'Photon200'
            ]
        }

        self._singleelectron_triggers = { #2017 and 2018 from monojet, applying dedicated trigger weights
            '2016': [
                'Ele27_WPTight_Gsf',
                'Ele105_CaloIdVT_GsfTrkIdT'
            ],
            '2017': [
                'Ele35_WPTight_Gsf',
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Photon200'
            ],
            '2018': [
                'Ele32_WPTight_Gsf',
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Photon200'
            ]
        }

        self._jec = {
        
            '2016': [
                'Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFPuppi',
                'Summer16_07Aug2017_V11_MC_L2L3Residual_AK4PFPuppi',
                'Summer16_07Aug2017_V11_MC_L2Relative_AK4PFPuppi',
                'Summer16_07Aug2017_V11_MC_L2Residual_AK4PFPuppi',
                'Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFPuppi'
            ],
            
            '2017':[
                'Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFPuppi',
                'Fall17_17Nov2017_V32_MC_L2L3Residual_AK4PFPuppi',
                'Fall17_17Nov2017_V32_MC_L2Relative_AK4PFPuppi',
                'Fall17_17Nov2017_V32_MC_L2Residual_AK4PFPuppi',
                'Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFPuppi'
            ],

            '2018':[
                'Autumn18_V19_MC_L1FastJet_AK4PFPuppi',
                'Autumn18_V19_MC_L2L3Residual_AK4PFPuppi',
                'Autumn18_V19_MC_L2Relative_AK4PFPuppi', #currently broken
                'Autumn18_V19_MC_L2Residual_AK4PFPuppi',  
                'Autumn18_V19_MC_L3Absolute_AK4PFPuppi'  
            ]
        }

        self._junc = {
    
            '2016':[
                'Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFPuppi'
            ],

            '2017':[
                'Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFPuppi'
            ],

            '2018':[
                'Autumn18_V19_MC_Uncertainty_AK4PFPuppi'
            ]
        }

        self._jr = {
        
            '2016': [
                'Summer16_25nsV1b_MC_PtResolution_AK4PFPuppi'
            ],
        
            '2017':[
                'Fall17_V3b_MC_PtResolution_AK4PFPuppi'
            ],

            '2018':[
                'Autumn18_V7b_MC_PtResolution_AK4PFPuppi'
            ]
        }

        self._jersf = {
    
            '2016':[
                'Summer16_25nsV1b_MC_SF_AK4PFPuppi'
            ],

            '2017':[
                'Fall17_V3b_MC_SF_AK4PFPuppi'
            ],

            '2018':[
                'Autumn18_V7b_MC_SF_AK4PFPuppi'
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
            'cutflow': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                #hist.Bin('cut', 'Cut index', 11, 0, 11),
                hist.Bin('cut', 'Cut index', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
            ),
            'template': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Cat('systematic', 'Systematic'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('recoil','Hadronic Recoil',[250,310,370,470,590,840,1020,1250,3000]),
                hist.Bin('fjmass','AK15 Jet Mass', [0,40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300]),#[0, 30, 60, 80, 120, 300]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'recoil': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('recoil','Hadronic Recoil',[250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0, 3000])
            ),
            'mindphirecoil': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('mindphirecoil','Min dPhi(Recoil,AK4s)',30,0,3.5)
            ),
            'minDphirecoil': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('minDphirecoil','Min dPhi(Recoil,AK15s)',30,0,3.5)
            ),
            'fjmass': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('fjmass','AK15 Jet Mass',30,0,300)
            ),
            'CaloMinusPfOverRecoil': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('CaloMinusPfOverRecoil','Calo - Pf / Recoil',35,0,1)
            ),
            'met': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('met','MET',30,0,600)
            ),
            'metphi': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('metphi','MET phi',35,-3.5,3.5)
            ),
            'mindphimet': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('mindphimet','Min dPhi(MET,AK4s)',30,0,3.5)
            ),
            'minDphimet': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('minDphimet','Min dPhi(MET,AK15s)',30,0,3.5)
            ),
            'j1pt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('j1pt','AK4 Leading Jet Pt',[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])
            ),
            'j1eta': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('j1eta','AK4 Leading Jet Eta',35,-3.5,3.5)
            ),
            'j1phi': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('j1phi','AK4 Leading Jet Phi',35,-3.5,3.5)
            ),
            'fj1pt': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('fj1pt','AK15 Leading SoftDrop Jet Pt',[160.0, 200.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])
            ),
            'fj1eta': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('fj1eta','AK15 Leading SoftDrop Jet Eta',35,-3.5,3.5)
            ),
            'fj1phi': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('fj1phi','AK15 Leading SoftDrop Jet Phi',35,-3.5,3.5)
            ),
            'njets': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('njets','AK4 Number of Jets',6,-0.5,5.5)
            ),
            'ndflvL': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('ndflvL','AK4 Number of deepFlavor Loose Jets',6,-0.5,5.5)
            ),
            'nfjclean': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('nfjclean','AK15 Number of Cleaned Jets',4,-0.5,3.5)
            ),
            'mT': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('mT','Transverse Mass',20,0,600)
            ),
            'l1pt': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('l1pt','Leading Lepton/Photon Pt',[0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])
            ),
            'l1eta': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('l1eta','Leading Lepton/Photon Eta',48,-2.4,2.4)
            ),
            'l1phi': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('l1phi','Leading Lepton/Photon Phi',64,-3.2,3.2)
            ),
            'dilepmass': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('dilepmass','Dilepton Mass',100,0,500)
            ),
            'dileppt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('dileppt','Dilepton Pt',150,0,800)
            ),
            'ZHbbvsQCD': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD',15,0,1)
            ),
            'TvsQCD': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
                hist.Bin('TvsQCD','TvsQCD',15,0,1)
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

        selected_regions = []
        for region, samples in self._samples.items():
            for sample in samples:
                if sample not in dataset: continue
                selected_regions.append(region)

        isData = 'genWeight' not in events.columns
        selection = processor.PackedSelection()
        hout = self.accumulator.identity()

        ###
        #Getting corrections, ids from .coffea files
        ###

        get_msd_weight          = self._corrections['get_msd_weight']
        get_ttbar_weight        = self._corrections['get_ttbar_weight']
        get_nnlo_nlo_weight     = self._corrections['get_nnlo_nlo_weight'][self._year]
        get_nlo_qcd_weight      = self._corrections['get_nlo_qcd_weight'][self._year]
        get_nlo_ewk_weight      = self._corrections['get_nlo_ewk_weight'][self._year]
        get_pu_weight           = self._corrections['get_pu_weight'][self._year]          
        get_met_trig_weight     = self._corrections['get_met_trig_weight'][self._year]    
        get_met_zmm_trig_weight = self._corrections['get_met_zmm_trig_weight'][self._year]
        get_ele_trig_weight     = self._corrections['get_ele_trig_weight'][self._year]    
        get_pho_trig_weight     = self._corrections['get_pho_trig_weight'][self._year]    
        get_ele_loose_id_sf     = self._corrections['get_ele_loose_id_sf'][self._year]
        get_ele_tight_id_sf     = self._corrections['get_ele_tight_id_sf'][self._year]
        get_pho_tight_id_sf     = self._corrections['get_pho_tight_id_sf'][self._year]
        get_pho_csev_sf         = self._corrections['get_pho_csev_sf'][self._year]
        get_mu_tight_id_sf      = self._corrections['get_mu_tight_id_sf'][self._year]
        get_mu_loose_id_sf      = self._corrections['get_mu_loose_id_sf'][self._year]
        get_ele_reco_sf         = self._corrections['get_ele_reco_sf'][self._year]
        get_ele_reco_lowet_sf   = self._corrections['get_ele_reco_lowet_sf']
        get_mu_tight_iso_sf     = self._corrections['get_mu_tight_iso_sf'][self._year]
        get_mu_loose_iso_sf     = self._corrections['get_mu_loose_iso_sf'][self._year]
        get_ecal_bad_calib      = self._corrections['get_ecal_bad_calib']
        get_deepflav_weight     = self._corrections['get_btag_weight']['deepflav'][self._year]
        Jetevaluator            = self._corrections['Jetevaluator']
        
        isLooseElectron = self._ids['isLooseElectron'] 
        isTightElectron = self._ids['isTightElectron'] 
        isLooseMuon     = self._ids['isLooseMuon']     
        isTightMuon     = self._ids['isTightMuon']     
        isLooseTau      = self._ids['isLooseTau']      
        isLoosePhoton   = self._ids['isLoosePhoton']   
        isTightPhoton   = self._ids['isTightPhoton']   
        isGoodJet       = self._ids['isGoodJet']       
        isGoodFatJet    = self._ids['isGoodFatJet']    
        isHEMJet        = self._ids['isHEMJet']        
        
        match = self._common['match']
        sigmoid = self._common['sigmoid'] #to calculate photon trigger efficiency
        deepflavWPs = self._common['btagWPs']['deepflav'][self._year]
        deepcsvWPs = self._common['btagWPs']['deepcsv'][self._year]

        ###
        # Derive jet corrector for JEC/JER
        ###
        
        JECcorrector = FactorizedJetCorrector(**{name: Jetevaluator[name] for name in self._jec[self._year]})
        JECuncertainties = JetCorrectionUncertainty(**{name:Jetevaluator[name] for name in self._junc[self._year]})
        JER = JetResolution(**{name:Jetevaluator[name] for name in self._jr[self._year]})
        JERsf = JetResolutionScaleFactor(**{name:Jetevaluator[name] for name in self._jersf[self._year]})
        Jet_transformer = JetTransformer(jec=JECcorrector,junc=JECuncertainties, jer = JER, jersf = JERsf)
        
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
        mu['isloose'] = isLooseMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.looseId,self._year)
        mu['istight'] = isTightMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.tightId,self._year)
        mu['T'] = TVector2Array.from_polar(mu.pt, mu.phi)
        mu_loose=mu[mu.isloose.astype(np.bool)]
        mu_tight=mu[mu.istight.astype(np.bool)]
        mu_ntot = mu.counts
        mu_nloose = mu_loose.counts
        mu_ntight = mu_tight.counts
        leading_mu = mu[mu.pt.argmax()]
        leading_mu = leading_mu[leading_mu.istight.astype(np.bool)]

        e = events.Electron
        e['isclean'] = ~match(e,mu_loose,0.3) 
        e['isloose'] = isLooseElectron(e.pt,e.eta+e.deltaEtaSC,e.dxy,e.dz,e.cutBased,self._year)
        e['istight'] = isTightElectron(e.pt,e.eta+e.deltaEtaSC,e.dxy,e.dz,e.cutBased,self._year)
        e['T'] = TVector2Array.from_polar(e.pt, e.phi)
        e_clean = e[e.isclean.astype(np.bool)]
        e_loose = e_clean[e_clean.isloose.astype(np.bool)]
        e_tight = e_clean[e_clean.istight.astype(np.bool)]
        e_ntot = e.counts
        e_nloose = e_loose.counts
        e_ntight = e_tight.counts
        leading_e = e[e.pt.argmax()]
        leading_e = leading_e[leading_e.isclean.astype(np.bool)]
        leading_e = leading_e[leading_e.istight.astype(np.bool)]

        tau = events.Tau
        tau['isclean']=~match(tau,mu_loose,0.4)&~match(tau,e_loose,0.4)
        tau['isloose']=isLooseTau(tau.pt,tau.eta,tau.idDecayMode,tau.idMVAoldDM2017v2,self._year)
        tau_clean=tau[tau.isclean.astype(np.bool)]
        tau_loose=tau_clean[tau_clean.isloose.astype(np.bool)]
        tau_ntot=tau.counts
        tau_nloose=tau_loose.counts

        pho = events.Photon
        pho['isclean']=~match(pho,mu_loose,0.5)&~match(pho,e_loose,0.5)
        _id = 'cutBasedBitmap'
        if self._year=='2016': 
            _id = 'cutBased'
        pho['isloose']=isLoosePhoton(pho.pt,pho.eta,pho[_id],self._year)&(pho.electronVeto) #added electron veto flag
        pho['istight']=isTightPhoton(pho.pt,pho[_id],self._year)&(pho.isScEtaEB)&(pho.electronVeto) #tight photons are barrel only
        pho['T'] = TVector2Array.from_polar(pho.pt, pho.phi)
        pho_clean=pho[pho.isclean.astype(np.bool)]
        pho_loose=pho_clean[pho_clean.isloose.astype(np.bool)]
        pho_tight=pho_clean[pho_clean.istight.astype(np.bool)]
        pho_ntot=pho.counts
        pho_nloose=pho_loose.counts
        pho_ntight=pho_tight.counts
        leading_pho = pho[pho.pt.argmax()]
        leading_pho = leading_pho[leading_pho.isclean.astype(np.bool)]
        leading_pho = leading_pho[leading_pho.istight.astype(np.bool)]

        fj = events.AK15Puppi
        fj['sd'] = fj.subjets.sum()
        fj['isclean'] =~match(fj.sd,pho_loose,1.5)&~match(fj.sd,mu_loose,1.5)&~match(fj.sd,e_loose,1.5)
        fj['isgood'] = isGoodFatJet(fj.sd.pt, fj.sd.eta, fj.jetId)
        fj['T'] = TVector2Array.from_polar(fj.pt, fj.phi)
        fj['msd_raw'] = (fj.subjets * (1 - fj.subjets.rawFactor)).sum().mass
        fj['msd_corr'] = fj.msd_raw * awkward.JaggedArray.fromoffsets(fj.array.offsets, np.maximum(1e-5, get_msd_weight(fj.sd.pt.flatten(),fj.sd.eta.flatten())))
        fj['rho'] = 2 * np.log(fj.msd_corr / fj.sd.pt)
        probQCD=fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
        probZHbb=fj.probZbb+fj.probHbb
        fj['ZHbbvsQCD'] = probZHbb/(probZHbb+probQCD)
        probT=fj.probTbcq+fj.probTbqq
        fj['TvsQCD'] = probT/(probT+probQCD)
        probV=fj.probWcq+fj.probWqq+fj.probZbb+fj.probZcc+fj.probZqq
        probX=probZHbb+probV
        fj['XvsQCD'] = probX/(probX+probQCD)
        fj_good = fj[fj.isgood.astype(np.bool)]
        fj_clean = fj_good[fj_good.isclean.astype(np.bool)]
        fj_ntot = fj.counts
        fj_ngood = fj_good.counts
        fj_nclean = fj_clean.counts

        j = events.Jet
        j['isgood'] = isGoodJet(j.pt, j.eta, j.jetId, j.puId, j.neHEF, j.chHEF)
        j['isHEM'] = isHEMJet(j.pt, j.eta, j.phi)
        j['isclean'] = ~match(j,e_loose,0.4)&~match(j,mu_loose,0.4)&~match(j,pho_loose,0.4)
        j['isiso'] = ~match(j,fj_clean[fj_clean.pt.argmax()],1.5)
        j['isdcsvL'] = (j.btagDeepB>deepcsvWPs['loose'])
        j['isdflvL'] = (j.btagDeepFlavB>deepflavWPs['loose'])
        j['T'] = TVector2Array.from_polar(j.pt, j.phi)
        j['p4'] = TLorentzVectorArray.from_ptetaphim(j.pt, j.eta, j.phi, j.mass)
        j['ptRaw'] =j.pt * (1-j.rawFactor)
        j['massRaw'] = j.mass * (1-j.rawFactor)
        j['rho'] = j.pt.ones_like()*events.fixedGridRhoFastjetAll.array
        j_good = j[j.isgood.astype(np.bool)]
        j_clean = j_good[j_good.isclean.astype(np.bool)]
        j_iso = j_clean[j_clean.isiso.astype(np.bool)]
        j_dcsvL = j_iso[j_iso.isdcsvL.astype(np.bool)]
        j_dflvL = j_iso[j_iso.isdflvL.astype(np.bool)]
        j_HEM = j[j.isHEM.astype(np.bool)]
        j_ntot=j.counts
        j_ngood=j_good.counts
        j_nclean=j_clean.counts
        j_niso=j_iso.counts
        j_ndcsvL=j_dcsvL.counts
        j_ndflvL=j_dflvL.counts
        j_nHEM = j_HEM.counts
        leading_j = j[j.pt.argmax()]
        leading_j = leading_j[leading_j.isgood.astype(np.bool)]
        leading_j = leading_j[leading_j.isclean.astype(np.bool)]

        ###
        #Calculating derivatives
        ###

        ele_pairs = e_loose.distincts()
        diele = ele_pairs.i0+ele_pairs.i1
        diele['T'] = TVector2Array.from_polar(diele.pt, diele.phi)
        leading_ele_pair = ele_pairs[diele.pt.argmax()]
        leading_diele = diele[diele.pt.argmax()]

        mu_pairs = mu_loose.distincts()
        dimu = mu_pairs.i0+mu_pairs.i1
        dimu['T'] = TVector2Array.from_polar(dimu.pt, dimu.phi)
        leading_mu_pair = mu_pairs[dimu.pt.argmax()]
        leading_dimu = dimu[dimu.pt.argmax()]

        ###
        # Calculate recoil and transverse mass
        ###

        u = {
            'sr'    : met.T,
            'wecr'  : met.T+leading_e.T.sum(),
            'tecr'  : met.T+leading_e.T.sum(),
            'wmcr'  : met.T+leading_mu.T.sum(),
            'tmcr'  : met.T+leading_mu.T.sum(),
            'zecr'  : met.T+leading_diele.T.sum(),
            'zmcr'  : met.T+leading_dimu.T.sum(),
            'gcr'   : met.T+leading_pho.T.sum()
        }

        mT = {
            'wecr'  : np.sqrt(2*leading_e.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_e.T.sum())))),
            'tecr'  : np.sqrt(2*leading_e.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_e.T.sum())))),
            'wmcr'  : np.sqrt(2*leading_mu.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_mu.T.sum())))),
            'tmcr'  : np.sqrt(2*leading_mu.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_mu.T.sum())))) 
        }

        ###
        #Calculating weights
        ###
        if not isData:
            
            ###
            # JEC/JER
            ###

            #j['ptGenJet'] = j.matched_gen.pt
            #Jet_transformer.transform(j)

            gen = events.GenPart

            ###
            # Fat-jet top matching at decay level
            ###

            qFromW = gen[
                (abs(gen.pdgId) < 4) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (abs(gen.distinctParent.pdgId) == 24)
            ]
            cFromW = gen[
                (abs(gen.pdgId) == 4) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (abs(gen.distinctParent.pdgId) == 24)
            ]

            def tbqqmatch(topid, dR=1.5):
                qFromWFromTop = qFromW[qFromW.distinctParent.distinctParent.pdgId == topid]
                bFromTop = gen[
                    (abs(gen.pdgId) == 5) &
                    gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                    (gen.distinctParent.pdgId == topid)
                ]
                jetgenWq = fj.sd.cross(qFromWFromTop, nested=True)
                jetgenb = fj.sd.cross(bFromTop, nested=True)
                qWmatch = ((jetgenWq.i0.delta_r(jetgenWq.i1) < dR).sum()==2) & (qFromWFromTop.counts>0)
                bmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < dR).sum()==1) & (bFromTop.counts>0)
                return qWmatch & bmatch
            fj['isTbqq'] = tbqqmatch(6)|tbqqmatch(-6)

            def tbcqmatch(topid, dR=1.5):
                qFromWFromTop = qFromW[qFromW.distinctParent.distinctParent.pdgId == topid]
                cFromWFromTop = cFromW[cFromW.distinctParent.distinctParent.pdgId == topid]
                bFromTop = gen[
                    (abs(gen.pdgId) == 5) &
                    gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                    (gen.distinctParent.pdgId == topid)
                ]
                jetgenWq = fj.sd.cross(qFromWFromTop, nested=True)
                jetgenWc = fj.sd.cross(cFromWFromTop, nested=True)
                jetgenb = fj.sd.cross(bFromTop, nested=True)
                qWmatch = ((jetgenWq.i0.delta_r(jetgenWq.i1) < dR).sum()==1) & (qFromWFromTop.counts>0)
                cWmatch = ((jetgenWc.i0.delta_r(jetgenWc.i1) < dR).sum()==1) & (cFromWFromTop.counts>0)
                bmatch =  ((jetgenb.i0.delta_r(jetgenb.i1) < dR).sum()==1)   & (bFromTop.counts>0)
                return cWmatch & qWmatch & bmatch
            fj['isTbcq'] = tbcqmatch(6)|tbcqmatch(-6)

            def tqqmatch(topid, dR=1.5):
                qFromWFromTop = qFromW[qFromW.distinctParent.distinctParent.pdgId == topid]
                jetgenWq = fj.sd.cross(qFromWFromTop, nested=True)
                qWmatch = ((jetgenWq.i0.delta_r(jetgenWq.i1) < dR).sum()==2) & (qFromWFromTop.counts>0)
                return qWmatch
            fj['isTqq'] = tqqmatch(6)|tqqmatch(-6)

            def tcqmatch(topid, dR=1.5):
                qFromWFromTop = qFromW[qFromW.distinctParent.distinctParent.pdgId == topid]
                cFromWFromTop = cFromW[cFromW.distinctParent.distinctParent.pdgId == topid]
                jetgenWq = fj.sd.cross(qFromWFromTop, nested=True)
                jetgenWc = fj.sd.cross(cFromWFromTop, nested=True)
                qWmatch = ((jetgenWq.i0.delta_r(jetgenWq.i1) < dR).sum()==1) & (qFromWFromTop.counts>0)
                cWmatch = ((jetgenWc.i0.delta_r(jetgenWc.i1) < dR).sum()==1) & (cFromWFromTop.counts>0)
                return cWmatch & qWmatch
            fj['isTcq'] = tcqmatch(6)|tcqmatch(-6)

            def tbqmatch(topid, dR=1.5):
                qFromWFromTop = qFromW[qFromW.distinctParent.distinctParent.pdgId == topid]
                bFromTop = gen[
                    (abs(gen.pdgId) == 5) &
                    gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                    (gen.distinctParent.pdgId == topid)
                ]
                jetgenWq = fj.sd.cross(qFromWFromTop, nested=True)
                jetgenb = fj.sd.cross(bFromTop, nested=True)
                qWmatch = ((jetgenWq.i0.delta_r(jetgenWq.i1) < dR).sum()==1) & (qFromWFromTop.counts>0)
                bmatch =  ((jetgenb.i0.delta_r(jetgenb.i1) < dR).sum()==1)   & (bFromTop.counts>0)
                return qWmatch & bmatch
            fj['isTbq'] = tbqmatch(6)|tbqmatch(-6)

            def tbcmatch(topid, dR=1.5):
                cFromWFromTop = cFromW[cFromW.distinctParent.distinctParent.pdgId == topid]
                bFromTop = gen[
                    (abs(gen.pdgId) == 5) &
                    gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                    (gen.distinctParent.pdgId == topid)
                ]
                jetgenWc = fj.sd.cross(cFromWFromTop, nested=True)
                jetgenb = fj.sd.cross(bFromTop, nested=True)
                cWmatch = ((jetgenWc.i0.delta_r(jetgenWc.i1) < dR).sum()==1) & (cFromWFromTop.counts>0)
                bmatch =  ((jetgenb.i0.delta_r(jetgenb.i1) < dR).sum()==1)   & (bFromTop.counts>0)
                return cWmatch & bmatch
            fj['isTbc'] = tbcmatch(6)|tbcmatch(-6)

            ###
            # Fat-jet W->qq matching at decay level
            ###
            def wqqmatch(wid, dR=1.5):
                qFromSameW = qFromW[qFromW.distinctParent.pdgId == wid]
                jetgenq = fj.sd.cross(qFromSameW, nested=True)
                qqmatch = ((jetgenq.i0.delta_r(jetgenq.i1) < dR).sum()==2) & (qFromSameW.counts>0)
                return qqmatch
            fj['isWqq']  = wqqmatch(24)|wqqmatch(-24)

            ###
            # Fat-jet W->cq matching at decay level
            ###
            def wcqmatch(wid, dR=1.5):
                qFromSameW = qFromW[qFromW.distinctParent.pdgId == wid]
                cFromSameW = cFromW[cFromW.distinctParent.pdgId == wid]
                jetgenq = fj.sd.cross(qFromSameW, nested=True)
                jetgenc = fj.sd.cross(cFromSameW, nested=True)
                qmatch = ((jetgenq.i0.delta_r(jetgenq.i1) < dR).sum()==1) & (qFromSameW.counts>0)
                cmatch = ((jetgenc.i0.delta_r(jetgenc.i1) < dR).sum()==1) & (cFromSameW.counts>0)
                return qmatch & cmatch
            fj['isWcq']  = wcqmatch(24)|wcqmatch(-24)

            ###
            # Fat-jet Z->bb matching at decay level
            ###
            bFromZ = gen[
                (abs(gen.pdgId) == 5) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (abs(gen.distinctParent.pdgId) == 23)
            ]
            def zbbmatch(zid, dR=1.5):
                bFromSameZ = bFromZ[bFromZ.distinctParent.pdgId == zid]
                jetgenb = fj.sd.cross(bFromSameZ, nested=True)
                bbmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < dR).sum()==2) & (bFromSameZ.counts>0)
                return bbmatch
            fj['isZbb']  = zbbmatch(23)|zbbmatch(-23)

            ###
            # Fat-jet Z->cc matching at decay level
            ###
            cFromZ = gen[
                (abs(gen.pdgId) == 4) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (abs(gen.distinctParent.pdgId) == 23)
            ]
            def zccmatch(zid, dR=1.5):
                cFromSameZ = cFromZ[cFromZ.distinctParent.pdgId == zid]
                jetgenc = fj.sd.cross(cFromSameZ, nested=True)
                ccmatch = ((jetgenc.i0.delta_r(jetgenc.i1) < dR).sum()==2) & (cFromSameZ.counts>0)
                return ccmatch
            fj['isZcc']  = zccmatch(23)|zccmatch(-23)

            ###
            # Fat-jet Z->qq matching at decay level
            ###
            qFromZ = gen[
                (abs(gen.pdgId) < 4) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (abs(gen.distinctParent.pdgId) == 23)
            ]
            def zqqmatch(zid, dR=1.5):
                qFromSameZ = qFromZ[qFromZ.distinctParent.pdgId == zid]
                jetgenq = fj.sd.cross(qFromSameZ, nested=True)
                qqmatch = ((jetgenq.i0.delta_r(jetgenq.i1) < dR).sum()==2) & (qFromSameZ.counts>0)
                return qqmatch
            fj['isZqq']  = zqqmatch(23)|zqqmatch(-23)

            ###
            # Fat-jet H->bb matching at decay level
            ###
            bFromH = gen[
                (abs(gen.pdgId) == 5) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (abs(gen.distinctParent.pdgId) == 25)
            ]
            def hbbmatch(hid, dR=1.5):
                bFromSameH = bFromH[bFromH.distinctParent.pdgId == hid]
                jetgenb = fj.sd.cross(bFromSameH, nested=True)
                bbmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < dR).sum()==2) & (bFromSameH.counts>0)
                return bbmatch
            fj['isHbb']  = hbbmatch(25)|hbbmatch(-25)

            ###
            # Fat-jet dark H->bb matching at decay level
            ###
            bFromHs = gen[
                (abs(gen.pdgId) == 5) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (abs(gen.distinctParent.pdgId) == 54)
            ]
            def hsbbmatch(hid, dR=1.5):
                bFromSameHs = bFromHs[bFromHs.distinctParent.pdgId == hid]
                jetgenb = fj.sd.cross(bFromSameHs, nested=True)
                bbmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < dR).sum()==2) & (bFromSameHs.counts>0)
                return bbmatch
            fj['isHsbb']  = hsbbmatch(54)|hsbbmatch(-54)

            gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            jetgenb = fj.sd.cross(gen[gen.isb], nested=True)
            bmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < 1.5).sum()==1)&(gen[gen.isb].counts>0)
            fj['isb']  = bmatch
        
            bmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < 1.5).sum()==2)&(gen[gen.isb].counts>0)
            fj['isbb']  = bmatch

            gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            jetgenc = fj.sd.cross(gen[gen.isc], nested=True)
            cmatch = ((jetgenc.i0.delta_r(jetgenc.i1) < 1.5).sum()==1)&(gen[gen.isc].counts>0)
            fj['isc']  = cmatch

            cmatch = ((jetgenc.i0.delta_r(jetgenc.i1) < 1.5).sum()==2)&(gen[gen.isc].counts>0)
            fj['iscc']  = cmatch

            gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            genTops = gen[gen.isTop]
            nlo = np.ones(events.size)
            if('TTJets' in dataset): 
                nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt.sum()) * get_ttbar_weight(genTops[:,1].pt.sum()))
                
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

            #genWs = gen[gen.isW&(gen.pt>=100)]
            genWs = gen[gen.isW] 
            genZs = gen[gen.isZ]
            genDYs = gen[gen.isZ&(gen.mass>30)]
            genIsoAs = gen[gen.isIsoA] 

            nnlo_nlo = {}
            nlo_qcd = np.ones(events.size)
            nlo_ewk = np.ones(events.size)
            if('GJets' in dataset): 
                if self._year=='2016':
                    nlo_qcd = get_nlo_qcd_weight['a'](genIsoAs.pt.max())
                    nlo_ewk = get_nlo_ewk_weight['a'](genIsoAs.pt.max())
                for systematic in get_nnlo_nlo_weight['a']:
                    nnlo_nlo[systematic] = get_nnlo_nlo_weight['a'][systematic](genIsoAs.pt.max())*((genIsoAs.counts>0)&(genIsoAs.pt.max()>=290)) + \
                                           get_nnlo_nlo_weight['a'][systematic](290)*((genIsoAs.counts>0)&~(genIsoAs.pt.max()>=290)&(genIsoAs.pt.max()>=100)) + \
                                           (~((genIsoAs.counts>0)&(genIsoAs.pt.max()>=100))).astype(np.int)
            elif('WJets' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['w'](genWs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['w'](genWs.pt.max())
                for systematic in get_nnlo_nlo_weight['w']:
                    nnlo_nlo[systematic] = get_nnlo_nlo_weight['w'][systematic](genWs.pt.max())*((genWs.counts>0)&(genWs.pt.max()>=100)) + \
                                           (~((genWs.counts>0)&(genWs.pt.max()>=100))).astype(np.int)
            elif('DY' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['dy'](genDYs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['dy'](genDYs.pt.max())
                for systematic in get_nnlo_nlo_weight['dy']:
                    nnlo_nlo[systematic] = get_nnlo_nlo_weight['dy'][systematic](genDYs.pt.max())*((genDYs.counts>0)&(genDYs.pt.max()>=100)) + \
                                           (~((genDYs.counts>0)&(genDYs.pt.max()>=100))).astype(np.int)
            elif('ZJets' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['z'](genZs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['z'](genZs.pt.max())
                for systematic in get_nnlo_nlo_weight['z']:
                    nnlo_nlo[systematic] = get_nnlo_nlo_weight['z'][systematic](genZs.pt.max())*((genZs.counts>0)&(genZs.pt.max()>=100)) + \
                                           (~((genZs.counts>0)&(genZs.pt.max()>=100))).astype(np.int)

            ###
            # Calculate PU weight and systematic variations
            ###

            pu = get_pu_weight(events.Pileup.nTrueInt)

            ###
            # Trigger efficiency weight
            ###

            e1sf = get_ele_trig_weight(leading_ele_pair.i0.eta.sum()+leading_ele_pair.i0.deltaEtaSC.sum(),leading_ele_pair.i0.pt.sum())*(leading_ele_pair.i0.pt.sum()>40).astype(np.int)
            e2sf = get_ele_trig_weight(leading_ele_pair.i1.eta.sum()+leading_ele_pair.i1.deltaEtaSC.sum(),leading_ele_pair.i1.pt.sum())*(leading_ele_pair.i1.pt.sum()>40).astype(np.int)

            if self._year == '2016':
                sf =  get_pho_trig_weight(leading_pho.pt.sum())
            elif self._year == '2017': #Sigmoid used for 2017 and 2018, values from monojet
                sf = sigmoid(leading_pho.pt.sum(),0.335,217.91,0.065,0.996) / sigmoid(leading_pho.pt.sum(),0.244,212.34,0.050,1.000)
                sf[np.isnan(sf) | np.isinf(sf)] == 1
            elif self._year == '2018':
                sf = sigmoid(leading_pho.pt.sum(),1.022, 218.39, 0.086, 0.999) / sigmoid(leading_pho.pt.sum(), 0.301,212.83,0.062,1.000)
                sf[np.isnan(sf) | np.isinf(sf)] == 1

            trig = {
                'sr':   get_met_trig_weight(met.pt),
                'wmcr': get_met_trig_weight(u['wmcr'].mag),
                'tmcr': get_met_trig_weight(u['tmcr'].mag),
                'zmcr': get_met_zmm_trig_weight(u['zmcr'].mag),
                'wecr': get_ele_trig_weight(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(), leading_e.pt.sum()),
                'tecr': get_ele_trig_weight(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(), leading_e.pt.sum()),
                'zecr': 1 - (1 - e1sf)*(1 - e2sf),
                'gcr':  sf
            }

            ### 
            # Calculating electron and muon ID weights
            ###

            mueta = abs(leading_mu.eta.sum())
            mu1eta=abs(leading_mu_pair.i0.eta.sum())
            mu2eta=abs(leading_mu_pair.i1.eta.sum())
            if self._year=='2016':
                mueta=leading_mu.eta.sum()
                mu1eta=leading_mu_pair.i0.eta.sum()
                mu2eta=leading_mu_pair.i1.eta.sum()
            if self._year=='2016':
                sf = get_pho_tight_id_sf(leading_pho.eta.sum(),leading_pho.pt.sum())
            else: #2017/2018 monojet measurement depends only on abs(eta)
                sf = get_pho_tight_id_sf(abs(leading_pho.eta.sum()))

            ids ={
                'sr':  np.ones(events.size),
                'wmcr': get_mu_tight_id_sf(mueta,leading_mu.pt.sum()),
                'tmcr': get_mu_tight_id_sf(mueta,leading_mu.pt.sum()),
                'zmcr': get_mu_loose_id_sf(mu1eta,leading_mu_pair.i0.pt.sum()) * get_mu_loose_id_sf(mu2eta,leading_mu_pair.i1.pt.sum()),
                'wecr': get_ele_tight_id_sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'tecr': get_ele_tight_id_sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'zecr': get_ele_loose_id_sf(leading_ele_pair.i0.eta.sum()+leading_ele_pair.i0.deltaEtaSC.sum(),leading_ele_pair.i0.pt.sum()) * get_ele_loose_id_sf(leading_ele_pair.i1.eta.sum()+leading_ele_pair.i1.deltaEtaSC.sum(),leading_ele_pair.i1.pt.sum()),
                'gcr':  sf
            }

            ###
            # Reconstruction weights for electrons
            ###

            def ele_reco_sf(pt, eta):#2017 has separate weights for low/high pT (threshold at 20 GeV)
                return get_ele_reco_sf(eta, pt)*(pt>20).astype(np.int) + get_ele_reco_lowet_sf(eta, pt)*(~(pt>20)).astype(np.int)

            if self._year == '2017':
                sf = ele_reco_sf
            else:
                sf = get_ele_reco_sf

            reco = {
                'sr': np.ones(events.size),
                'wmcr': np.ones(events.size),
                'tmcr': np.ones(events.size),
                'zmcr': np.ones(events.size),
                'wecr': sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'tecr': sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'zecr': sf(leading_ele_pair.i0.eta.sum()+leading_ele_pair.i0.deltaEtaSC.sum(),leading_ele_pair.i0.pt.sum()) * sf(leading_ele_pair.i1.eta.sum()+leading_ele_pair.i1.deltaEtaSC.sum(),leading_ele_pair.i1.pt.sum()),
                'gcr': np.ones(events.size)
            }

            ###
            # Isolation weights for muons
            ###

            isolation = {
                'sr'  : np.ones(events.size),
                'wmcr': get_mu_tight_iso_sf(mueta,leading_mu.pt.sum()),
                'tmcr': get_mu_tight_iso_sf(mueta,leading_mu.pt.sum()),
                'zmcr': get_mu_loose_iso_sf(mu1eta,leading_mu_pair.i0.pt.sum()) * get_mu_loose_iso_sf(mu2eta,leading_mu_pair.i1.pt.sum()),
                'wecr': np.ones(events.size),
                'tecr': np.ones(events.size),
                'zecr': np.ones(events.size),
                'gcr':  np.ones(events.size)
            }

            ###
            # CSEV weight for photons: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_Veto_CSEV_or_pixel_seed
            ###

            if self._year == '2016':
                csev_weight = get_pho_csev_sf(abs(leading_pho.eta.sum()), leading_pho.pt.sum())
            elif self._year == '2017':
                csev_sf_index = 0.5*(leading_pho.isScEtaEB.sum()).astype(np.int)+3.5*(~(leading_pho.isScEtaEB.sum())).astype(np.int)+1*(leading_pho.r9.sum()>0.94).astype(np.int)+2*(leading_pho.r9.sum()<=0.94).astype(np.int)
                csev_weight = get_pho_csev_sf(csev_sf_index)
            elif self._year == '2018':
                csev_weight = get_pho_csev_sf(leading_pho.pt.sum(), abs(leading_pho.eta.sum()))
            csev_weight[csev_weight==0] = 1

            csev = {
                'sr'  : np.ones(events.size),
                'wmcr': np.ones(events.size),
                'tmcr': np.ones(events.size),
                'zmcr': np.ones(events.size),
                'wecr': np.ones(events.size),
                'tecr': np.ones(events.size),
                'zecr': np.ones(events.size),
                'gcr':  csev_weight
            }

            ###
            # AK4 b-tagging weights
            ###

            btag = {}
            btagUp = {}
            btagDown = {}
            btag['sr'],   btagUp['sr'],   btagDown['sr']   = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'0')
            btag['wmcr'], btagUp['wmcr'], btagDown['wmcr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'0')
            btag['tmcr'], btagUp['tmcr'], btagDown['tmcr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'-1')
            btag['wecr'], btagUp['wecr'], btagDown['wecr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'0')
            btag['tecr'], btagUp['tecr'], btagDown['tecr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'-1')
            btag['zmcr'], btagUp['zmcr'], btagDown['zmcr'] = np.ones(events.size), np.ones(events.size), np.ones(events.size)
            btag['zecr'], btagUp['zecr'], btagDown['zecr'] = np.ones(events.size), np.ones(events.size), np.ones(events.size)
            btag['gcr'],  btagUp['gcr'],  btagDown['gcr']  = np.ones(events.size), np.ones(events.size), np.ones(events.size)

        ###
        # Selections
        ###

        met_filters =  np.ones(events.size, dtype=np.bool)
        if isData: met_filters = met_filters & events.Flag['eeBadScFilter']#this filter is recommended for data only
        for flag in AnalysisProcessor.met_filter_flags[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters',met_filters)

        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._met_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('met_triggers', triggers)

        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._singleelectron_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('singleelectron_triggers', triggers)

        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._singlephoton_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('singlephoton_triggers', triggers)

        noHEMj = np.ones(events.size, dtype=np.bool)
        if self._year=='2018': noHEMj = (j_nHEM==0)

        noHEMmet = np.ones(events.size, dtype=np.bool)
        if self._year=='2018': noHEMmet = (met.pt>470)|(met.phi>-0.62)|(met.phi<-1.62)

        leading_fj = fj[fj.sd.pt.argmax()]
        leading_fj = leading_fj[leading_fj.isgood.astype(np.bool)]
        leading_fj = leading_fj[leading_fj.isclean.astype(np.bool)]
        selection.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('isoneE', (e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('istwoM', (e_nloose==0)&(mu_nloose==2)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('istwoE',(e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('isoneA', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_ntight==1)&(pho_nloose==1))
        selection.add('leading_e_pt',(e_loose.pt.max()>40))
        selection.add('dimu_mass',(leading_dimu.mass.sum()>60)&(leading_dimu.mass.sum()<120))
        selection.add('diele_mass',(leading_diele.mass.sum()>60)&(leading_diele.mass.sum()<120))
        selection.add('noextrab', (j_ndflvL==0))
        selection.add('extrab', (j_ndflvL>0))
        selection.add('fatjet', (fj_nclean>0))
        selection.add('noHEMj', noHEMj)
        selection.add('noHEMmet', noHEMmet)
        selection.add('met120',(met.pt<120))
        selection.add('met100',(met.pt>100))
        #selection.add('mindphimet',(abs(met.T.delta_phi(j_clean.T)).min())>0.7)

        regions = {
            #'sr': ['iszeroL','fatjet','noextrab','noHEMmet','met_filters','met_triggers','noHEMj'],
            'sr': ['fatjet', 'noHEMj', 'iszeroL', 'noextrab','met_filters','met_triggers','noHEMmet'],
            'wmcr': ['isoneM','fatjet','noextrab','noHEMj','met_filters','met_triggers'],
            'tmcr': ['isoneM','fatjet','extrab','noHEMj','met_filters','met_triggers'],
            'wecr': ['isoneE','fatjet','noextrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'tecr': ['isoneE','fatjet','extrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'zmcr': ['istwoM','fatjet','noHEMj','met_filters','met_triggers', 'dimu_mass','met120'],
            'zecr': ['istwoE','fatjet','noHEMj','met_filters','singleelectron_triggers', 'diele_mass','met120','leading_e_pt'],
            'gcr': ['isoneA','fatjet','noHEMj','met_filters','singlephoton_triggers']
        }

        isFilled = False

        #for region in selected_regions: 
        for region, cuts in regions.items():
            if region not in selected_regions: continue
            print('Considering region:', region)

            ###
            # Adding recoil and minDPhi requirements
            ###

            selection.add('recoil_'+region, (u[region].mag>250))
            selection.add('mindphi_'+region, (abs(u[region].delta_phi(j_clean.T)).min()>0.5))
            selection.add('minDphi_'+region, (abs(u[region].delta_phi(fj_clean.T)).min()>1.5))
            selection.add('calo_'+region, ( (abs(calomet.pt - met.pt) / u[region].mag) < 0.5))
            #regions[region].update({'recoil_'+region,'mindphi_'+region})
            regions[region].insert(0, 'recoil_'+region)
            regions[region].insert(3, 'mindphi_'+region)
            regions[region].insert(4, 'minDphi_'+region)
            regions[region].insert(5, 'calo_'+region)
            print('Selection:',regions[region])
            variables = {
                'recoil':                 u[region].mag,
                'mindphirecoil':          abs(u[region].delta_phi(j_clean.T)).min(),
                'minDphirecoil':          abs(u[region].delta_phi(fj_clean.T)).min(),
                'CaloMinusPfOverRecoil':  abs(calomet.pt - met.pt) / u[region].mag,
                'met':                    met.pt,
                'metphi':                 met.phi,
                'mindphimet':             abs(met.T.delta_phi(j_clean.T)).min(),
                'minDphimet':             abs(met.T.delta_phi(fj_clean.T)).min(),
                'j1pt':                   leading_j.pt,
                'j1eta':                  leading_j.eta,
                'j1phi':                  leading_j.phi,
                'fj1pt':                  leading_fj.sd.pt,
                'fj1eta':                 leading_fj.sd.eta,
                'fj1phi':                 leading_fj.sd.phi,
                'fjmass':                 leading_fj.msd_corr,
                'njets':                  j_nclean,
                'ndflvL':                 j_ndflvL,
                'nfjclean':               fj_nclean,
                'ZHbbvsQCD':              leading_fj.ZHbbvsQCD,
                'TvsQCD':                 leading_fj.TvsQCD,
            }
            if region in mT:
                variables['mT']           = mT[region]
            if 'e' in region:
                variables['l1pt']      = leading_e.pt
                variables['l1phi']     = leading_e.phi
                variables['l1eta']     = leading_e.eta
                if 'z' in region:
                    variables['dilepmass']  = leading_diele.mass
                    variables['dileppt']    = leading_diele.pt
            if 'm' in region:
                variables['l1pt']      = leading_mu.pt
                variables['l1phi']     = leading_mu.phi
                variables['l1eta']     = leading_mu.eta
                if 'z' in region:
                    variables['dilepmass']  = leading_dimu.mass
                    variables['dileppt']    = leading_dimu.pt
            if 'g' in region:
                variables['l1pt']      = leading_pho.pt
                variables['l1phi']     = leading_pho.phi
                variables['l1eta']     = leading_pho.eta
            print('Variables:',variables.keys())

            def fill(dataset, gentype, weight, cut):

                flat_variables = {k: v[cut].flatten() for k, v in variables.items()}
                flat_gentype = {k: (~np.isnan(v[cut])*gentype[cut]).flatten() for k, v in variables.items()}
                flat_weight = {k: (~np.isnan(v[cut])*weight[cut]).flatten() for k, v in variables.items()}
            
                for histname, h in hout.items():
                    if not isinstance(h, hist.Hist):
                        continue
                    if histname not in variables:
                        continue
                    elif histname == 'sumw':
                        continue
                    elif histname == 'template':
                        continue
                    else:
                        flat_variable = {histname: flat_variables[histname]}
                        h.fill(dataset=dataset, 
                               region=region, 
                               gentype=flat_gentype[histname], 
                               **flat_variable, 
                               weight=flat_weight[histname])

            if isData:
                if not isFilled:
                    hout['sumw'].fill(dataset=dataset, sumw=1, weight=1)
                    isFilled=True
                cut = selection.all(*regions[region])
                hout['template'].fill(dataset=dataset,
                                      region=region,
                                      systematic='nominal',
                                      gentype=np.zeros(events.size, dtype=np.int),
                                      recoil=u[region].mag,
                                      fjmass=leading_fj.msd_corr.sum(),
                                      ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                      weight=np.ones(events.size)*cut)
                fill(dataset, np.zeros(events.size, dtype=np.int), np.ones(events.size), cut)
            else:
                weights = processor.Weights(len(events))
                if 'L1PreFiringWeight' in events.columns: weights.add('prefiring',events.L1PreFiringWeight.Nom)
                weights.add('genw',events.genWeight)
                weights.add('nlo_qcd',nlo_qcd)
                weights.add('nlo_ewk',nlo_ewk)
                if 'cen' in nnlo_nlo:
                    #weights.add('nnlo_nlo',nnlo_nlo['cen'])
                    weights.add('qcd1',np.ones(events.size), nnlo_nlo['qcd1up']/nnlo_nlo['cen'], nnlo_nlo['qcd1do']/nnlo_nlo['cen'])
                    weights.add('qcd2',np.ones(events.size), nnlo_nlo['qcd2up']/nnlo_nlo['cen'], nnlo_nlo['qcd2do']/nnlo_nlo['cen'])
                    weights.add('qcd3',np.ones(events.size), nnlo_nlo['qcd3up']/nnlo_nlo['cen'], nnlo_nlo['qcd3do']/nnlo_nlo['cen'])
                    weights.add('ew1',np.ones(events.size), nnlo_nlo['ew1up']/nnlo_nlo['cen'], nnlo_nlo['ew1do']/nnlo_nlo['cen'])
                    weights.add('ew2G',np.ones(events.size), nnlo_nlo['ew2Gup']/nnlo_nlo['cen'], nnlo_nlo['ew2Gdo']/nnlo_nlo['cen'])
                    weights.add('ew3G',np.ones(events.size), nnlo_nlo['ew3Gup']/nnlo_nlo['cen'], nnlo_nlo['ew3Gdo']/nnlo_nlo['cen'])
                    weights.add('ew2W',np.ones(events.size), nnlo_nlo['ew2Wup']/nnlo_nlo['cen'], nnlo_nlo['ew2Wdo']/nnlo_nlo['cen'])
                    weights.add('ew3W',np.ones(events.size), nnlo_nlo['ew3Wup']/nnlo_nlo['cen'], nnlo_nlo['ew3Wdo']/nnlo_nlo['cen'])
                    weights.add('ew2Z',np.ones(events.size), nnlo_nlo['ew2Zup']/nnlo_nlo['cen'], nnlo_nlo['ew2Zdo']/nnlo_nlo['cen'])
                    weights.add('ew3Z',np.ones(events.size), nnlo_nlo['ew3Zup']/nnlo_nlo['cen'], nnlo_nlo['ew3Zdo']/nnlo_nlo['cen'])
                    weights.add('mix',np.ones(events.size), nnlo_nlo['mixup']/nnlo_nlo['cen'], nnlo_nlo['mixdo']/nnlo_nlo['cen'])
                    weights.add('muF',np.ones(events.size), nnlo_nlo['muFup']/nnlo_nlo['cen'], nnlo_nlo['muFdo']/nnlo_nlo['cen'])
                    weights.add('muR',np.ones(events.size), nnlo_nlo['muRup']/nnlo_nlo['cen'], nnlo_nlo['muRdo']/nnlo_nlo['cen'])
                weights.add('pileup',pu)
                weights.add('trig', trig[region])
                weights.add('ids', ids[region])
                weights.add('reco', reco[region])
                weights.add('isolation', isolation[region])
                weights.add('csev', csev[region])
                weights.add('btag',btag[region], btagUp[region], btagDown[region])

                wgentype = {
                    'xbb' : (
                        (leading_fj.isHsbb | leading_fj.isHbb | leading_fj.isZbb)
                    ).sum(),
                    'tbcq' : (
                        ~(leading_fj.isHsbb |leading_fj.isHbb | leading_fj.isZbb)&
                        leading_fj.isTbcq 
                    ).sum(),
                    'tbqq' : (
                        ~(leading_fj.isHsbb |leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        leading_fj.isTbqq 
                    ).sum(),
                    'zcc' : (
                        ~(leading_fj.isHsbb |leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        leading_fj.isZcc
                    ).sum(),
                    'wcq' : (
                        ~(leading_fj.isHsbb |leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        (leading_fj.isTcq | leading_fj.isWcq)
                    ).sum(),
                    'vqq' : (
                        ~(leading_fj.isHsbb |leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        ~(leading_fj.isTcq | leading_fj.isWcq) &
                        (leading_fj.isWqq | leading_fj.isZqq | leading_fj.isTqq)
                    ).sum(),
                    'bb' : (
                        ~(leading_fj.isHsbb | leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        ~(leading_fj.isTcq | leading_fj.isWcq) &
                        ~(leading_fj.isWqq | leading_fj.isZqq | leading_fj.isTqq) &
                        leading_fj.isbb
                    ).sum(),
                    'bc' : (
                        ~(leading_fj.isHsbb | leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        ~(leading_fj.isTcq | leading_fj.isWcq) &
                        ~(leading_fj.isWqq | leading_fj.isZqq | leading_fj.isTqq) &
                        ~leading_fj.isbb &
                        (leading_fj.isTbc | (leading_fj.isb & leading_fj.isc))
                    ).sum(),
                    'b' : (
                        ~(leading_fj.isHsbb | leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        ~(leading_fj.isTcq | leading_fj.isWcq) &
                        ~(leading_fj.isWqq | leading_fj.isZqq | leading_fj.isTqq) &
                        ~leading_fj.isbb &
                        ~(leading_fj.isTbc | (leading_fj.isb & leading_fj.isc)) &
                        (leading_fj.isTbq | leading_fj.isb)
                    ).sum(),
                    'cc' : (
                        ~(leading_fj.isHsbb | leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        ~(leading_fj.isTcq | leading_fj.isWcq) &
                        ~(leading_fj.isWqq | leading_fj.isZqq | leading_fj.isTqq) &
                        ~leading_fj.isbb &
                        ~(leading_fj.isTbc | (leading_fj.isb & leading_fj.isc)) &
                        ~(leading_fj.isTbq | leading_fj.isb) &
                        leading_fj.iscc
                    ).sum(),
                    'c' : (
                        ~(leading_fj.isHsbb | leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        ~(leading_fj.isTcq | leading_fj.isWcq) &
                        ~(leading_fj.isWqq | leading_fj.isZqq | leading_fj.isTqq) &
                        ~leading_fj.isbb &
                        ~(leading_fj.isTbc | (leading_fj.isb & leading_fj.isc)) &
                        ~(leading_fj.isTbq | leading_fj.isb) &
                        ~leading_fj.iscc &
                        leading_fj.isc
                    ).sum(),
                    'other' : (
                        ~(leading_fj.isHsbb | leading_fj.isHbb | leading_fj.isZbb)&
                        ~leading_fj.isTbcq &
                        ~leading_fj.isTbqq &
                        ~leading_fj.isZcc &
                        ~(leading_fj.isTcq | leading_fj.isWcq) &
                        ~(leading_fj.isWqq | leading_fj.isZqq | leading_fj.isTqq) &
                        ~leading_fj.isbb &
                        ~(leading_fj.isTbc | (leading_fj.isb & leading_fj.isc)) &
                        ~(leading_fj.isTbq | leading_fj.isb) &
                        ~leading_fj.iscc &
                        ~leading_fj.isc
                    ).sum(),
                }
                vgentype=np.zeros(events.size, dtype=np.int)
                for gentype in self._gentype_map.keys():
                    vgentype += self._gentype_map[gentype]*wgentype[gentype]

                if 'WJets' in dataset or 'ZJets' in dataset or 'DY' in dataset or 'GJets' in dataset:
                    if not isFilled:
                        hout['sumw'].fill(dataset='HF--'+dataset, sumw=1, weight=events.genWeight.sum())
                        hout['sumw'].fill(dataset='LF--'+dataset, sumw=1, weight=events.genWeight.sum())
                        isFilled=True
                    whf = ((gen[gen.isb].counts>0)|(gen[gen.isc].counts>0)).astype(np.int)
                    wlf = (~(whf.astype(np.bool))).astype(np.int)
                    cut = selection.all(*regions[region])
                    systematics = [None,
                                   'btagUp',
                                   'btagDown',
                                   'qcd1Up',
                                   'qcd1Down',
                                   'qcd2Up',
                                   'qcd2Down',
                                   'qcd3Up',
                                   'qcd3Down',
                                   'muFUp',
                                   'muFDown',
                                   'muRUp',
                                   'muRDown',
                                   'ew1Up',
                                   'ew1Down',
                                   'ew2GUp',
                                   'ew2GDown',
                                   'ew2WUp',
                                   'ew2WDown',
                                   'ew2ZUp',
                                   'ew2ZDown',
                                   'ew3GUp',
                                   'ew3GDown',
                                   'ew3WUp',
                                   'ew3WDown',
                                   'ew3ZUp',
                                   'ew3ZDown',
                                   'mixUp',
                                   'mixDown']
                    for systematic in systematics:
                        sname = 'nominal' if systematic is None else systematic
                        hout['template'].fill(dataset='HF--'+dataset,
                                              region=region,
                                              systematic=sname,
                                              gentype=vgentype,
                                              recoil=u[region].mag,
                                              fjmass=leading_fj.msd_corr.sum(),
                                              ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                              weight=weights.weight(modifier=systematic)*whf*cut)
                        hout['template'].fill(dataset='LF--'+dataset,
                                              region=region,
                                              systematic=sname,
                                              gentype=vgentype,
                                              recoil=u[region].mag,
                                              fjmass=leading_fj.msd_corr.sum(),
                                              ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                              weight=weights.weight(modifier=systematic)*wlf*cut)

                    ## Cutflow loop
                    vcut=np.zeros(events.size, dtype=np.int)
                    hout['cutflow'].fill(dataset='HF--'+dataset, region=region, cut=vcut, weight=weights.weight()*whf)
                    hout['cutflow'].fill(dataset='LF--'+dataset, region=region, cut=vcut, weight=weights.weight()*wlf)
                    allcuts = set()
                    for i, icut in enumerate(cuts):
                        allcuts.add(icut)
                        jcut = selection.all(*allcuts)
                        vcut = (i+1)*jcut
                        hout['cutflow'].fill(dataset='HF--'+dataset, region=region, cut=vcut, weight=weights.weight()*jcut*whf)
                        hout['cutflow'].fill(dataset='LF--'+dataset, region=region, cut=vcut, weight=weights.weight()*jcut*wlf)

                    fill('HF--'+dataset, vgentype, weights.weight()*whf, cut)
                    fill('LF--'+dataset, vgentype, weights.weight()*wlf, cut)
                else:
                    if not isFilled:
                        hout['sumw'].fill(dataset=dataset, sumw=1, weight=events.genWeight.sum())
                        isFilled=True
                    cut = selection.all(*regions[region])
                    for systematic in [None, 'btagUp', 'btagDown']:
                        sname = 'nominal' if systematic is None else systematic
                        hout['template'].fill(dataset=dataset,
                                              region=region,
                                              systematic=sname,
                                              gentype=vgentype,
                                              recoil=u[region].mag,
                                              fjmass=leading_fj.msd_corr.sum(),
                                              ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                              #TvsQCD=leading_fj.TvsQCD.sum(),
                                              #XvsQCD=leading_fj.XvsQCD.sum(),
                                              weight=weights.weight(modifier=systematic)*cut)
                    ## Cutflow loop
                    vcut=np.zeros(events.size, dtype=np.int)
                    hout['cutflow'].fill(dataset=dataset, region=region, cut=vcut, weight=weights.weight())
                    allcuts = set()
                    for i, icut in enumerate(cuts):
                        allcuts.add(icut)
                        jcut = selection.all(*allcuts)
                        vcut = (i+1)*jcut
                        hout['cutflow'].fill(dataset=dataset, region=region, cut=vcut, weight=weights.weight()*jcut)

                    fill(dataset, vgentype, weights.weight(), cut)

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

    save(processor_instance, 'data/darkhiggs'+options.year+'.processor')
