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

        self._samples = {
            'sr':('ZJets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET','mhs'),
            'wmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'tmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'wecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'tecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'qcdcr':('ZJets','WJets','TT','ST','WW','WZ','ZZ','QCD','MET'),
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
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'template': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Cat('systematic', 'Systematic'),
                hist.Bin('recoil','Hadronic Recoil',[250,310,370,470,590,840,1020,1250,3000]),
                hist.Bin('fjmass','AK15 Jet Mass', [40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300]),#[0, 30, 60, 80, 120, 300]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'ZHbbvsQCD': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD',15,0,1)
            ),
            'mindphirecoil': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('mindphirecoil','Min dPhi(Recoil,AK4s)',30,0,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'minDphirecoil': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('minDphirecoil','Min dPhi(Recoil,AK15s)',30,0,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'CaloMinusPfOverRecoil': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('CaloMinusPfOverRecoil','Calo - Pf / Recoil',35,0,1),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'met': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('met','MET',30,0,600),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'metphi': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('metphi','MET phi',35,-3.5,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'mindphimet': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('mindphimet','Min dPhi(MET,AK4s)',30,0,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'minDphimet': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('minDphimet','Min dPhi(MET,AK15s)',30,0,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'j1pt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('j1pt','AK4 Leading Jet Pt',[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'j1eta': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('j1eta','AK4 Leading Jet Eta',35,-3.5,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'j1phi': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('j1phi','AK4 Leading Jet Phi',35,-3.5,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'fj1pt': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('fj1pt','AK15 Leading SoftDrop Jet Pt',[160.0, 200.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'fj1eta': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('fj1eta','AK15 Leading SoftDrop Jet Eta',35,-3.5,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'fj1phi': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('fj1phi','AK15 Leading SoftDrop Jet Phi',35,-3.5,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'njets': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('njets','AK4 Number of Jets',6,-0.5,5.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'ndflvL': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('ndflvL','AK4 Number of deepFlavor Loose Jets',6,-0.5,5.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'nfjclean': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('nfjclean','AK15 Number of Cleaned Jets',4,-0.5,3.5),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'mT': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('mT','Transverse Mass',20,0,600),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'l1pt': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('l1pt','Leading Lepton/Photon Pt',[0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'l1eta': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('l1eta','Leading Lepton/Photon Eta',48,-2.4,2.4),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'l1phi': hist.Hist(
                'Events', 
                hist.Cat('dataset', 'Dataset'), 
                hist.Cat('region', 'Region'), 
                hist.Bin('l1phi','Leading Lepton/Photon Phi',64,-3.2,3.2),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
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
        get_ele_loose_id_sf     = self._corrections['get_ele_loose_id_sf'][self._year]
        get_ele_tight_id_sf     = self._corrections['get_ele_tight_id_sf'][self._year]
        get_mu_tight_id_sf      = self._corrections['get_mu_tight_id_sf'][self._year]
        get_mu_loose_id_sf      = self._corrections['get_mu_loose_id_sf'][self._year]
        get_ele_reco_sf         = self._corrections['get_ele_reco_sf'][self._year]
        get_ele_reco_lowet_sf   = self._corrections['get_ele_reco_lowet_sf']
        get_mu_tight_iso_sf     = self._corrections['get_mu_tight_iso_sf'][self._year]
        get_mu_loose_iso_sf     = self._corrections['get_mu_loose_iso_sf'][self._year]
        get_ecal_bad_calib      = self._corrections['get_ecal_bad_calib']
        get_deepflav_weight     = self._corrections['get_btag_weight']['deepflav'][self._year]
        get_doublebtag_weight   = self._corrections['get_doublebtag_weight'][self._year]
        get_mu_rochester_sf     = self._corrections['get_mu_rochester_sf'][self._year]
        
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
        deepflavWPs = self._common['btagWPs']['deepflav'][self._year]
        deepcsvWPs = self._common['btagWPs']['deepcsv'][self._year]

        ###
        #Initialize global quantities (MET ecc.)
        ###

        met = events.MET
        if self._year == '2017': met = events.METFixEE2017#Recommended for 2017
        met['T']  = TVector2Array.from_polar(met.pt, met.phi)
        calomet = events.CaloMET
        puppimet = events.PuppiMET

        ###
        #Initialize physics objects
        ###

        mu = events.Muon
        rochester = get_mu_rochester_sf
        _muon_offsets = mu.pt.offsets
        _charge = mu.charge
        _pt = mu.pt
        _eta = mu.eta
        _phi = mu.phi
        if isData:
            _k = rochester.kScaleDT(_charge, _pt, _eta, _phi)
        else:
            # for default if gen present
            _gpt = mu.matched_gen.pt
            # for backup w/o gen
            _nl = mu.nTrackerLayers
            _u = awkward.JaggedArray.fromoffsets(_muon_offsets, np.random.rand(*_pt.flatten().shape))
            _hasgen = (_gpt.fillna(-1) > 0)
            _kspread = rochester.kSpreadMC(_charge[_hasgen], _pt[_hasgen], _eta[_hasgen], _phi[_hasgen],
                                           _gpt[_hasgen])
            _ksmear = rochester.kSmearMC(_charge[~_hasgen], _pt[~_hasgen], _eta[~_hasgen], _phi[~_hasgen],
                                         _nl[~_hasgen], _u[~_hasgen])
            _k = _pt.ones_like()
            _k[_hasgen] = _kspread
            _k[~_hasgen] = _ksmear
        mask = _pt < 200
        rochester_pt = _pt.ones_like()
        rochester_pt[~mask] = _pt[~mask]
        rochester_pt[mask] = (_k * _pt)[mask]
        mu['pt'] = rochester_pt
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
        tau['isloose']=isLooseTau(tau.pt,tau.eta,tau.idDecayMode,tau.idDecayModeNewDMs,tau.idDeepTau2017v2p1VSe,tau.idDeepTau2017v2p1VSjet,tau.idDeepTau2017v2p1VSmu,self._year)
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
        # Calculate recoil and transverse mass
        ###

        u = {
            'sr'    : met.T,
            'wecr'  : met.T+leading_e.T.sum(),
            'tecr'  : met.T+leading_e.T.sum(),
            'wmcr'  : met.T+leading_mu.T.sum(),
            'tmcr'  : met.T+leading_mu.T.sum(),
            'qcdcr' : met.T,
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
            
            gen = events.GenPart

            gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            genTops = gen[gen.isTop]
            nlo = np.ones(events.size)
            if('TTJets' in dataset): 
                nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt.sum()) * get_ttbar_weight(genTops[:,1].pt.sum()))
                
            gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            
            genWs = gen[gen.isW] 
            genZs = gen[gen.isZ]
            genDYs = gen[gen.isZ&(gen.mass>30)]
            
            nnlo_nlo = {}
            nlo_qcd = np.ones(events.size)
            nlo_ewk = np.ones(events.size)
            if('WJets' in dataset): 
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

            trig = {
                'sr':   get_met_trig_weight(met.pt),
                'wmcr': get_met_trig_weight(u['wmcr'].mag),
                'tmcr': get_met_trig_weight(u['tmcr'].mag),
                'wecr': get_ele_trig_weight(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(), leading_e.pt.sum()),
                'tecr': get_ele_trig_weight(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(), leading_e.pt.sum()),
                'qcdcr':   get_met_trig_weight(met.pt),
            }

            ### 
            # Calculating electron and muon ID weights
            ###

            mueta = abs(leading_mu.eta.sum())
            if self._year=='2016':
                mueta=leading_mu.eta.sum()
            ids ={
                'sr':  np.ones(events.size),
                'wmcr': get_mu_tight_id_sf(mueta,leading_mu.pt.sum()),
                'tmcr': get_mu_tight_id_sf(mueta,leading_mu.pt.sum()),
                'wecr': get_ele_tight_id_sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'tecr': get_ele_tight_id_sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'qcdcr':  np.ones(events.size),
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
                'wecr': sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'tecr': sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'qcdcr': np.ones(events.size),
            }

            ###
            # Isolation weights for muons
            ###

            isolation = {
                'sr'  : np.ones(events.size),
                'wmcr': get_mu_tight_iso_sf(mueta,leading_mu.pt.sum()),
                'tmcr': get_mu_tight_iso_sf(mueta,leading_mu.pt.sum()),
                'wecr': np.ones(events.size),
                'tecr': np.ones(events.size),
                'qcdcr'  : np.ones(events.size),
            }

            ###
            # AK4 b-tagging weights
            ###

            btagSF = {}
            btagSFbc_correlatedUp = {}
            btagSFbc_correlatedDown = {}
            btagSFbc_uncorrelatedUp = {}
            btagSFbc_uncorrelatedDown = {}
            btagSFlight_correlatedUp = {}
            btagSFlight_correlatedDown = {}
            btagSFlight_uncorrelatedUp = {}
            btagSFlight_uncorrelatedDown = {}
            btagSF['sr'], btagSFbc_correlatedUp['sr'], btagSFbc_correlatedDown['sr'], btagSFbc_uncorrelatedUp['sr'], btagSFbc_uncorrelatedDown['sr'], btagSFlight_correlatedUp['sr'], btagSFlight_correlatedDown['sr'], btagSFlight_uncorrelatedUp['sr'], btagSFlight_uncorrelatedDown['sr']   = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'0')
            btagSF['wmcr'], btagSFbc_correlatedUp['wmcr'], btagSFbc_correlatedDown['wmcr'], btagSFbc_uncorrelatedUp['wmcr'], btagSFbc_uncorrelatedDown['wmcr'], btagSFlight_correlatedUp['wmcr'], btagSFlight_correlatedDown['wmcr'], btagSFlight_uncorrelatedUp['wmcr'], btagSFlight_uncorrelatedDown['wmcr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'0')
            btagSF['wecr'], btagSFbc_correlatedUp['wecr'], btagSFbc_correlatedDown['wecr'], btagSFbc_uncorrelatedUp['wecr'], btagSFbc_uncorrelatedDown['wecr'], btagSFlight_correlatedUp['wecr'], btagSFlight_correlatedDown['wecr'], btagSFlight_uncorrelatedUp['wecr'], btagSFlight_uncorrelatedDown['wecr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'0')
            btagSF['tmcr'], btagSFbc_correlatedUp['tmcr'], btagSFbc_correlatedDown['tmcr'], btagSFbc_uncorrelatedUp['tmcr'], btagSFbc_uncorrelatedDown['tmcr'], btagSFlight_correlatedUp['tmcr'], btagSFlight_correlatedDown['tmcr'], btagSFlight_uncorrelatedUp['tmcr'], btagSFlight_uncorrelatedDown['tmcr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'-1')
            btagSF['tecr'], btagSFbc_correlatedUp['tecr'], btagSFbc_correlatedDown['tecr'], btagSFbc_uncorrelatedUp['tecr'], btagSFbc_uncorrelatedDown['tecr'], btagSFlight_correlatedUp['tecr'], btagSFlight_correlatedDown['tecr'], btagSFlight_uncorrelatedUp['tecr'], btagSFlight_uncorrelatedDown['tecr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'-1')
            btagSF['qcdcr'], btagSFbc_correlatedUp['qcdcr'], btagSFbc_correlatedDown['qcdcr'], btagSFbc_uncorrelatedUp['qcdcr'], btagSFbc_uncorrelatedDown['qcdcr'], btagSFlight_correlatedUp['qcdcr'], btagSFlight_correlatedDown['qcdcr'], btagSFlight_uncorrelatedUp['qcdcr'], btagSFlight_uncorrelatedDown['qcdcr'] = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,'0')

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
        selection.add('leading_e_pt',(e_loose.pt.max()>40))
        selection.add('noextrab', (j_ndflvL==0))
        selection.add('extrab', (j_ndflvL>0))
        selection.add('fatjet', (fj_nclean>0))
        selection.add('noHEMj', noHEMj)
        selection.add('noHEMmet', noHEMmet)
        selection.add('met120',(met.pt<120))
        selection.add('met100',(met.pt>100))
        selection.add('msd40',(leading_fj.msd_corr.sum()>40))
        selection.add('recoil_qcdcr', (u['qcdcr'].mag>250))
        selection.add('mindphi_qcdcr', (abs(u['qcdcr'].delta_phi(j_clean.T)).min()<0.1))
        selection.add('minDphi_qcdcr', (abs(u['qcdcr'].delta_phi(fj_clean.T)).min()>1.5))
        selection.add('calo_qcdcr', ( (abs(calomet.pt - met.pt) / u['qcdcr'].mag)<0.5))
            
        #selection.add('mindphimet',(abs(met.T.delta_phi(j_clean.T)).min())>0.7)

        regions = {
            #'sr': ['iszeroL','fatjet','noextrab','noHEMmet','met_filters','met_triggers','noHEMj'],
            'sr': ['msd40','fatjet', 'noHEMj','iszeroL','noextrab','met_filters','met_triggers','noHEMmet'],
            'wmcr': ['msd40','isoneM','fatjet','noextrab','noHEMj','met_filters','met_triggers'],
            'tmcr': ['msd40','isoneM','fatjet','extrab','noHEMj','met_filters','met_triggers'],
            'wecr': ['msd40','isoneE','fatjet','noextrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'tecr': ['msd40','isoneE','fatjet','extrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'qcdcr': ['recoil_qcdcr','mindphi_qcdcr','minDphi_qcdcr','calo_qcdcr','msd40','fatjet', 'noHEMj','iszeroL','noextrab','met_filters','met_triggers','noHEMmet'],
        }

        isFilled = False

        #for region in selected_regions: 
        for region, cuts in regions.items():
            if region not in selected_regions: continue

            ###
            # Adding recoil and minDPhi requirements
            ###

            selection.add('recoil_'+region, (u[region].mag>250))
            selection.add('mindphi_'+region, (abs(u[region].delta_phi(j_clean.T)).min()>0.5))
            selection.add('minDphi_'+region, (abs(u[region].delta_phi(fj_clean.T)).min()>1.5))
            selection.add('calo_'+region, ( (abs(calomet.pt - met.pt) / u[region].mag) < 0.5))
            #regions[region].update({'recoil_'+region,'mindphi_'+region})
            if 'qcd' not in region:
                regions[region].insert(0, 'recoil_'+region)
                regions[region].insert(3, 'mindphi_'+region)
                regions[region].insert(4, 'minDphi_'+region)
                regions[region].insert(5, 'calo_'+region)
            variables = {
                'mindphirecoil':          abs(u[region].delta_phi(j_clean.T)).min(),
                'minDphirecoil':          abs(u[region].delta_phi(fj_clean.T)).min(),
                'CaloMinusPfOverRecoil':  abs(calomet.pt - met.pt) / u[region].mag,
                'met':                    met.pt.flatten(),
                'metphi':                 met.phi.flatten(),
                'mindphimet':             abs(met.T.delta_phi(j_clean.T)).min(),
                'minDphimet':             abs(met.T.delta_phi(fj_clean.T)).min(),
                'j1pt':                   leading_j.pt.sum(),
                'j1eta':                  leading_j.eta.sum(),
                'j1phi':                  leading_j.phi.sum(),
                'fj1pt':                  leading_fj.sd.pt.sum(),
                'fj1eta':                 leading_fj.sd.eta.sum(),
                'fj1phi':                 leading_fj.sd.phi.sum(),
                'njets':                  j_nclean,
                'ndflvL':                 j_ndflvL,
                'nfjclean':               fj_nclean,
            }
            if region in mT:
                variables['mT']           = mT[region]
            if 'e' in region:
                variables['l1pt']      = leading_e.pt.sum()
                variables['l1phi']     = leading_e.phi.sum()
                variables['l1eta']     = leading_e.eta.sum()
            if 'm' in region:
                variables['l1pt']      = leading_mu.pt.sum()
                variables['l1phi']     = leading_mu.phi.sum()
                variables['l1eta']     = leading_mu.eta.sum()

            def fill(dataset, weight, cut):
                for histname, h in hout.items():
                    if not isinstance(h, hist.Hist):
                        continue
                    if histname not in variables:
                        continue
                    flat_variable = {histname: variables[histname]}
                    h.fill(dataset=dataset, 
                           region=region, 
                           **flat_variable, 
                           ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                           weight=weight*cut)

            if isData:
                if not isFilled:
                    hout['sumw'].fill(dataset=dataset, sumw=1, weight=1)
                    isFilled=True
                cut = selection.all(*regions[region])
                hout['template'].fill(dataset=dataset,
                                      region=region,
                                      systematic='nominal',
                                      recoil=u[region].mag,
                                      fjmass=leading_fj.msd_corr.sum(),
                                      ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                      weight=np.ones(events.size)*cut)
                hout['ZHbbvsQCD'].fill(dataset=dataset,
                                      region=region,
                                      ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                      weight=np.ones(events.size)*cut)
                fill(dataset, np.ones(events.size), cut)
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
                weights.add('btagSF',btagSF[region])
                weights.add('btagSFbc_correlated',np.ones(events.size), btagSFbc_correlatedUp[region]/btagSF[region], btagSFbc_correlatedDown[region]/btagSF[region])
                weights.add('btagSFbc_uncorrelated',np.ones(events.size), btagSFbc_uncorrelatedUp[region]/btagSF[region], btagSFbc_uncorrelatedDown[region]/btagSF[region])
                weights.add('btagSFlight_correlated',np.ones(events.size), btagSFlight_correlatedUp[region]/btagSF[region], btagSFlight_correlatedDown[region]/btagSF[region])
                weights.add('btagSFlight_uncorrelated',np.ones(events.size), btagSFlight_uncorrelatedUp[region]/btagSF[region], btagSFlight_uncorrelatedDown[region]/btagSF[region])

                ###
                # AK15 doubleb-tagging weights
                ###
                
                if('mhs' in dataset):
                    for k in get_doublebtag_weight(leading_fj.sd.pt.sum())[0]:
                        doublebtag = get_doublebtag_weight(leading_fj.sd.pt.sum())[0][k]
                        doublebtagUp = get_doublebtag_weight(leading_fj.sd.pt.sum())[1][k]
                        doublebtagDown = get_doublebtag_weight(leading_fj.sd.pt.sum())[2][k]
                        weights.add('doublebtag'+k,doublebtag, doublebtagUp, doublebtagDown)

                if 'WJets' in dataset or 'ZJets' in dataset or 'DY' in dataset:
                    if not isFilled:
                        hout['sumw'].fill(dataset='HF--'+dataset, sumw=1, weight=events.genWeight.sum())
                        hout['sumw'].fill(dataset='LF--'+dataset, sumw=1, weight=events.genWeight.sum())
                        isFilled=True
                    whf = ((gen[gen.isb].counts>0)|(gen[gen.isc].counts>0)).astype(np.int)
                    wlf = (~(whf.astype(np.bool))).astype(np.int)
                    cut = selection.all(*regions[region])
                    systematics = [None,
                                   'btagSFbc_correlatedUp',
                                   'btagSFbc_correlatedDown',
                                   'btagSFbc_uncorrelatedUp',
                                   'btagSFbc_uncorrelatedDown',
                                   'btagSFlight_correlatedUp',
                                   'btagSFlight_correlatedDown',
                                   'btagSFlight_uncorrelatedUp',
                                   'btagSFlight_uncorrelatedDown',
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
                                              recoil=u[region].mag,
                                              fjmass=leading_fj.msd_corr.sum(),
                                              ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                              weight=weights.weight(modifier=systematic)*whf*cut)
                        hout['template'].fill(dataset='LF--'+dataset,
                                              region=region,
                                              systematic=sname,
                                              recoil=u[region].mag,
                                              fjmass=leading_fj.msd_corr.sum(),
                                              ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                              weight=weights.weight(modifier=systematic)*wlf*cut)
                    ## Cutflow loop
                    vcut=np.zeros(events.size, dtype=np.int)
                    hout['cutflow'].fill(dataset='HF--'+dataset, region=region, cut=vcut, ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(), weight=weights.weight()*whf)
                    hout['cutflow'].fill(dataset='LF--'+dataset, region=region, cut=vcut, ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(), weight=weights.weight()*wlf)
                    allcuts = set()
                    for i, icut in enumerate(cuts):
                        allcuts.add(icut)
                        jcut = selection.all(*allcuts)
                        vcut = (i+1)*jcut
                        hout['cutflow'].fill(dataset='HF--'+dataset, region=region, cut=vcut, ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(), weight=weights.weight()*jcut*whf)
                        hout['cutflow'].fill(dataset='LF--'+dataset, region=region, cut=vcut, ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(), weight=weights.weight()*jcut*wlf)

                    hout['ZHbbvsQCD'].fill(dataset='HF--'+dataset,
                                           region=region,
                                           ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                           weight=weights.weight()*whf*cut)
                    hout['ZHbbvsQCD'].fill(dataset='LF--'+dataset,
                                           region=region,
                                           ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                           weight=weights.weight()*wlf*cut)
                    fill('HF--'+dataset, weights.weight()*whf, cut)
                    fill('LF--'+dataset, weights.weight()*wlf, cut)
                else:
                    if not isFilled:
                        hout['sumw'].fill(dataset=dataset, sumw=1, weight=events.genWeight.sum())
                        isFilled=True
                    cut = selection.all(*regions[region])
                    systematics = [None, 
                                   'btagSFbc_correlatedUp', 
                                   'btagSFbc_correlatedDown', 
                                   'btagSFbc_uncorrelatedUp', 
                                   'btagSFbc_uncorrelatedDown',
                                   'btagSFlight_correlatedUp', 
                                   'btagSFlight_correlatedDown', 
                                   'btagSFlight_uncorrelatedUp', 
                                   'btagSFlight_uncorrelatedDown',
                               ]
                    if('mhs' in dataset):
                        for k in get_doublebtag_weight(leading_fj.sd.pt.sum())[0]:
                            systematics.append('doublebtag'+k+'Up')
                            systematics.append('doublebtag'+k+'Down')
                    for systematic in systematics:
                        sname = 'nominal' if systematic is None else systematic
                        hout['template'].fill(dataset=dataset,
                                              region=region,
                                              systematic=sname,
                                              recoil=u[region].mag,
                                              fjmass=leading_fj.msd_corr.sum(),
                                              ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                              weight=weights.weight(modifier=systematic)*cut)
                    ## Cutflow loop
                    vcut=np.zeros(events.size, dtype=np.int)
                    hout['cutflow'].fill(dataset=dataset, region=region, cut=vcut, ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(), weight=weights.weight())
                    allcuts = set()
                    for i, icut in enumerate(cuts):
                        allcuts.add(icut)
                        jcut = selection.all(*allcuts)
                        vcut = (i+1)*jcut
                        hout['cutflow'].fill(dataset=dataset, region=region, cut=vcut, ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(), weight=weights.weight()*jcut)

                    hout['ZHbbvsQCD'].fill(dataset=dataset,
                                           region=region,
                                           ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                           weight=weights.weight()*cut)
                    fill(dataset, weights.weight(), cut)

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

    save(processor_instance, 'data/darkhiggs'+options.name+'.processor')
