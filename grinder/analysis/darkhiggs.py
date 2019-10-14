#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import pprint
import numpy as np
import awkward
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea.arrays import Initialize
from coffea import hist, processor
from utils.triggers import met_trigger_paths, singleele_trigger_paths, singlepho_trigger_paths
from utils.corrections import get_ttbar_weight, get_nlo_weight, get_pu_weight
from utils.corrections import get_met_trig_weight, get_met_zmm_trig_weight, get_ele_trig_weight, get_pho_trig_weight
from utils.corrections import get_ecal_bad_calib
from utils.ids import e_id, isLooseElectron, isTightElectron
from utils.ids import mu_id, isLooseMuon, isTightMuon
from utils.ids import tau_id, isLooseTau
from utils.ids import pho_id, isLoosePhoton, isTightPhoton
from utils.ids import j_id, fj_id, isGoodJet, isGoodFatJet, isHEMJet
from utils.metfilters import met_filter_flags
from utils.deep import deep

samples = {
    "iszeroL":('ZJets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','MET','Mhs_50','Mhs_70','Mhs_90','MonoJet','MonoW','MonoZ','Monojet'),
    "isoneM":('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','MET'),
    "isoneE":('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','SingleElectron','EGamma'),
    "istwoM":('WJets','DY','TT','ST','WW','WZ','ZZ','HToBB','MET'),
    "istwoE":('WJets','DY','TT','ST','WW','WZ','ZZ','HToBB','SingleElectron','EGamma'),
    "isoneA":('GJets','QCD','SinglePhoton','EGamma')
}

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, selected_regions, year, xsec, lumi):
        self._selected_regions = selected_regions
        self._year = year
        self._xsec = xsec
        self._lumi = lumi
        
        self._accumulator = processor.dict_accumulator({
            'sumw': hist.Hist("sumw", hist.Cat("dataset", "Primary dataset"), hist.Bin("sumw", "Weight value", [0.])),
            'CaloMinusPfOverRecoil': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("CaloMinusPfOverRecoil","Calo - Pf / Recoil",35,0,1)),
            'recoil': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("recoil","Hadronic Recoil",[250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
            'mindphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mindphi","Min dPhi(MET,AK4s)",30,0,3.5)),
            'diledphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("diledphi","Min dPhi(Dileptons,AK4s)",30,0,3.5)),
            'ledphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("ledphi","Min dPhi(lepton,AK4s)",30,0,3.5)),
            'j1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("j1pt","AK4 Leading Jet Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
            'j1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("j1eta","AK4 Leading Jet Eta",35,-3.5,3.5)),
            'j1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("j1phi","AK4 Leading Jet Phi",35,-3.5,3.5)),
            'fj1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fj1pt","AK15 Leading Jet Pt",[200.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
            'fj1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fj1eta","AK15 Leading Jet Eta",35,-3.5,3.5)),
            'fj1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fj1phi","AK15 Leading Jet Phi",35,-3.5,3.5)),
            'njets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("njets","AK4 Number of Jets",6,0,5)),
            'ndcsvL': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("ndcsvL","AK4 Number of deepCSV Loose Jets",6,0,5)),
            'ndflvL': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("ndflvL","AK4 Number of deepFlavor Loose Jets",6,0,5)),
            'ndcsvM': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("ndcsvM","AK4 Number of deepCSV Medium Jets",6,0,5)),
            'ndflvM': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("ndflvM","AK4 Number of deepFlavor Medium Jets",6,0,5)),
            'ndcsvT': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("ndcsvT","AK4 Number of deepCSV Tight Jets",6,0,5)),
            'ndflvT': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("ndflvT","AK4 Number of deepFlavor Tight Jets",6,0,5)),
            'nfjtot': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("nfjtot","AK15 Number of Jets",4,0,3)),
            'nfjgood': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("nfjgood","AK15 Number of Good Jets",4,0,3)),
            'nfjclean': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("nfjclean","AK15 Number of cleaned Jets",4,0,3)),
            'fjmass': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fjmass","AK15 Jet Mass",[0.0, 40.0, 60.0, 75.0, 85.0, 115.0, 135.0, 225.0, 500.0])),
            'e1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("e1pt","Leading Electron Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
            'e1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("e1eta","Leading Electron Eta",48,-2.4,2.4)),
            'e1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("e1phi","Leading Electron Phi",64,-3.2,3.2)),
            'dielemass': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("dielemass","Dielectron mass",100,0,500)),
            'mu1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mu1pt","Leading Muon Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
            'mu1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mu1eta","Leading Muon Eta",48,-2.4,2.4)),
            'mu1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mu1phi","Leading Muon Phi",64,-3.2,3.2)),
            'dimumass': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("dimumass","Dimuon mass",100,0,500)),
            'TopTagger': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("TopTagger","TopTagger",15,0,1)),
            'DarkHiggsTagger': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("DarkHiggsTagger","DarkHiggsTagger",15,0,1)),
            'VvsQCDTagger': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("VvsQCDTagger","VvsQCDTagger",15,0,1)),
            'probTbcq': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probTbcq","probTbcq",15,0,1)),
            'probTbqq': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probTbqq","probTbqq",15,0,1)),
            'probTbc': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probTbc","probTbc",15,0,1)),
            'probTbq': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probTbq","probTbq",15,0,1)),
            'probWcq': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probWcq","probWcq",15,0,1)),
            'probWqq': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probWqq","probWqq",15,0,1)),
            'probZbb': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probZbb","probZbb",15,0,1)),
            'probZcc': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probZcc","probZcc",15,0,1)),
            'probZqq': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probZqq","probZqq",15,0,1)),
            'probHbb': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probHbb","probHbb",15,0,1)),
            'probHcc': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probHcc","probHcc",15,0,1)),
            'probHqqqq': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probHqqqq","probHqqqq",15,0,1)),
            'probQCDbb': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probQCDbb","probQCDbb",15,0,1)),
            'probQCDcc': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probQCDcc","probQCDcc",15,0,1)),
            'probQCDb': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probQCDb","probQCDb",15,0,1)),
            'probQCDc': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probQCDc","probQCDc",15,0,1)),
            'probQCDothers': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("probQCDothers","probQCDothers",15,0,1)),
            'recoilVSmindphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("recoil","Hadronic Recoil",[250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]), hist.Bin("mindphi","Min dPhi(MET,AK4s)",30,0,3.5)),

        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):

            dataset = df['dataset']

            ###
            #Initialize global quantities (MET ecc.)
            ###

            met = Initialize({'pt':df['MET_pt'],
                              'eta':0,
                              'phi':df['MET_phi'],
                              'mass':0})

            calomet = Initialize({'pt':df['CaloMET_pt'],
                                  'eta':0,
                                  'phi':df['CaloMET_phi'],
                                  'mass':0})

            ###
            #Initialize physics objects
            ###

            #Define first and empty object that will use as protection against arrays with size 0
            #Will use MET to set the correct size for the arrays
            #Not used at the moment

            #empty_jagged = awkward.JaggedArray.fromcounts(np.ones_like(met.pt, dtype=int),np.zeros_like(met.pt))
            #empty_obj = Initialize({'pt':empty_jagged,
            #                        'eta':empty_jagged,
            #                        'phi':empty_jagged,
            #                        'mass':empty_jagged})

            e = Initialize({'pt':df['Electron_pt'],
                            'eta':df['Electron_eta'],
                            'phi':df['Electron_phi'],
                            'mass':df['Electron_mass']})

            for key in e_id[self._year]:
                e[key] = e.pt.zeros_like()
                if e_id[self._year][key] in df:
                    e[key] = df[e_id[self._year][key]]

            e['isloose'] = isLooseElectron(e.pt,e.eta,e.dxy,e.dz,e.iso,e.loose_id,self._year)
            e['istight'] = isTightElectron(e.pt,e.eta,e.dxy,e.dz,e.iso,e.tight_id,self._year)

            leading_e = e[e.pt.argmax()]
            leading_e = leading_e[leading_e.istight]

            e_loose = e[e.isloose]
            e_tight = e[e.istight]

            e_ntot = e.counts
            e_nloose = e_loose.counts
            e_ntight = e_tight.counts

            mu = Initialize({'pt':df['Muon_pt'],
                             'eta':df['Muon_eta'],
                             'phi':df['Muon_phi'],
                             'mass':df['Muon_mass']})

            for key in mu_id[self._year]:
                mu[key] = mu.pt.zeros_like()
                if mu_id[self._year][key] in df:
                    mu[key] = df[mu_id[self._year][key]]

            mu['isloose'] = isLooseMuon(mu.pt,mu.eta,mu.dxy,mu.dz,mu.iso,mu.med_id,self._year)
            mu['istight'] = isTightMuon(mu.pt,mu.eta,mu.dxy,mu.dz,mu.iso,mu.tight_id,self._year)

            #print ("muon pt maximum argument:",mu.pt.argmax())
            #print ("subleading mu:",mu.pt.argsort()[:,1:2])
            #print("leading mu:", mu[:,:1])
            #print("subleading mu:", mu[:,:2])
            #print("leading mu with arg:", mu[mu.pt.argmax()])
            leading_mu = mu[mu.pt.argmax()]
            leading_mu = leading_mu[leading_mu.istight]
            #subleading_mu = mu[mu.pt.argsort()[:,1:2]]

            mu_loose=mu[mu.isloose]
            mu_tight=mu[mu.istight]

            mu_ntot = mu.counts
            mu_nloose = mu_loose.counts
            mu_ntight = mu_tight.counts

            tau = Initialize({'pt':df['Tau_pt'],
                              'eta':df['Tau_eta'],
                              'phi':df['Tau_phi'],
                              'mass':df['Tau_mass']})

            for key in tau_id[self._year]:
                tau[key] = tau.pt.zeros_like()
                if tau_id[self._year][key] in df:
                    tau[key] = df[tau_id[self._year][key]]


            tau['isclean'] =~tau.match(mu_loose,0.3)&~tau.match(e_loose,0.3)
            tau['isloose']=isLooseTau(tau.pt,tau.eta,tau.decayMode,tau.id,self._year)&tau.isclean
            tau_loose=tau[tau.isloose]

            tau_ntot=tau.counts
            tau_nloose=tau_loose.counts

            pho = Initialize({'pt':df['Photon_pt'],
                              'eta':df['Photon_eta'],
                              'phi':df['Photon_phi'],
                              'mass':df['Photon_mass']})

            for key in pho_id[self._year]:
                pho[key] = pho.pt.zeros_like()
                if pho_id[self._year][key] in df:
                    pho[key] = df[pho_id[self._year][key]]

            pho['isclean'] =~pho.match(e_loose,0.4)
            pho['isloose']=isLoosePhoton(pho.pt,pho.eta,pho.loose_id,pho.eleveto,self._year)&pho.isclean
            pho['istight']=isTightPhoton(pho.pt,pho.eta,pho.tight_id,pho.eleveto,self._year)&pho.isclean

            leading_pho = pho[pho.pt.argmax()]
            leading_pho = leading_pho[leading_pho.istight]

            pho_loose=pho[pho.isloose]
            pho_tight=pho[pho.istight]

            pho_ntot=pho.counts
            pho_nloose=pho_loose.counts
            pho_ntight=pho_tight.counts

            fj = Initialize({'pt':df['AK15Puppi_pt'],
                             'eta':df['AK15Puppi_eta'],
                             'phi':df['AK15Puppi_phi'],
                             'mass':df['AK15Puppi_mass']})

            for key in fj_id[self._year]:
                fj[key] = fj.pt.zeros_like()
                if fj_id[self._year][key] in df:
                    fj[key] = df[fj_id[self._year][key]]

            fj['isgood'] = isGoodFatJet(fj.pt, fj.eta, fj.id)
            fj['isclean'] =~fj.match(pho_loose,1.5)&~fj.match(mu_loose,1.5)&~fj.match(e_loose,1.5)&fj.isgood
            #fj['isclean'] =~fj.match(pho_tight,1.5)&~fj.match(mu_tight,1.5)&~fj.match(e_tight,1.5)&fj.isgood

            for key in deep[self._year]:
                fj[key] = fj.pt.zeros_like()
                if deep[self._year][key] in df:
                    fj[key] = df[deep[self._year][key]]

            #fj['probQCD'] = fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
            fj['TopTagger'] = fj.probTbcq+fj.probTbqq
            fj['DarkHiggsTagger'] = fj.probZbb + fj.probHbb #/ (fj.probZbb+fj.probHbb+fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq+fj.probHcc+fj.probHqqqq+fj.probQCD)
            fj['VvsQCDTagger'] = (fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq) / (fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq+fj.probQCDothers+fj.probQCDcc)

            leading_fj = fj[fj.pt.argmax()]
            leading_fj = leading_fj[leading_fj.isclean]

            fj_good = fj[fj.isgood]
            fj_clean=fj[fj.isclean]

            fj_ntot=fj.counts
            fj_ngood=fj_good.counts
            fj_nclean=fj_clean.counts

            j = Initialize({'pt':df['Jet_pt'],
                            'eta':df['Jet_eta'],
                            'phi':df['Jet_phi'],
                            'mass':df['Jet_mass']})

            #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
            j['deepcsv'] = df['Jet_btagDeepB']
            j['deepflv'] = df['Jet_btagDeepFlavB']

            for key in j_id[self._year]:
                j[key] = j.pt.zeros_like()
                if j_id[self._year][key] in df:
                    j[key] = df[j_id[self._year][key]]

            j['isgood'] = isGoodJet(j.pt, j.eta, j.id, j.nhf, j.nef, j.chf, j.cef)
            j['isHEM'] = isHEMJet(j.pt, j.eta, j.phi)
            j['isclean'] = ~j.match(e_loose,0.4)&~j.match(mu_loose,0.4)&~j.match(pho_loose,0.4)&j.isgood
            #j['isclean'] = ~j.match(e_tight,0.4)&~j.match(mu_tight,0.4)&~j.match(pho_tight,0.4)&j.isgood
            j['isiso'] =  ~(j.match(fj_clean,1.5))&j.isclean
            j['isdcsvL'] = (j.deepcsv>0.1241)&j.isiso
            j['isdflvL'] = (j.deepflv>0.0494)&j.isiso
            j['isdcsvM'] = (j.deepcsv>0.4184)&j.isiso
            j['isdflvM'] = (j.deepflv>0.2770)&j.isiso
            j['isdcsvT'] = (j.deepcsv>0.7527)&j.isiso
            j['isdflvT'] = (j.deepflv>0.7264)&j.isiso

            leading_j = j[j.pt.argmax()]
            leading_j = leading_j[leading_j.isclean]

            j_good = j[j.isgood]
            j_clean = j[j.isclean]
            j_iso = j[j.isiso]
            j_dcsvL = j[j.isdcsvL]
            j_dflvL = j[j.isdflvL]
            j_dcsvM = j[j.isdcsvM]
            j_dflvM = j[j.isdflvM]
            j_dcsvT = j[j.isdcsvT]
            j_dflvT = j[j.isdflvT]
            j_HEM = j[j.isHEM]

            j_ntot=j.counts
            j_ngood=j_good.counts
            j_nclean=j_clean.counts
            j_niso=j_iso.counts
            j_ndcsvL=j_dcsvL.counts
            j_ndflvL=j_dflvL.counts
            j_ndcsvM=j_dcsvM.counts
            j_ndflvM=j_dflvM.counts
            j_ndcsvT=j_dcsvT.counts
            j_ndflvT=j_dflvT.counts
            j_nHEM = j_HEM.counts

            ###
            #Calculating derivatives
            ###
            ele_pairs = e_loose.distincts()
            diele = leading_e
            leading_diele = leading_e
            if ele_pairs.i0.content.size>0:
                diele = ele_pairs.i0+ele_pairs.i1
                leading_diele = diele[diele.pt.argmax()]

            mu_pairs = mu_loose.distincts()
            dimu = leading_mu
            leading_dimu = leading_mu
            if mu_pairs.i0.content.size>0:
                dimu = mu_pairs.i0+mu_pairs.i1
                leading_dimu = dimu[dimu.pt.argmax()]

            u={}
            u["iszeroL"] = met
            u["isoneM"] = met+leading_mu.sum()
            u["isoneE"] = met+leading_e.sum()
            u["istwoM"] = met+leading_dimu.sum()
            u["istwoE"] = met+leading_diele.sum()
            u["isoneA"] = met+leading_pho.sum()

            lepSys={}
            lepSys["iszeroL"] = met
            lepSys["isoneM"] = leading_mu.sum()
            lepSys["isoneE"] = leading_e.sum()
            lepSys["istwoM"] = leading_dimu.sum()
            lepSys["istwoE"] = leading_diele.sum()
            lepSys["isoneA"] = leading_pho.sum()

            leadlepton={}
            leadlepton["iszeroL"] = met
            leadlepton["isoneM"] = leading_mu.sum()
            leadlepton["isoneE"] = leading_e.sum()
            leadlepton["istwoM"] = leading_mu.sum()
            leadlepton["istwoE"] = leading_e.sum()
            leadlepton["isoneA"] = leading_pho.sum()

            ###
            #Calculating weights
            ###

            ###
            # For MC, retrieve the LHE weights, to take into account NLO destructive interference, and their sum
            ###

            genw = np.ones_like(df['MET_pt'])
            sumw = 1.
            wnlo = np.ones_like(df['MET_pt'])
            if self._xsec[dataset] != -1:
                genw = df['genWeight']
                sumw = genw.sum()

                gen = Initialize({'pt':df['GenPart_pt'],
                                  'eta':df['GenPart_eta'],
                                  'phi':df['GenPart_phi'],
                                  'mass':df['GenPart_mass'],
                                  'pdgid':df['GenPart_pdgId'],
                                  'status':df['GenPart_status'], 
                                  'flags':df['GenPart_statusFlags'],
                                  'motherid':df['GenPart_genPartIdxMother']})

                genLastCopy = gen[gen.flags&(1 << 13)==0]
                genTops = genLastCopy[abs(genLastCopy.pdgid)==6 ]
                genWs   = genLastCopy[abs(genLastCopy.pdgid)==24]
                genZs   = genLastCopy[abs(genLastCopy.pdgid)==23]
                genAs   = genLastCopy[abs(genLastCopy.pdgid)==22]
                genHs   = genLastCopy[abs(genLastCopy.pdgid)==25]

                #isTT = (genTops.counts==2)
                #isW  = (genTops.counts==0)&(genWs.counts==1)&(genZs.counts==0)&(genAs.counts==0)&(genDs.counts==0)&(genHs.counts==0)
                #isZ  = (genTops.counts==0)&(genWs.counts==0)&(genZs.counts==1)&(genAs.counts==0)&(genDs.counts==0)&(genHs.counts==0)
                #isA  = (genTops.counts==0)&(genWs.counts==0)&(genZs.counts==0)&(genAs.counts==1)&(genDs.counts==0)&(genHs.counts==0)
                #isDY = (genTops.counts==0)&(genWs.counts==0)&(genZs.counts==0)&(genAs.counts==0)&(genDs.counts==1)&(genHs.counts==0)

                if  ('TTJets'   in dataset): wnlo = np.sqrt(get_ttbar_weight(genTops[0].pt.sum()) * get_ttbar_weight(genTops[1].pt.sum()))
                elif('WJets'    in dataset): wnlo = get_nlo_weight('w',genWs[0].pt.sum())
                elif('DY' in dataset or 'ZJets' in dataset): wnlo = get_nlo_weight('z',genZs[0].pt.sum())
                elif('GJets' in dataset): wnlo = get_nlo_weight('a',genAs[0].pt.sum())    
                

            ###
            # Calculate PU weight and systematic variations
            ###

            nvtx = df['PV_npvs']
            pu,puUp,puDown = get_pu_weight(nvtx,self._year)

            ###
            #Importing the MET filters per year from metfilters.py and constructing the filter boolean
            ###

            met_filters = {}
            for flag in met_filter_flags[self._year]:
                if flag in df:
                    met_filters[flag] = df[flag]

            ###
            #Importing the trigger paths per year from trigger.py and constructing the trigger boolean
            ###

            pass_trig = {}
            met_trigger = {}
            for path in met_trigger_paths[self._year]:
                if path in df:
                    met_trigger[path] = df[path]
            passMetTrig = False
            for path in met_trigger:
                passMetTrig |= met_trigger[path]

            singleele_trigger = {}
            for path in singleele_trigger_paths[self._year]:
                if path in df:
                    singleele_trigger[path] = df[path]
            passSingleEleTrig = False
            for path in singleele_trigger:
                passSingleEleTrig |= singleele_trigger[path]

            singlepho_trigger = {}
            for path in singlepho_trigger_paths[self._year]:
                if path in df:
                    singlepho_trigger[path] = df[path]
            passSinglePhoTrig = False
            for path in singlepho_trigger:
                passSinglePhoTrig |= singlepho_trigger[path]

            pass_trig['iszeroL'] = passMetTrig
            pass_trig['isoneM'] = passMetTrig
            pass_trig['istwoM'] = passMetTrig
            pass_trig['isoneE'] = passSingleEleTrig
            pass_trig['istwoE'] = passSingleEleTrig
            pass_trig['isoneA'] =passSinglePhoTrig

            ###
            # Trigger efficiency weight
            ###

            trig = {}
            trig['iszeroL'] = get_met_trig_weight(u["iszeroL"].pt,self._year)
            trig['isoneM'] = get_met_trig_weight(u["isoneM"].pt,self._year)
            trig['istwoM'] = get_met_zmm_trig_weight(u["istwoM"].pt,self._year)
            trig['isoneE'] = get_ele_trig_weight(leading_e.eta.sum(), leading_e.pt.sum(),
                                                      np.full_like(leading_e.eta.sum(),-99),
                                                      np.full_like(leading_e.pt.sum(),-99),self._year)
            trig['istwoE'] = trig['isoneE']
            if ele_pairs.i0.content.size>0:
                trig['istwoE'] =get_ele_trig_weight(ele_pairs[diele.pt.argmax()].i0.eta.sum(),ele_pairs[diele.pt.argmax()].i0.pt.sum(),
                                                         ele_pairs[diele.pt.argmax()].i1.eta.sum(),ele_pairs[diele.pt.argmax()].i1.pt.sum(),self._year)
            trig['isoneA'] = get_pho_trig_weight(leading_pho.pt.sum(),self._year)

            ###
            #Event selection
            ###

            selections = processor.PackedSelection()

            selections.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
            selections.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(tau_nloose==0)&(pho_nloose==0))
            selections.add('isoneE', (e_ntight==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(met.pt>50))
            selections.add('istwoM', (e_nloose==0) & (mu_ntight>=1) & (mu_nloose==2) & (tau_nloose==0)&(pho_nloose==0)&(leading_dimu.mass.sum()>60) & (leading_dimu.mass.sum()<120))
            selections.add('istwoE', (e_ntight>=1) & (e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(leading_diele.mass.sum()>60)&(leading_diele.mass.sum()<120))
            selections.add('isoneA', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_ntight==1))
            selections.add('topveto', (leading_fj.TopTagger.sum()<0.25))
            selections.add('noextrab', (j_ndflvL==0))
            selections.add('extrab', (j_ndflvL>0))
            selections.add('ismonohs', (leading_fj.DarkHiggsTagger.sum()>0.2))
            selections.add('ismonoV', ~(leading_fj.DarkHiggsTagger.sum()>0.2)&(leading_fj.VvsQCDTagger.sum()>0.8))
            selections.add('ismonojet', ~(leading_fj.DarkHiggsTagger.sum()>0.2)&~(leading_fj.VvsQCDTagger.sum()>0.8))
            selections.add('noHEMj', (j_nHEM==0))

            #selections.add('istwoM', (e_nloose==0) & (mu_ntight==1) & (mu_nloose==2) & (tau_nloose==0)&(pho_nloose==0)&(leading_dimu.mass.sum()>60) & (leading_dimu.mass.sum()<120))
            #selections.add('istwoE', (e_ntight==1)&(e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(leading_diele.mass.sum()>60)&(leading_diele.mass.sum()<120))

            ###
            #Adding weights and selections
            ###

            weights = {}
            regions = {}
            for k in self._selected_regions[dataset]:
                weights[k] = processor.Weights(df.size)
                weights[k].add('nlo',wnlo)
                weights[k].add('genw',genw)
                weights[k].add('pileup',pu,puUp,puDown)
                weights[k].add('passMetFilters',np.prod([met_filters[key] for key in met_filters], axis=0))
                weights[k].add('trig', trig[k])
                weights[k].add('pass_trig', pass_trig[k])


                #baggy = (fj_nclean>0)&(fj_clean.pt.max()>160)&(abs(u[k].delta_phi(j_clean)).min()>0.8)&(u[k].pt>250)
                #skinny = (j_nclean>0) & (j_clean.pt.max()>100) & (abs(u[k].delta_phi(j_clean)).min()>0.5) & (u[k].pt>250)
                #skinny_no_baggy = ~baggy&skinny
                selections.add(k+'baggy', (fj_nclean>0)&(fj_clean.pt.max()>160)&(abs(u[k].delta_phi(j_clean)).min()>0.8)&(u[k].pt>250))
                #selections.add(k+'baggy', (fj_nclean>0)&(fj_clean.pt.max()>160)&(u[k].pt>250))

                regions[k+'_baggy'] =  {k,k+'baggy','noHEMj'}
                regions[k+'_topveto'] =  {k,k+'baggy','topveto','noextrab','noHEMj'}
                for s in ['ismonohs','ismonoV','ismonojet']:
                    regions[k+'_'+s] = {k,k+'baggy','topveto','noextrab','noHEMj',s}
                regions[k+'_ismonohs'+'_extrab'] = {k,k+'baggy','topveto','extrab','ismonohs','noHEMj'}

            variables = {}
            variables['j1pt'] = leading_j.pt
            variables['j1eta'] = leading_j.eta
            variables['j1phi'] = leading_j.phi
            variables['fj1pt'] = leading_fj.pt
            variables['fj1eta'] = leading_fj.eta
            variables['fj1phi'] = leading_fj.phi
            variables['e1pt'] = leading_e.pt
            variables['e1phi'] = leading_e.phi
            variables['e1eta'] = leading_e.eta
            variables['dielemass'] = leading_diele.mass
            variables['mu1pt'] = leading_mu.pt
            variables['mu1phi'] = leading_mu.phi
            variables['mu1eta'] = leading_mu.eta
            variables['dimumass'] = leading_dimu.mass
            variables['njets'] = j_nclean
            variables['ndcsvL'] = j_ndcsvL
            variables['ndflvL'] = j_ndflvL
            variables['ndcsvM'] = j_ndcsvM
            variables['ndflvM'] = j_ndflvM
            variables['ndcsvT'] = j_ndcsvT
            variables['ndflvT'] = j_ndflvT
            variables['nfjtot'] = fj_ntot
            variables['nfjgood'] = fj_ngood
            variables['nfjclean'] = fj_nclean
            variables['fjmass'] = leading_fj.mass
            variables['TopTagger'] = leading_fj.TopTagger
            variables['DarkHiggsTagger'] = leading_fj.DarkHiggsTagger
            variables['VvsQCDTagger'] = leading_fj.VvsQCDTagger
            variables['probTbcq']      = leading_fj.probTbcq
            variables['probTbqq']      = leading_fj.probTbqq
            variables['probTbc']       = leading_fj.probTbc
            variables['probTbq']       = leading_fj.probTbq
            variables['probWcq']       = leading_fj.probWcq
            variables['probWqq']       = leading_fj.probWqq
            variables['probZbb']       = leading_fj.probZbb
            variables['probZcc']       = leading_fj.probZcc
            variables['probZqq']       = leading_fj.probZqq
            variables['probHbb']       = leading_fj.probHbb
            variables['probHcc']       = leading_fj.probHcc
            variables['probHqqqq']     = leading_fj.probHqqqq
            variables['probQCDbb']     = leading_fj.probQCDbb
            variables['probQCDcc']     = leading_fj.probQCDcc
            variables['probQCDb']      = leading_fj.probQCDb
            variables['probQCDc']      = leading_fj.probQCDc
            variables['probQCDothers'] = leading_fj.probQCDothers

            hout = self.accumulator.identity()
            i = 0
            while i < len(self._selected_regions[dataset]):
                r = self._selected_regions[dataset][i]
                weight = weights[r].weight()
                for s in ['ismonohs','ismonoV','ismonojet','baggy','topveto','ismonohs_extrab']:
                    cut = selections.all(*regions[r+'_'+s])
                    flat_variables = {k: v[cut].flatten() for k, v in variables.items()}
                    flat_weights = {k: (~np.isnan(v[cut])*weight[cut]).flatten() for k, v in variables.items()}
                    for histname, h in hout.items():
                        if not isinstance(h, hist.Hist):
                            continue
                        if histname == 'sumw':
                            h.fill(dataset=dataset, sumw=1, weight=sumw)
                        elif histname == 'recoil':
                            h.fill(dataset=dataset, region=r, jet_selection=s, recoil=u[r].pt, weight=weight*cut)
                        elif histname == 'CaloMinusPfOverRecoil':
                            h.fill(dataset=dataset, region=r, jet_selection=s, CaloMinusPfOverRecoil= abs(calomet.pt - met.pt) / u[r].pt, weight=weight*cut)
                        elif histname == 'mindphi':
                            h.fill(dataset=dataset, region=r, jet_selection=s, mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=weight*cut)
                        elif histname == 'diledphi':
                            h.fill(dataset=dataset, region=r, jet_selection=s, diledphi=abs(lepSys[r].delta_phi(j_clean)).min(), weight=weight*cut)
                        elif histname == 'ledphi':
                            h.fill(dataset=dataset, region=r, jet_selection=s, ledphi=abs(leadlepton[r].delta_phi(j_clean)).min(), weight=weight*cut)
                        elif histname == 'recoilVSmindphi':
                            h.fill(dataset=dataset, region=r, jet_selection=s, recoil=u[r].pt, mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=weight*cut)
                        else:
                            h.fill(dataset=dataset, region=r, jet_selection=s, **flat_variables, weight=flat_weights[histname])
                i += 1
            return hout

    def postprocess(self, accumulator):
            scale = {}
            for d in accumulator['sumw'].identifiers('dataset'):
                dataset = d.name
                if self._xsec[dataset]!= -1: scale[dataset] = self._lumi*self._xsec[dataset]
                else: scale[dataset] = 1

            for histname, h in accumulator.items():
                if histname == 'sumw': continue
                if isinstance(h, hist.Hist):
                    h.scale(scale, axis="dataset")

            return accumulator

