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
from utils.corrections import get_bad_ecal_weight
from utils.ids import e_id, isLooseElectron, isTightElectron
from utils.ids import mu_id, isLooseMuon, isTightMuon
from utils.ids import tau_id, isLooseTau
from utils.ids import pho_id, isLoosePhoton, isTightPhoton
from utils.ids import j_id, fj_id, isGoodJet, isGoodFatJet
from utils.metfilters import met_filter_flags
from utils.deep import deep

samples = {
    "iszeroL":('ZJets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','MET'),
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
                'mindphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mindphi","Min dPhi(MET,AK4s)",15,0,6.28)),
                'j1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("j1pt","AK4 Leading Jet Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
                'j1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("j1eta","AK4 Leading Jet Eta",35,-3.5,3.5)),
                'j1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("j1phi","AK4 Leading Jet Phi",35,-3.5,3.5)),
                'fj1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fj1pt","AK15 Leading Jet Pt",[200.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
                'fj1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fj1eta","AK15 Leading Jet Eta",35,-3.5,3.5)),
                'fj1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fj1phi","AK15 Leading Jet Phi",35,-3.5,3.5)),
                'njets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("njets","AK4 Number of Jets",6,0,5)),
                'nfjets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("nfjets","AK15 Number of Jets",4,0,3)),
                'fjmass': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("fjmass","AK15 Jet Mass",50,20,250)),
                'e1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("e1pt","Leading Electron Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
                'e1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("e1eta","Leading Electron Eta",48,-2.4,2.4)),
                'e1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("e1phi","Leading Electron Phi",64,-3.2,3.2)),
                'mu1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mu1pt","Leading Muon Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
                'mu1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mu1eta","Leading Muon Eta",48,-2.4,2.4)),
                'mu1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("mu1phi","Leading Muon Phi",64,-3.2,3.2)),
                'TvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("TvsQCD","TvsQCD",15,0,1)),
                'hSvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("hSvsQCD","hSvsQCD",15,0,1)),
                'VvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("jet_selection", "JetSelection"), hist.Bin("VvsQCD","VvsQCD",15,0,1)),
                })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):

            dataset = df['dataset']
            genw = np.ones_like(df['MET_pt'])
            sumw = 1.

            ###
            # For MC, retrieve the LHE weights, to take into account NLO destructive interference, and their sum
            ###
            if self._xsec[dataset] != -1:
                genw = df['genWeight']
                sumw = genw.sum()

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

            empty_jagged = awkward.JaggedArray.fromcounts(np.ones_like(met.pt, dtype=int),np.zeros_like(met.pt))
            empty_obj = Initialize({'pt':empty_jagged,
                                    'eta':empty_jagged,
                                    'phi':empty_jagged,
                                    'mass':empty_jagged})

            e = Initialize({'pt':df['Electron_pt'],
                            'eta':df['Electron_eta'],
                            'phi':df['Electron_phi'],
                            'mass':df['Electron_mass']})

            #checking content size, if zero, initialize to empty object
            if not (e.content.size > 0): e = empty_obj
            
            for key in e_id[self._year]:
                e[key] = e.pt.zeros_like()
                if e_id[self._year][key] in df:
                    e[key] = df[e_id[self._year][key]]
            e['isloose'] = isLooseElectron(e.pt,e.eta,e.dxy,e.dz,e.iso,e.loose_id,self._year)
            e['istight'] = isTightElectron(e.pt,e.eta,e.dxy,e.dz,e.iso,e.loose_id,self._year)

            e_loose = e
            if e[e.isloose].content.size > 0: e_loose = e[e.isloose]

            e_tight = e
            if e[e.istight].content.size > 0: e_tight = e[e.istight]
            
            e_ntot = e.counts
            e_nloose = e_loose.counts
            e_ntight = e_tight.counts

            mu = Initialize({'pt':df['Muon_pt'],
                             'eta':df['Muon_eta'],
                             'phi':df['Muon_phi'],
                             'mass':df['Muon_mass']})

            #checking content size, if zero, initialize to empty object
            if not (mu.content.size > 0): mu = empty_obj
            
            for key in mu_id[self._year]:
                mu[key] = mu.pt.zeros_like()
                if mu_id[self._year][key] in df:
                    mu[key] = df[mu_id[self._year][key]]
                    
            mu['isloose'] = isLooseMuon(mu.pt,mu.eta,mu.dxy,mu.dz,mu.iso,self._year)
            mu['istight'] = isTightMuon(mu.pt,mu.eta,mu.dxy,mu.dz,mu.iso,mu.tight_id,self._year)

            mu_loose = mu
            if mu[mu.isloose].content.size > 0:
                mu_loose=mu[mu.isloose]

            mu_tight = mu
            if mu[mu.istight].content.size > 0:
                mu_tight=mu[mu.istight]

            mu_ntot = mu.counts
            mu_nloose = mu_loose.counts
            mu_ntight = mu_tight.counts

            tau = Initialize({'pt':df['Tau_pt'],
                              'eta':df['Tau_eta'],
                              'phi':df['Tau_phi'],
                              'mass':df['Tau_mass']})

            #checking content size, if zero, initialize to empty object
            if not (tau.content.size > 0): tau = empty_obj

            for key in tau_id[self._year]:
                tau[key] = tau.pt.zeros_like()
                if tau_id[self._year][key] in df:
                    tau[key] = df[tau_id[self._year][key]]

            
            tau['isloose']=isLooseTau(tau.pt,tau.eta,tau.decayMode,tau.id,self._year)
            tau_loose = tau

            if tau[tau.isloose].content.size>0:
                tau_loose=tau[tau.isloose]
                
            tau_ntot=tau.counts
            tau_nloose=tau_loose.counts

            pho = Initialize({'pt':df['Photon_pt'],
                              'eta':df['Photon_eta'],
                              'phi':df['Photon_phi'],
                              'mass':df['Photon_mass']})

            #checking content size, if zero, initialize to empty object
            if not (pho.content.size > 0): pho = empty_obj

            for key in pho_id[self._year]:
                pho[key] = pho.pt.zeros_like()
                if pho_id[self._year][key] in df:
                    pho[key] = df[pho_id[self._year][key]]
                    
            pho['isloose']=isLoosePhoton(pho.pt,pho.eta,pho.loose_id,pho.eleveto,self._year)
            pho['istight']=isTightPhoton(pho.pt,pho.eta,pho.tight_id,pho.eleveto,self._year)

            pho_loose = pho
            if pho[pho.isloose].content.size>0:
                pho_loose=pho[pho.isloose]

            pho_tight = pho
            if pho[pho.istight].content.size>0:
                pho_tight=pho[pho.istight]
                
            pho_ntot=pho.counts
            pho_nloose=pho_loose.counts
            pho_ntight=pho_tight.counts

            fj = Initialize({'pt':df['AK15Puppi_pt'],
                             'eta':df['AK15Puppi_eta'],
                             'phi':df['AK15Puppi_phi'],
                             'mass':df['AK15Puppi_mass']})

            #checking content size, if zero, initialize to empty object
            if not (fj.content.size > 0): fj = empty_obj

            for key in fj_id[self._year]:
                fj[key] = fj.pt.zeros_like()
                if fj_id[self._year][key] in df:
                    fj[key] = df[fj_id[self._year][key]]

            fj['isgood'] = isGoodFatJet(fj.pt, fj.eta, fj.id)
            fj['isclean'] =~fj.match(pho,1.5)&~fj.match(mu,1.5)&~fj.match(e,1.5)&fj.isgood

            for key in deep[self._year]:
                fj[key] = fj.pt.zeros_like()
                if deep[self._year][key] in df:
                    fj[key] = df[deep[self._year][key]]
            
            fj['probQCD'] = fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
            fj['TvsQCD'] = fj.probTbcq+fj.probTbqq+fj.probTbc+fj.probTbq
            fj['hSvsQCD'] = (fj.probZbb + fj.probHbb) / (fj.probZbb+fj.probHbb+fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq+fj.probHcc+fj.probHqqqq+fj.probQCD)
            fj['VvsQCD'] = (fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq) / (fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq+fj.probHcc+fj.probHqqqq+fj.probQCD)

            fj_good=fj
            if fj[fj.isgood].content.size > 0: fj_good = fj[fj.isgood]
            fj_clean=fj
            if fj[fj.isclean].content.size > 0: fj_clean=fj[fj.isclean]
            
            fj_ntot=fj.counts
            fj_ngood=fj_good.counts
            fj_nclean=fj_clean.counts

            j = Initialize({'pt':df['Jet_pt'],
                            'eta':df['Jet_eta'],
                            'phi':df['Jet_phi'],
                            'mass':df['Jet_mass']})

            #checking content size, if zero, initialize to empty object
            if not (j.content.size > 0): j = empty_obj

            for key in j_id[self._year]:
                j[key] = j.pt.zeros_like()
                if j_id[self._year][key] in df:
                    j[key] = df[j_id[self._year][key]]

            j['isgood'] = isGoodJet(j.pt, j.eta, j.id, j.nhf, j.nef, j.chf, j.cef)
            j['isclean'] = ~j.match(e,0.4)&~j.match(mu,0.4)&~j.match(pho,0.4)&j.isgood
            j['isiso'] =  ~(j.match(fj,1.5))&j.isclean

            j_good = j
            if j[j.isgood].content.size > 0: j_good = j[j.isgood]
            j_clean = j
            if j[j.isclean].content.size > 0: j_clean = j[j.isclean]
            j_iso = j
            if j[j.isiso].content.size > 0: j_iso = j[j.isiso]

            j_ntot=j.counts
            j_ngood=j_good.counts
            j_nclean=j_clean.counts
            j_niso=j_iso.counts

            ###
            #Getting leading pT objects
            ###
            leading_mu = mu_tight[mu_tight.pt.argmax()]
            leading_e = e_tight[e_tight.pt.argmax()]
            leading_pho = pho_tight[pho_tight.pt.argmax()]
            leading_j = j_clean[j_clean.pt.argmax()]
            leading_fj = fj_clean[fj_clean.pt.argmax()]

            ###
            #Calculating derivatives
            ###
            ele_pairs = e_loose.distincts()
            diele = e_loose[e_loose.pt.argmax()]
            if ele_pairs.i0.content.size>0:
                diele = ele_pairs.i0+ele_pairs.i1
            leading_diele = diele[diele.pt.argmax()]
  
            mu_pairs = mu_loose.distincts()
            dimu = mu_loose[mu_loose.pt.argmax()]
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
            
            ###
            #Calculating weights
            ###

            weights = {}
            for k in u.keys():
                weights[k] = processor.Weights(df.size)
                
            weights['iszeroL'].add('trig', get_met_trig_weight(u["iszeroL"].pt,self._year))
            weights['iszeroL'].add('pass_trig',passMetTrig)
            weights['isoneM'].add('trig', get_met_trig_weight(u["isoneM"].pt,self._year))
            weights['isoneM'].add('pass_trig',passMetTrig)
            weights['istwoM'].add('trig', get_met_zmm_trig_weight(u["istwoM"].pt,self._year))
            weights['istwoM'].add('pass_trig',passMetTrig)
            weights['isoneE'].add('trig',get_ele_trig_weight(leading_e.eta.sum(), leading_e.pt.sum(), np.full_like(leading_e.eta.sum(),-99),np.full_like(leading_e.pt.sum(),-99),self._year))
            weights['isoneE'].add('passSingleEleTrig',passSingleEleTrig)
            weights['istwoE'].add('trig', get_ele_trig_weight(ele_pairs[diele.pt.argmax()].i0.eta.sum(),ele_pairs[diele.pt.argmax()].i0.pt.sum(),
                                                               ele_pairs[diele.pt.argmax()].i1.eta.sum(),ele_pairs[diele.pt.argmax()].i1.pt.sum(),self._year))
            weights['istwoE'].add('passSingleEleTrig',passSingleEleTrig)
            #weights.add('trig_istwoE', get_ele_trig_weight(leading_e.eta.sum(), leading_e.pt.sum(), 
            #                                                   np.full_like(leading_e.eta.sum(),-99),np.full_like(leading_e.pt.sum(),-99),self._year))
            weights['isoneA'].add('trig', get_pho_trig_weight(leading_pho.pt.sum(),self._year))
            weights['isoneA'].add('passSinglePhoTrig',passSinglePhoTrig)

            wnlo = np.ones_like(df['MET_pt'])
            if self._xsec[dataset] != -1:
                gen = Initialize({'pt':df['GenPart_pt'],
                                  'eta':df['GenPart_eta'],
                                  'phi':df['GenPart_phi'],
                                  'mass':df['GenPart_mass'],
                                  'pdgid':df['GenPart_pdgId'],
                                  'status':df['GenPart_status'], 
                                  'flags':df['GenPart_statusFlags'],
                                  'motherid':df['GenPart_genPartIdxMother']})
        
                genLastCopy = gen[gen.flags&(1 << 13)==0]
                genTops = genLastCopy[abs(genLastCopy.pdgid)==6]
                genWs = genLastCopy[abs(genLastCopy.pdgid)==24]
                genZs = genLastCopy[abs(genLastCopy.pdgid)==23]
                genAs = genLastCopy[abs(genLastCopy.pdgid)==22]
                genHs = genLastCopy[abs(genLastCopy.pdgid)==25]

                isTT = (genTops.counts==2)
                isW  = (genTops.counts==0)&(genWs.counts==1)&(genZs.counts==0)&(genAs.counts==0)&(genHs.counts==0)
                isZ  = (genTops.counts==0)&(genWs.counts==0)&(genZs.counts==1)&(genAs.counts==0)&(genHs.counts==0)
                isA  = (genTops.counts==0)&(genWs.counts==0)&(genZs.counts==0)&(genAs.counts==1)&(genHs.counts==0)

                if('TTJets' in dataset): wnlo = np.sqrt(get_ttbar_weight(genTops[0].pt.sum()) * get_ttbar_weight(genTops[1].pt.sum()))
                elif('WJets' in dataset): wnlo = get_nlo_weight('w',genWs[0].pt.sum())
                elif('DY' in dataset or 'ZJets' in dataset): wnlo = get_nlo_weight('z',genZs[0].pt.sum())
                elif('GJets' in dataset): wnlo = get_nlo_weight('a',genAs[0].pt.sum())

            for k in u.keys():
                weights[k].add('nlo',wnlo)
                weights[k].add('genw',genw)
                weights[k].add('pileup',pu,puUp,puDown)
                weights[k].add('passMetFilters',np.prod([met_filters[key] for key in met_filters], axis=0))


            ###
            #Event selection
            ###

            selections = processor.PackedSelection()

            selections.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
            selections.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(tau_nloose==0)&(pho_nloose==0))
            selections.add('isoneE', (e_ntight==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
            selections.add('istwoM', (e_nloose==0) & (mu_ntight==1) & (mu_nloose==2) & (tau_nloose==0)&(pho_nloose==0)&(leading_dimu.mass.sum()>60) & (leading_dimu.mass.sum()<120))
            selections.add('istwoE', (e_ntight==1)&(e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(leading_diele.mass.sum()>60)&(leading_diele.mass.sum()<120))
            selections.add('isoneA', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_ntight==1))

            for k in u.keys():
                selections.add(k+'baggy', (fj_nclean>0)&(fj_clean.pt.max()>200)&(abs(u[k].delta_phi(j_clean)).min()>0.8)&(u[k].pt>250))
                selections.add(k+'skinny', ~((fj_nclean>0)&(fj_clean.pt.max()>200))&(j_nclean>0)&(j_clean.pt.max()>100)&(abs(u[k].delta_phi(j_clean)).min()>0.5)&(u[k].pt>250))
                selections.add(k+'inclusive', (~((fj_nclean>0)&(fj_clean.pt.max()>200))&(j_nclean>0)&(j_clean.pt.max()>100)&(abs(u[k].delta_phi(j_clean)).min()>0.5)&(u[k].pt>250)) | ((fj_nclean>0)&(fj_clean.pt.max()>200)&(abs(u[k].delta_phi(j_clean)).min()>0.8)&(u[k].pt>250)))

            regions = {}
            for k in u.keys():
                for j in ["baggy","skinny","inclusive"]:
                    regions[k+'_'+j] = {k,k+j}
                    
            variables = {}
            variables['j1pt'] = leading_j.pt.sum()
            variables['j1eta'] = leading_j.eta.sum()
            variables['j1phi'] = leading_j.phi.sum()
            variables['fj1pt'] = leading_fj.pt.sum()
            variables['fj1eta'] = leading_fj.eta.sum()
            variables['fj1phi'] = leading_fj.phi.sum()
            variables['e1pt'] = leading_e.pt.sum()
            variables['e1phi'] = leading_e.phi.sum()
            variables['e1eta'] = leading_e.eta.sum()
            variables['mu1pt'] = leading_mu.pt.sum()
            variables['mu1phi'] = leading_mu.phi.sum()
            variables['mu1eta'] = leading_mu.eta.sum()
            variables['njets'] = j_nclean
            variables['nfjets'] = fj_nclean
            variables['fjmass'] = leading_fj.mass.sum()
            variables['TvsQCD'] = leading_fj.TvsQCD.sum()
            variables['hSvsQCD'] = leading_fj.hSvsQCD.sum()
            variables['VvsQCD'] = leading_fj.VvsQCD.sum()
            #variables['recoil'] = u[r].pt
            #variables['CaloMinusPfOverRecoil'] = abs(calomet.pt - met.pt) / u[r].pt
            #variables['mindphi'] = abs(u[r].delta_phi(j_clean)).min()

            hout = self.accumulator.identity()
            for histname, h in hout.items():
                if not isinstance(h, hist.Hist):
                    continue
                i = 0
                if histname == 'sumw':
                    h.fill(dataset=dataset, sumw=1, weight=sumw)                
                else:
                    while i < len(self._selected_regions):
                        r = self._selected_regions[i]
                        for s in ["baggy","skinny","inclusive"]:
                            weight = weights[r].weight()
                            cut = selections.all(*regions[r+'_'+s])
                            if histname == 'recoil':
                                h.fill(dataset=dataset, region=r, jet_selection=s, recoil=u[r].pt, weight=weight*cut)
                            elif histname == 'CaloMinusPfOverRecoil':
                                h.fill(dataset=dataset, region=r, jet_selection=s, CaloMinusPfOverRecoil= abs(calomet.pt - met.pt) / u[r].pt, weight=weight*cut)
                            elif histname == 'mindphi':
                                h.fill(dataset=dataset, region=r, jet_selection=s, mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=weight*cut)
                            else:
                                h.fill(dataset=dataset, region=r, jet_selection=s, **variables, weight=weight*cut)
                        i += 1
            return hout

    def postprocess(self, accumulator):

            scale = {}
            for d in accumulator['sumw'].identifiers('dataset'):
                dataset = d.name
                if self._xsec[dataset]!= -1: scale[dataset] = self._lumi*self._xsec[dataset]
                else: scale[dataset] = 1

            for h in accumulator.values():
                if isinstance(h, hist.Hist):
                    h.scale(scale, axis="dataset")

            return accumulator
                    
            
