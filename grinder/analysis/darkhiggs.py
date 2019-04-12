#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from Builder import Initialize
from fnal_column_analysis_tools import hist
from utils.triggers import met_trigger_paths, singleele_trigger_paths, singlepho_trigger_paths
from utils.corrections import get_ttbar_weight, get_nlo_weight, get_pu_weight
from utils.corrections import get_met_trig_weight, get_met_zmm_trig_weight, get_ele_trig_weight, get_pho_trig_weight
from utils.ids import e_id, isLooseElectron, isTightElectron
from utils.ids import mu_id, isLooseMuon, isTightMuon
from utils.ids import tau_id, isLooseTau
from utils.ids import pho_id, isLoosePhoton, isTightPhoton
from utils.metfilters import met_filter_flags

hists = {
    'sumw': hist.Hist("sumw", hist.Cat("dataset", "Primary dataset"), hist.Bin("sumw", "Weight value", [0.])),
    'CaloMinusPfOverRecoil': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("CaloMinusPfOverRecoil","Calo - Pf / Recoil",35,0,1)),
    'recoil': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("recoil","Hadronic Recoil",[250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
    'mindphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("mindphi","Min dPhi(MET,AK4s)",15,0,6.28)),
    'j1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("j1pt","AK4 Leading Jet Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
    'j1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("j1eta","AK4 Leading Jet Eta",35,-3.5,3.5)),
    'j1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("j1phi","AK4 Leading Jet Phi",35,-3.5,3.5)),
    'fj1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("fj1pt","AK15 Leading Jet Pt",[200.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
    'fj1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("fj1eta","AK15 Leading Jet Eta",35,-3.5,3.5)),
    'fj1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("fj1phi","AK15 Leading Jet Phi",35,-3.5,3.5)),
    'njets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("njets","AK4 Number of Jets",6,0,5)),
    'nfjets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("nfjets","AK15 Number of Jets",4,0,3)),
    'fjmass': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("fjmass","AK15 Jet Mass",50,20,250)),
    #'l1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("l1pt","Leading Lepton Pt",[30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
    #'l1eta': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("l1eta","Leading Lepton Eta",35,-3.5,3.5)),
    #'l1phi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("l1phi","Leading Lepton Phi",35,-3.5,3.5)),
    'TvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("TvsQCD","TvsQCD",15,0,1)),
    'hSvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("hSvsQCD","hSvsQCD",15,0,1)),
    'VvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Cat("selection", "Selection"), hist.Bin("VvsQCD","VvsQCD",15,0,1)),

}


samples = {
    "iszeroL":('ZJets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','MET'),
    "isoneM":('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','MET'),
    "isoneE":('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','SingleElectron','EGamma'),
    "istwoM":('WJets','DY','TT','ST','WW','WZ','ZZ','HToBB','MET'),
    "istwoE":('WJets','DY','TT','ST','WW','WZ','ZZ','HToBB','SingleElectron','EGamma'),
    "isoneA":('GJets','QCD','SinglePhoton','EGamma')
}

def analysis(selection, year, xsec, dataset, file):
    weight = {}
    tree = uproot.open(file)["Events"]
    genw = 1
    sumw = 1

    ###
    # For MC, retrieve the LHE weights, to take into account NLO destructive interference, and their sum
    ###
    if xsec != -1:
        genw = tree.array("genWeight")
        run_tree = uproot.open(file)["Runs"]
        sumw = run_tree.array("genEventSumw")[0]
        #nScaleWeights = run_tree.array("nLHEScaleSumw")

    ###
    # Calculate PU weight and systematic variations
    ###

    nvtx = tree.array("PV_npvs")
    weight["pu"],weight["puUp"],weight["puDown"] = get_pu_weight(nvtx,year)
    #print(weight["pu"],weight["puUp"],weight["puDown"])

    ###
    #Importing the MET filters per year from metfilters.py and constructing the filter boolean
    ###

    met_filters = {}
    for flag in met_filter_flags[year]:
        try:
            met_filters[flag] = tree.array(flag)
        except KeyError:
            pass
    passMetFilters = np.prod([met_filters[key] for key in met_filters], axis=0)

    ###
    #Importing the trigger paths per year from trigger.py and constructing the trigger boolean
    ###

    met_trigger = {}
    for path in met_trigger_paths[year]:
        try:
            met_trigger[path] = tree.array(path)
        except KeyError:
            #print("No trigger bit in file for path ",path)
            pass
    passMetTrig = False
    for path in met_trigger:
        passMetTrig |= met_trigger[path]

    singleele_trigger = {}
    for path in singleele_trigger_paths[year]:
        try:
            singleele_trigger[path] = tree.array(path)
        except KeyError:
            #print("No trigger bit in file for path ",path)
            pass
    passSingleEleTrig = False
    for path in singleele_trigger:
        passSingleEleTrig |= singleele_trigger[path] 

    singlepho_trigger = {}
    for path in singlepho_trigger_paths[year]:
        try:
            singlepho_trigger[path] = tree.array(path)
        except KeyError:
            #print("No trigger bit in file for path ",path)
            pass
    passSinglePhoTrig = False
    for path in singlepho_trigger:
        passSinglePhoTrig |= singlepho_trigger[path] 

    ###
    #Initialize physics objects
    ###
    e = Initialize({'pt':tree.array("Electron_pt"),
                    'eta':tree.array("Electron_eta"),
                    'phi':tree.array("Electron_phi"),
                    'mass':tree.array("Electron_mass"),
                    'dxy':tree.array('Electron_dxy'),
                    'dz':tree.array('Electron_dz')})
    for key in e_id[year]:
        try:
            e[key] = tree.array(e_id[year][key])
        except KeyError:
            e[key] = e.pt.zeros_like()
    e['isloose'] = isLooseElectron(e.pt,e.eta,e.dxy,e.dz,e.iso,e.loose_id,year)
    e['istight'] = isTightElectron(e.pt,e.eta,e.dxy,e.dz,e.iso,e.loose_id,year)
    e_loose = e[e.isloose]
    e_tight = e[e.istight]
    e_ntot = e.counts
    e_nloose = e_loose.counts
    e_ntight = e[e.istight].counts

    mu = Initialize({'pt':tree.array("Muon_pt"),
                     'eta':tree.array("Muon_eta"),
                     'phi':tree.array("Muon_phi"),
                     'mass':tree.array("Muon_mass"),
                     'dxy':tree.array('Muon_dxy'),
                     'dz':tree.array('Muon_dz')})
    for key in mu_id[year]:
        try:
            mu[key] = tree.array(mu_id[year][key])
        except KeyError:
            mu[key] = mu.pt.zeros_like()
    mu['isloose'] = isLooseMuon(mu.pt,mu.eta,mu.dxy,mu.dz,mu.iso,year)
    mu['istight'] = isTightMuon(mu.pt,mu.eta,mu.dxy,mu.dz,mu.iso,mu.tight_id,year)
    mu_loose=mu[mu.isloose]
    mu_tight=mu[mu.istight]
    mu_ntot = mu.counts
    mu_nloose = mu_loose.counts
    mu_ntight = mu_tight.counts

    tau = Initialize({'pt':tree.array('Tau_pt'),
                      'eta':tree.array('Tau_eta'),
                      'phi':tree.array('Tau_phi'),
                      'mass':tree.array('Tau_mass')})
    for key in tau_id[year]:
        try:
            tau[key] = tree.array(tau_id[year][key])
        except KeyError:
            tau[key] = tau.pt.zeros_like()
    tau['isloose']=isLooseTau(tau.pt,tau.eta,tau.decayMode,tau.id,year)
    tau_loose=tau[tau.isloose]
    tau_ntot=tau.counts
    tau_nloose=tau_loose.counts

    pho = Initialize({'pt':tree.array('Photon_pt'),
                      'eta':tree.array('Photon_eta'),
                      'phi':tree.array('Photon_phi'),
                      'mass':tree.array('Photon_mass')})
    for key in pho_id[year]:
        try:
            pho[key] = tree.array(pho_id[year][key])
        except KeyError:
            pho[key] = pho.pt.zeros_like()
    pho['isloose']=isLoosePhoton(pho.pt,pho.eta,pho.loose_id,pho.eleveto,year)
    pho['istight']=isTightPhoton(pho.pt,pho.eta,pho.tight_id,pho.eleveto,year)
    pho_loose=pho[pho.isloose]
    pho_tight=pho[pho.istight]
    pho_ntot=pho.counts
    pho_nloose=pho_loose.counts
    pho_ntight=pho_tight.counts

    fj = Initialize({'pt':tree.array('AK15Puppi_pt'),
                     'eta':tree.array('AK15Puppi_eta'),
                     'phi':tree.array('AK15Puppi_phi'),
                     'mass':tree.array('AK15Puppi_mass'),
                     'id':tree.array('AK15Puppi_jetId'),
                     'probTbcq':tree.array('AK15Puppi_probTbcq'),
                     'probTbqq':tree.array('AK15Puppi_probTbqq'),
                     'probTbc':tree.array('AK15Puppi_probTbc'),
                     'probTbq':tree.array('AK15Puppi_probTbq'),
                     'probWcq':tree.array('AK15Puppi_probWcq'),
                     'probWqq':tree.array('AK15Puppi_probWqq'),
                     'probZbb':tree.array('AK15Puppi_probZbb'),
                     'probZcc':tree.array('AK15Puppi_probZcc'),
                     'probZqq':tree.array('AK15Puppi_probZqq'),
                     'probHbb':tree.array('AK15Puppi_probHbb'),
                     'probHcc':tree.array('AK15Puppi_probHcc'),
                     'probHqqqq':tree.array('AK15Puppi_probHqqqq'),
                     'probQCDbb':tree.array('AK15Puppi_probQCDbb'),
                     'probQCDcc':tree.array('AK15Puppi_probQCDcc'),
                     'probQCDb':tree.array('AK15Puppi_probQCDb'),
                     'probQCDc':tree.array('AK15Puppi_probQCDc'),
                     'probQCDothers':tree.array('AK15Puppi_probQCDothers')})
    fj['probQCD'] = fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
    fj['TvsQCD'] = fj.probTbcq+fj.probTbqq+fj.probTbc+fj.probTbq
    fj['hSvsQCD'] = (fj.probZbb + fj.probHbb) / (fj.probZbb+fj.probHbb+fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq+fj.probHcc+fj.probHqqqq+fj.probQCD)
    fj['VvsQCD'] = (fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq) / (fj.probWcq+fj.probWqq+fj.probZcc+fj.probZqq+fj.probHcc+fj.probHqqqq+fj.probQCD)
    fj['isgood'] = (fj.pt > 200)&(abs(fj.eta)<2.4)&((fj.id&2)!=0)
    fj['isclean'] =~fj.match(pho,1.5)&~fj.match(mu,1.5)&~fj.match(e,1.5)&fj.isgood
    fj_good=fj[fj.isgood]
    fj_clean=fj[fj.isclean]
    fj_ntot=fj.counts
    fj_ngood=fj_good.counts
    fj_nclean=fj_clean.counts

    j = Initialize({'pt':tree.array('Jet_pt'),
                    'eta':tree.array('Jet_eta'),
                    'phi':tree.array('Jet_phi'),
                    'mass':tree.array('Jet_mass'),
                    'id':tree.array('Jet_jetId'),
                    'nhf':tree.array('Jet_neHEF'),
                    'nef':tree.array('Jet_neEmEF'),
                    'chf':tree.array('Jet_chHEF'),
                    'cef':tree.array('Jet_chEmEF')})
    j['isgood'] = (j.pt>25)&(abs(j.eta)<2.4)&((j.id&2)!=0)&(j.nhf<0.8)&(j.nef<0.99)&(j.chf>0.1)&(j.cef<0.99)
    j['isclean'] = ~j.match(e,0.4)&~j.match(mu,0.4)&~j.match(pho,0.4)&j.isgood
    j['isiso'] =  ~(j.match(fj,1.5))&j.isclean
    j_good = j[j.isgood]
    j_clean = j[j.isclean]
    j_iso = j[j.isiso]
    j_ntot=j.counts
    j_ngood=j_good.counts
    j_nclean=j_clean.counts
    j_niso=j_iso.counts

    met = Initialize({'pt':tree.array("MET_pt"),
                      'eta':0,
                      'phi':tree.array("MET_phi"),
                      'mass':0})

    calomet = Initialize({'pt':tree.array("CaloMET_pt"),
                      'eta':0,
                      'phi':tree.array("CaloMET_phi"),
                      'mass':0})

    weight["nlo"] = 1
    if xsec != -1:
        gen = Initialize({'pt':tree.array('GenPart_pt'),
                          'eta':tree.array('GenPart_eta'),
                          'phi':tree.array('GenPart_phi'),
                          'mass':tree.array('GenPart_mass'),
                          'pdgid':tree.array('GenPart_pdgId'),
                          'status':tree.array('GenPart_status'), 
                          'flags':tree.array('GenPart_statusFlags'),
                          'motherid':tree.array('GenPart_genPartIdxMother')})
        
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

        if('TTJets' in dataset): weight["nlo"] = np.sqrt(get_ttbar_weight(genTops[0].pt.sum()) * get_ttbar_weight(genTops[1].pt.sum()))
        elif('WJets' in dataset): weight["nlo"] = get_nlo_weight('w',genWs[0].pt.sum())
        elif('DY' in dataset or 'ZJets' in dataset): weight["nlo"] = get_nlo_weight('z',genZs[0].pt.sum())
        elif('GJets' in dataset): weight["nlo"] = get_nlo_weight('a',genAs[0].pt.sum())




    ###
    #Calculating derivatives
    ###
    ele_pairs = e_loose.distincts()
    diele = ele_pairs.i0+ele_pairs.i1
    mu_pairs = mu_loose.distincts()
    dimu = mu_pairs.i0+mu_pairs.i1
    
    u={}
    u["iszeroL"] = met

    if mu_tight.content.size>0:
        u["isoneM"] = met+mu_tight[mu_tight.pt.argmax()].sum()
    else:
        u["isoneM"] = met

    if e_tight.content.size>0:
        u["isoneE"] = met+e_tight[e_tight.pt.argmax()].sum()
    else:
        u["isoneE"] = met

    if dimu.content.size>0:
        u["istwoM"] = met+dimu[dimu.pt.argmax()].sum()
    else:
        u["istwoM"] = met

    if diele.content.size>0:
        u["istwoE"] = met+diele[diele.pt.argmax()].sum()
    else:
        u["istwoE"] = met

    if pho_tight.content.size>0:
        u["isoneA"] = met+pho_tight[pho_tight.pt.argmax()].sum()
    else:
        u["isoneA"] = met

    weight["trig"] = {}
    weight["trig"]["iszeroL"] = get_met_trig_weight(u["iszeroL"].pt,year)
    weight["trig"]["isoneM"] = get_met_trig_weight(u["isoneM"].pt,year)
    weight["trig"]["istwoM"] = get_met_zmm_trig_weight(u["istwoM"].pt,year)
    weight["trig"]["isoneE"] = 1
    if e_tight.content.size>0:
        weight["trig"]["isoneE"] = get_ele_trig_weight(e_tight[e_tight.pt.argmax()].eta.sum(), e_tight[e_tight.pt.argmax()].pt.sum(), np.full_like(e_tight[e_tight.pt.argmax()].eta.sum(),-99),np.full_like(e_tight[e_tight.pt.argmax()].pt.sum(),-99),year)
    weight["trig"]["istwoE"] = 1
    if diele.content.size>0:
        weight["trig"]["istwoE"] = get_ele_trig_weight(ele_pairs[diele.pt.argmax()].i0.eta.sum(),ele_pairs[diele.pt.argmax()].i0.pt.sum(),ele_pairs[diele.pt.argmax()].i1.eta.sum(),ele_pairs[diele.pt.argmax()].i1.pt.sum(),year)
    weight["trig"]["isoneA"] = 1
    if pho_tight.content.size>0:
        weight["trig"]["isoneA"] = get_pho_trig_weight(pho_tight[pho_tight.pt.argmax()].pt.sum(),year)
    #print(weight["trig"]["iszeroL"],weight["trig"]["isoneM"],weight["trig"]["istwoM"])
    #print(weight["trig"]["isoneE"],weight["trig"]["istwoE"])
    #print(weight["trig"]["isoneA"])

    ###
    #Event selection
    ###
    selections={}
    for k in u.keys():
        selections[k]={} 
        selections[k]["baggy"] = (fj_nclean>0)&(fj_clean.pt.max()>200)&(abs(u[k].delta_phi(j_clean)).min()>0.8)&(u[k].pt>250)&(passMetFilters)
        selections[k]["skinny"] = ~((fj_nclean>0)&(fj_clean.pt.max()>200))&(j_nclean>0)&(j_clean.pt.max()>100)&(abs(u[k].delta_phi(j_clean)).min()>0.5)&(u[k].pt>250)&(passMetFilters)
        selections[k]["inclusive"] = selections[k]["skinny"]|selections[k]["baggy"]
    for s in ["baggy","skinny","inclusive"]:
        selections["iszeroL"][s] = selections["iszeroL"][s]&(e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(passMetTrig)
        selections["isoneM"][s] = selections["isoneM"][s]&(e_nloose==0)&(mu_ntight==1)&(tau_nloose==0)&(pho_nloose==0)&(passMetTrig)
        selections["isoneE"][s] = selections["isoneE"][s]&(e_ntight==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(passSingleEleTrig)
        selections["istwoM"][s] = selections["istwoM"][s]&(e_nloose==0)&(mu_ntight==1)&(mu_nloose==2)&(tau_nloose==0)&(pho_nloose==0)&(passMetTrig)
        if dimu.content.size > 0:
            selections["istwoM"][s] = selections["istwoM"][s]&(e_nloose==0)&(mu_ntight==1)&(mu_nloose==2)&(tau_nloose==0)&(pho_nloose==0)&(dimu[dimu.pt.argmax()].mass.sum()>60)&(dimu[dimu.pt.argmax()].mass.sum()<120)&(passMetTrig)
        selections["istwoE"][s] = selections["istwoE"][s]&(e_ntight==1)&(e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(passSingleEleTrig)
        if diele.content.size > 0:
            selections["istwoE"][s] = selections["istwoE"][s]&(e_ntight==1)&(e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(diele[diele.pt.argmax()].mass.sum()>60)&(diele[diele.pt.argmax()].mass.sum()<120)&(passSingleEleTrig)
        selections["isoneA"][s] = selections["isoneA"][s]&(e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_ntight==1)&(passSinglePhoTrig)
 
    variables = {}
    variables['j1pt'] = j_clean.pt.max()
    variables['j1eta'] = j_clean.eta.max()
    variables['j1phi'] = j_clean.phi.max()
    variables['fj1pt'] = fj_clean.pt.max()
    variables['fj1eta'] = fj_clean.eta.max()
    variables['fj1phi'] = fj_clean.phi.max()
    variables['njets'] = j_nclean
    variables['nfjets'] = fj_nclean
    if fj_clean.content.size > 0:
        variables['fjmass'] = fj_clean[fj_clean.pt.argmax()].mass.sum()
        variables['TvsQCD'] = fj_clean[fj_clean.pt.argmax()].TvsQCD.sum()
        variables['hSvsQCD'] = fj_clean[fj_clean.pt.argmax()].hSvsQCD.sum()
        variables['VvsQCD'] = fj_clean[fj_clean.pt.argmax()].VvsQCD.sum()
    # Filler; does not matter anyway since fj_clean is empty
    #For a proper fix, need to make sure we are not using max on an empty numpy array
    else:
        variables['fjmass'] = -1
        variables['TvsQCD'] = -1
        variables['hSvsQCD'] = -1
        variables['VvsQCD'] = -1
    hout = {}
    for k in hists.keys():
        h = hists[k].copy(content=False)
        i = 0
        if k == 'sumw':
            h.fill(dataset=dataset, sumw=1, weight=sumw)
        else:
            while i < len(selection):
                r = selection[i]
                if k == 'recoil':
                    h.fill(dataset=dataset, region=r, selection='inclusive', recoil=u[r].pt, weight=genw*weight['nlo']*selections[r]['inclusive'])
                    h.fill(dataset=dataset, region=r, selection='baggy', recoil=u[r].pt, weight=genw*weight['nlo']*selections[r]['baggy'])
                    h.fill(dataset=dataset, region=r, selection='skinny', recoil=u[r].pt, weight=genw*weight['nlo']*selections[r]['skinny'])
                elif k == 'CaloMinusPfOverRecoil':
                    h.fill(dataset=dataset, region=r, selection='inclusive', CaloMinusPfOverRecoil= abs(calomet.pt - met.pt) / u[r].pt, weight=genw*weight['nlo']*selections[r]['inclusive'])
                    h.fill(dataset=dataset, region=r, selection='baggy', CaloMinusPfOverRecoil= abs(calomet.pt - met.pt) / u[r].pt, weight=genw*weight['nlo']*selections[r]['baggy'])
                    h.fill(dataset=dataset, region=r, selection='skinny', CaloMinusPfOverRecoil= abs(calomet.pt - met.pt) / u[r].pt, weight=genw*weight['nlo']*selections[r]['skinny'])

                elif k == 'mindphi':
                    h.fill(dataset=dataset, region=r, selection='inclusive', mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=genw*weight['nlo']*selections[r]['inclusive'])
                    h.fill(dataset=dataset, region=r, selection='baggy', mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=genw*weight['nlo']*selections[r]['baggy'])
                    h.fill(dataset=dataset, region=r, selection='skinny', mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=genw*weight['nlo']*selections[r]['skinny'])
                else:
                    h.fill(dataset=dataset, region=r, selection='inclusive', **variables, weight=genw*weight['nlo']*selections[r]['inclusive'])
                    h.fill(dataset=dataset, region=r, selection='baggy', **variables, weight=genw*weight['nlo']*selections[r]['baggy'])
                    h.fill(dataset=dataset, region=r, selection='skinny', **variables, weight=genw*weight['nlo']*selections[r]['skinny'])
                i += 1
        hout[k] = h
    
    return dataset, sumw, tree.numentries, hout
