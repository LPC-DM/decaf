#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from Builder import Initialize
from fnal_column_analysis_tools import hist
from analysis.triggers import met_trigger_paths, singleele_trigger_paths, singlepho_trigger_paths
from analysis.corrections import get_ttbar_weight, get_nlo_weight
from analysis.ids import e_id, isLooseElectron, isTightElectron
from analysis.ids import mu_id, isLooseMuon
from analysis.ids import tau_id, isLooseTau
from analysis.ids import pho_id, isLoosePhoton

hists = {
    'sumw': hist.Hist("sumw", hist.Cat("dataset", "Primary dataset"), hist.Bin("sumw", "Weight value", [0.])),
    'recoil': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("recoil","Hadronic Recoil",[250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
    'mindphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("mindphi","Min dPhi(MET,AK4s)",15,0,6.28)),
    'j1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("j1pt","AK4 Leading Jet Pt",50,30,500)),
    'fj1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("fj1pt","AK15 Leading Jet Pt",50,200,700)),
    'njets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("njets","AK4 Number of Jets",6,0,5)),
    'nfjets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("nfjets","AK15 Number of Jets",4,0,3)),
    'fjmass': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("fjmass","AK15 Jet Mass",50,20,250)),
    'TvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("TvsQCD","TvsQCD",15,0,1)),
    'WvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("WvsQCD","WvsQCD",15,0,1)),
    'ZvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("ZvsQCD","ZvsQCD",15,0,1)),
    'VvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("VvsQCD","VvsQCD",15,0,1)),
    'ZHbbvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("ZHbbvsQCD","ZHbbvsQCD",15,0,1)),
    'ZHccvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("ZHccvsQCD","ZHccvsQCD",15,0,1)),
    'WcqvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("WcqvsQCD","WcqvsQCD",15,0,1)),
    'WqqvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("WqqvsQCD","WqqvsQCD",15,0,1)),
    'ZbbvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("ZbbvsQCD","ZbbvsQCD",15,0,1)),
    'ZccvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("ZccvsQCD","ZccvsQCD",15,0,1)),
    'ZqqvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Cat("region", "Region"), hist.Bin("ZqqvsQCD","ZqqvsQCD",15,0,1))
}


samples = {
    "iszeroL":('ZJets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','MET'),
    "isoneM":('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','MET'),
    "isoneE":('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','SingleElectron'),
    "istwoM":('WJets','DY','TT','ST','WW','WZ','ZZ','HToBB','MET'),
    "istwoE":('WJets','DY','TT','ST','WW','WZ','ZZ','HToBB','SingleElectron'),
    "isoneA":('GJets','QCD','SinglePhoton')
}

def analysis(selection, year, xsec, dataset, file):
    tree = uproot.open(file)["Events"]
    genw = 1
    sumw = 1

    ###
    #For MC, retrieve the LHE weights, to take into account NLO destructive interference, and their sum
    ###
    if xsec != -1:
        genw = tree.array("genWeight")
        run_tree = uproot.open(file)["Runs"]
        sumw = run_tree.array("genEventSumw")[0]
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
    passMetTrig = np.prod([met_trigger[key] for key in met_trigger], axis=0)

    singleele_trigger = {}
    for path in singleele_trigger_paths[year]:
        try:
            singleele_trigger[path] = tree.array(path)
        except KeyError:
            #print("No trigger bit in file for path ",path)
            pass
    passSingleEleTrig = np.prod([singleele_trigger[key] for key in singleele_trigger], axis=0)

    singlepho_trigger = {}
    for path in singlepho_trigger_paths[year]:
        try:
            singlepho_trigger[path] = tree.array(path)
        except KeyError:
            #print("No trigger bit in file for path ",path)
            pass
    passSinglePhoTrig = np.prod([singlepho_trigger[key] for key in singlepho_trigger], axis=0)

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
    mu_loose=mu[mu.isloose]
    mu_ntot = mu.counts
    mu_nloose = mu_loose.counts

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
    pho['isloose']=isLoosePhoton(pho.pt,pho.eta,year)
    pho_loose=pho[pho.isloose]
    pho_ntot=pho.counts
    pho_nloose=pho_loose.counts

    fj = Initialize({'pt':tree.array('AK15Puppi_pt'),
                     'eta':tree.array('AK15Puppi_eta'),
                     'phi':tree.array('AK15Puppi_phi'),
                     'mass':tree.array('AK15Puppi_mass'),
                     'jetId':tree.array('AK15Puppi_jetId'),
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
    fj['TvsQCD'] = (fj.probTbcq+fj.probTbqq)/fj.probQCD
    fj['WvsQCD'] = (fj.probWcq+fj.probWqq)/fj.probQCD
    fj['ZvsQCD'] = (fj.probZbb+fj.probZcc+fj.probZqq)/fj.probQCD
    fj['VvsQCD'] = (fj.probZbb+fj.probZcc+fj.probZqq+fj.probWcq+fj.probWqq)/fj.probQCD
    fj['ZHbbvsQCD'] = (fj.probZbb+fj.probHbb)/fj.probQCD
    fj['ZHccvsQCD'] = (fj.probZcc+fj.probHcc)/fj.probQCD
    fj['WcqvsQCD'] = fj.probWcq/fj.probQCD
    fj['WqqvsQCD'] = fj.probWqq/fj.probQCD
    fj['ZbbvsQCD'] = fj.probZbb/fj.probQCD
    fj['ZccvsQCD'] = fj.probZcc/fj.probQCD
    fj['ZqqvsQCD'] = fj.probZqq/fj.probQCD
    fj['isgood'] = (fj.pt > 200)&(abs(fj.eta)<2.4)&(fj.jetId > 0)
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
                    'id':tree.array('Jet_jetId')})
    j['isgood'] = (j.pt>25)&(abs(j.eta)<4.5)&((j.id&2)!=0)
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


    weight = {}
    weight["nlo"] = 1
    if('TTJets' in dataset): weight["nlo"] = np.sqrt(get_ttbar_weight(genTops[0].pt.sum()) * get_ttbar_weight(genTops[1].pt.sum()))
    elif('WJets' in dataset): weight["nlo"] = get_nlo_weight('w',genWs[0].pt.sum())
    elif('DY' in dataset or 'ZJets' in dataset): weight["nlo"] = get_nlo_weight('z',genZs[0].pt.sum())
    elif('GJets' in dataset): weight["nlo"] = get_nlo_weight('a',genAs[0].pt.sum())


    ###
    #Calculating derivatives
    ###
    diele = e_loose.distincts().i0+e_loose.distincts().i1
    dimu = mu_loose.distincts().i0+mu_loose.distincts().i1
    
    u={}
    u["iszeroL"] = met

    if mu_loose.content.size>0:
        u["isoneM"] = met+mu_loose[mu_loose.pt.argmax()].sum()
    else:
        u["isoneM"] = met

    if e_loose.content.size>0:
        u["isoneE"] = met+e_loose[e_loose.pt.argmax()].sum()
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

    if pho_loose.content.size>0:
        u["isoneA"] = met+pho_loose[pho_loose.pt.argmax()].sum()
    else:
        u["isoneA"] = met

    ###
    #Event selection
    ###
    skinny={}
    loose={}
    inclusive={}
    for k in u.keys():
#        if selection in k:
        skinny[k] = (j_nclean>0)&(j_clean.pt.max()>100)&(abs(u[k].delta_phi(j_clean)).min()>0.5)
        loose[k] = (fj_nclean>0)&(fj_clean.pt.max()>200)&(abs(u[k].delta_phi(j_clean)).min()>0.8)
        inclusive[k] = skinny[k]|loose[k]
 
    selections={}
    selections["iszeroL"] = (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(passMetTrig)
    selections["isoneM"] = (e_nloose==0)&(mu_nloose==1)&(tau_nloose==0)&(pho_nloose==0)&(passMetTrig)
    selections["isoneE"] = (e_nloose==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(passSingleEleTrig)
    if dimu.content.size > 0:
        selections["istwoM"] = (e_nloose==0)&(mu_nloose==2)&(tau_nloose==0)&(pho_nloose==0)&(dimu[dimu.pt.argmax()].mass.sum()>60)&(dimu[dimu.pt.argmax()].mass.sum()<120)&(passMetTrig)
    else:
        selections["istwoM"] = (e_nloose==0)&(mu_nloose==2)&(tau_nloose==0)&(pho_nloose==0)&(passMetTrig)
    if diele.content.size > 0:
        selections["istwoE"] = (e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(diele[diele.pt.argmax()].mass.sum()>60)&(diele[diele.pt.argmax()].mass.sum()<120)&(passSingleEleTrig)
    else:
        selections["istwoE"] = (e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(passSingleEleTrig)
    selections["isoneA"] = (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==1)&(passSinglePhoTrig)

    for k in u.keys():
#        if selection in k:
        skinny[k] = skinny[k]&selections[k]&(u[k].pt>200)
        loose[k] = loose[k]&selections[k]&(u[k].pt>200)
        inclusive[k] = inclusive[k]&selections[k]&(u[k].pt>200)

    variables = {}
    variables['j1pt'] = j_clean.pt.max()
    variables['fj1pt'] = fj_clean.pt.max()
    variables['njets'] = j_nclean
    variables['nfjets'] = fj_nclean
    if fj_clean.content.size > 0:
        variables['fjmass'] = fj_clean[fj_clean.pt.argmax()].mass.sum()
        variables['TvsQCD'] = fj_clean[fj_clean.pt.argmax()].TvsQCD.sum()
        variables['WvsQCD'] = fj_clean[fj_clean.pt.argmax()].WvsQCD.sum()
        variables['ZvsQCD'] = fj_clean[fj_clean.pt.argmax()].ZvsQCD.sum()
        variables['VvsQCD'] = fj_clean[fj_clean.pt.argmax()].VvsQCD.sum()
        variables['ZHbbvsQCD'] = fj_clean[fj_clean.pt.argmax()].ZHbbvsQCD.sum()
        variables['ZHccvsQCD'] = fj_clean[fj_clean.pt.argmax()].ZHccvsQCD.sum()
        variables['WcqvsQCD'] = fj_clean[fj_clean.pt.argmax()].WcqvsQCD.sum()
        variables['WqqvsQCD'] = fj_clean[fj_clean.pt.argmax()].WqqvsQCD.sum()
        variables['ZbbvsQCD'] = fj_clean[fj_clean.pt.argmax()].ZbbvsQCD.sum()
        variables['ZccvsQCD'] = fj_clean[fj_clean.pt.argmax()].ZccvsQCD.sum()
        variables['ZqqvsQCD'] = fj_clean[fj_clean.pt.argmax()].ZqqvsQCD.sum()
    # Filler; does not matter anyway since fj_clean is empty
    #For a proper fix, need to make sure we are not using max on an empty numpy array
    else:
        variables['fjmass'] = -1
        variables['TvsQCD'] = -1
        variables['WvsQCD'] = -1
        variables['ZvsQCD'] = -1
        variables['VvsQCD'] = -1
        variables['ZHbbvsQCD'] = -1
        variables['ZHccvsQCD'] = -1
        variables['WcqvsQCD'] = -1
        variables['WqqvsQCD'] = -1
        variables['ZbbvsQCD'] = -1
        variables['ZccvsQCD'] = -1
        variables['ZqqvsQCD'] = -1
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
                    h.fill(dataset=dataset, region=r, recoil=u[r].pt, weight=genw*weight['nlo']*inclusive[r])
                elif k == 'mindphi':
                    h.fill(dataset=dataset, region=r, mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=genw*weight['nlo']*inclusive[r])
                else:
                    h.fill(dataset=dataset, region=r, **variables, weight=genw*weight['nlo']*inclusive[r])
                i += 1
        hout[k] = h
    
    return dataset, sumw, tree.numentries, hout
