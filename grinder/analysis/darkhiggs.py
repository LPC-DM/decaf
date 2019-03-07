#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from Builder import Initialize
from fnal_column_analysis_tools import hist

hists = {
    'recoil': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("recoil","Hadronic Recoil",[250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])),
    'mindphi': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("mindphi","Min dPhi(MET,AK4s)",15,0,6.28)),
    'j1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("j1pt","AK4 Leading Jet Pt",50,30,500)),
    'fj1pt': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("fj1pt","AK15 Leading Jet Pt",50,200,700)),
    'njets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("njets","AK4 Number of Jets",6,0,5)),
    'nfjets': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("nfjets","AK15 Number of Jets",4,0,3)),
    'fjmass': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("fjmass","AK15 Jet Mass",50,20,250)),
    'TvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("TvsQCD","TvsQCD",15,0,1)),
    'WvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("WvsQCD","WvsQCD",15,0,1)),
    'ZvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("ZvsQCD","ZvsQCD",15,0,1)),
    'VvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("VvsQCD","VvsQCD",15,0,1)),
    'ZHbbvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("ZHbbvsQCD","ZHbbvsQCD",15,0,1)),
    'ZHccvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("ZHccvsQCD","ZHccvsQCD",15,0,1)),
    'WcqvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("WcqvsQCD","WcqvsQCD",15,0,1)),
    'WqqvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("WqqvsQCD","WqqvsQCD",15,0,1)),
    'ZbbvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("ZbbvsQCD","ZbbvsQCD",15,0,1)),
    'ZccvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("ZccvsQCD","ZccvsQCD",15,0,1)),
    'ZqqvsQCD': hist.Hist("Events", hist.Cat("dataset", "Primary dataset"), hist.Bin("ZqqvsQCD","ZqqvsQCD",15,0,1))
}


samples = {
    "iszeroL":('ZJets','WJets','DY','TT_TuneCUETP8M2T4','ST_t-channel','ST_tW','WW_TuneCUETP8M1','WZ_TuneCUETP8M1','ZZ_TuneCUETP8M1','QCD','VH_HToBB','WminusH','WplusH','ttHTobb','GluGluHToBB','VBFHToBB','MET'),
    "isoneM":('WJets','DYJetsToLL','TT_TuneCUETP8M2T4','ST_t-channel','ST_tW','WW_TuneCUETP8M1','WZ_TuneCUETP8M1','ZZ_TuneCUETP8M1','QCD','VH_HToBB','WminusH','WplusH','ttHTobb','GluGluHToBB','VBFHToBB','MET'),
    "isoneE":('WJets','DYJetsToLL','TT_TuneCUETP8M2T4','ST_t-channel','ST_tW','WW_TuneCUETP8M1','WZ_TuneCUETP8M1','ZZ_TuneCUETP8M1','QCD','VH_HToBB','WminusH','WplusH','ttHTobb','GluGluHToBB','VBFHToBB','SingleElectron'),
    "istwoM":('WJets','DYJetsToLL','TT_TuneCUETP8M2T4','ST_t-channel','ST_tW','WW_TuneCUETP8M1','WZ_TuneCUETP8M1','ZZ_TuneCUETP8M1','VH_HToBB','WminusH','WplusH','ttHTobb','MET'),
    "istwoE":('WJets','DYJetsToLL','TT_TuneCUETP8M2T4','ST_t-channel','ST_tW','WW_TuneCUETP8M1','WZ_TuneCUETP8M1','ZZ_TuneCUETP8M1','VH_HToBB','WminusH','WplusH','ttHTobb','SingleElectron'),
    "isoneA":('GJets','QCD','SinglePhoton')
}

def analysis(selection, year, xsec, dataset, file):
    tree = uproot.open(file)["Events"]
    genw = 1
    sumw = 1

    if xsec != -1:
        genw = tree.array("genWeight")
        run_tree = uproot.open(file)["Runs"]
        sumw = run_tree.array("genEventSumw")[0]

    e = Initialize({'pt':tree.array("Electron_pt"),
                    'eta':tree.array("Electron_eta"),
                    'phi':tree.array("Electron_phi"),
                    'mass':tree.array("Electron_mass"),
                    'dxy':tree.array('Electron_dxy'),
                    'dz':tree.array('Electron_dz')})

    if '2016' in year:
        e['loose_id'] = tree.array('Electron_mvaSpring16GP_WP90')
        e['tight_id']  = tree.array('Electron_mvaSpring16GP_WP80')
        e['iso'] = tree.array('Electron_pfRelIso03_all')
        e['isloose'] = (e.pt>7)&(abs(e.eta)<2.4)&(abs(e.dxy)<0.05)&(abs(e.dz)<0.2)&(e.iso<0.4)&(e.loose_id)
        e['istight'] = (e.pt>30)&(abs(e.eta)<2.4)&(abs(e.dxy)<0.05)&(abs(e.dz)<0.2)&(e.tight_id)&(e.iso<0.06)

    elif '2017' in year:
        e['loose_id'] = tree.array('Electron_mvaFall17Iso_WP90')
        e['tight_id'] = tree.array('Electron_mvaFall17Iso_WP80')
        e['isloose'] = (e.pt>7)&(abs(e.eta)<2.4)&(abs(e.dxy)<0.05)&(abs(e.dz)<0.2)&(e.loose_id)                    
        e['istight'] = (e.pt>30)&(abs(e.eta)<2.4)&(abs(e.dxy)<0.05)&(abs(e.dz)<0.2)&(e.tight_id)

    e_loose = e[e.isloose]
    e_tight = e[e.istight]
    e_ntot = e.counts
    e_nloose = e_loose.counts
    e_ntight = e[e.istight].counts

    mu = Initialize({'pt':tree.array("Muon_pt"),
                     'eta':tree.array("Muon_eta"),
                     'phi':tree.array("Muon_phi"),
                     'mass':tree.array("Muon_mass"),
                     'iso':tree.array('Muon_pfRelIso04_all'),
                     'dxy':tree.array('Muon_dxy'),
                     'dz':tree.array('Muon_dz')})
    mu['isloose']=(mu.counts>0)&(mu.pt>5)&(abs(mu.eta)<2.4)&(abs(mu.dxy)<0.5)&(abs(mu.dz)<1.0)&(mu.iso<0.4)
    m_loose=mu[mu.isloose]
    mu_ntot = mu.counts
    mu_nloose = m_loose.counts
    tau = Initialize({'pt':tree.array('Tau_pt'),
                      'eta':tree.array('Tau_eta'),
                      'phi':tree.array('Tau_phi'),
                      'mass':tree.array('Tau_mass'),
                      'decayMode':tree.array('Tau_idDecayMode')})
    if '2016' in year:
        tau['id'] = tree.array('Tau_idMVAnew')
        tau['isloose']=(tau.counts>0)&(tau.pt>18)&(abs(tau.eta)<2.3)&(tau.decayMode)&((tau.id&2)!=0)
    elif '2017' in year:
        #Need to find equivalent for 2017
        tau['isloose']=(tau.counts>0)&(tau.pt>18)&(abs(tau.eta)<2.3)&(tau.decayMode)
    tau_loose=tau[tau.isloose]
    tau_ntot=tau.counts
    tau_nloose=tau_loose.counts

    pho = Initialize({'pt':tree.array('Photon_pt'),
                      'eta':tree.array('Photon_eta'),
                      'phi':tree.array('Photon_phi'),
                      'mass':tree.array('Photon_mass')})
    pho['isloose'] = (pho.counts>0)&(pho.pt>15)*(abs(pho.eta)<2.5)
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

    diele = e_loose.distincts().i0+e_loose.distincts().i1
    dimu = m_loose.distincts().i0+m_loose.distincts().i1

    u={}
    u["iszeroL"] = met
    u["isoneM"] = met+m_loose
    u["isoneE"] = met+e_loose
    u["istwoM"] = met+dimu
    u["istwoE"] = met+diele
    u["isoneA"] = met+pho_loose

    skinny={}
    loose={}
    inclusive={}
    for k in u.keys():
        if selection in k:
            skinny[k] = (j_nclean>0)&(j_clean.pt.max()>100)&(abs(u[k].delta_phi(j_clean)).min()>0.5)
            loose[k] = (fj_nclean>0)&(fj_clean.pt.max()>200)&(abs(u[k].delta_phi(j_clean)).min()>0.8)
            inclusive[k] = skinny[k]|loose[k]
 
    selections={}
    selections["iszeroL"] = (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)
    selections["isoneM"] = (e_nloose==0)&(mu_nloose==1)&(tau_nloose==0)&(pho_nloose==0)
    selections["isoneE"] = (e_nloose==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)
    selections["istwoM"] = (e_nloose==0)&(mu_nloose==2)&(tau_nloose==0)&(pho_nloose==0)&(dimu.mass>60)&(dimu.mass<120)
    selections["istwoE"] = (e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)&(diele.mass>60)&(diele.mass<120)
    selections["isoneA"] = (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==1)

    for k in u.keys():
        if selection in k:
            skinny[k] = skinny[k]&selections[k]&(u[k].pt>200)
            loose[k] = loose[k]&selections[k]&(u[k].pt>200)
            inclusive[k] = inclusive[k]&selections[k]&(u[k].pt>200)

    variables = {}
    variables['j1pt'] = j_clean.pt.max()
    variables['fj1pt'] = fj_clean.pt.max()
    variables['njets'] = j_nclean
    variables['nfjets'] = fj_nclean
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
    

    hout = {}
    for k in hists.keys():
        h = hists[k].copy(content=False)
        for r in u.keys():
            if selection in r:
                if k == 'recoil':
                    h.fill(dataset=dataset, recoil=u[r].pt, weight=genw*inclusive[r])
                elif k == 'mindphi':
                    h.fill(dataset=dataset, mindphi=abs(u[r].delta_phi(j_clean)).min(), weight=genw*inclusive[r])
                else:
                    h.fill(dataset=dataset, **variables, weight=genw*inclusive[r])
                hout[k+'_'+r] = h
            else:
                continue
    
    return dataset, sumw, tree.numentries, hout
