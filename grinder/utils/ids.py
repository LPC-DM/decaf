import awkward
import uproot, uproot_methods
import numpy as np

e_id = {}
e_id['2016'] = {}
e_id['2016']['loose_id'] = 'Electron_mvaFall17V2Iso_WP90'
e_id['2016']['tight_id'] = 'Electron_mvaFall17V2Iso_WP80'
e_id['2016']['iso'] = 'Null'
e_id['2016']['dxy'] = 'Electron_dxy'
e_id['2016']['dz'] = 'Electron_dz'

e_id['2017'] = e_id['2016']
e_id['2017']['loose_id'] = 'Electron_mvaFall17V2Iso_WP90'
e_id['2017']['tight_id'] = 'Electron_mvaFall17V2Iso_WP80'
e_id['2017']['iso'] = 'Null'

e_id['2018'] = e_id['2016']
e_id['2018']['loose_id'] = 'Electron_mvaFall17V2Iso_WP90'
e_id['2018']['tight_id'] = 'Electron_mvaFall17V2Iso_WP80'
e_id['2018']['iso'] = 'Null'

def isLooseElectron(pt,eta,dxy,dz,iso,loose_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = ((pt>10)&(abs(eta)<2.5)&(loose_id))
    elif year=='2017':
        mask = ((pt>10)&(abs(eta)<2.5)&(loose_id))
    elif year=='2018':
        mask = ((pt>10)&(abs(eta)<2.5)&(loose_id))
    return mask

def isTightElectron(pt,eta,dxy,dz,iso,tight_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = ((pt>30)&(abs(eta)<2.5)&(tight_id)) # Trigger: HLT_Ele27_WPTight_Gsf_v
    elif year=='2017':
        mask = ((pt>38)&(abs(eta)<2.5)&(tight_id)) # Trigger: HLT_Ele35_WPTight_Gsf_v
    elif year=='2018':
        mask = ((pt>35)&(abs(eta)<2.5)&(tight_id)) # Trigger: HLT_Ele32_WPTight_Gsf_v
    return mask

mu_id = {}
mu_id['2016'] = {}
mu_id['2016']['iso'] = 'Muon_pfRelIso04_all'
mu_id['2016']['tight_id'] = 'Muon_tightId'
mu_id['2016']['dxy'] = 'Muon_dxy'
mu_id['2016']['dz'] = 'Muon_dz'

mu_id['2017'] = mu_id['2016']
mu_id['2017']['iso'] = 'Muon_pfRelIso04_all'
mu_id['2017']['tight_id'] = 'Muon_tightId'

mu_id['2018'] = mu_id['2016']
mu_id['2018']['iso'] = 'Muon_pfRelIso04_all'
mu_id['2018']['tight_id'] = 'Muon_tightId'

def isLooseMuon(pt,eta,dxy,dz,iso,year):
    #dxy and dz cuts are missing from loose_id; very loose isolation is 0.4
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>20)&(abs(eta)<2.4)&(abs(dxy)<0.5)&(abs(dz)<1.0)&(iso<0.4)
    elif year=='2017':
        mask = (pt>20)&(abs(eta)<2.4)&(abs(dxy)<0.5)&(abs(dz)<1.0)&(iso<0.4)
    elif year=='2018':
        mask = (pt>20)&(abs(eta)<2.4)&(abs(dxy)<0.5)&(abs(dz)<1.0)&(iso<0.4)
    return mask

def isTightMuon(pt,eta,dxy,dz,iso,tight_id,year):
    #dxy and dz cuts are baked on tight_id; tight isolation is 0.15
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>20)&(abs(eta)<2.4)&(tight_id)&(iso<0.15)
    elif year=='2017':
        mask = (pt>20)&(abs(eta)<2.4)&(tight_id)&(iso<0.15)
    elif year=='2018':
        mask = (pt>20)&(abs(eta)<2.4)&(tight_id)&(iso<0.15)
    return mask

tau_id = {}
tau_id['2016'] = {}
tau_id['2016']['id'] = 'Tau_idMVAnew'
tau_id['2016']['decayMode'] = 'Tau_idDecayMode'

tau_id['2017'] = tau_id['2016']
tau_id['2017']['id'] = 'Null'
tau_id['2017']['decayMode'] = 'Tau_idDecayMode'

tau_id['2018'] = tau_id['2016']
tau_id['2018']['id'] = 'Tau_idMVAoldDM2017v2'
tau_id['2018']['decayMode'] = 'Tau_idDecayMode'

def isLooseTau(pt,eta,decayMode,_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)&((_id&2)!=0)
    elif year=='2017':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)&((_id&2)!=0)
    elif year=='2018':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)&((_id&2)!=0)
    return mask

pho_id = {}
pho_id['2016'] = {}
pho_id['2016']['loose_id'] = 'Photon_mvaID_WP90'
pho_id['2016']['tight_id'] = 'Photon_mvaID_WP80'
pho_id['2016']['eleveto']  = 'Photon_electronVeto'

pho_id['2017'] = pho_id['2016']
pho_id['2017']['loose_id'] = 'Photon_mvaID_WP90'
pho_id['2017']['tight_id'] = 'Photon_mvaID_WP80'
pho_id['2017']['eleveto']  = 'Photon_electronVeto'

pho_id['2018'] = pho_id['2016']
pho_id['2018']['loose_id'] = 'Photon_mvaID_WP90'
pho_id['2018']['tight_id'] = 'Photon_mvaID_WP80'
pho_id['2018']['eleveto']  = 'Photon_electronVeto'


def isLoosePhoton(pt,eta,loose_id,eleveto,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>20)&(abs(eta)<2.5)&(loose_id)&(eleveto)
    elif year=='2017':
        mask = (pt>20)&(abs(eta)<2.5)&(loose_id)&(eleveto)
    elif year=='2018':
        mask = (pt>20)&(abs(eta)<2.5)&(loose_id)&(eleveto)
    return mask

def isTightPhoton(pt,eta,tight_id,eleveto,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>185)&(abs(eta)<2.5)&(tight_id)&(eleveto) # Trigger threshold is at 175
    elif year=='2017':
        mask = (pt>210)&(abs(eta)<2.5)&(tight_id)&(eleveto) # Trigger threshold is at 200
    elif year=='2018':
        mask = (pt>210)&(abs(eta)<2.5)&(tight_id)&(eleveto) # Trigger threshold is at 200
    return mask

fj_id = {}
fj_id['2016'] = {}
fj_id['2016']['id'] = 'AK15Puppi_jetId'
fj_id['2017'] = fj_id['2016']
fj_id['2018'] = fj_id['2016']

def isGoodFatJet(pt,eta, jet_id):
    mask = (pt > 160)&(abs(eta)<2.4)&((jet_id&2)!=0)
    return mask

j_id = {}
j_id['2016'] = {}
j_id['2016']['id'] = 'Jet_jetId'
j_id['2016']['nhf'] = 'Jet_neHEF'
j_id['2016']['nef'] = 'Jet_neEmEF'
j_id['2016']['chf'] = 'Jet_chHEF'
j_id['2016']['cef'] = 'Jet_chEmEF'
j_id['2017'] =	j_id['2016']
j_id['2018'] =	j_id['2016']

def isGoodJet(pt, eta, jet_id, nhf, nef, chf, cef):
    mask = (pt>25) & (abs(eta)<2.4) & ((jet_id&2)!=0) & (nhf<0.8) & (nef<0.99) & (chf>0.1) & (cef<0.99)
    return mask
