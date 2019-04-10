import awkward
import uproot, uproot_methods
import numpy as np

e_id = {}
e_id['2016'] = {}
e_id['2017'] = {}
e_id['2018'] = {}

e_id['2016']['loose_id'] = 'Electron_mvaSpring16GP_WP90'
e_id['2016']['tight_id'] = 'Electron_mvaSpring16GP_WP80'
e_id['2016']['iso'] = 'Electron_pfRelIso03_all'
e_id['2017']['loose_id'] = 'Electron_mvaFall17Iso_WP90'
e_id['2017']['tight_id'] = 'Electron_mvaFall17Iso_WP80'
e_id['2017']['iso'] = 'Null'

def isLooseElectron(pt,eta,dxy,dz,iso,loose_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>7)&(abs(eta)<2.4)&(abs(dxy)<0.05)&(abs(dz)<0.2)&(iso<0.4)&(loose_id)
    elif year=='2017':
        mask = (pt>7)&(abs(eta)<2.4)&(abs(dxy)<0.05)&(abs(dz)<0.2)&(loose_id)
    elif year=='2018':
        mask = (pt>10)&(abs(eta)<2.5)&(loose_id)
    return mask

def isTightElectron(pt,eta,dxy,dz,iso,tight_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        return ((pt>30)&(abs(eta)<2.4)&(abs(dxy)<0.05)&(abs(dz)<0.2)&(tight_id)&(iso<0.06))
    elif year=='2017':
        return ((pt>30)&(abs(eta)<2.4)&(abs(dxy)<0.05)&(abs(dz)<0.2)&(tight_id))
    elif year=='2018':
        return ((pt>40)&(abs(eta)<2.5)&(tight_id))
    return mask

mu_id = {}
mu_id['2016'] = {}
mu_id['2017'] = {}
mu_id['2018'] = {}

#mu_id['2016']['loose_id'] = 'Electron_mvaSpring16GP_WP90'
#mu_id['2016']['tight_id'] = 'Electron_mvaSpring16GP_WP80'
mu_id['2016']['iso'] = 'Muon_pfRelIso04_all'
mu_id['2016']['tight_id'] = 'Muon_tightId'
#mu_id['2017']['loose_id'] = 'Electron_mvaFall17Iso_WP90'
#mu_id['2017']['tight_id'] = 'Electron_mvaFall17Iso_WP80'
mu_id['2017']['iso'] = 'Muon_pfRelIso04_all'
mu_id['2017']['tight_id'] = 'Muon_tightId'
mu_id['2018']['iso'] = 'Null'
mu_id['2018']['tight_id'] = 'Muon_tightId'

def isLooseMuon(pt,eta,dxy,dz,iso,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>5)&(abs(eta)<2.4)&(abs(dxy)<0.5)&(abs(dz)<1.0)&(iso<0.4)
    elif year=='2017':
        mask = (pt>5)&(abs(eta)<2.4)&(abs(dxy)<0.5)&(abs(dz)<1.0)&(iso<0.4)
    elif year=='2018':
        mask = (pt>10)&(abs(eta)<2.4)
    return mask

def isTightMuon(pt,eta,dxy,dz,iso,tight_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>5)&(abs(eta)<2.4)&(abs(dxy)<0.5)&(abs(dz)<1.0)&(iso<0.4)
    elif year=='2017':
        mask = (pt>5)&(abs(eta)<2.4)&(abs(dxy)<0.5)&(abs(dz)<1.0)&(iso<0.4)
    elif year=='2018':
        mask = (pt>10)&(abs(eta)<2.4)&(tight_id)
    return mask

tau_id = {}
tau_id['2016'] = {}
tau_id['2017'] = {}
tau_id['2018'] = {}

tau_id['2016']['id'] = 'Tau_idMVAnew'
tau_id['2016']['decayMode'] = 'Tau_idDecayMode'
tau_id['2017']['id'] = 'Null'
tau_id['2017']['decayMode'] = 'Tau_idDecayMode'
tau_id['2018']['id'] = 'Tau_idMVAoldDM2017v2'
tau_id['2018']['decayMode'] = 'Tau_idDecayMode'

def isLooseTau(pt,eta,decayMode,_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>18)&(abs(eta)<2.3)&(decayMode)#&((_id&2)!=0)
    elif year=='2017':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)
    elif year=='2018':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)&((id&2)!=0)
    return mask

pho_id = {}
pho_id['2016'] = {}
pho_id['2017'] = {}
pho_id['2018'] = {}

#pho_id['2016']['loose_id'] = 'Electron_mvaSpring16GP_WP90'
#pho_id['2017']['loose_id'] = 'Electron_mvaFall17Iso_WP90'

def isLoosePhoton(pt,eta,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>15)*(abs(eta)<2.5)
    elif year=='2017':
        mask = (pt>15)*(abs(eta)<2.5)
    return mask

def isTightPhoton(pt,eta,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>15)*(abs(eta)<2.5)
    elif year=='2017':
        mask = (pt>15)*(abs(eta)<2.5)
    return mask
