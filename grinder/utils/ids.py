import awkward
import uproot, uproot_methods
import numpy as np

e_id = {}
e_id['2016'] = {}
e_id['2016']['loose_id'] = 'Electron_cutBased'
e_id['2016']['tight_id'] = 'Electron_cutBased'
e_id['2016']['dxy'] = 'Electron_dxy'
e_id['2016']['dz'] = 'Electron_dz'
e_id['2017'] = e_id['2016']
e_id['2018'] = e_id['2016']

#POG  Tight - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Offline_selection_criteria
#Electron_cutBased  Int_t   cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
def isLooseElectron(pt,eta,dxy,dz,iso,loose_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = ((pt>10)&(abs(eta)<2.5)&(loose_id>=1))
    elif year=='2017':
        mask = ((pt>10)&(abs(eta)<2.5)&(loose_id>=1))
    elif year=='2018':
        mask = ((pt>10)&(abs(eta)<2.5)&(loose_id>=1))
    return mask

def isTightElectron(pt,eta,dxy,dz,iso,tight_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = ((pt>40)&(abs(eta)<2.5)&(tight_id==4)) # Trigger: HLT_Ele27_WPTight_Gsf_v
    elif year=='2017':
        mask = ((pt>40)&(abs(eta)<2.5)&(tight_id==4)) # Trigger: HLT_Ele35_WPTight_Gsf_v
    elif year=='2018':
        mask = ((pt>40)&(abs(eta)<2.5)&(tight_id==4)) # Trigger: HLT_Ele32_WPTight_Gsf_v
    return mask

mu_id = {}
mu_id['2016'] = {}
mu_id['2016']['iso'] = 'Muon_pfRelIso04_all'
mu_id['2016']['tight_id'] = 'Muon_tightId'
mu_id['2016']['med_id'] = 'Muon_mediumId'
mu_id['2016']['dxy'] = 'Muon_dxy'
mu_id['2016']['dz'] = 'Muon_dz'
mu_id['2017'] = mu_id['2016']
mu_id['2018'] = mu_id['2016']


def isLooseMuon(pt,eta,dxy,dz,iso,med_id,year):
    #dxy and dz cuts are missing from med_id; very loose isolation is 0.4
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>10)&(abs(eta)<2.4)&(med_id>0)&(iso<0.25)
    elif year=='2017':
        mask = (pt>10)&(abs(eta)<2.4)&(med_id>0)&(iso<0.25)
    elif year=='2018':
        mask = (pt>10)&(abs(eta)<2.4)&(med_id>0)&(iso<0.25)
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
tau_id['2016']['id'] = 'Tau_idMVAoldDM2017v2'
tau_id['2016']['decayMode'] = 'Tau_idDecayMode'
#tau_id['2016']['clean']= 'Tau_cleanmask' 
tau_id['2017'] = tau_id['2016']
tau_id['2018'] = tau_id['2016']
#bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight

def isLooseTau(pt,eta,decayMode,_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)&((_id&2)==2)
    elif year=='2017':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)&((_id&2)==2)
    elif year=='2018':
        mask = (pt>20)&(abs(eta)<2.3)&(decayMode)&((_id&2)==2)
    return mask

pho_id = {}
pho_id['2016'] = {}
pho_id['2016']['loose_id'] = 'Photon_cutBasedBitmap'
pho_id['2016']['tight_id'] = 'Photon_cutBasedBitmap'
pho_id['2016']['eleveto']  = 'Photon_electronVeto'
pho_id['2016']['phoeta']   = 'Photon_eta'
pho_id['2017'] = pho_id['2016']
pho_id['2018'] = pho_id['2016']

#Photon_cutBasedBitmap  Int_t   cut-based ID bitmap, 2^(0:loose, 1:medium, 2:tight)
#Photon IDs:  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=36
def isLoosePhoton(pt,eta,loose_id,eleveto,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>15)&(abs(eta)<2.5)&((loose_id&2)>=1)&(eleveto)
    elif year=='2017':
        mask = (pt>15)&(abs(eta)<2.5)&((loose_id&2)>=1)&(eleveto)
    elif year=='2018':
        mask = (pt>15)&(abs(eta)<2.5)&((loose_id&2)>=1)&(eleveto)
    return mask

def isTightPhoton(pt,eta,tight_id,eleveto,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>215)&(abs(eta)<1.4442)&((tight_id&2)==2)&(eleveto) # Trigger threshold is at 175
    elif year=='2017':
        mask = (pt>215)&(abs(eta)<1.4442)&((tight_id&2)==2)&(eleveto) # Trigger threshold is at 200
    elif year=='2018':
        mask = (pt>215)&(abs(eta)<1.4442)&((tight_id&2)==2)&(eleveto) # Trigger threshold is at 200
    return mask

fj_id = {}
fj_id['2016'] = {}
fj_id['2016']['id'] = 'AK15Puppi_jetId'
fj_id['2017'] = fj_id['2016']
fj_id['2018'] = fj_id['2016']

def isGoodFatJet(pt,eta, jet_id):
    mask = (pt > 160)&(abs(eta)<2.4)&((jet_id&2)==2)
    return mask

#Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto
#POG use tight jetID as a standart JetID 
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
    mask = (pt>25) & (abs(eta)<2.4) & ((jet_id&2)==2) & (nhf<0.8) & (nef<0.99) & (chf>0.1) & (cef<0.99)
    return mask