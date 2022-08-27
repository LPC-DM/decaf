import awkward
import uproot, uproot_methods
import numpy as np
from coffea.util import save

#Electron_cutBased  Int_t   cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
def isLooseElectron(pt,eta,dxy,dz,veto_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = ((pt>10)&(abs(eta)<1.4442)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(veto_id>=1)) | ((pt>10)&(abs(eta)>1.5660)&(abs(eta)<2.5)&(abs(dxy)<0.1)&(abs(dz)<0.2)&(veto_id>=1))
    elif year=='2017':
        mask = ((pt>10)&(abs(eta)<1.4442)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(veto_id>=1)) | ((pt>10)&(abs(eta)>1.5660)&(abs(eta)<2.5)&(abs(dxy)<0.1)&(abs(dz)<0.2)&(veto_id>=1))
    elif year=='2018':
        mask = ((pt>10)&(abs(eta)<1.4442)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(veto_id>=1)) | ((pt>10)&(abs(eta)>1.5660)&(abs(eta)<2.5)&(abs(dxy)<0.1)&(abs(dz)<0.2)&(veto_id>=1))
    return mask

#2017/18 pT thresholds adjusted to match monojet, using dedicated ID SFs
def isTightElectron(pt,eta,dxy,dz,tight_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = ((pt>29)&(abs(eta)<1.4442)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(tight_id==4)) | ((pt>29)&(abs(eta)>1.5660)&(abs(eta)<2.5)&(abs(dxy)<0.1)&(abs(dz)<0.2)&(tight_id==4)) # Trigger: HLT_Ele27_WPTight_Gsf_v
    elif year=='2017':
        mask = ((pt>40)&(abs(eta)<1.4442)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(tight_id==4)) | ((pt>40)&(abs(eta)>1.5660)&(abs(eta)<2.5)&(abs(dxy)<0.1)&(abs(dz)<0.2)&(tight_id==4)) # Trigger: HLT_Ele35_WPTight_Gsf_v
    elif year=='2018':
        mask = ((pt>40)&(abs(eta)<1.4442)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(tight_id==4)) | ((pt>40)&(abs(eta)>1.5660)&(abs(eta)<2.5)&(abs(dxy)<0.1)&(abs(dz)<0.2)&(tight_id==4)) # Trigger: HLT_Ele32_WPTight_Gsf_v
    return mask


def isLooseMuon(pt,eta,iso,loose_id,year):
    #dxy and dz cuts are missing from med_id; loose isolation is 0.25
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>20)&(abs(eta)<2.4)&(loose_id>0)&(iso<0.25)
    elif year=='2017':
        mask = (pt>20)&(abs(eta)<2.4)&(loose_id>0)&(iso<0.25)
    elif year=='2018':
        mask = (pt>15)&(abs(eta)<2.4)&(loose_id>0)&(iso<0.25)
    return mask

def isTightMuon(pt,eta,iso,tight_id,year):
    #dxy and dz cuts are baked on tight_id; tight isolation is 0.15
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>30)&(abs(eta)<2.4)&(tight_id)&(iso<0.15)
    elif year=='2017':
        mask = (pt>30)&(abs(eta)<2.4)&(tight_id)&(iso<0.15)
    elif year=='2018':
        mask = (pt>30)&(abs(eta)<2.4)&(tight_id)&(iso<0.15)
    return mask

#bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight
def isLooseTau(pt,eta,decayMode,_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>18)&(abs(eta)<2.3)&(decayMode)&((_id&2)==2)
    elif year=='2017':
        mask = (pt>18)&(abs(eta)<2.3)&(decayMode)&((_id&2)==2)
    elif year=='2018':
        mask = (pt>18)&(abs(eta)<2.3)&(decayMode)&((_id&2)==2)
    return mask

#Photon_cutBased Int_t "cut-based spring16-V2p2 ID (0:fail, 1:loose, 2:medium, 3:tight" for 2016 NanoAOD
#Photon_cutBasedBitmap  Int_t   cut-based ID bitmap, 2^(0:loose, 1:medium, 2:tight)
#Photon IDs:  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=36
def isLoosePhoton(pt,eta,loose_id,year):
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>15)&~((abs(eta)>1.4442)&(abs(eta)<1.5660))&(abs(eta)<2.5)&(loose_id>=1)
    elif year=='2017':
        mask = (pt>15)&~((abs(eta)>1.4442)&(abs(eta)<1.5660))&(abs(eta)<2.5)&((loose_id&1)==1)
    elif year=='2018':
        mask = (pt>15)&~((abs(eta)>1.4442)&(abs(eta)<1.5660))&(abs(eta)<2.5)&((loose_id&1)==1)
    return mask

def isTightPhoton(pt,tight_id,year):
    #isScEtaEB is used (barrel only), so no eta requirement
    #2017/18 pT requirement adjusted to match monojet, using dedicated ID SFs
    #Tight photon use medium ID, as in monojet
    mask = ~(pt==np.nan)#just a complicated way to initialize a jagged array with the needed shape to True
    if year=='2016':
        mask = (pt>200)&(tight_id>=2) # Trigger threshold is at 175
    elif year=='2017':
        mask = (pt>230)&((tight_id&2)==2) # Trigger threshold is at 200
    elif year=='2018':
        mask = (pt>230)&((tight_id&2)==2) # Trigger threshold is at 200
    return mask

def isGoodFatJet(pt, eta, jet_id):
    mask = (pt > 160) & (abs(eta)<2.4) & ((jet_id&2)==2) 
    return mask

#Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto
#POG use tight jetID as a standart JetID 

def isGoodJet(pt, eta, jet_id, pu_id, nhf, chf):
    mask = (pt>30) & (abs(eta)<2.4) & ((jet_id&2)==2) & (nhf<0.8) & (chf>0.1)# & (nef<0.99) & (cef<0.99)
    mask = ((pt>=50)&mask) | ((pt<50)&mask&((pu_id&1)==1)) #https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID, using loose wp
    return mask

def isHEMJet(pt, eta, phi):
    mask = (pt>30) & ((eta>-3.0)&(eta<-1.3)) & ((phi>-1.57)&(phi<-0.87))
    return mask

ids = {}
ids['isLooseElectron'] = isLooseElectron
ids['isTightElectron'] = isTightElectron
ids['isLooseMuon']     = isLooseMuon
ids['isTightMuon']     = isTightMuon
ids['isLooseTau']      = isLooseTau
ids['isLoosePhoton']   = isLoosePhoton
ids['isTightPhoton']   = isTightPhoton
ids['isGoodJet']       = isGoodJet
ids['isGoodFatJet']    = isGoodFatJet
ids['isHEMJet']        = isHEMJet
save(ids, 'data/ids.coffea')
