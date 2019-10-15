#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
from coffea.arrays import Initialize
from coffea import hist, lookup_tools

'''
def get_pu_weight(year, nvtx):
    mc["2016"] = "../data/pileup/pileup_profile_Summer16.root"
    mc["2017"] = "../data/pileup/pileup_profile_Summer16.root"

    fmc = uproot.open(mc[year])
    mc_pu = fmc["pileup"].values

    data["2016"] = "../data/pileup/PileupData_GoldenJSON_Full2016.root"
    data["2017"] = "../data/pileup/pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root"

    fdata = uproot.open(data[year])
    data_pu = fdata["pileup"].values

    corr = (data_pu / np.maximum(mc_pu, 1)) / (data_pu.sum() / mc_pu.sum())
    corr[mc_pu==0.] = 1.
    return lookup_tools.dense_lookup.dense_lookup(corr, fdata["pileup"].edges)
'''
def get_pu_weight(nvtx, year):
    pu = {}
    pu["2018"] = "data/pileup/puWeights_10x_56ifb.root"
    pu["2017"] = "data/pileup/puWeights_90x_41ifb.root"
    pu["2016"] = "data/pileup/puWeights_80x_37ifb.root"
    fpu = uproot.open(pu[year])
    pu_cen = fpu["puWeights"].values
    pu_up = fpu["puWeightsUp"].values
    pu_down= fpu["puWeightsDown"].values
    sf_pu_cen = lookup_tools.dense_lookup.dense_lookup(pu_cen, fpu["puWeights"].edges)
    sf_pu_up = lookup_tools.dense_lookup.dense_lookup(pu_up, fpu["puWeightsUp"].edges)    
    sf_pu_down = lookup_tools.dense_lookup.dense_lookup(pu_down, fpu["puWeightsDown"].edges)

    return sf_pu_cen(nvtx), sf_pu_up(nvtx), sf_pu_down(nvtx)

def get_met_trig_weight(metnomu, year):
    met_trig = {}
    met_trig["2016"] = "data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root"
    met_trig["2017"] = "data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root"
    met_trig["2018"] = "data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root"
    fmet_trig = uproot.open(met_trig[year])
    met_trig_corr = fmet_trig["hden_monojet_recoil_clone_passed"].values
    sf_met_trig = lookup_tools.dense_lookup.dense_lookup(met_trig_corr, fmet_trig["hden_monojet_recoil_clone_passed"].edges)

    return sf_met_trig(metnomu)

def get_met_zmm_trig_weight(metnomu, year):
    met_zmm_trig = {}
    met_zmm_trig["2016"] = "data/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root"
    met_zmm_trig["2017"] = "data/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root"
    met_zmm_trig["2018"] = "data/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root"
    fmet_zmm_trig = uproot.open(met_zmm_trig[year])
    met_zmm_trig_corr = fmet_zmm_trig["hden_monojet_recoil_clone_passed"].values
    sf_met_zmm_trig = lookup_tools.dense_lookup.dense_lookup(met_zmm_trig_corr, fmet_zmm_trig["hden_monojet_recoil_clone_passed"].edges)
    return sf_met_zmm_trig(metnomu)


def get_ele_trig_weight(eta1,pt1,eta2,pt2,year):
    eff1=np.zeros_like(pt1)
    eff2=np.zeros_like(pt2)
    ele_trig = {}
    ele_trig["2016"] = "data/trigger_eff/eleTrig.root"
    ele_trig["2017"] = "data/trigger_eff/eleTrig.root"
    ele_trig["2018"] = "data/trigger_eff/eleTrig.root"
    fele_trig = uproot.open(ele_trig[year])
    ele_trig_corr = fele_trig["hEffEtaPt"].values
    sf_ele_trig = lookup_tools.dense_lookup.dense_lookup(ele_trig_corr, fele_trig["hEffEtaPt"].edges)
    eff1=sf_ele_trig(eta1,pt1)
    if (eta2!=-99).all() and (pt2!=-99).all():
        eff2=sf_ele_trig(eta2,pt2)
    return 1 - (1-eff1)*(1-eff2)

def get_pho_trig_weight(pt, year):
    pho_trig = {}
    pho_trig["2016"] = "data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"
    pho_trig["2017"] = "data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"
    pho_trig["2018"] = "data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"
    fpho_trig = uproot.open(pho_trig[year])
    pho_trig_corr = fpho_trig["hden_photonpt_clone_passed"].values
    sf_pho_trig = lookup_tools.dense_lookup.dense_lookup(pho_trig_corr, fpho_trig["hden_photonpt_clone_passed"].edges)

    return sf_pho_trig(pt)

def get_ttbar_weight(pt):
    return np.exp(0.0615 - 0.0005 * np.clip(pt, 0, 800))

def get_nlo_weight(type,year, pt):
    #print('The pT is:',pt)
    kfactor = uproot.open("data/nlo/kfactors.root")

    sf_qcd = 1
    sf_ewk = 1
    #sf_adhoc = 1
    #sf_qcd2j = 1

    lo = {}
    lo['z'] = "ZJets_LO/inv_pt"    
    lo['w'] = "WJets_LO/inv_pt"
    lo['a'] = "GJets_LO/inv_pt_G"
    #lo['d'] = "DYJets_LO/inv_pt"

    nlo = {}
    nlo['z'] = "ZJets_012j_NLO/nominal"
    nlo['w'] = "WJets_012j_NLO/nominal"
    nlo['a'] = "GJets_1j_NLO/nominal_G"
    #nlo['d'] = "DYJets_012j_NLO/nominal"

    ewk = {}
    ewk['z'] = "EWKcorr/Z"
    ewk['w'] = "EWKcorr/W"
    ewk['a'] = "EWKcorr/photon"
    #ewk['d'] = "EWKcorr/DY"

    LO = kfactor[lo[type]].values
    NLO = kfactor[nlo[type]].values
    EWK = kfactor[ewk[type]].values

    sf_qcd = NLO / LO
    sf_ewk = EWK / LO

    correction=lookup_tools.dense_lookup.dense_lookup(sf_qcd*sf_ewk, kfactor[nlo[type]].edges)
    if (year != '2016' and type != 'a'):
        adhoc = uproot.open("data/nlo/2017_gen_v_pt_stat1_qcd_sf.root")
        nlo_lo = {}
        nlo_lo['z'] = "dy_monojet"
        nlo_lo['w'] = "wjet_monojet"
        sf_qcd = adhoc[nlo_lo[type]].values
        
        correction=lookup_tools.dense_lookup.dense_lookup(sf_ewk, kfactor[nlo[type]].edges)
        correction_qcd=lookup_tools.dense_lookup.dense_lookup(sf_qcd, adhoc[nlo_lo[type]].edges)
        return correction(pt)*correction_qcd(pt)

    return correction(pt)

### Obsolete
#def get_bad_ecal_weight(eta,phi):
#    badecal = "data/badecal/hotjets-runBCDEFGH.root"
#    fbadecal = uproot.open(badecal)
#    badecal_corr = fbadecal["h2jet"].values
#    correction=lookup_tools.dense_lookup.dense_lookup(badecal_corr, fbadecal["h2jet"].edges)
#    return correction(eta,phi)

def get_ecal_bad_calib(run_number, lumi_number, event_number, year, dataset):
    bad = {}
    bad["2016"] = {}
    bad["2017"] = {}
    bad["2018"] = {}
    bad["2016"]["MET"]            = "data/ecalBadCalib/Run2016_MET.root"
    bad["2016"]["SinglePhoton"]   = "data/ecalBadCalib/Run2016_SinglePhoton.root"
    bad["2016"]["SingleElectron"] = "data/ecalBadCalib/Run2016_SingleElectron.root"
    bad["2017"]["MET"]            = "data/ecalBadCalib/Run2017_MET.root"
    bad["2017"]["SinglePhoton"]   = "data/ecalBadCalib/Run2017_SinglePhoton.root"
    bad["2017"]["SingleElectron"] = "data/ecalBadCalib/Run2017_SingleElectron.root"
    bad["2018"]["MET"]            = "data/ecalBadCalib/Run2018_MET.root"
    bad["2018"]["EGamma"]         = "data/ecalBadCalib/Run2018_EGamma.root"
    
    regular_dataset = ""
    regular_dataset = [name for name in ["MET","SinglePhoton","SingleElectron","EGamma"] if (name in dataset)]
    fbad = uproot.open(bad[year][regular_dataset[0]])
    bad_tree = fbad["vetoEvents"]
    runs_to_veto = bad_tree.array("Run")
    lumis_to_veto = bad_tree.array("LS")
    events_to_veto = bad_tree.array("Event")

    # We want events that do NOT have (a vetoed run AND a vetoed LS and a vetoed event number)
    return np.logical_not(np.isin(run_number, runs_to_veto) * np.isin(lumi_number, lumis_to_veto) * np.isin(event_number, events_to_veto))
