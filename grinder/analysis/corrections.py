#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
from Builder import Initialize
from fnal_column_analysis_tools import hist, lookup_tools

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

def get_ttbar_weight(pt):
    return np.exp(0.0615 - 0.0005 * np.clip(pt, 0, 800))

def get_nlo_weight(type, pt):
    #print('The pT is:',pt)
    kfactor = uproot.open("data/nlo/kfactors.root")
    sf_qcd = 1
    sf_ewk = 1
    #sf_qcd2j = 1

    lo = {}
    lo['z'] = "ZJets_LO/inv_pt"    
    lo['w'] = "WJets_LO/inv_pt"
    lo['a'] = "GJets_LO/inv_pt_G"

    nlo = {}
    nlo['z'] = "ZJets_012j_NLO/nominal"
    nlo['w'] = "WJets_012j_NLO/nominal"
    nlo['a'] = "GJets_1j_NLO/nominal_G"

    ewk = {}
    ewk['z'] = "EWKcorr/Z"
    ewk['w'] = "EWKcorr/W"
    ewk['a'] = "EWKcorr/photon"

    LO = kfactor[lo[type]].values
    NLO = kfactor[nlo[type]].values
    EWK = kfactor[ewk[type]].values

    sf_qcd = NLO / LO
    sf_ewk = EWK / NLO

    correction=lookup_tools.dense_lookup.dense_lookup(sf_qcd*sf_ewk, kfactor[nlo[type]].edges)
    return correction(pt)
