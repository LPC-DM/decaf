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

def get_nlo_weight(gen):
    
