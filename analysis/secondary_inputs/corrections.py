#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
from coffea.arrays import Initialize
from coffea import hist, lookup_tools
from coffea.util import save

get_pu_weight = {}
get_pu_weight['2017'] = {}
get_pu_weight['2018'] = {}

pu = {}
pu["2018"] = uproot.open("secondary_inputs/pileup/puWeights_10x_56ifb.root")
pu["2017"] = uproot.open("secondary_inputs/pileup/puWeights_90x_41ifb.root")
pu["2016"] = uproot.open("secondary_inputs/pileup/puWeights_80x_37ifb.root")
for year in ['2016','2017','2018']:
    fpu = pu[year]
    pu_cen = fpu["puWeights"].values
    pu_up = fpu["puWeightsUp"].values
    pu_down= fpu["puWeightsDown"].values
    get_pu_weight[year] = {}
    get_pu_weight[year]['cen'] = lookup_tools.dense_lookup.dense_lookup(pu_cen, fpu["puWeights"].edges)
    get_pu_weight[year]['up'] = lookup_tools.dense_lookup.dense_lookup(pu_up, fpu["puWeightsUp"].edges)    
    get_pu_weight[year]['down'] = lookup_tools.dense_lookup.dense_lookup(pu_down, fpu["puWeightsDown"].edges)

get_met_trig_weight = {}

met_trig = {}
met_trig["2016"] = uproot.open("secondary_inputs/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")
met_trig["2017"] = uproot.open("secondary_inputs/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")
met_trig["2018"] = uproot.open("secondary_inputs/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")
for year in ['2016','2017','2018']:
    fmet_trig = met_trig[year]
    met_trig_corr = fmet_trig["hden_monojet_recoil_clone_passed"].values
    get_met_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(met_trig_corr, fmet_trig["hden_monojet_recoil_clone_passed"].edges)

get_met_zmm_trig_weight = {}

met_zmm_trig = {}
met_zmm_trig["2016"] = uproot.open("secondary_inputs/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root")
met_zmm_trig["2017"] = uproot.open("secondary_inputs/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root")
met_zmm_trig["2018"] = uproot.open("secondary_inputs/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root")
for year in ['2016','2017','2018']:
    fmet_zmm_trig = met_zmm_trig[year]
    met_zmm_trig_corr = fmet_zmm_trig["hden_monojet_recoil_clone_passed"].values
    get_met_zmm_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(met_zmm_trig_corr, fmet_zmm_trig["hden_monojet_recoil_clone_passed"].edges)


get_ele_trig_weight = {}

ele_trig = {}
ele_trig["2016"] = uproot.open("secondary_inputs/trigger_eff/eleTrig.root")
ele_trig["2017"] = uproot.open("secondary_inputs/trigger_eff/eleTrig.root")
ele_trig["2018"] = uproot.open("secondary_inputs/trigger_eff/eleTrig.root")

for year in ['2016','2017','2018']:
    fele_trig = ele_trig[year]
    ele_trig_corr = fele_trig["hEffEtaPt"].values
    get_ele_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(ele_trig_corr, fele_trig["hEffEtaPt"].edges)

get_pho_trig_weight = {}

pho_trig = {}
pho_trig["2016"] = uproot.open("secondary_inputs/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
pho_trig["2017"] = uproot.open("secondary_inputs/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
pho_trig["2018"] = uproot.open("secondary_inputs/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
for year in ['2016','2017','2018']:
    fpho_trig = pho_trig[year]
    pho_trig_corr = fpho_trig["hden_photonpt_clone_passed"].values
    get_pho_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(pho_trig_corr, fpho_trig["hden_photonpt_clone_passed"].edges)

def get_ttbar_weight(pt):
    return np.exp(0.0615 - 0.0005 * np.clip(pt, 0, 800))

get_nlo_weight = {}

kfactor = uproot.open("secondary_inputs/nlo/kfactors.root")
for year in ['2016','2017','2018']:

    get_nlo_weight[year] = {}

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

    for type in ['z','w','a']:
        LO = kfactor[lo[type]].values
        NLO = kfactor[nlo[type]].values
        EWK = kfactor[ewk[type]].values

        sf_qcd = NLO / LO
        sf_ewk = EWK / LO

        get_nlo_weight[year][type]=lookup_tools.dense_lookup.dense_lookup(sf_qcd*sf_ewk, kfactor[nlo[type]].edges)
        if (year != '2016' and type != 'a'): get_nlo_weight[year][type]=lookup_tools.dense_lookup.dense_lookup(sf_ewk, kfactor[nlo[type]].edges)

get_adhoc_weight = {}                                       
kfactor = uproot.open("secondary_inputs/nlo/2017_gen_v_pt_stat1_qcd_sf.root")
get_adhoc_weight['z']=lookup_tools.dense_lookup.dense_lookup(kfactor["dy_monojet"].values, kfactor["dy_monojet"].edges)
get_adhoc_weight['w']=lookup_tools.dense_lookup.dense_lookup(kfactor["wjet_monojet"].values, kfactor["wjet_monojet"].edges)


def get_ecal_bad_calib(run_number, lumi_number, event_number, year, dataset):
    bad = {}
    bad["2016"] = {}
    bad["2017"] = {}
    bad["2018"] = {}
    bad["2016"]["MET"]            = "secondary_inputs/ecalBadCalib/Run2016_MET.root"
    bad["2016"]["SinglePhoton"]   = "secondary_inputs/ecalBadCalib/Run2016_SinglePhoton.root"
    bad["2016"]["SingleElectron"] = "secondary_inputs/ecalBadCalib/Run2016_SingleElectron.root"
    bad["2017"]["MET"]            = "secondary_inputs/ecalBadCalib/Run2017_MET.root"
    bad["2017"]["SinglePhoton"]   = "secondary_inputs/ecalBadCalib/Run2017_SinglePhoton.root"
    bad["2017"]["SingleElectron"] = "secondary_inputs/ecalBadCalib/Run2017_SingleElectron.root"
    bad["2018"]["MET"]            = "secondary_inputs/ecalBadCalib/Run2018_MET.root"
    bad["2018"]["EGamma"]         = "secondary_inputs/ecalBadCalib/Run2018_EGamma.root"
    
    regular_dataset = ""
    regular_dataset = [name for name in ["MET","SinglePhoton","SingleElectron","EGamma"] if (name in dataset)]
    fbad = uproot.open(bad[year][regular_dataset[0]])
    bad_tree = fbad["vetoEvents"]
    runs_to_veto = bad_tree.array("Run")
    lumis_to_veto = bad_tree.array("LS")
    events_to_veto = bad_tree.array("Event")

    # We want events that do NOT have (a vetoed run AND a vetoed LS and a vetoed event number)
    return np.logical_not(np.isin(run_number, runs_to_veto) * np.isin(lumi_number, lumis_to_veto) * np.isin(event_number, events_to_veto))

corrections = {}
corrections['get_ttbar_weight']        = get_ttbar_weight
corrections['get_nlo_weight']          = get_nlo_weight
corrections['get_adhoc_weight']        = get_adhoc_weight
corrections['get_pu_weight']           = get_pu_weight
corrections['get_met_trig_weight']     = get_met_trig_weight
corrections['get_met_zmm_trig_weight'] = get_met_zmm_trig_weight
corrections['get_ele_trig_weight']     = get_ele_trig_weight
corrections['get_pho_trig_weight']     = get_pho_trig_weight
corrections['get_ecal_bad_calib']      = get_ecal_bad_calib
save(corrections, 'secondary_inputs/corrections.coffea')
