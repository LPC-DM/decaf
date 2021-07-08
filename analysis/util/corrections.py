#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
import os
from coffea import hist, lookup_tools
from coffea.lookup_tools import extractor, dense_lookup
from coffea.util import save, load
from coffea.btag_tools import BTagScaleFactor

###
# Pile-up weight
###

pu_files = {
    '2018': uproot.open("data/pileup/PileupHistograms_2018_69mb_pm5.root"),
    '2017': uproot.open("data/pileup/PileupHistograms_2017_69mb_pm5.root"),
    '2016': uproot.open("data/pileup/PileupHistograms_2016_69mb_pm5.root")
}
get_pu_weight = {}
for year in ['2016','2017','2018']:
    pu_hist=pu_files[year]['pu_weights_central']
    get_pu_weight[year] = lookup_tools.dense_lookup.dense_lookup(pu_hist.values, pu_hist.edges)

###
# MET trigger efficiency SFs, 2017/18 from monojet. Depends on recoil.
###

met_trig_hists = {
    '2016': uproot.open("data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")['hden_monojet_recoil_clone_passed'],
    '2017': uproot.open("data/trigger_eff/met_trigger_sf.root")['120pfht_hltmu_1m_2017'],
    '2018': uproot.open("data/trigger_eff/met_trigger_sf.root")['120pfht_hltmu_1m_2018']
}
get_met_trig_weight = {}
for year in ['2016','2017','2018']:
    met_trig_hist=met_trig_hists[year]
    get_met_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(met_trig_hist.values, met_trig_hist.edges)

###
# MET z->mumu efficiency SF. 2017/18 using 1m as done in monojet, 2m used only for systematics. Depends on recoil.
###

zmm_trig_hists ={
    '2016': uproot.open("data/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root")['hden_monojet_recoil_clone_passed'],
    '2017': uproot.open("data/trigger_eff/met_trigger_sf.root")['120pfht_hltmu_1m_2017'],
    '2018': uproot.open("data/trigger_eff/met_trigger_sf.root")['120pfht_hltmu_1m_2018']
}
get_met_zmm_trig_weight = {}
for year in ['2016','2017','2018']:
    zmm_trig_hist = zmm_trig_hists[year]
    get_met_zmm_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(zmm_trig_hist.values, zmm_trig_hist.edges)

###
# Electron trigger efficiency SFs. depends on supercluster eta and pt:
###

ele_trig_hists = {
    '2016': uproot.open("data/trigger_eff/eleTrig.root")['hEffEtaPt'],
    '2017': uproot.open("data/trigger_eff/electron_trigger_sf_2017.root")['EGamma_SF2D'],#monojet measurement for the combined trigger path
    '2018': uproot.open("data/trigger_eff/electron_trigger_sf_2018.root")['EGamma_SF2D'] #approved by egamma group: https://indico.cern.ch/event/924522/
}
get_ele_trig_weight = {}
for year in ['2016','2017','2018']:
    ele_trig_hist = ele_trig_hists[year]
    get_ele_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(ele_trig_hist.values, ele_trig_hist.edges)

###
# Photon trigger efficiency SFs. 2017/18 not actually used, sigmoid is used instead.
###

pho_trig_files = {
    '2016': uproot.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"),
    "2017": uproot.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root"),
    "2018": uproot.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
}
get_pho_trig_weight = {}
for year in ['2016','2017','2018']:
    pho_trig_hist = pho_trig_files[year]["hden_photonpt_clone_passed"]
    get_pho_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(pho_trig_hist.values, pho_trig_hist.edges)

###
# Electron id SFs. 2017/18 used dedicated weights from monojet. depends on supercluster eta and pt.
###

ele_loose_files = {
    '2016': uproot.open("data/ScaleFactor/2016_ElectronWPVeto_Fall17V2.root"),
    '2017': uproot.open("data/ScaleFactor/2017_ElectronWPVeto_Fall17V2_BU.root"),
    '2018': uproot.open("data/ScaleFactor/2018_ElectronWPVeto_Fall17V2_BU.root")
}
ele_tight_files = {
    '2016': uproot.open("data/ScaleFactor/2016LegacyReReco_ElectronTight_Fall17V2.root"),
    '2017': uproot.open("data/ScaleFactor/2017_ElectronTight_Fall17V2_BU.root"),
    '2018': uproot.open("data/ScaleFactor/2018_ElectronTight_Fall17V2_BU.root")
}
get_ele_loose_id_sf = {}
get_ele_tight_id_sf = {}
for year in ['2016','2017','2018']:
    ele_loose_sf_hist = ele_loose_files[year]["EGamma_SF2D"]
    get_ele_loose_id_sf[year]  = lookup_tools.dense_lookup.dense_lookup(ele_loose_sf_hist.values, ele_loose_sf_hist.edges)
    ele_tight_sf_hist =ele_tight_files[year]["EGamma_SF2D"]
    get_ele_tight_id_sf[year]  = lookup_tools.dense_lookup.dense_lookup(ele_tight_sf_hist.values, ele_tight_sf_hist.edges)

###
# Electron reconstruction SFs. Depends on supercluster eta and pt.    
###

ele_reco_files = {
    '2016': uproot.open("data/ScaleFactor/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root"),
    '2017': uproot.open("data/ScaleFactor/2017_egammaEffi_txt_EGM2D_runBCDEF_passingRECO.root"),
    '2018': uproot.open("data/ScaleFactor/2018_egammaEffi_txt_EGM2D_updatedAll.root")
}
get_ele_reco_sf = {}
for year in ['2016','2017','2018']:
    ele_reco_hist = ele_reco_files[year]["EGamma_SF2D"]
    get_ele_reco_sf[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.values, ele_reco_hist.edges)
#2017 has a separate set of weights for low pt electrons (pt<20).
ele_reco_lowet_hist = uproot.open("data/ScaleFactor/2017_egammaEffi_txt_EGM2D_runBCDEF_passingRECO_lowEt.root")['EGamma_SF2D']
get_ele_reco_lowet_sf=lookup_tools.dense_lookup.dense_lookup(ele_reco_lowet_hist.values, ele_reco_lowet_hist.edges)

###
# Photon ID SFs. Tight photons use medium id. 2017/18 use dedicated measurement from monojet, depends only on abs(eta): https://indico.cern.ch/event/879924/
###

pho_tight_hists = {
    '2016': uproot.open("data/ScaleFactor/Fall17V2_2016_Medium_photons.root")['EGamma_SF2D'],
    '2017': uproot.open("data/ScaleFactor/photon_medium_id_sf_v0.root")['photon_medium_id_sf_2017'],
    '2018': uproot.open("data/ScaleFactor/photon_medium_id_sf_v0.root")['photon_medium_id_sf_2018']
}
get_pho_tight_id_sf = {}
for year in ['2016','2017','2018']:
    pho_tight_hist=pho_tight_hists[year]
    get_pho_tight_id_sf[year] = lookup_tools.dense_lookup.dense_lookup(pho_tight_hist.values, pho_tight_hist.edges)

###
# Photon CSEV weight: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_Veto_CSEV_or_pixel_seed
###

pho_csev_hists = {
    '2016': uproot.open("data/ScaleFactor/ScalingFactors_80X_Summer16_rename.root")['Scaling_Factors_CSEV_R9_Inclusive'],
    '2017': uproot.open("data/ScaleFactor/CSEV_ScaleFactors_2017.root")['Medium_ID'],
    '2018': uproot.open("data/ScaleFactor/CSEV_2018.root")['eleVeto_SF'],
}
get_pho_csev_sf = {}
for year in ['2016','2017','2018']:
    pho_csev_hist=pho_csev_hists[year]
    get_pho_csev_sf[year] = lookup_tools.dense_lookup.dense_lookup(pho_csev_hist.values, pho_csev_hist.edges)

###
# Muon ID SFs
###

mu_files = {
    '2016': uproot.open("data/ScaleFactor/2016LegacyReReco_Muon_SF_ID.root"),
    '2017': uproot.open("data/ScaleFactor/2017_Muon_RunBCDEF_SF_ID.root"),
    '2018': uproot.open("data/ScaleFactor/2018_Muon_RunABCD_SF_ID.root")
}
mu_tight_hist = {
    '2016': mu_files['2016']["NUM_TightID_DEN_genTracks_eta_pt"],
    '2017': mu_files['2017']["NUM_TightID_DEN_genTracks_pt_abseta"],
    '2018': mu_files['2018']["NUM_TightID_DEN_TrackerMuons_pt_abseta"]
}
mu_loose_hist = {
    '2016': mu_files['2016']["NUM_LooseID_DEN_genTracks_eta_pt"],
    '2017': mu_files['2017']["NUM_LooseID_DEN_genTracks_pt_abseta"],
    '2018': mu_files['2018']["NUM_LooseID_DEN_TrackerMuons_pt_abseta"]
}
get_mu_tight_id_sf = {}
get_mu_loose_id_sf = {}
for year in ['2016','2017','2018']:
    get_mu_tight_id_sf[year] = lookup_tools.dense_lookup.dense_lookup(mu_tight_hist[year].values, mu_tight_hist[year].edges)
    get_mu_loose_id_sf[year] = lookup_tools.dense_lookup.dense_lookup(mu_loose_hist[year].values, mu_loose_hist[year].edges)

###
# Muon isolation SFs
###

mu_iso_files = {
    '2016': uproot.open("data/ScaleFactor/Merged_SF_ISO.root"),
    '2017': uproot.open("data/ScaleFactor/RunBCDEF_SF_ISO_syst.root"),
    '2018': uproot.open("data/ScaleFactor/RunABCD_SF_ISO.root")
}
mu_iso_tight_hist = {
    '2016': mu_iso_files['2016']["NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"],
    '2017': mu_iso_files['2017']["NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"],
    '2018': mu_iso_files['2018']["NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"]
}
mu_iso_loose_hist = {
    '2016': mu_iso_files['2016']["NUM_LooseRelIso_DEN_LooseID_eta_pt"],
    '2017': mu_iso_files['2017']["NUM_LooseRelIso_DEN_LooseID_pt_abseta"],
    '2018': mu_iso_files['2018']["NUM_LooseRelIso_DEN_LooseID_pt_abseta"]
}
get_mu_tight_iso_sf = {}
get_mu_loose_iso_sf = {}
for year in ['2016','2017','2018']:
    get_mu_tight_iso_sf[year] = lookup_tools.dense_lookup.dense_lookup(mu_iso_tight_hist[year].values, mu_iso_tight_hist[year].edges)
    get_mu_loose_iso_sf[year] = lookup_tools.dense_lookup.dense_lookup(mu_iso_loose_hist[year].values, mu_iso_loose_hist[year].edges)
###
# V+jets NLO k-factors
###

nlo_qcd_hists = {
    '2016':{
        'dy': uproot.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_qcd"],
        'w': uproot.open("data/vjets_SFs/merged_kfactors_wjets.root")["kfactor_monojet_qcd"],
        'z': uproot.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_qcd"],
        'a': uproot.open("data/vjets_SFs/merged_kfactors_gjets.root")["kfactor_monojet_qcd"]
    },
    '2017':{
        'z': uproot.open("data/vjets_SFs/SF_QCD_NLO_ZJetsToNuNu.root")["kfac_znn_filter"],
        'w': uproot.open("data/vjets_SFs/SF_QCD_NLO_WJetsToLNu.root")["wjet_dress_monojet"],
        'dy': uproot.open("data/vjets_SFs/SF_QCD_NLO_DYJetsToLL.root")["kfac_dy_filter"],
        'a': uproot.open("data/vjets_SFs/SF_QCD_NLO_GJets.root")["gjets_stat1_monojet"]
    },
    '2018':{
        'z': uproot.open("data/vjets_SFs/SF_QCD_NLO_ZJetsToNuNu.root")["kfac_znn_filter"],
        'w': uproot.open("data/vjets_SFs/SF_QCD_NLO_WJetsToLNu.root")["wjet_dress_monojet"],
        'dy': uproot.open("data/vjets_SFs/SF_QCD_NLO_DYJetsToLL.root")["kfac_dy_filter"],
        'a': uproot.open("data/vjets_SFs/SF_QCD_NLO_GJets.root")["gjets_stat1_monojet"]
    }
}
nlo_ewk_hists = {
    'dy': uproot.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_ewk"],
    'w': uproot.open("data/vjets_SFs/merged_kfactors_wjets.root")["kfactor_monojet_ewk"],
    'z': uproot.open("data/vjets_SFs/merged_kfactors_zjets.root")["kfactor_monojet_ewk"],
    'a': uproot.open("data/vjets_SFs/merged_kfactors_gjets.root")["kfactor_monojet_ewk"]
}    
get_nlo_qcd_weight = {}
get_nlo_ewk_weight = {}
for year in ['2016','2017','2018']:
    get_nlo_qcd_weight[year] = {}
    get_nlo_ewk_weight[year] = {}
    for p in ['dy','w','z','a']:
        get_nlo_qcd_weight[year][p] = lookup_tools.dense_lookup.dense_lookup(nlo_qcd_hists[year][p].values, nlo_qcd_hists[year][p].edges)
        get_nlo_ewk_weight[year][p] = lookup_tools.dense_lookup.dense_lookup(nlo_ewk_hists[p].values, nlo_ewk_hists[p].edges)

###
# V+jets NNLO weights
# The schema is process_NNLO_NLO_QCD1QCD2QCD3_EW1EW2EW3_MIX, where 'n' stands for 'nominal', 'u' for 'up', and 'd' for 'down'
###

histname={
    'dy': 'eej_NNLO_NLO_',
    'w':  'evj_NNLO_NLO_',
    'z': 'vvj_NNLO_NLO_',
    'a': 'aj_NNLO_NLO_'
}
correlated_variations = {
    'cen':    'nnn_nnn_n',
    'qcd1up': 'unn_nnn_n',
    'qcd1do': 'dnn_nnn_n',
    'qcd2up': 'nun_nnn_n',
    'qcd2do': 'ndn_nnn_n',
    'qcd3up': 'nnu_nnn_n',
    'qcd3do': 'nnd_nnn_n',
    'ew1up' : 'nnn_unn_n',
    'ew1do' : 'nnn_dnn_n',
    'mixup' : 'nnn_nnn_u',
    'mixdo' : 'nnn_nnn_d',
    'muFup' : 'nnn_nnn_n_Weight_scale_variation_muR_1p0_muF_2p0',
    'muFdo' : 'nnn_nnn_n_Weight_scale_variation_muR_1p0_muF_0p5',
    'muRup' : 'nnn_nnn_n_Weight_scale_variation_muR_2p0_muF_1p0',
    'muRdo' : 'nnn_nnn_n_Weight_scale_variation_muR_0p5_muF_1p0'
}
uncorrelated_variations = {
    'dy': {
        'ew2Gup': 'nnn_nnn_n',
        'ew2Gdo': 'nnn_nnn_n',
        'ew2Wup': 'nnn_nnn_n',
        'ew2Wdo': 'nnn_nnn_n',
        'ew2Zup': 'nnn_nun_n',
        'ew2Zdo': 'nnn_ndn_n',
        'ew3Gup': 'nnn_nnn_n',
        'ew3Gdo': 'nnn_nnn_n',
        'ew3Wup': 'nnn_nnn_n',
        'ew3Wdo': 'nnn_nnn_n',
        'ew3Zup': 'nnn_nnu_n',
        'ew3Zdo': 'nnn_nnd_n'
    },
    'w': {
        'ew2Gup': 'nnn_nnn_n',
        'ew2Gdo': 'nnn_nnn_n',
        'ew2Wup': 'nnn_nun_n',
        'ew2Wdo': 'nnn_ndn_n',
        'ew2Zup': 'nnn_nnn_n',
        'ew2Zdo': 'nnn_nnn_n',
        'ew3Gup': 'nnn_nnn_n',
        'ew3Gdo': 'nnn_nnn_n',
        'ew3Wup': 'nnn_nnu_n',
        'ew3Wdo': 'nnn_nnd_n',
        'ew3Zup': 'nnn_nnn_n',
        'ew3Zdo': 'nnn_nnn_n'
    },
    'z': {
        'ew2Gup': 'nnn_nnn_n',
        'ew2Gdo': 'nnn_nnn_n',
        'ew2Wup': 'nnn_nnn_n',
        'ew2Wdo': 'nnn_nnn_n',
        'ew2Zup': 'nnn_nun_n',
        'ew2Zdo': 'nnn_ndn_n',
        'ew3Gup': 'nnn_nnn_n',
        'ew3Gdo': 'nnn_nnn_n',
        'ew3Wup': 'nnn_nnn_n',
        'ew3Wdo': 'nnn_nnn_n',
        'ew3Zup': 'nnn_nnu_n',
        'ew3Zdo': 'nnn_nnd_n'
    },
    'a': {
        'ew2Gup': 'nnn_nun_n',
        'ew2Gdo': 'nnn_ndn_n',
        'ew2Wup': 'nnn_nnn_n',
        'ew2Wdo': 'nnn_nnn_n',
        'ew2Zup': 'nnn_nnn_n',
        'ew2Zdo': 'nnn_nnn_n',
        'ew3Gup': 'nnn_nnu_n',
        'ew3Gdo': 'nnn_nnd_n',
        'ew3Wup': 'nnn_nnn_n',
        'ew3Wdo': 'nnn_nnn_n',
        'ew3Zup': 'nnn_nnn_n',
        'ew3Zdo': 'nnn_nnn_n'
    }
}
get_nnlo_nlo_weight = {}
for year in ['2016','2017','2018']:
    get_nnlo_nlo_weight[year] = {}
    nnlo_file = {
        'dy': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_eej_madgraph_"+year+".root"),
        'w': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_evj_madgraph_"+year+".root"),
        'z': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_vvj_madgraph_"+year+".root"),
        'a': uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_aj_madgraph_"+year+".root")
    }
    for p in ['dy','w','z','a']:
        get_nnlo_nlo_weight[year][p] = {}
        for cv in correlated_variations:
            hist=nnlo_file[p][histname[p]+correlated_variations[cv]]
            get_nnlo_nlo_weight[year][p][cv]=lookup_tools.dense_lookup.dense_lookup(hist.values, hist.edges)
        for uv in uncorrelated_variations[p]:
            hist=nnlo_file[p][histname[p]+uncorrelated_variations[p][uv]]
            get_nnlo_nlo_weight[year][p][uv]=lookup_tools.dense_lookup.dense_lookup(hist.values, hist.edges)

def get_ttbar_weight(pt):
    return np.exp(0.0615 - 0.0005 * np.clip(pt, 0, 800))

def get_msd_weight(pt, eta):
    gpar = np.array([1.00626, -1.06161, 0.0799900, 1.20454])
    cpar = np.array([1.09302, -0.000150068, 3.44866e-07, -2.68100e-10, 8.67440e-14, -1.00114e-17])
    fpar = np.array([1.27212, -0.000571640, 8.37289e-07, -5.20433e-10, 1.45375e-13, -1.50389e-17])
    genw = gpar[0] + gpar[1]*np.power(pt*gpar[2], -gpar[3])
    ptpow = np.power.outer(pt, np.arange(cpar.size))
    cenweight = np.dot(ptpow, cpar)
    forweight = np.dot(ptpow, fpar)
    weight = np.where(np.abs(eta)<1.3, cenweight, forweight)
    return genw*weight


def get_ecal_bad_calib(run_number, lumi_number, event_number, year, dataset):
    bad = {}
    bad["2016"] = {}
    bad["2017"] = {}
    bad["2018"] = {}
    bad["2016"]["MET"]            = "ecalBadCalib/Run2016_MET.root"
    bad["2016"]["SinglePhoton"]   = "ecalBadCalib/Run2016_SinglePhoton.root"
    bad["2016"]["SingleElectron"] = "ecalBadCalib/Run2016_SingleElectron.root"
    bad["2017"]["MET"]            = "ecalBadCalib/Run2017_MET.root"
    bad["2017"]["SinglePhoton"]   = "ecalBadCalib/Run2017_SinglePhoton.root"
    bad["2017"]["SingleElectron"] = "ecalBadCalib/Run2017_SingleElectron.root"
    bad["2018"]["MET"]            = "ecalBadCalib/Run2018_MET.root"
    bad["2018"]["EGamma"]         = "ecalBadCalib/Run2018_EGamma.root"
    
    regular_dataset = ""
    regular_dataset = [name for name in ["MET","SinglePhoton","SingleElectron","EGamma"] if (name in dataset)]
    fbad = uproot.open(bad[year][regular_dataset[0]])
    bad_tree = fbad["vetoEvents"]
    runs_to_veto = bad_tree.array("Run")
    lumis_to_veto = bad_tree.array("LS")
    events_to_veto = bad_tree.array("Event")

    # We want events that do NOT have (a vetoed run AND a vetoed LS and a vetoed event number)
    return np.logical_not(np.isin(run_number, runs_to_veto) * np.isin(lumi_number, lumis_to_veto) * np.isin(event_number, events_to_veto))

class BTagCorrector:

    def __init__(self, tagger, year, workingpoint):
        self._year = year
        common = load('data/common.coffea')
        self._wp = common['btagWPs'][tagger][year][workingpoint]
        files = {
            'deepflav': {
                '2016': 'DeepJet_2016LegacySF_V1.csv',
                '2017': 'DeepFlavour_94XSF_V4_B_F.csv',
                '2018': 'DeepJet_102XSF_V1.csv'
                },
            'deepcsv': {
                '2016': 'DeepCSV_2016LegacySF_V1.csv',
                '2017': 'DeepCSV_94XSF_V5_B_F.csv',
                '2018': 'DeepCSV_102XSF_V1.csv'      
                }
        }
        filename = 'data/'+files[tagger][year]
        self.sf = BTagScaleFactor(filename, workingpoint)
        files = {
            '2016': 'btag2016.merged',
            '2017': 'btag2017.merged',
            '2018': 'btag2018.merged',
        }
        filename = 'hists/'+files[year]
        btag = load(filename)
        bpass = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag', 'pass').values()[()]
        ball = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag').values()[()]
        nom = bpass / np.maximum(ball, 1.)
        self.eff = lookup_tools.dense_lookup.dense_lookup(nom, [ax.edges() for ax in btag[tagger].axes()[3:]])

    def btag_weight(self, pt, eta, flavor, tag):
        abseta = abs(eta)
        
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1b_Event_reweighting_using_scale
        def zerotag(eff):
            return (1 - eff).prod()

        eff = self.eff(flavor, pt, abseta)
        sf_nom = self.sf.eval('central', flavor, abseta, pt)
        sf_up = self.sf.eval('up', flavor, abseta, pt)
        sf_down = self.sf.eval('down', flavor, abseta, pt)

        eff_data_nom  = np.minimum(1., sf_nom*eff)
        eff_data_up   = np.minimum(1., sf_up*eff)
        eff_data_down = np.minimum(1., sf_down*eff)

        nom = zerotag(eff_data_nom)/zerotag(eff)
        up = zerotag(eff_data_up)/zerotag(eff)
        down = zerotag(eff_data_down)/zerotag(eff)

        if '-1' in tag: 
            nom = (1 - zerotag(eff_data_nom)) / (1 - zerotag(eff))
            up = (1 - zerotag(eff_data_up)) / (1 - zerotag(eff))
            down = (1 - zerotag(eff_data_down)) / (1 - zerotag(eff))

        return np.nan_to_num(nom), np.nan_to_num(up), np.nan_to_num(down)

get_btag_weight = {
    'deepflav': {
        '2016': {
            'loose'  : BTagCorrector('deepflav','2016','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2016','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2016','tight').btag_weight
        },
        '2017': {
            'loose'  : BTagCorrector('deepflav','2017','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2017','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2017','tight').btag_weight
        },
        '2018': {
            'loose'  : BTagCorrector('deepflav','2018','loose').btag_weight,
            'medium' : BTagCorrector('deepflav','2018','medium').btag_weight,
            'tight'  : BTagCorrector('deepflav','2018','tight').btag_weight
        }
    },
    'deepcsv' : {
        '2016': {
            'loose'  : BTagCorrector('deepcsv','2016','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2016','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2016','tight').btag_weight
        },
        '2017': {
            'loose'  : BTagCorrector('deepcsv','2017','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2017','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2017','tight').btag_weight
        },
        '2018': {
            'loose'  : BTagCorrector('deepcsv','2018','loose').btag_weight,
            'medium' : BTagCorrector('deepcsv','2018','medium').btag_weight,
            'tight'  : BTagCorrector('deepcsv','2018','tight').btag_weight
        }
    }
}

Jetext = extractor()
for directory in ['jec', 'jersf', 'jr', 'junc']:
    directory='data/'+directory
    print('Loading files in:',directory)
    for filename in os.listdir(directory):
        if '~' in filename: continue
        if 'AK4PFPuppi' not in filename: continue
        if 'DATA' in filename: continue
        filename=directory+'/'+filename
        print('Loading file:',filename)
        Jetext.add_weight_sets(['* * '+filename])
    print('All files in',directory,'loaded')
Jetext.finalize()
Jetevaluator = Jetext.make_evaluator()

corrections = {
    'get_msd_weight':           get_msd_weight,
    'get_ttbar_weight':         get_ttbar_weight,
    'get_nnlo_nlo_weight':      get_nnlo_nlo_weight,
    'get_nlo_qcd_weight':       get_nlo_qcd_weight,
    'get_nlo_ewk_weight':       get_nlo_ewk_weight,
    'get_pu_weight':            get_pu_weight,
    'get_met_trig_weight':      get_met_trig_weight,
    'get_met_zmm_trig_weight':  get_met_zmm_trig_weight,
    'get_ele_trig_weight':      get_ele_trig_weight,
    'get_pho_trig_weight':      get_pho_trig_weight,
    'get_ele_loose_id_sf':      get_ele_loose_id_sf,
    'get_ele_tight_id_sf':      get_ele_tight_id_sf,
    'get_pho_tight_id_sf':      get_pho_tight_id_sf,
    'get_pho_csev_sf':          get_pho_csev_sf,
    'get_mu_tight_id_sf':       get_mu_tight_id_sf,
    'get_mu_loose_id_sf':       get_mu_loose_id_sf,
    'get_ele_reco_sf':          get_ele_reco_sf,
    'get_ele_reco_lowet_sf':    get_ele_reco_lowet_sf,
    'get_mu_tight_iso_sf':      get_mu_tight_iso_sf,
    'get_mu_loose_iso_sf':      get_mu_loose_iso_sf,
    'get_ecal_bad_calib':       get_ecal_bad_calib,
    'get_btag_weight':          get_btag_weight,
    'Jetevaluator':             Jetevaluator,
}
save(corrections, 'data/corrections.coffea')



 
