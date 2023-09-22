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
# Muon scale and resolution (i.e. Rochester)
###

tag = 'roccor.Run2.v5'
get_mu_rochester_sf = {}
for year in ['2016','2017','2018']:
    fname = f'data/{tag}/RoccoR{year}.txt'
    sfs = lookup_tools.txt_converters.convert_rochester_file(fname,loaduncs=True)
    get_mu_rochester_sf[year] = lookup_tools.rochester_lookup.rochester_lookup(sfs)


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
                '2016': 'DeepJet_2016LegacySF_V1_YearCorrelation-V1.csv',
                '2017': 'DeepFlavour_94XSF_V3_B_F_comb_YearCorrelation-V1.csv',
                '2018': 'DeepJet_102XSF_V1_YearCorrelation-V1.csv'
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
            '2016': 'btageff2016.merged',
            '2017': 'btageff2017.merged',
            '2018': 'btageff2018.merged',
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

        bc = flavor > 0
        light = ~bc
        
        eff = self.eff(flavor, pt, abseta)
        
        sf_nom = self.sf.eval('central', flavor, abseta, pt)
        
        bc_sf_up_correlated = pt.ones_like()
        bc_sf_up_correlated[~bc] = sf_nom[~bc]
        bc_sf_up_correlated[bc] = self.sf.eval('up_correlated', flavor, eta, pt)[bc]
        
        bc_sf_down_correlated = pt.ones_like()
        bc_sf_down_correlated[~bc] = sf_nom[~bc]
        bc_sf_down_correlated[bc] = self.sf.eval('down_correlated', flavor, eta, pt)[bc]

        bc_sf_up_uncorrelated = pt.ones_like()
        bc_sf_up_uncorrelated[~bc] = sf_nom[~bc]
        bc_sf_up_uncorrelated[bc] = self.sf.eval('up_uncorrelated', flavor, eta, pt)[bc]

        bc_sf_down_uncorrelated = pt.ones_like()
        bc_sf_down_uncorrelated[~bc] = sf_nom[~bc]
        bc_sf_down_uncorrelated[bc] = self.sf.eval('down_uncorrelated', flavor, eta, pt)[bc]

        light_sf_up_correlated = pt.ones_like()
        light_sf_up_correlated[~light] = sf_nom[~light]
        light_sf_up_correlated[light] = self.sf.eval('up_correlated', flavor, abseta, pt)[light]

        light_sf_down_correlated = pt.ones_like()
        light_sf_down_correlated[~light] = sf_nom[~light]
        light_sf_down_correlated[light] = self.sf.eval('down_correlated', flavor, abseta, pt)[light]

        light_sf_up_uncorrelated = pt.ones_like()
        light_sf_up_uncorrelated[~light] = sf_nom[~light]
        light_sf_up_uncorrelated[light] = self.sf.eval('up_uncorrelated', flavor, abseta, pt)[light]

        light_sf_down_uncorrelated = pt.ones_like()
        light_sf_down_uncorrelated[~light] = sf_nom[~light]
        light_sf_down_uncorrelated[light] = self.sf.eval('down_uncorrelated', flavor, abseta, pt)[light]

        eff_data_nom  = np.minimum(1., sf_nom*eff)
        bc_eff_data_up_correlated   = np.minimum(1., bc_sf_up_correlated*eff)
        bc_eff_data_down_correlated = np.minimum(1., bc_sf_down_correlated*eff)
        bc_eff_data_up_uncorrelated   = np.minimum(1., bc_sf_up_uncorrelated*eff)
        bc_eff_data_down_uncorrelated = np.minimum(1., bc_sf_down_uncorrelated*eff)
        light_eff_data_up_correlated   = np.minimum(1., light_sf_up_correlated*eff)
        light_eff_data_down_correlated = np.minimum(1., light_sf_down_correlated*eff)
        light_eff_data_up_uncorrelated   = np.minimum(1., light_sf_up_uncorrelated*eff)
        light_eff_data_down_uncorrelated = np.minimum(1., light_sf_down_uncorrelated*eff)
       
        nom = zerotag(eff_data_nom)/zerotag(eff)
        bc_up_correlated = zerotag(bc_eff_data_up_correlated)/zerotag(eff)
        bc_down_correlated = zerotag(bc_eff_data_down_correlated)/zerotag(eff)
        bc_up_uncorrelated = zerotag(bc_eff_data_up_uncorrelated)/zerotag(eff)
        bc_down_uncorrelated = zerotag(bc_eff_data_down_uncorrelated)/zerotag(eff)
        light_up_correlated = zerotag(light_eff_data_up_correlated)/zerotag(eff)
        light_down_correlated = zerotag(light_eff_data_down_correlated)/zerotag(eff)
        light_up_uncorrelated = zerotag(light_eff_data_up_uncorrelated)/zerotag(eff)
        light_down_uncorrelated = zerotag(light_eff_data_down_uncorrelated)/zerotag(eff)
        
        if '-1' in tag: 
            nom = (1 - zerotag(eff_data_nom)) / (1 - zerotag(eff))
            bc_up_correlated = (1 - zerotag(bc_eff_data_up_correlated)) / (1 - zerotag(eff))
            bc_down_correlated = (1 - zerotag(bc_eff_data_down_correlated)) / (1 - zerotag(eff))
            bc_up_uncorrelated = (1 - zerotag(bc_eff_data_up_uncorrelated)) / (1 - zerotag(eff))
            bc_down_uncorrelated = (1 - zerotag(bc_eff_data_down_uncorrelated)) / (1 - zerotag(eff))
            light_up_correlated = (1 - zerotag(light_eff_data_up_correlated)) / (1 - zerotag(eff))
            light_down_correlated = (1 - zerotag(light_eff_data_down_correlated)) / (1 - zerotag(eff))
            light_up_uncorrelated = (1 - zerotag(light_eff_data_up_uncorrelated)) / (1 - zerotag(eff))
            light_down_uncorrelated = (1 - zerotag(light_eff_data_down_uncorrelated)) / (1 - zerotag(eff))

        return np.nan_to_num(nom, nan=1.), np.nan_to_num(bc_up_correlated, nan=1.), np.nan_to_num(bc_down_correlated, nan=1.), np.nan_to_num(bc_up_uncorrelated, nan=1.), np.nan_to_num(bc_down_uncorrelated, nan=1.), np.nan_to_num(light_up_correlated, nan=1.), np.nan_to_num(light_down_correlated, nan=1.), np.nan_to_num(light_up_uncorrelated, nan=1.), np.nan_to_num(light_down_uncorrelated, nan=1.)

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

class DoubleBTagCorrector:

    def __init__(self, year):
        self._year = year
        sf = {
            '2018': {
                'value': np.array([0.82, 0.82, 0.75, 0.81]),
                'unc': np.array([2*np.sqrt(0.07**2 + 0.11**2), np.sqrt(0.07**2 + 0.11**2), np.sqrt(0.06**2 + 0.06**2), np.sqrt(0.05**2 + 0.01**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
            '2017': {
                'value': np.array([0.84, 0.84, 0.98, 0.86]),
                'unc': np.array([2*np.sqrt(0.05**2 + 0.13**2), np.sqrt(0.05**2 + 0.13**2), np.sqrt(0.05**2 + 0.12**2), np.sqrt(0.05**2 + 0.05**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
            '2016': {
                'value': np.array([1.01, 1.01, 0.95, 0.99]),
                'unc': np.array([2*np.sqrt(0.06**2 + 0.02**2), np.sqrt(0.06**2 + 0.02**2), np.sqrt(0.05**2 + 0.09**2), np.sqrt(0.06**2 + 0.00**2)]),
                'edges': np.array([160, 350, 400, 500, 2500])
            },
        }
        self.sf_nom={}
        self.sf_up={} 
        self.sf_down={}
        
        for i in range(len(sf[year]['value'])):
            
            nom = np.ones_like(sf[year]['value'])
            nom[i] = sf[year]['value'][i]
            unc = np.zeros_like(sf[year]['value'])
            unc[i] = sf[year]['unc'][i]
            print(i,nom)
            self.sf_nom['Pt'+str(i)] = lookup_tools.dense_lookup.dense_lookup(nom, sf[year]['edges'])
            self.sf_up['Pt'+str(i)] = lookup_tools.dense_lookup.dense_lookup(nom+unc, sf[year]['edges'])
            self.sf_down['Pt'+str(i)] = lookup_tools.dense_lookup.dense_lookup(nom-unc, sf[year]['edges'])
        

    def doublebtag_weight(self, pt):
        sf_nom={} 
        sf_up={} 
        sf_down={}
        for k in self.sf_nom:
            sf_nom[k], sf_up[k], sf_down[k] = self.sf_nom[k](pt), self.sf_up[k](pt), self.sf_down[k](pt)
        return sf_nom, sf_up, sf_down

get_doublebtag_weight = {
    '2016': DoubleBTagCorrector('2016').doublebtag_weight,
    '2017': DoubleBTagCorrector('2017').doublebtag_weight,
    '2018': DoubleBTagCorrector('2018').doublebtag_weight,
}

class Reweighting:

    def __init__(self, year):
        self._year = year
        files = {
            '2016': 'reweighting2016.scaled',
            '2017': 'reweighting2017.scaled',
            '2018': 'reweighting2018.scaled',
        }
        filename = 'hists/'+files[year]
        reweighting = load(filename)
        hists = load(filename)
        if year == '2018':
            hists['data']['reweighting'].scale(2.)
        num=hists['data']['reweighting'].integrate('process').values()[()]
        den=hists['bkg']['reweighting'].integrate('process').values()[()]
        reweighting = num / np.maximum(den, 1.)
        self.reweighting = lookup_tools.dense_lookup.dense_lookup(reweighting, [ax.edges() for ax in hists['data']['reweighting'].axes()[1:]])

    def weight(self, tau21, pt, eta):
        weight = self.reweighting(tau21, pt, eta)
        return weight

get_reweighting = {
    '2016': Reweighting('2016').weight,
    '2017': Reweighting('2017').weight,
    '2018': Reweighting('2018').weight,
}
'''
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
'''
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
    'get_doublebtag_weight':    get_doublebtag_weight,
    'get_reweighting':          get_reweighting,
    'get_mu_rochester_sf':      get_mu_rochester_sf,
    #'Jetevaluator':             Jetevaluator,
}
save(corrections, 'data/corrections.coffea')



 
