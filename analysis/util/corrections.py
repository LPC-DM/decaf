#!/usr/bin/env python
import uproot, uproot_methods
import numpy as np
import os
from coffea.arrays import Initialize
from coffea import hist, lookup_tools
from coffea.lookup_tools import extractor, dense_lookup
from coffea.util import save, load
from coffea.btag_tools import BTagScaleFactor


get_pu_weight = {}
get_pu_weight['2017'] = {}
get_pu_weight['2018'] = {}

pu = {}
pu["2018"] = uproot.open("data/pileup/puWeights_10x_56ifb.root")
pu["2017"] = uproot.open("data/pileup/puWeights_90x_41ifb.root")
pu["2016"] = uproot.open("data/pileup/puWeights_80x_37ifb.root")
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
met_trig["2016"] = uproot.open("data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")
met_trig["2017"] = uproot.open("data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")
met_trig["2018"] = uproot.open("data/trigger_eff/metTriggerEfficiency_recoil_monojet_TH1F.root")
for year in ['2016','2017','2018']:
    fmet_trig = met_trig[year]
    met_trig_corr = fmet_trig["hden_monojet_recoil_clone_passed"].values
    get_met_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(met_trig_corr, fmet_trig["hden_monojet_recoil_clone_passed"].edges)

get_met_zmm_trig_weight = {}

met_zmm_trig = {}
met_zmm_trig["2016"] = uproot.open("data/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root")
met_zmm_trig["2017"] = uproot.open("data/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root")
met_zmm_trig["2018"] = uproot.open("data/trigger_eff/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root")
for year in ['2016','2017','2018']:
    fmet_zmm_trig = met_zmm_trig[year]
    met_zmm_trig_corr = fmet_zmm_trig["hden_monojet_recoil_clone_passed"].values
    get_met_zmm_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(met_zmm_trig_corr, fmet_zmm_trig["hden_monojet_recoil_clone_passed"].edges)


get_ele_trig_weight = {}

ele_trig = {}
ele_trig["2016"] = uproot.open("data/trigger_eff/eleTrig.root")
ele_trig["2017"] = uproot.open("data/trigger_eff/eleTrig.root")
ele_trig["2018"] = uproot.open("data/trigger_eff/eleTrig.root")

for year in ['2016','2017','2018']:
    fele_trig = ele_trig[year]
    ele_trig_corr = fele_trig["hEffEtaPt"].values
    get_ele_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(ele_trig_corr, fele_trig["hEffEtaPt"].edges)

get_pho_trig_weight = {}

pho_trig = {}
pho_trig["2016"] = uproot.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
pho_trig["2017"] = uproot.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
pho_trig["2018"] = uproot.open("data/trigger_eff/photonTriggerEfficiency_photon_TH1F.root")
for year in ['2016','2017','2018']:
    fpho_trig = pho_trig[year]
    pho_trig_corr = fpho_trig["hden_photonpt_clone_passed"].values
    get_pho_trig_weight[year] = lookup_tools.dense_lookup.dense_lookup(pho_trig_corr, fpho_trig["hden_photonpt_clone_passed"].edges)

get_ele_loose_id_sf = {}
get_ele_tight_id_sf = {}
get_ele_loose_id_eff = {}
get_ele_tight_id_eff = {}

ele_loose = {}
ele_loose['2016'] = uproot.open("data/ScaleFactor/2016LegacyReReco_ElectronLoose_Fall17V2.root")
ele_loose['2017'] = uproot.open("data/ScaleFactor/2017_ElectronLoose.root")
ele_loose['2018'] = uproot.open("data/ScaleFactor/2018_ElectronLoose.root")
ele_tight = {}
ele_tight['2016'] = uproot.open("data/ScaleFactor/2016LegacyReReco_ElectronTight_Fall17V2.root")
ele_tight['2017'] = uproot.open("data/ScaleFactor/2017_ElectronTight.root")
ele_tight['2018'] = uproot.open("data/ScaleFactor/2018_ElectronTight.root")
for year in ['2016','2017','2018']:
    fele_loose = ele_loose[year]
    get_ele_loose_id_sf[year]  = lookup_tools.dense_lookup.dense_lookup(fele_loose["EGamma_SF2D"].values, fele_loose["EGamma_SF2D"].edges)
    get_ele_loose_id_eff[year] = lookup_tools.dense_lookup.dense_lookup(fele_loose["EGamma_EffMC2D"].values, fele_loose["EGamma_EffMC2D"].edges)
    fele_tight  = ele_tight[year]
    get_ele_tight_id_sf[year]  = lookup_tools.dense_lookup.dense_lookup(fele_tight["EGamma_SF2D"].values, fele_tight["EGamma_SF2D"].edges)
    get_ele_tight_id_eff[year] = lookup_tools.dense_lookup.dense_lookup(fele_tight["EGamma_EffMC2D"].values, fele_tight["EGamma_EffMC2D"].edges)

get_pho_tight_id_sf = {}

pho_tight = {}
pho_tight['2016'] = uproot.open("data/ScaleFactor/Fall17V2_2016_Tight_photons.root")
pho_tight['2017'] = uproot.open("data/ScaleFactor/2017_PhotonsTight.root")
pho_tight['2018'] = uproot.open("data/ScaleFactor/2018_PhotonsTight.root")
for year in ['2016','2017','2018']:
    fpho_tight=pho_tight[year]
    get_pho_tight_id_sf[year] = lookup_tools.dense_lookup.dense_lookup(fpho_tight["EGamma_SF2D"].values, fpho_tight["EGamma_SF2D"].edges)

get_mu_tight_id_sf = {}
get_mu_loose_id_sf = {}

mu_id2016 = uproot.open("data/ScaleFactor/2016LegacyReReco_Muon_SF_ID.root")
get_mu_tight_id_sf['2016'] = lookup_tools.dense_lookup.dense_lookup(mu_id2016["NUM_TightID_DEN_genTracks_eta_pt"].values, mu_id2016["NUM_TightID_DEN_genTracks_eta_pt"].edges)
get_mu_loose_id_sf['2016'] = lookup_tools.dense_lookup.dense_lookup(mu_id2016["NUM_LooseID_DEN_genTracks_eta_pt"].values, mu_id2016["NUM_LooseID_DEN_genTracks_eta_pt"].edges)
mu_id2017 = uproot.open("data/ScaleFactor/2017_Muon_RunBCDEF_SF_ID.root")
get_mu_tight_id_sf['2017'] = lookup_tools.dense_lookup.dense_lookup(mu_id2017["NUM_TightID_DEN_genTracks_pt_abseta"].values, mu_id2017["NUM_TightID_DEN_genTracks_pt_abseta"].edges)
get_mu_loose_id_sf['2017'] = lookup_tools.dense_lookup.dense_lookup(mu_id2017["NUM_LooseID_DEN_genTracks_pt_abseta"].values, mu_id2017["NUM_LooseID_DEN_genTracks_pt_abseta"].edges)
mu_id2018 = uproot.open("data/ScaleFactor/2018_Muon_RunABCD_SF_ID.root")
get_mu_tight_id_sf['2018'] = lookup_tools.dense_lookup.dense_lookup(mu_id2018["NUM_TightID_DEN_TrackerMuons_pt_abseta"].values, mu_id2018["NUM_TightID_DEN_TrackerMuons_pt_abseta"].edges)
get_mu_loose_id_sf['2018'] = lookup_tools.dense_lookup.dense_lookup(mu_id2018["NUM_LooseID_DEN_TrackerMuons_pt_abseta"].values, mu_id2018["NUM_LooseID_DEN_TrackerMuons_pt_abseta"].edges)

get_ele_reco_sf = {}

ele_reco = {}
ele_reco['2016']=uproot.open("data/ScaleFactor/2016_ElectronReco.root")    
ele_reco['2017']=uproot.open("data/ScaleFactor/2017_ElectronReco.root")
ele_reco['2018']=uproot.open("data/ScaleFactor/2018_ElectronReco.root")
for year in ['2016','2017','2018']:
    fele_reco = ele_reco[year]
    get_ele_reco_sf[year]=lookup_tools.dense_lookup.dense_lookup(fele_reco["EGamma_SF2D"].values, fele_reco["EGamma_SF2D"].edges)

get_mu_tight_iso_sf = {}
get_mu_loose_iso_sf = {}

mu_iso2016=uproot.open("data/ScaleFactor/Merged_SF_ISO.root")
mu_iso2017=uproot.open("data/ScaleFactor/RunBCDEF_SF_ISO_syst.root")
mu_iso2018=uproot.open("data/ScaleFactor/RunABCD_SF_ISO.root")
get_mu_tight_iso_sf['2016'] = lookup_tools.dense_lookup.dense_lookup(mu_iso2016["NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"].values,mu_iso2016["NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"].edges)
get_mu_loose_iso_sf['2016'] = lookup_tools.dense_lookup.dense_lookup(mu_iso2016["NUM_LooseRelIso_DEN_LooseID_eta_pt"].values,mu_iso2016["NUM_LooseRelIso_DEN_LooseID_eta_pt"].edges)
get_mu_tight_iso_sf['2017'] = lookup_tools.dense_lookup.dense_lookup(mu_iso2017["NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"].values,mu_iso2017["NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"].edges)
get_mu_loose_iso_sf['2017'] = lookup_tools.dense_lookup.dense_lookup(mu_iso2017["NUM_LooseRelIso_DEN_LooseID_pt_abseta"].values,mu_iso2017["NUM_LooseRelIso_DEN_LooseID_pt_abseta"].edges)
get_mu_tight_iso_sf['2018'] = lookup_tools.dense_lookup.dense_lookup(mu_iso2018["NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"].values,mu_iso2018["NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"].edges)
get_mu_loose_iso_sf['2018'] = lookup_tools.dense_lookup.dense_lookup(mu_iso2018["NUM_LooseRelIso_DEN_LooseID_pt_abseta"].values,mu_iso2018["NUM_LooseRelIso_DEN_LooseID_pt_abseta"].edges)

get_nlo_weight = {}
kfactor = uproot.open("data/nlo/kfactors.root")
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
kfactor = uproot.open("data/nlo/2017_gen_v_pt_stat1_qcd_sf.root")
get_adhoc_weight['z']=lookup_tools.dense_lookup.dense_lookup(kfactor["dy_monojet"].values, kfactor["dy_monojet"].edges)
get_adhoc_weight['w']=lookup_tools.dense_lookup.dense_lookup(kfactor["wjet_monojet"].values, kfactor["wjet_monojet"].edges)

get_nnlo_weight = {}
kfactor = uproot.open("data/nnlo/lindert_qcd_nnlo_sf.root")
get_nnlo_weight['dy'] = lookup_tools.dense_lookup.dense_lookup(kfactor["eej"].values, kfactor["eej"].edges)
get_nnlo_weight['w'] = lookup_tools.dense_lookup.dense_lookup(kfactor["evj"].values, kfactor["evj"].edges)
get_nnlo_weight['z'] = lookup_tools.dense_lookup.dense_lookup(kfactor["vvj"].values, kfactor["vvj"].edges)


get_nnlo_nlo_weight = {}
kfactor_eej = uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_eej_madgraph_"+year+".root")
kfactor_evj = uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_evj_madgraph_"+year+".root")
kfactor_vvj = uproot.open("data/Vboson_Pt_Reweighting/"+year+"/TheoryXS_vvj_madgraph_"+year+".root")
get_nnlo_nlo_weight['dy']=lookup_tools.dense_lookup.dense_lookup(kfactor_eej['eej_NNLO_NLO_nnn_nnn_n'].values, kfactor_eej['eej_NNLO_NLO_nnn_nnn_n'].edges)
get_nnlo_nlo_weight['w']=lookup_tools.dense_lookup.dense_lookup(kfactor_evj['evj_NNLO_NLO_nnn_nnn_n'].values, kfactor_evj['evj_NNLO_NLO_nnn_nnn_n'].edges)
get_nnlo_nlo_weight['z']=lookup_tools.dense_lookup.dense_lookup(kfactor_vvj['vvj_NNLO_NLO_nnn_nnn_n'].values, kfactor_vvj['vvj_NNLO_NLO_nnn_nnn_n'].edges)

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
                '2016': 'DeepJet_102XSF_V1.csv',
                '2017': 'DeepJet_102XSF_V1.csv',
                '2018': 'DeepJet_102XSF_V1.csv',
                },
            'deepcsv': {
                '2016': 'DeepCSV_102XSF_V1.csv',
                '2017': 'DeepCSV_102XSF_V1.csv',
                '2018': 'DeepCSV_102XSF_V1.csv',
                }
        }
        filename = 'data/'+files[tagger][year]
        self.sf = BTagScaleFactor(filename, workingpoint)
        files = {
            '2016': 'btag2018.merged',
            '2017': 'btag2018.merged',
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
        def zerotag(eff, sf):
            eff_data = np.minimum(1, sf*eff)
            return ((1 - eff_data) / (1 - eff)).prod()

        eff = self.eff(flavor, pt, abseta)
        sf_nom = self.sf.eval('central', flavor, abseta, pt)
        sf_up = self.sf.eval('up', flavor, abseta, pt)
        sf_down = self.sf.eval('down', flavor, abseta, pt)

        nom = zerotag(eff, sf_nom)
        up = zerotag(eff, sf_up)
        down = zerotag(eff, sf_down)
        if '-1' in tag: 
            nom = np.maximum(1 - zerotag(eff, sf_nom), 0)
            up = np.maximum(1 - zerotag(eff, sf_up), 0)
            down = np.maximum(1 - zerotag(eff, sf_down), 0)
        return nom, up, down

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

corrections = {}
corrections['get_msd_weight']          = get_msd_weight
corrections['get_ttbar_weight']        = get_ttbar_weight
corrections['get_nlo_weight']          = get_nlo_weight
corrections['get_nnlo_weight']         = get_nnlo_weight
corrections['get_nnlo_nlo_weight']     = get_nnlo_nlo_weight
corrections['get_adhoc_weight']        = get_adhoc_weight
corrections['get_pu_weight']           = get_pu_weight
corrections['get_met_trig_weight']     = get_met_trig_weight
corrections['get_met_zmm_trig_weight'] = get_met_zmm_trig_weight
corrections['get_ele_trig_weight']     = get_ele_trig_weight
corrections['get_pho_trig_weight']     = get_pho_trig_weight
corrections['get_ele_loose_id_sf']     = get_ele_loose_id_sf
corrections['get_ele_tight_id_sf']     = get_ele_tight_id_sf
corrections['get_ele_loose_id_eff']    = get_ele_loose_id_eff
corrections['get_ele_tight_id_eff']    = get_ele_tight_id_eff
corrections['get_pho_tight_id_sf']     = get_pho_tight_id_sf 
corrections['get_mu_tight_id_sf']      = get_mu_tight_id_sf
corrections['get_mu_loose_id_sf']      = get_mu_loose_id_sf
corrections['get_ele_reco_sf']         = get_ele_reco_sf
corrections['get_mu_tight_iso_sf']     = get_mu_tight_iso_sf
corrections['get_mu_loose_iso_sf']     = get_mu_loose_iso_sf
corrections['get_ecal_bad_calib']      = get_ecal_bad_calib
corrections['get_btag_weight']         = get_btag_weight
corrections['Jetevaluator']            = Jetevaluator

save(corrections, 'data/corrections.coffea')



 
