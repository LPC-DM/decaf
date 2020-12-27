from __future__ import print_function, division
from optparse import OptionParser
from collections import defaultdict, OrderedDict
import concurrent.futures
import sys
import os
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import gzip
import json
from coffea import hist, processor
from coffea.util import load, save
import ROOT

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

mass_binning = [
    0,
    40,
    50,
    60,
    70,
    80,
    90,
    100,
    120,
    150,
    180,
    240,
    300,
]
# recoil_binning=[250,310,370,470,590,840,1020,1250,3000]
recoil_binning = [250, 310, 370, 470, 590, 3000]

category_map = {"pass": 1, "fail": 0}

with open("data/hf_systematic.json") as fin:
    hf_systematic = json.load(fin)

def template(dictionary, process, systematic, recoil, region, category):
    histogram = dictionary[region].integrate("process", process)
    nominal = histogram.integrate("systematic", "nominal").values()[()][
        recoil, :, category_map[category]
    ]
    output = nominal
    if "nominal" not in systematic and "data" not in systematic:
        # print('Normalizing',systematic,'histogram of',process,'in region',region)
        output = np.nan_to_num(
            histogram.integrate("systematic", systematic).values()[()][
                recoil, :, category_map[category]
            ]
            / nominal.sum()
        )
    if "data" not in systematic:
        # print('Remiving zeros from',systematic,'histogram of',process,'in region',region)
        output[output <= 0] = 1e-7
    binning = (
        dictionary[region]
        .integrate("process", process)
        .integrate("systematic", systematic)
        .axis("fjmass")
        .edges()
    )
    return (output, binning, "fjmass")

def remap_histograms(hists):
    data_hists = {}
    bkg_hists = {}
    signal_hists = {}

    process = hist.Cat("process", "Process", sorting="placement")
    cats = ("process",)
    sig_map = OrderedDict()
    bkg_map = OrderedDict()
    data_map = OrderedDict()
    bkg_map["Hbb"] = ("Hbb*",)
    bkg_map["DY+jets"] = ("DY+*",)
    bkg_map["VV"] = (["WW", "WZ", "ZZ"],)
    bkg_map["ST"] = ("ST*",)
    bkg_map["TT"] = ("TT*",)
    bkg_map["W+jets"] = ("W+*",)
    bkg_map["Z+jets"] = ("Z+*",)
    bkg_map["G+jets"] = ("G+*",)
    bkg_map["QCD"] = ("QCD*",)
    sig_map["Mhs_50"] = ("*Mhs_50*",)  ## signals
    sig_map["Mhs_70"] = ("*Mhs_70*",)
    sig_map["Mhs_90"] = ("*Mhs_90*",)
    sig_map["MonoJet"] = ("MonoJet*",)  ## signals
    sig_map["MonoW"] = ("MonoW*",)  ## signals
    sig_map["MonoZ"] = ("MonoZ*",)  ## signals
    data_map["MET"] = ("MET*",)
    data_map["SingleElectron"] = ("SingleElectron*",)
    data_map["SinglePhoton"] = ("SinglePhoton*",)
    data_map["EGamma"] = ("EGamma*",)

    for key in hists["data"].keys():
        bkg_hists[key] = hists["bkg"][key].group(cats, process, bkg_map)
        signal_hists[key] = hists["sig"][key].group(cats, process, sig_map)
        data_hists[key] = hists["data"][key].group(cats, process, data_map)

    bkg_hists["template"] = bkg_hists["template"].rebin(
        "fjmass", hist.Bin("fjmass", "Mass", mass_binning)
    )
    signal_hists["template"] = signal_hists["template"].rebin(
        "fjmass", hist.Bin("fjmass", "Mass", mass_binning)
    )
    data_hists["template"] = data_hists["template"].rebin(
        "fjmass", hist.Bin("fjmass", "Mass", mass_binning)
    )

    bkg_hists["template"] = bkg_hists["template"].rebin(
        "recoil", hist.Bin("recoil", "Recoil", recoil_binning)
    )
    signal_hists["template"] = signal_hists["template"].rebin(
        "recoil", hist.Bin("recoil", "Recoil", recoil_binning)
    )
    data_hists["template"] = data_hists["template"].rebin(
        "recoil", hist.Bin("recoil", "Recoil", recoil_binning)
    )

    hists = {"bkg": bkg_hists, "sig": signal_hists, "data": data_hists}

    return hists


def initialize_nuisances(hists, year):

    ###
    # Let's start from TT
    ###

    ###
    # First, tagging efficiency and SF, DEPRECATED
    ###

    '''
    tt_efficiency = {"2018": 0.5}A
    sf_tt = rl.IndependentParameter(
        "sf_tt" + year, 1.0, 0.01, 1.0 / tt_efficiency[year]
    )

    tt_weight = {
        "pass": sf_tt * tt_efficiency[year],
        "fail": 1 - (sf_tt * tt_efficiency[year]),
    }
    '''
    sr_tt = (
        hists["bkg"]["template"]
        .integrate("region", "sr")
        .integrate("process", "TT")
        .integrate("systematic", "nominal")
    )

    sr_ttNuisances={}
    recoilbins = np.array(recoil_binning)
    nrecoil = len(recoilbins) - 1
    for recoilbin in range(nrecoil):
        sr_ttPass = sr_tt.sum("gentype").values()[()][
            recoilbin, :, 1
                  ]
        print(sr_ttPass)
        sr_ttNuisances[recoilbin] = np.array(  # one nuisance per mass shape bin in pass
            [
                rl.IndependentParameter(
                    "sr" + year + "_tt_pass_recoil"+str(recoilbin)+"_mass%d" % i,
                    b,
                    0,
                    sr_ttPass.max() * 2,
                    )
                for i, b in enumerate(sr_ttPass)
                ]
            ),
            

    ###
    # Let's move to V+jets
    ###

    sr_zjets = (
        hists["bkg"]["template"]
        .integrate("region", "sr")
        .integrate("process", "Z+jets")
        .integrate("systematic", "nominal")
    )

    ###
    # First, the mass shape
    ###

    sr_zjetsMass = sr_zjets.sum("gentype", "recoil")
    sr_zjetsMassFail = sr_zjetsMass.values()[()][
        :, 0
    ]  # get the fail histogram, inclusive in recoil
    sr_zjetsMassFail = (
        sr_zjetsMassFail / sr_zjetsMassFail.sum()
    )  # normalize to the integral to get the shape in fail
    sr_zjetsShape = np.array(
        [
            rl.IndependentParameter(
                "sr" + year + "_zjshape_fail_mass%d" % i,
                b,
                0,
                sr_zjetsMassFail.max() * 2,
            )
            for i, b in enumerate(sr_zjetsMassFail)
        ]
    )

    ###
    # Then, recoil rate
    ###

    sr_zjetsRecoil = sr_zjets.sum("gentype", "fjmass")
    sr_zjetsRecoilFail = sr_zjetsRecoil.values()[()][:, 0]
    sr_zjetsRate = np.array(
        [
            rl.IndependentParameter(
                "sr" + year + "_zj_fail_recoil%d" % i,
                b,
                0,
                sr_zjetsRecoilFail.max() * 2,
            )
            for i, b in enumerate(sr_zjetsRecoilFail)
        ]
    )
    return sr_zjetsShape, sr_zjetsRate, sr_ttNuisances


def computeTFs(hists, year, recoil, category):

    model_id = year + category + "recoil" + str(recoil)

    bkg_hists = hists["bkg"]
    background = {}
    for r in bkg_hists["template"].identifiers("region"):
        background[str(r)] = bkg_hists["template"].integrate("region", r).sum("gentype")

    def addBtagSyst(process, region, templ):
        btagUp = template(background, process, "btagUp", recoil, region, category)[0]
        btagDown = template(background, process, "btagDown", recoil, region, category)[0]
        templ.setParamEffect(btag, btagUp, btagDown)

    def addVJetsSyst(process, region, templ):
        ew1Up = template(background, process, "ew1Up", recoil, region, category)[0]
        ew1Down = template(background, process, "ew1Down", recoil, region, category)[0]
        templ.setParamEffect(ew1, ew1Up, ew1Down)
        ew2GUp = template(background, process, "ew2GUp", recoil, region, category)[0]
        ew2GDown = template(background, process, "ew2GDown", recoil, region, category)[0]
        templ.setParamEffect(ew2G, ew2GUp, ew2GDown)
        ew2WUp = template(background, process, "ew2WUp", recoil, region, category)[0]
        ew2WDown = template(background, process, "ew2WDown", recoil, region, category)[0]
        templ.setParamEffect(ew2W, ew2WUp, ew2WDown)
        ew2ZUp = template(background, process, "ew2ZUp", recoil, region, category)[0]
        ew2ZDown = template(background, process, "ew2ZDown", recoil, region, category)[0]
        templ.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
        ew3GUp = template(background, process, "ew3GUp", recoil, region, category)[0]
        ew3GDown = template(background, process, "ew3GDown", recoil, region, category)[0]
        templ.setParamEffect(ew3G, ew3GUp, ew3GDown)
        ew3WUp = template(background, process, "ew3WUp", recoil, region, category)[0]
        ew3WDown = template(background, process, "ew3WDown", recoil, region, category)[0]
        templ.setParamEffect(ew3W, ew3WUp, ew3WDown)
        ew3ZUp = template(background, process, "ew3ZUp", recoil, region, category)[0]
        ew3ZDown = template(background, process, "ew3ZDown", recoil, region, category)[0]
        templ.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
        mixUp = template(background, process, "mixUp", recoil, region, category)[0]
        mixDown = template(background, process, "mixDown", recoil, region, category)[0]
        templ.setParamEffect(mix, mixUp, mixDown)
        muFUp = template(background, process, "muFUp", recoil, region, category)[0]
        muFDown = template(background, process, "muFDown", recoil, region, category)[0]
        templ.setParamEffect(muF, muFUp, muFDown)
        muRUp = template(background, process, "muRUp", recoil, region, category)[0]
        muRDown = template(background, process, "muRDown", recoil, region, category)[0]
        templ.setParamEffect(muR, muRUp, muRDown)
        qcd1Up = template(background, process, "qcd1Up", recoil, region, category)[0]
        qcd1Down = template(background, process, "qcd1Down", recoil, region, category)[0]
        templ.setParamEffect(qcd1, qcd1Up, qcd1Down)
        qcd2Up = template(background, process, "qcd2Up", recoil, region, category)[0]
        qcd2Down = template(background, process, "qcd2Down", recoil, region, category)[0]
        templ.setParamEffect(qcd2, qcd2Up, qcd2Down)
        qcd3Up = template(background, process, "qcd3Up", recoil, region, category)[0]
        qcd3Down = template(background, process, "qcd3Down", recoil, region, category)[0]
        templ.setParamEffect(qcd3, qcd3Up, qcd3Down)


    ###
    # Z+jets template
    ###

    sr_zjets = rl.TemplateSample(
        "sr" + model_id + "_zjets",
        rl.Sample.BACKGROUND,
        template(background, "Z+jets", "nominal", recoil, "sr", category),
    )
    sr_zjets.setParamEffect(lumi, 1.027)
    sr_zjets.setParamEffect(zjets_norm, 1.4)
    sr_zjets.setParamEffect(trig_met, 1.01)
    sr_zjets.setParamEffect(veto_tau, 1.03)
    sr_zjets.setParamEffect(jec, 1.05)
    sr_zjets.setParamEffect(zhf_fraction, hf_systematic["Z+jets"]["sr"][category])
    addBtagSyst("Z+jets", "sr", sr_zjets)
    addVJetsSyst("Z+jets", "sr", sr_zjets)

    ###
    # W+jets templates
    ###
    sr_wjets = rl.TemplateSample(
        "sr" + model_id + "_wjets",
        rl.Sample.BACKGROUND,
        template(background, "W+jets", "nominal", recoil, "sr", category),
    )
    sr_wjets.setParamEffect(lumi, 1.027)
    sr_wjets.setParamEffect(wjets_norm, 1.4)
    sr_wjets.setParamEffect(trig_met, 1.01)
    sr_wjets.setParamEffect(veto_tau, 1.03)
    sr_wjets.setParamEffect(jec, 1.05)
    sr_wjets.setParamEffect(whf_fraction, hf_systematic["W+jets"]["sr"][category])
    addBtagSyst("W+jets", "sr", sr_wjets)
    addVJetsSyst("W+jets", "sr", sr_wjets)

    wmcr_wjets = rl.TemplateSample(
        "wmcr" + model_id + "_wjets",
        rl.Sample.BACKGROUND,
        template(background, "W+jets", "nominal", recoil, "wmcr", category),
    )
    wmcr_wjets.setParamEffect(lumi, 1.027)
    wmcr_wjets.setParamEffect(trig_met, 1.01)
    wmcr_wjets.setParamEffect(veto_tau, 1.03)
    wmcr_wjets.setParamEffect(wjets_norm, 1.4)
    wmcr_wjets.setParamEffect(jec, 1.05)
    wmcr_wjets.setParamEffect(id_mu, 1.02)
    wmcr_wjets.setParamEffect(iso_mu, 1.02)
    wmcr_wjets.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["wmcr"][category]
    )
    addBtagSyst("W+jets", "wmcr", wmcr_wjets)
    addVJetsSyst("W+jets", "wmcr", wmcr_wjets)

    wecr_wjets = rl.TemplateSample(
        "wecr" + model_id + "_wjets",
        rl.Sample.BACKGROUND,
        template(background, "W+jets", "nominal", recoil, "wecr", category),
    )
    wecr_wjets.setParamEffect(lumi, 1.027)
    wecr_wjets.setParamEffect(trig_e, 1.01)
    wecr_wjets.setParamEffect(veto_tau, 1.03)
    wecr_wjets.setParamEffect(wjets_norm, 1.4)
    wecr_wjets.setParamEffect(jec, 1.05)
    wecr_wjets.setParamEffect(id_e, 1.02)
    wecr_wjets.setParamEffect(reco_e, 1.02)
    wecr_wjets.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["wecr"][category]
    )
    addBtagSyst("W+jets", "wecr", wecr_wjets)
    addVJetsSyst("W+jets", "wecr", wecr_wjets)

    tmcr_wjets = rl.TemplateSample(
        "tmcr" + model_id + "_wjets",
        rl.Sample.BACKGROUND,
        template(background, "W+jets", "nominal", recoil, "tmcr", category),
    )
    tmcr_wjets.setParamEffect(lumi, 1.027)
    tmcr_wjets.setParamEffect(trig_met, 1.01)
    tmcr_wjets.setParamEffect(veto_tau, 1.03)
    tmcr_wjets.setParamEffect(wjets_norm, 1.4)
    tmcr_wjets.setParamEffect(jec, 1.05)
    tmcr_wjets.setParamEffect(id_mu, 1.02)
    tmcr_wjets.setParamEffect(iso_mu, 1.02)
    tmcr_wjets.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["tmcr"][category]
    )
    addBtagSyst("W+jets", "tmcr", tmcr_wjets)
    addVJetsSyst("W+jets", "tmcr", tmcr_wjets)

    tecr_wjets = rl.TemplateSample(
        "tecr" + model_id + "_wjets",
        rl.Sample.BACKGROUND,
        template(background, "W+jets", "nominal", recoil, "tecr", category),
    )
    tecr_wjets.setParamEffect(lumi, 1.027)
    tecr_wjets.setParamEffect(trig_e, 1.01)
    tecr_wjets.setParamEffect(veto_tau, 1.03)
    tecr_wjets.setParamEffect(wjets_norm, 1.4)
    tecr_wjets.setParamEffect(jec, 1.05)
    tecr_wjets.setParamEffect(id_e, 1.02)
    tecr_wjets.setParamEffect(reco_e, 1.02)
    tecr_wjets.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["tecr"][category]
    )
    addBtagSyst("W+jets", "tecr", tecr_wjets)
    addVJetsSyst("W+jets", "tecr", tecr_wjets)

    ###
    # TT templates
    ###

    sr_tt = rl.TemplateSample(
        "sr" + model_id + "_tt",
        rl.Sample.BACKGROUND,
        template(background, "TT", "nominal", recoil, "sr", category),
    )
    sr_tt.setParamEffect(lumi, 1.027)
    sr_tt.setParamEffect(tt_norm, 1.2)
    sr_tt.setParamEffect(trig_met, 1.01)
    sr_tt.setParamEffect(veto_tau, 1.03)
    sr_tt.setParamEffect(jec, 1.05)
    addBtagSyst("TT", "sr", sr_tt)

    wmcr_tt = rl.TemplateSample(
        "wmcr" + model_id + "_tt",
        rl.Sample.BACKGROUND,
        template(background, "TT", "nominal", recoil, "wmcr", category),
    )
    wmcr_tt.setParamEffect(lumi, 1.027)
    wmcr_tt.setParamEffect(trig_met, 1.01)
    wmcr_tt.setParamEffect(veto_tau, 1.03)
    wmcr_tt.setParamEffect(tt_norm, 1.2)
    wmcr_tt.setParamEffect(jec, 1.05)
    wmcr_tt.setParamEffect(id_mu, 1.02)
    wmcr_tt.setParamEffect(iso_mu, 1.02)
    addBtagSyst("TT", "wmcr", wmcr_tt)

    wecr_tt = rl.TemplateSample(
        "wecr" + model_id + "_tt",
        rl.Sample.BACKGROUND,
        template(background, "TT", "nominal", recoil, "wecr", category),
    )
    wecr_tt.setParamEffect(lumi, 1.027)
    wecr_tt.setParamEffect(trig_e, 1.01)
    wecr_tt.setParamEffect(veto_tau, 1.03)
    wecr_tt.setParamEffect(tt_norm, 1.2)
    wecr_tt.setParamEffect(jec, 1.05)
    wecr_tt.setParamEffect(id_e, 1.02)
    wecr_tt.setParamEffect(reco_e, 1.02)
    addBtagSyst("TT", "wecr", wecr_tt)

    tmcr_tt = rl.TemplateSample(
        "tmcr" + model_id + "_tt",
        rl.Sample.BACKGROUND,
        template(background, "TT", "nominal", recoil, "tmcr", category),
    )
    tmcr_tt.setParamEffect(lumi, 1.027)
    tmcr_tt.setParamEffect(trig_met, 1.01)
    tmcr_tt.setParamEffect(veto_tau, 1.03)
    tmcr_tt.setParamEffect(tt_norm, 1.2)
    tmcr_tt.setParamEffect(jec, 1.05)
    tmcr_tt.setParamEffect(id_mu, 1.02)
    tmcr_tt.setParamEffect(iso_mu, 1.02)
    addBtagSyst("TT", "tmcr", tmcr_tt) 

    tecr_tt = rl.TemplateSample(
        "tecr" + model_id + "_tt",
        rl.Sample.BACKGROUND,
        template(background, "TT", "nominal", recoil, "tecr", category),
    )
    tecr_tt.setParamEffect(lumi, 1.027)
    tecr_tt.setParamEffect(trig_e, 1.01)
    tecr_tt.setParamEffect(veto_tau, 1.03)
    tecr_tt.setParamEffect(tt_norm, 1.2)
    tecr_tt.setParamEffect(jec, 1.05)
    tecr_tt.setParamEffect(id_e, 1.02)
    tecr_tt.setParamEffect(reco_e, 1.02)
    addBtagSyst("TT", "tecr", tecr_tt)

    ###
    # DY+jets templates
    ###

    zmcr_dyjets = rl.TemplateSample(
        "zmcr" + model_id + "_dyjets",
        rl.Sample.BACKGROUND,
        template(background, "DY+jets", "nominal", recoil, "zmcr", category),
    )
    zmcr_dyjets.setParamEffect(lumi, 1.027)
    zmcr_dyjets.setParamEffect(trig_met, 1.01)
    zmcr_dyjets.setParamEffect(veto_tau, 1.03)
    zmcr_dyjets.setParamEffect(zjets_norm, 1.4)
    zmcr_dyjets.setParamEffect(jec, 1.05)
    zmcr_dyjets.setParamEffect(id_mu, 1.02)
    zmcr_dyjets.setParamEffect(iso_mu, 1.02)
    zmcr_dyjets.setParamEffect(
        zhf_fraction, hf_systematic["DY+jets"]["zmcr"][category]
    )
    addVJetsSyst("DY+jets", "zmcr", zmcr_dyjets)

    zecr_dyjets = rl.TemplateSample(
        "zecr" + model_id + "_dyjets",
        rl.Sample.BACKGROUND,
        template(background, "DY+jets", "nominal", recoil, "zecr", category),
    )
    zecr_dyjets.setParamEffect(lumi, 1.027)
    zecr_dyjets.setParamEffect(trig_e, 1.01)
    zecr_dyjets.setParamEffect(veto_tau, 1.03)
    zecr_dyjets.setParamEffect(zjets_norm, 1.4)
    zecr_dyjets.setParamEffect(jec, 1.05)
    zecr_dyjets.setParamEffect(id_e, 1.02)
    zecr_dyjets.setParamEffect(reco_e, 1.02)
    zecr_dyjets.setParamEffect(
        zhf_fraction, hf_systematic["DY+jets"]["zecr"][category]
    )
    addVJetsSyst("DY+jets", "zecr", zecr_dyjets)
    ###
    # G+jets templates
    ###

    gcr_gjets = rl.TemplateSample(
        "gcr" + model_id + "_gjets",
        rl.Sample.BACKGROUND,
        template(background, "G+jets", "nominal", recoil, "gcr", category),
    )
    gcr_gjets.setParamEffect(lumi, 1.027)
    gcr_gjets.setParamEffect(trig_pho, 1.01)
    gcr_gjets.setParamEffect(veto_tau, 1.03)
    gcr_gjets.setParamEffect(gjets_norm, 1.4)
    gcr_gjets.setParamEffect(jec, 1.05)
    gcr_gjets.setParamEffect(id_pho, 1.02)
    gcr_gjets.setParamEffect(ghf_fraction, hf_systematic["G+jets"]["gcr"][category])
    addVJetsSyst("G+jets", "gcr", gcr_gjets)
    ###
    # Compute TFs
    ###

    sr_wjetsTransferFactor = sr_wjets.getExpectation() / sr_zjets.getExpectation()
    wmcr_wjetsTransferFactor = wmcr_wjets.getExpectation() / sr_wjets.getExpectation()
    wmcr_ttTransferFactor = wmcr_tt.getExpectation() / sr_tt.getExpectation()
    tmcr_wjetsTransferFactor = tmcr_wjets.getExpectation() / sr_wjets.getExpectation()
    tmcr_ttTransferFactor = tmcr_tt.getExpectation() / sr_tt.getExpectation()
    wecr_wjetsTransferFactor = wecr_wjets.getExpectation() / sr_wjets.getExpectation()
    wecr_ttTransferFactor = wecr_tt.getExpectation() / sr_tt.getExpectation()
    tecr_wjetsTransferFactor = tecr_wjets.getExpectation() / sr_wjets.getExpectation()
    tecr_ttTransferFactor = tecr_tt.getExpectation() / sr_tt.getExpectation()
    zmcr_dyjetsTransferFactor = zmcr_dyjets.getExpectation() / sr_zjets.getExpectation()
    zecr_dyjetsTransferFactor = zecr_dyjets.getExpectation() / sr_zjets.getExpectation()
    gcr_gjetsTransferFactor = gcr_gjets.getExpectation() / sr_zjets.getExpectation()

    return (
        sr_wjetsTransferFactor,
        wmcr_wjetsTransferFactor,
        wmcr_ttTransferFactor,
        tmcr_wjetsTransferFactor,
        tmcr_ttTransferFactor,
        wecr_wjetsTransferFactor,
        wecr_ttTransferFactor,
        tecr_wjetsTransferFactor,
        tecr_ttTransferFactor,
        zmcr_dyjetsTransferFactor,
        zecr_dyjetsTransferFactor,
        gcr_gjetsTransferFactor,
    )


def rhalphabeth(msdbins):

    process = hist.Cat("process", "Process", sorting="placement")
    cats = ("process",)
    bkg_map = OrderedDict()
    # bkg_map['V+jets'] = (['Z+jets','W+jets'],)
    bkg_map["V+jets"] = (["Z+jets"],)
    vjets_hists = {}
    for key in hists["data"].keys():
        vjets_hists[key] = hists["bkg"][key].group(cats, process, bkg_map)

    # Build qcd MC pass+fail model and fit to polynomial
    qcdmodel = rl.Model("qcdmodel")
    qcdpass, qcdfail = 0.0, 0.0
    msds = np.meshgrid(msdbins[:-1] + 0.5 * np.diff(msdbins), indexing="ij")[0]
    msds = np.sqrt(msds) * np.sqrt(msds)
    print(msds)
    msdscaled = msds / 300.0
    msd = rl.Observable("fjmass", msdbins)
    failCh = rl.Channel("fail")
    passCh = rl.Channel("pass")
    qcdmodel.addChannel(failCh)
    qcdmodel.addChannel(passCh)
    # mock template
    ptnorm = 1
    vjetsHistFail = (
        vjets_hists["template"]
        .integrate("region", "sr")
        .sum("gentype", "recoil")
        .integrate("process", "V+jets")
        .integrate("systematic", "nominal")
        .values()[()][:, 0]
    )
    vjetsHistFail[vjetsHistFail <= 0] = 1e-7
    failTempl = (
        vjetsHistFail,
        vjets_hists["template"]
        .integrate("region", "sr")
        .sum("gentype", "recoil")
        .integrate("process", "V+jets")
        .integrate("systematic", "nominal")
        .axis("fjmass")
        .edges(),
        "fjmass",
    )
    vjetsHistPass = (
        vjets_hists["template"]
        .integrate("region", "sr")
        .sum("gentype", "recoil")
        .integrate("process", "V+jets")
        .integrate("systematic", "nominal")
        .values()[()][:, 1]
    )
    vjetsHistPass[vjetsHistPass <= 0] = 1e-7
    passTempl = (
        vjetsHistPass,
        vjets_hists["template"]
        .integrate("region", "sr")
        .sum("gentype", "recoil")
        .integrate("process", "V+jets")
        .integrate("systematic", "nominal")
        .axis("fjmass")
        .edges(),
        "fjmass",
    )
    failCh.setObservation(failTempl)
    passCh.setObservation(passTempl)
    qcdfail += failCh.getObservation().sum()
    qcdpass += passCh.getObservation().sum()

    qcdeff = qcdpass / qcdfail
    tf_MCtempl = rl.BernsteinPoly("tf_MCtempl", (2,), ["fjmass"])
    tf_MCtempl_params = qcdeff * tf_MCtempl(msdscaled)
    failCh = qcdmodel["fail"]
    passCh = qcdmodel["pass"]
    failObs = failCh.getObservation()
    qcdparams = np.array(
        [rl.IndependentParameter("qcdparam_msdbin%d" % i, 0) for i in range(msd.nbins)]
    )
    sigmascale = 10.0
    scaledparams = (
        failObs * (1 + sigmascale / np.maximum(1.0, np.sqrt(failObs))) ** qcdparams
    )
    fail_qcd = rl.ParametericSample("fail_qcd", rl.Sample.BACKGROUND, msd, scaledparams)
    failCh.addSample(fail_qcd)
    print(tf_MCtempl_params)
    pass_qcd = rl.TransferFactorSample(
        "pass_qcd", rl.Sample.BACKGROUND, tf_MCtempl_params, fail_qcd
    )
    passCh.addSample(pass_qcd)

    qcdfit_ws = ROOT.RooWorkspace("qcdfit_ws")
    simpdf, obs = qcdmodel.renderRoofit(qcdfit_ws)
    qcdfit = simpdf.fitTo(
        obs,
        ROOT.RooFit.Extended(True),
        ROOT.RooFit.SumW2Error(True),
        ROOT.RooFit.Strategy(2),
        ROOT.RooFit.Save(),
        ROOT.RooFit.Minimizer("Minuit2", "migrad"),
        ROOT.RooFit.PrintLevel(-1),
    )
    qcdfit_ws.add(qcdfit)
    if "pytest" not in sys.modules:
        qcdfit_ws.writeToFile(os.path.join(str("models"), "testModel_qcdfit.root"))
    if qcdfit.status() != 0:
        raise RuntimeError("Could not fit qcd")

    param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
    decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(
        tf_MCtempl.name + "_deco", qcdfit, param_names
    )
    tf_MCtempl.parameters = decoVector.correlated_params.reshape(
        tf_MCtempl.parameters.shape
    )
    tf_MCtempl_params_final = tf_MCtempl(msdscaled)
    tf_dataResidual = rl.BernsteinPoly(
        "tf_dataResidual", (2,), ["fjmass"], limits=(0, 10)
    )
    tf_dataResidual_params = tf_dataResidual(msdscaled)
    tf_params = qcdeff * tf_MCtempl_params_final * tf_dataResidual_params
    return tf_params


def model(year, recoil, category):

    model_id = year + category + "recoil" + str(recoil)
    print(model_id)
    model = rl.Model("darkhiggs" + model_id)

    data_hists = hists["data"]
    bkg_hists = hists["bkg"]
    signal_hists = hists["sig"]

    ###
    # Preparing histograms for fit
    ##

    data = {}
    for r in data_hists["template"].identifiers("region"):
        data[str(r)] = data_hists["template"].integrate("region", r).sum("gentype")

    background = {}
    for r in bkg_hists["template"].identifiers("region"):
        background[str(r)] = bkg_hists["template"].integrate("region", r).sum("gentype")

    signal = {}
    for r in bkg_hists["template"].identifiers("region"):
        signal[str(r)] = signal_hists["template"].integrate("region", r).sum("gentype")

    ###
    ###
    # Signal region
    ###
    ###

    ch_name = "sr" + model_id
    sr = rl.Channel(ch_name)
    model.addChannel(sr)

    ###
    # Add data distribution to the channel
    ###

    sr.setObservation(template(data, "MET", "data", recoil, "sr", category))

    ###
    # Z(->nunu)+jets data-driven model
    ###
    sr_zjetsTemplate = template(background, "Z+jets", "nominal", recoil, "sr", category)
    sr_zjetsObservable = rl.Observable("fjmass", sr_zjetsTemplate[1])
    if category == "pass":
        sr_zjets = rl.ParametericSample(
            ch_name + "_zjets",
            rl.Sample.BACKGROUND,
            sr_zjetsObservable,
            sr_zjetsBinYields * tf_params,
        )
    else:
        sr_zjets = rl.ParametericSample(
            ch_name + "_zjets",
            rl.Sample.BACKGROUND,
            sr_zjetsObservable,
            sr_zjetsBinYields * 1.0,
        )
    sr.addSample(sr_zjets)

    ###
    # W(->lnu)+jets data-driven model
    ###

    # Adding W-Z link
    sr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, sr_wjetsTransferFactor, sr_zjets
    )
    sr.addSample(sr_wjets)

    ###
    # top-antitop model
    ###

    sr_ttTemplate = template(background, "TT", "nominal", recoil, "sr", category)
    if category == "pass":
        sr_ttObservable = rl.Observable("fjmass", sr_ttTemplate[1])
        print(sr_ttTemplate[1])
        sr_tt = rl.ParametericSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields
            )
    else:
        sr_tt = rl.TemplateSample(ch_name + "_ttMC", rl.Sample.BACKGROUND, sr_ttTemplate)
        sr_tt.setParamEffect(lumi, 1.027)
        sr_tt.setParamEffect(trig_met, 1.01)
        sr_tt.setParamEffect(veto_tau, 1.03)
        sr_tt.setParamEffect(st_norm, 1.2)
        sr_tt.setParamEffect(jec, 1.05)
        btagUp = template(background, "TT", "btagUp", recoil, "sr", category)[0]
        btagDown = template(background, "TT", "btagDown", recoil, "sr", category)[0]
        sr_tt.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_tt)
        

    ###
    # Other MC-driven processes
    ###

    sr_stTemplate = template(background, "ST", "nominal", recoil, "sr", category)
    sr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, sr_stTemplate)
    sr_st.setParamEffect(lumi, 1.027)
    sr_st.setParamEffect(trig_met, 1.01)
    sr_st.setParamEffect(veto_tau, 1.03)
    sr_st.setParamEffect(st_norm, 1.2)
    sr_st.setParamEffect(jec, 1.05)
    btagUp = template(background, "ST", "btagUp", recoil, "sr", category)[0]
    btagDown = template(background, "ST", "btagDown", recoil, "sr", category)[0]
    sr_st.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_st)

    sr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "sr", category)
    sr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, sr_dyjetsTemplate
    )
    sr_dyjets.setParamEffect(lumi, 1.027)
    sr_dyjets.setParamEffect(trig_met, 1.01)
    sr_dyjets.setParamEffect(veto_tau, 1.03)
    sr_dyjets.setParamEffect(zjets_norm, 1.4)
    sr_dyjets.setParamEffect(jec, 1.05)
    btagUp = template(background, "DY+jets", "btagUp", recoil, "sr", category)[0]
    btagDown = template(background, "DY+jets", "btagDown", recoil, "sr", category)[0]
    btagDown[btagDown <= 0] = 1e-7
    sr_dyjets.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_dyjets)

    sr_vvTemplate = template(background, "VV", "nominal", recoil, "sr", category)
    sr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, sr_vvTemplate)
    sr_vv.setParamEffect(lumi, 1.027)
    sr_vv.setParamEffect(trig_met, 1.01)
    sr_vv.setParamEffect(veto_tau, 1.03)
    sr_vv.setParamEffect(vv_norm, 1.2)
    sr_vv.setParamEffect(jec, 1.05)
    btagUp = template(background, "VV", "btagUp", recoil, "sr", category)[0]
    btagDown = template(background, "VV", "btagDown", recoil, "sr", category)[0]
    btagDown[btagDown <= 0] = 1e-7
    sr_vv.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_vv)

    sr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "sr", category)
    sr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, sr_hbbTemplate)
    sr_hbb.setParamEffect(lumi, 1.027)
    sr_hbb.setParamEffect(trig_met, 1.01)
    sr_hbb.setParamEffect(veto_tau, 1.03)
    sr_hbb.setParamEffect(hbb_norm, 1.2)
    sr_hbb.setParamEffect(jec, 1.05)
    btagUp = template(background, "Hbb", "btagUp", recoil, "sr", category)[0]
    btagDown = template(background, "Hbb", "btagDown", recoil, "sr", category)[0]
    sr_hbb.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_hbb)

    sr_qcdTemplate = template(background, "QCD", "nominal", recoil, "sr", category)
    sr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, sr_qcdTemplate)
    sr_qcd.setParamEffect(lumi, 1.027)
    sr_qcd.setParamEffect(trig_met, 1.01)
    sr_qcd.setParamEffect(veto_tau, 1.03)
    sr_qcd.setParamEffect(qcdsig_norm, 2.0)
    sr_qcd.setParamEffect(jec, 1.05)
    btagUp = template(background, "QCD", "btagUp", recoil, "sr", category)[0]
    btagDown = template(background, "QCD", "btagDown", recoil, "sr", category)[0]
    sr_qcd.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_qcd)

    for s in signal["sr"].identifiers("process"):
        if "Mhs_50" not in str(s):
            continue
        sr_signalTemplate = template(signal, s, "nominal", recoil, "sr", category)
        sr_signal = rl.TemplateSample(
            ch_name + "_" + str(s), rl.Sample.SIGNAL, sr_signalTemplate
        )
        sr_signal.setParamEffect(lumi, 1.027)
        sr_signal.setParamEffect(trig_met, 1.01)
        sr_signal.setParamEffect(veto_tau, 1.03)
        sr_signal.setParamEffect(jec, 1.05)
        btagUp = template(signal, s, "btagUp", recoil, "sr", category)[0]
        btagDown = template(signal, s, "btagDown", recoil, "sr", category)[0]
        sr_signal.setParamEffect(btag, btagUp, btagDown)
        sr.addSample(sr_signal)

    ###
    # End of SR
    ###

    ###
    ###
    # Single muon W control region
    ###
    ###

    ch_name = "wmcr" + model_id
    wmcr = rl.Channel(ch_name)
    model.addChannel(wmcr)

    ###
    # Add data distribution to the channel
    ###

    wmcr.setObservation(template(data, "MET", "data", recoil, "wmcr", category))

    ###
    # W(->lnu)+jets data-driven model
    ###

    wmcr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, wmcr_wjetsTransferFactor, sr_wjets
    )
    wmcr.addSample(wmcr_wjets)

    ###
    # top-antitop model
    ###

    if category == "pass":
        wmcr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt
            )
    else:
        wmcr_ttTemplate = template(background, "TT", "nominal", recoil, "wmcr", category)
        wmcr_tt = rl.TemplateSample(
            ch_name + "_ttMC", rl.Sample.BACKGROUND, wmcr_ttTemplate
            )
        wmcr_tt.setParamEffect(lumi, 1.027)
        wmcr_tt.setParamEffect(trig_met, 1.01)
        wmcr_tt.setParamEffect(veto_tau, 1.03)
        wmcr_tt.setParamEffect(tt_norm, 1.2)
        wmcr_tt.setParamEffect(jec, 1.05)
        wmcr_tt.setParamEffect(id_mu, 1.02)
        wmcr_tt.setParamEffect(iso_mu, 1.02)
        btagUp = template(background, "TT", "btagUp", recoil, "wmcr", category)[0]
        btagDown = template(background, "TT", "btagDown", recoil, "wmcr", category)[0]
    wmcr_tt.setParamEffect(btag, btagUp, btagDown)
        
    wmcr.addSample(wmcr_tt)

    ###
    # Other MC-driven processes
    ###

    wmcr_stTemplate = template(background, "ST", "nominal", recoil, "wmcr", category)
    wmcr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, wmcr_stTemplate
    )
    wmcr_st.setParamEffect(lumi, 1.027)
    wmcr_st.setParamEffect(trig_met, 1.01)
    wmcr_st.setParamEffect(veto_tau, 1.03)
    wmcr_st.setParamEffect(st_norm, 1.2)
    wmcr_st.setParamEffect(jec, 1.05)
    wmcr_st.setParamEffect(id_mu, 1.02)
    wmcr_st.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "ST", "btagUp", recoil, "wmcr", category)[0]
    btagDown = template(background, "ST", "btagDown", recoil, "wmcr", category)[0]
    wmcr_st.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_st)

    wmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wmcr", category)
    wmcr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wmcr_dyjetsTemplate
    )
    wmcr_dyjets.setParamEffect(lumi, 1.027)
    wmcr_dyjets.setParamEffect(trig_met, 1.01)
    wmcr_dyjets.setParamEffect(veto_tau, 1.03)
    wmcr_dyjets.setParamEffect(zjets_norm, 1.4)
    wmcr_dyjets.setParamEffect(jec, 1.05)
    wmcr_dyjets.setParamEffect(id_mu, 1.02)
    wmcr_dyjets.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "DY+jets", "btagUp", recoil, "wmcr", category)[0]
    btagDown = template(background, "DY+jets", "btagDown", recoil, "wmcr", category)[0]
    wmcr_dyjets.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_dyjets)

    wmcr_vvTemplate = template(background, "VV", "nominal", recoil, "wmcr", category)
    wmcr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, wmcr_vvTemplate
    )
    wmcr_vv.setParamEffect(lumi, 1.027)
    wmcr_vv.setParamEffect(trig_met, 1.01)
    wmcr_vv.setParamEffect(veto_tau, 1.03)
    wmcr_vv.setParamEffect(vv_norm, 1.2)
    wmcr_vv.setParamEffect(jec, 1.05)
    wmcr_vv.setParamEffect(id_mu, 1.02)
    wmcr_vv.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "VV", "btagUp", recoil, "wmcr", category)[0]
    btagDown = template(background, "VV", "btagDown", recoil, "wmcr", category)[0]
    wmcr_vv.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_vv)

    wmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wmcr", category)
    wmcr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, wmcr_hbbTemplate
    )
    wmcr_hbb.setParamEffect(lumi, 1.027)
    wmcr_hbb.setParamEffect(trig_met, 1.01)
    wmcr_hbb.setParamEffect(veto_tau, 1.03)
    wmcr_hbb.setParamEffect(hbb_norm, 1.2)
    wmcr_hbb.setParamEffect(jec, 1.05)
    wmcr_hbb.setParamEffect(id_mu, 1.02)
    wmcr_hbb.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "Hbb", "btagUp", recoil, "wmcr", category)[0]
    btagDown = template(background, "Hbb", "btagDown", recoil, "wmcr", category)[0]
    wmcr_hbb.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_hbb)

    wmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wmcr", category)
    wmcr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, wmcr_qcdTemplate
    )
    wmcr_qcd.setParamEffect(lumi, 1.027)
    wmcr_qcd.setParamEffect(trig_met, 1.01)
    wmcr_qcd.setParamEffect(veto_tau, 1.03)
    wmcr_qcd.setParamEffect(qcdmu_norm, 2.0)
    wmcr_qcd.setParamEffect(jec, 1.05)
    wmcr_qcd.setParamEffect(id_mu, 1.02)
    wmcr_qcd.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "QCD", "btagUp", recoil, "wmcr", category)[0]
    btagDown = template(background, "QCD", "btagDown", recoil, "wmcr", category)[0]
    wmcr_qcd.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_qcd)

    ###
    # End of single muon W control region
    ###

    ###
    ###
    # Single muon top control region
    ###
    ###

    ch_name = "tmcr" + model_id
    tmcr = rl.Channel(ch_name)
    model.addChannel(tmcr)

    ###
    # Add data distribution to the channel
    ###

    tmcr.setObservation(template(data, "MET", "data", recoil, "tmcr", category))

    ###
    # W(->lnu)+jets data-driven model
    ###

    tmcr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, tmcr_wjetsTransferFactor, sr_wjets
    )
    tmcr.addSample(tmcr_wjets)

    ###
    # top-antitop model
    ###

    if category == "pass":
        tmcr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt
            )
    else:
        tmcr_ttTemplate = template(background, "TT", "nominal", recoil, "tmcr", category)
        tmcr_tt = rl.TemplateSample(
            ch_name + "_ttMC", rl.Sample.BACKGROUND, tmcr_ttTemplate
            )
        tmcr_tt.setParamEffect(lumi, 1.027)
        tmcr_tt.setParamEffect(trig_met, 1.01)
        tmcr_tt.setParamEffect(veto_tau, 1.03)
        tmcr_tt.setParamEffect(tt_norm, 1.2)
        tmcr_tt.setParamEffect(jec, 1.05)
        tmcr_tt.setParamEffect(id_mu, 1.02)
        tmcr_tt.setParamEffect(iso_mu, 1.02)
        btagUp = template(background, "TT", "btagUp", recoil, "tmcr", category)[0]
        btagDown = template(background, "TT", "btagDown", recoil, "tmcr", category)[0]
        tmcr_tt.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_tt)

    ###
    # Other MC-driven processes
    ###

    tmcr_stTemplate = template(background, "ST", "nominal", recoil, "tmcr", category)
    tmcr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, tmcr_stTemplate
    )
    tmcr_st.setParamEffect(lumi, 1.027)
    tmcr_st.setParamEffect(trig_met, 1.01)
    tmcr_st.setParamEffect(veto_tau, 1.03)
    tmcr_st.setParamEffect(st_norm, 1.2)
    tmcr_st.setParamEffect(jec, 1.05)
    tmcr_st.setParamEffect(id_mu, 1.02)
    tmcr_st.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "ST", "btagUp", recoil, "tmcr", category)[0]
    btagDown = template(background, "ST", "btagDown", recoil, "tmcr", category)[0]
    tmcr_st.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_st)

    tmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tmcr", category)
    tmcr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tmcr_dyjetsTemplate
    )
    tmcr_dyjets.setParamEffect(lumi, 1.027)
    tmcr_dyjets.setParamEffect(trig_met, 1.01)
    tmcr_dyjets.setParamEffect(veto_tau, 1.03)
    tmcr_dyjets.setParamEffect(zjets_norm, 1.4)
    tmcr_dyjets.setParamEffect(jec, 1.05)
    tmcr_dyjets.setParamEffect(id_mu, 1.02)
    tmcr_dyjets.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "DY+jets", "btagUp", recoil, "tmcr", category)[0]
    btagDown = template(background, "DY+jets", "btagDown", recoil, "tmcr", category)[0]
    tmcr_dyjets.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_dyjets)

    tmcr_vvTemplate = template(background, "VV", "nominal", recoil, "tmcr", category)
    tmcr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, tmcr_vvTemplate
    )
    tmcr_vv.setParamEffect(lumi, 1.027)
    tmcr_vv.setParamEffect(trig_met, 1.01)
    tmcr_vv.setParamEffect(veto_tau, 1.03)
    tmcr_vv.setParamEffect(vv_norm, 1.2)
    tmcr_vv.setParamEffect(jec, 1.05)
    tmcr_vv.setParamEffect(id_mu, 1.02)
    tmcr_vv.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "VV", "btagUp", recoil, "tmcr", category)[0]
    btagDown = template(background, "VV", "btagDown", recoil, "tmcr", category)[0]
    tmcr_vv.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_vv)

    tmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tmcr", category)
    tmcr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, tmcr_hbbTemplate
    )
    tmcr_hbb.setParamEffect(lumi, 1.027)
    tmcr_hbb.setParamEffect(trig_met, 1.01)
    tmcr_hbb.setParamEffect(veto_tau, 1.03)
    tmcr_hbb.setParamEffect(hbb_norm, 1.2)
    tmcr_hbb.setParamEffect(jec, 1.05)
    tmcr_hbb.setParamEffect(id_mu, 1.02)
    tmcr_hbb.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "Hbb", "btagUp", recoil, "tmcr", category)[0]
    btagDown = template(background, "Hbb", "btagDown", recoil, "tmcr", category)[0]
    tmcr_hbb.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_hbb)

    tmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tmcr", category)
    tmcr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, tmcr_qcdTemplate
    )
    tmcr_qcd.setParamEffect(lumi, 1.027)
    tmcr_qcd.setParamEffect(trig_met, 1.01)
    tmcr_qcd.setParamEffect(veto_tau, 1.03)
    tmcr_qcd.setParamEffect(qcdmu_norm, 2.0)
    tmcr_qcd.setParamEffect(jec, 1.05)
    tmcr_qcd.setParamEffect(id_mu, 1.02)
    tmcr_qcd.setParamEffect(iso_mu, 1.02)
    btagUp = template(background, "QCD", "btagUp", recoil, "tmcr", category)[0]
    btagDown = template(background, "QCD", "btagDown", recoil, "tmcr", category)[0]
    tmcr_qcd.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_qcd)

    ###
    # End of single muon top control region
    ###

    ###
    ###
    # Single electron W control region
    ###
    ###

    ch_name = "wecr" + model_id
    wecr = rl.Channel(ch_name)
    model.addChannel(wecr)

    ###
    # Add data distribution to the channel
    ###

    if year == "2018":
        wecr.setObservation(template(data, "EGamma", "data", recoil, "wecr", category))
    else:
        wecr.setObservation(template(data, "SingleElectron", "data", recoil, "wecr", category))

    ###
    # W(->lnu)+jets data-driven model
    ###

    wecr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, wecr_wjetsTransferFactor, sr_wjets
    )
    wecr.addSample(wecr_wjets)

    ###
    # top-antitop model
    ###

    if category == "pass":
        wecr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt
            )
    else:
        wecr_ttTemplate = template(background, "TT", "nominal", recoil, "wecr", category)
        wecr_tt = rl.TemplateSample(
            ch_name + "_ttMC", rl.Sample.BACKGROUND, wecr_ttTemplate
            )
        wecr_tt.setParamEffect(lumi, 1.027)
        wecr_tt.setParamEffect(trig_e, 1.01)
        wecr_tt.setParamEffect(veto_tau, 1.03)
        wecr_tt.setParamEffect(tt_norm, 1.2)
        wecr_tt.setParamEffect(jec, 1.05)
        wecr_tt.setParamEffect(id_e, 1.02)
        wecr_tt.setParamEffect(reco_e, 1.02)
        btagUp = template(background, "TT", "btagUp", recoil, "wecr", category)[0]
        btagDown = template(background, "TT", "btagDown", recoil, "wecr", category)[0]
        wecr_tt.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_tt)

    ###
    # Other MC-driven processes
    ###

    wecr_stTemplate = template(background, "ST", "nominal", recoil, "wecr", category)
    wecr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, wecr_stTemplate
    )
    wecr_st.setParamEffect(lumi, 1.027)
    wecr_st.setParamEffect(trig_e, 1.01)
    wecr_st.setParamEffect(veto_tau, 1.03)
    wecr_st.setParamEffect(st_norm, 1.2)
    wecr_st.setParamEffect(jec, 1.05)
    wecr_st.setParamEffect(id_e, 1.02)
    wecr_st.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "ST", "btagUp", recoil, "wecr", category)[0]
    btagDown = template(background, "ST", "btagDown", recoil, "wecr", category)[0]
    wecr_st.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_st)

    wecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wecr", category)
    wecr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wecr_dyjetsTemplate
    )
    wecr_dyjets.setParamEffect(lumi, 1.027)
    wecr_dyjets.setParamEffect(trig_e, 1.01)
    wecr_dyjets.setParamEffect(veto_tau, 1.03)
    wecr_dyjets.setParamEffect(zjets_norm, 1.4)
    wecr_dyjets.setParamEffect(jec, 1.05)
    wecr_dyjets.setParamEffect(id_e, 1.02)
    wecr_dyjets.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "DY+jets", "btagUp", recoil, "wecr", category)[0]
    btagDown = template(background, "DY+jets", "btagDown", recoil, "wecr", category)[0]
    wecr_dyjets.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_dyjets)

    wecr_vvTemplate = template(background, "VV", "nominal", recoil, "wecr", category)
    wecr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, wecr_vvTemplate
    )
    wecr_vv.setParamEffect(lumi, 1.027)
    wecr_vv.setParamEffect(trig_e, 1.01)
    wecr_vv.setParamEffect(veto_tau, 1.03)
    wecr_vv.setParamEffect(vv_norm, 1.2)
    wecr_vv.setParamEffect(jec, 1.05)
    wecr_vv.setParamEffect(id_e, 1.02)
    wecr_vv.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "VV", "btagUp", recoil, "wecr", category)[0]
    btagDown = template(background, "VV", "btagDown", recoil, "wecr", category)[0]
    wecr_vv.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_vv)

    wecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wecr", category)
    wecr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, wecr_hbbTemplate
    )
    wecr_hbb.setParamEffect(lumi, 1.027)
    wecr_hbb.setParamEffect(trig_e, 1.01)
    wecr_hbb.setParamEffect(veto_tau, 1.03)
    wecr_hbb.setParamEffect(hbb_norm, 1.2)
    wecr_hbb.setParamEffect(jec, 1.05)
    wecr_hbb.setParamEffect(id_e, 1.02)
    wecr_hbb.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "Hbb", "btagUp", recoil, "wecr", category)[0]
    btagDown = template(background, "Hbb", "btagDown", recoil, "wecr", category)[0]
    wecr_hbb.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_hbb)

    wecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wecr", category)
    wecr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, wecr_qcdTemplate
    )
    wecr_qcd.setParamEffect(lumi, 1.027)
    wecr_qcd.setParamEffect(trig_e, 1.01)
    wecr_qcd.setParamEffect(veto_tau, 1.03)
    wecr_qcd.setParamEffect(qcde_norm, 2.0)
    wecr_qcd.setParamEffect(jec, 1.05)
    wecr_qcd.setParamEffect(id_e, 1.02)
    wecr_qcd.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "QCD", "btagUp", recoil, "wecr", category)[0]
    btagDown = template(background, "QCD", "btagDown", recoil, "wecr", category)[0]
    wecr_qcd.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_qcd)

    ###
    # End of single electron W control region
    ###

    ###
    ###
    # Single electron top control region
    ###
    ###

    ch_name = "tecr" + model_id
    tecr = rl.Channel(ch_name)
    model.addChannel(tecr)

    ###
    # Add data distribution to the channel
    ###

    if year == "2018":
        tecr.setObservation(template(data, "EGamma", "data", recoil, "tecr", category))
    else:
        tecr.setObservation(template(data, "SingleElectron", "data", recoil, "tecr", category))

    ###
    # W(->lnu)+jets data-driven model
    ###

    tecr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, tecr_wjetsTransferFactor, sr_wjets
    )
    tecr.addSample(tecr_wjets)

    ###
    # top-antitop model
    ###

    if category == "pass":
        tecr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt
            )
    else:
        tecr_ttTemplate = template(background, "TT", "nominal", recoil, "tecr", category)
        tecr_tt = rl.TemplateSample(
            ch_name + "_ttMC", rl.Sample.BACKGROUND, tecr_ttTemplate
            )
        tecr_tt.setParamEffect(lumi, 1.027)
        tecr_tt.setParamEffect(trig_e, 1.01)
        tecr_tt.setParamEffect(veto_tau, 1.03)
        tecr_tt.setParamEffect(tt_norm, 1.2)
        tecr_tt.setParamEffect(jec, 1.05)
        tecr_tt.setParamEffect(id_e, 1.02)
        tecr_tt.setParamEffect(reco_e, 1.02)
        btagUp = template(background, "TT", "btagUp", recoil, "tecr", category)[0]
        btagDown = template(background, "TT", "btagDown", recoil, "tecr", category)[0]
        tecr_tt.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_tt)

    ###
    # Other MC-driven processes
    ###

    tecr_stTemplate = template(background, "ST", "nominal", recoil, "tecr", category)
    tecr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, tecr_stTemplate
    )
    tecr_st.setParamEffect(lumi, 1.027)
    tecr_st.setParamEffect(trig_e, 1.01)
    tecr_st.setParamEffect(veto_tau, 1.03)
    tecr_st.setParamEffect(st_norm, 1.2)
    tecr_st.setParamEffect(jec, 1.05)
    tecr_st.setParamEffect(id_e, 1.02)
    tecr_st.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "ST", "btagUp", recoil, "tecr", category)[0]
    btagDown = template(background, "ST", "btagDown", recoil, "tecr", category)[0]
    tecr_st.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_st)

    tecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tecr", category)
    tecr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tecr_dyjetsTemplate
    )
    tecr_dyjets.setParamEffect(lumi, 1.027)
    tecr_dyjets.setParamEffect(trig_e, 1.01)
    tecr_dyjets.setParamEffect(veto_tau, 1.03)
    tecr_dyjets.setParamEffect(zjets_norm, 1.4)
    tecr_dyjets.setParamEffect(jec, 1.05)
    tecr_dyjets.setParamEffect(id_e, 1.02)
    tecr_dyjets.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "DY+jets", "btagUp", recoil, "tecr", category)[0]
    btagDown = template(background, "DY+jets", "btagDown", recoil, "tecr", category)[0]
    tecr_dyjets.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_dyjets)

    tecr_vvTemplate = template(background, "VV", "nominal", recoil, "tecr", category)
    tecr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, tecr_vvTemplate
    )
    tecr_vv.setParamEffect(lumi, 1.027)
    tecr_vv.setParamEffect(trig_e, 1.01)
    tecr_vv.setParamEffect(veto_tau, 1.03)
    tecr_vv.setParamEffect(vv_norm, 1.2)
    tecr_vv.setParamEffect(jec, 1.05)
    tecr_vv.setParamEffect(id_e, 1.02)
    tecr_vv.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "VV", "btagUp", recoil, "tecr", category)[0]
    btagDown = template(background, "VV", "btagDown", recoil, "tecr", category)[0]
    tecr_vv.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_vv)

    tecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tecr", category)
    tecr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, tecr_hbbTemplate
    )
    tecr_hbb.setParamEffect(lumi, 1.027)
    tecr_hbb.setParamEffect(trig_e, 1.01)
    tecr_hbb.setParamEffect(veto_tau, 1.03)
    tecr_hbb.setParamEffect(hbb_norm, 1.2)
    tecr_hbb.setParamEffect(jec, 1.05)
    tecr_hbb.setParamEffect(id_e, 1.02)
    tecr_hbb.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "Hbb", "btagUp", recoil, "tecr", category)[0]
    btagDown = template(background, "Hbb", "btagDown", recoil, "tecr", category)[0]
    tecr_hbb.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_hbb)

    tecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tecr", category)
    tecr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, tecr_qcdTemplate
    )
    tecr_qcd.setParamEffect(lumi, 1.027)
    tecr_qcd.setParamEffect(trig_e, 1.01)
    tecr_qcd.setParamEffect(veto_tau, 1.03)
    tecr_qcd.setParamEffect(qcde_norm, 2.0)
    tecr_qcd.setParamEffect(jec, 1.05)
    tecr_qcd.setParamEffect(id_e, 1.02)
    tecr_qcd.setParamEffect(reco_e, 1.02)
    btagUp = template(background, "QCD", "btagUp", recoil, "tecr", category)[0]
    btagDown = template(background, "QCD", "btagDown", recoil, "tecr", category)[0]
    tecr_qcd.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_qcd)

    ###
    # End of single electron top control region
    ###

    ###
    ###
    # Double muon control region
    ###
    ###

    ch_name = "zmcr" + model_id
    zmcr = rl.Channel(ch_name)
    model.addChannel(zmcr)

    ###
    # Add data distribution to the channel
    ###

    zmcr.setObservation(template(data, "MET", "data", recoil, "zmcr", category))

    zmcr_dyjets = rl.TransferFactorSample(
        ch_name + "_dyjets", rl.Sample.BACKGROUND, zmcr_dyjetsTransferFactor, sr_zjets
    )
    zmcr.addSample(zmcr_dyjets)

    ###
    # Other MC-driven processes
    ###

    zmcr_ttTemplate = template(background, "TT", "nominal", recoil, "zmcr", category)
    zmcr_tt = rl.TemplateSample(
        ch_name + "_ttMC", rl.Sample.BACKGROUND, zmcr_ttTemplate
    )
    zmcr_tt.setParamEffect(lumi, 1.027)
    zmcr_tt.setParamEffect(trig_met, 1.01)
    zmcr_tt.setParamEffect(veto_tau, 1.03)
    zmcr_tt.setParamEffect(tt_norm, 1.2)
    zmcr_tt.setParamEffect(jec, 1.05)
    zmcr_tt.setParamEffect(id_mu, 1.02)
    zmcr_tt.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_tt)

    zmcr_stTemplate = template(background, "ST", "nominal", recoil, "zmcr", category)
    zmcr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, zmcr_stTemplate
    )
    zmcr_st.setParamEffect(lumi, 1.027)
    zmcr_st.setParamEffect(trig_met, 1.01)
    zmcr_st.setParamEffect(veto_tau, 1.03)
    zmcr_st.setParamEffect(st_norm, 1.2)
    zmcr_st.setParamEffect(jec, 1.05)
    zmcr_st.setParamEffect(id_mu, 1.02)
    zmcr_st.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_st)

    zmcr_vvTemplate = template(background, "VV", "nominal", recoil, "zmcr", category)
    zmcr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, zmcr_vvTemplate
    )
    zmcr_vv.setParamEffect(lumi, 1.027)
    zmcr_vv.setParamEffect(trig_met, 1.01)
    zmcr_vv.setParamEffect(veto_tau, 1.03)
    zmcr_vv.setParamEffect(vv_norm, 1.2)
    zmcr_vv.setParamEffect(jec, 1.05)
    zmcr_vv.setParamEffect(id_mu, 1.02)
    zmcr_vv.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_vv)

    zmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "zmcr", category)
    zmcr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, zmcr_hbbTemplate
    )
    zmcr_hbb.setParamEffect(lumi, 1.027)
    zmcr_hbb.setParamEffect(trig_met, 1.01)
    zmcr_hbb.setParamEffect(veto_tau, 1.03)
    zmcr_hbb.setParamEffect(hbb_norm, 1.2)
    zmcr_hbb.setParamEffect(jec, 1.05)
    zmcr_hbb.setParamEffect(id_mu, 1.02)
    zmcr_hbb.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_hbb)

    ###
    # End of double muon control region
    ###

    ###
    ###
    # Double electron control region
    ###
    ###

    ch_name = "zecr" + model_id
    zecr = rl.Channel(ch_name)
    model.addChannel(zecr)

    ###
    # Add data distribution to the channel
    ###

    if year == "2018":
        zecr.setObservation(template(data, "EGamma", "data", recoil, "zecr", category))
    else:
        zecr.setObservation(template(data, "SingleElectron", "data", recoil, "zecr", category))

    zecr_dyjets = rl.TransferFactorSample(
        ch_name + "_dyjets", rl.Sample.BACKGROUND, zecr_dyjetsTransferFactor, sr_zjets
    )
    zecr.addSample(zecr_dyjets)

    ###
    # Other MC-driven processes
    ###

    zecr_ttTemplate = template(background, "TT", "nominal", recoil, "zecr", category)
    zecr_tt = rl.TemplateSample(
        ch_name + "_ttMC", rl.Sample.BACKGROUND, zecr_ttTemplate
    )
    zecr_tt.setParamEffect(lumi, 1.027)
    zecr_tt.setParamEffect(trig_e, 1.01)
    zecr_tt.setParamEffect(veto_tau, 1.03)
    zecr_tt.setParamEffect(tt_norm, 1.2)
    zecr_tt.setParamEffect(jec, 1.05)
    zecr_tt.setParamEffect(id_e, 1.02)
    zecr_tt.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_tt)

    zecr_stTemplate = template(background, "ST", "nominal", recoil, "zecr", category)
    zecr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, zecr_stTemplate
    )
    zecr_st.setParamEffect(lumi, 1.027)
    zecr_st.setParamEffect(trig_e, 1.01)
    zecr_st.setParamEffect(veto_tau, 1.03)
    zecr_st.setParamEffect(st_norm, 1.2)
    zecr_st.setParamEffect(jec, 1.05)
    zecr_st.setParamEffect(id_e, 1.02)
    zecr_st.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_st)

    zecr_vvTemplate = template(background, "VV", "nominal", recoil, "zecr", category)
    zecr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, zecr_vvTemplate
    )
    zecr_vv.setParamEffect(lumi, 1.027)
    zecr_vv.setParamEffect(trig_e, 1.01)
    zecr_vv.setParamEffect(veto_tau, 1.03)
    zecr_vv.setParamEffect(vv_norm, 1.2)
    zecr_vv.setParamEffect(jec, 1.05)
    zecr_vv.setParamEffect(id_e, 1.02)
    zecr_vv.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_vv)

    zecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "zecr", category)
    zecr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, zecr_hbbTemplate
    )
    zecr_hbb.setParamEffect(lumi, 1.027)
    zecr_hbb.setParamEffect(trig_e, 1.01)
    zecr_hbb.setParamEffect(veto_tau, 1.03)
    zecr_hbb.setParamEffect(hbb_norm, 1.2)
    zecr_hbb.setParamEffect(jec, 1.05)
    zecr_hbb.setParamEffect(id_e, 1.02)
    zecr_hbb.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_hbb)

    ###
    # End of double electron control region
    ###

    ###
    ###
    # Single photon control region
    ###
    ###

    ch_name = "gcr" + model_id
    gcr = rl.Channel(ch_name)
    model.addChannel(gcr)

    ###
    # Add data distribution to the channel
    ###

    if year == "2018":
        gcr.setObservation(template(data, "EGamma", "data", recoil, "gcr", category))
    else:
        gcr.setObservation(template(data, "SinglePhoton", "data", recoil, "gcr", category))

    gcr_gjets = rl.TransferFactorSample(
        ch_name + "_gjets", rl.Sample.BACKGROUND, gcr_gjetsTransferFactor, sr_zjets
    )
    gcr.addSample(gcr_gjets)

    gcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "gcr", category)
    gcr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, gcr_qcdTemplate
    )
    gcr_qcd.setParamEffect(lumi, 1.027)
    gcr_qcd.setParamEffect(trig_pho, 1.01)
    gcr_qcd.setParamEffect(veto_tau, 1.03)
    gcr_qcd.setParamEffect(qcdpho_norm, 2.0)
    gcr_qcd.setParamEffect(jec, 1.05)
    gcr_qcd.setParamEffect(id_pho, 1.02)
    gcr.addSample(gcr_qcd)

    return model


if __name__ == "__main__":
    if not os.path.exists("datacards"):
        os.mkdir("datacards")
    parser = OptionParser()
    parser.add_option("-y", "--year", help="year", dest="year", default="")
    (options, args) = parser.parse_args()
    year = options.year

    ###
    # Extract histograms from input file
    ###

    print("Grouping histograms")
    hists = load("hists/darkhiggs" + year + ".scaled")
    hists = remap_histograms(hists)
    (
        sr_zjetsShape,
        sr_zjetsRate,
        sr_ttNuisances,
    ) = initialize_nuisances(hists, year)
    tf_params = rhalphabeth(mass_binning)
    #tf_params = 0.05

    ###
    ###
    # Setting up systematics
    ###
    ###
    lumi = rl.NuisanceParameter("lumi" + year, "lnN")
    qcdpho_norm = rl.NuisanceParameter("qcdpho_norm", "lnN")
    qcde_norm = rl.NuisanceParameter("qcde_norm", "lnN")
    qcdmu_norm = rl.NuisanceParameter("qcdmu_norm", "lnN")
    qcdsig_norm = rl.NuisanceParameter("qcdsig_norm", "lnN")
    st_norm = rl.NuisanceParameter("st_norm", "lnN")
    tt_norm = rl.NuisanceParameter("tt_norm", "lnN")
    vv_norm = rl.NuisanceParameter("vv_norm", "lnN")
    hbb_norm = rl.NuisanceParameter("hbb_norm", "lnN")
    zjets_norm = rl.NuisanceParameter("zjets_norm", "lnN")
    wjets_norm = rl.NuisanceParameter("wjets_norm", "lnN")
    gjets_norm = rl.NuisanceParameter("gjets_norm", "lnN")
    whf_fraction = rl.NuisanceParameter("whf_fraction", "lnN")
    zhf_fraction = rl.NuisanceParameter("zhf_fraction", "lnN")
    ghf_fraction = rl.NuisanceParameter("ghf_fraction", "lnN")
    id_e = rl.NuisanceParameter("id_e" + year, "lnN")
    id_mu = rl.NuisanceParameter("id_mu" + year, "lnN")
    id_pho = rl.NuisanceParameter("id_pho" + year, "lnN")
    reco_e = rl.NuisanceParameter("reco_e" + year, "lnN")
    iso_mu = rl.NuisanceParameter("iso_mu" + year, "lnN")
    trig_e = rl.NuisanceParameter("trig_e" + year, "lnN")
    trig_met = rl.NuisanceParameter("trig_met" + year, "lnN")
    trig_pho = rl.NuisanceParameter("trig_pho" + year, "lnN")
    veto_tau = rl.NuisanceParameter("veto_tau" + year, "lnN")
    jec = rl.NuisanceParameter("jec" + year, "lnN")
    btag = rl.NuisanceParameter("btag" + year, "shape")  # AK4 btag
    ew1 = rl.NuisanceParameter("ew1", "shape")
    ew2G = rl.NuisanceParameter("ew2G", "shape")
    ew2W = rl.NuisanceParameter("ew2W", "shape")
    ew2Z = rl.NuisanceParameter("ew2Z", "shape")
    ew3G = rl.NuisanceParameter("ew3G", "shape")
    ew3W = rl.NuisanceParameter("ew3W", "shape")
    ew3Z = rl.NuisanceParameter("ew3Z", "shape")
    mix = rl.NuisanceParameter("mix", "shape")
    muF = rl.NuisanceParameter("muF", "shape")
    muR = rl.NuisanceParameter("muR", "shape")
    qcd1 = rl.NuisanceParameter("qcd1", "shape")
    qcd2 = rl.NuisanceParameter("qcd2", "shape")
    qcd3 = rl.NuisanceParameter("qcd3", "shape")

    model_dict = {}
    recoilbins = np.array(recoil_binning)
    nrecoil = len(recoilbins) - 1
    for recoilbin in range(nrecoil):
        sr_zjetsBinYields = sr_zjetsShape * sr_zjetsRate[recoilbin]
        for category in ["pass", "fail"]:
            sr_ttBinYields = sr_ttNuisances[recoilbin]
            (
                sr_wjetsTransferFactor,
                wmcr_wjetsTransferFactor,
                wmcr_ttTransferFactor,
                tmcr_wjetsTransferFactor,
                tmcr_ttTransferFactor,
                wecr_wjetsTransferFactor,
                wecr_ttTransferFactor,
                tecr_wjetsTransferFactor,
                tecr_ttTransferFactor,
                zmcr_dyjetsTransferFactor,
                zecr_dyjetsTransferFactor,
                gcr_gjetsTransferFactor,
            ) = computeTFs(hists, year, recoilbin, category)
            with open(
                "data/darkhiggs-"
                + year
                + "-"
                + category
                + "-recoil"
                + str(recoilbin)
                + ".model",
                "wb",
            ) as fout:
                pickle.dump(model(year, recoilbin, category), fout, protocol=2)
