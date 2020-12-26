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

def template(dictionary, process, systematic, region, category):
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
    # First, tagging efficiency and SF
    ###

    tt_efficiency = {"2018": 0.5}
    sf_tt = rl.IndependentParameter(
        "sf_tt" + year, 1.0, 0.01, 1.0 / tt_efficiency[year]
    )

    tt_weight = {
        "pass": sf_tt * tt_efficiency[year],
        "fail": 1 - (sf_tt * tt_efficiency[year]),
    }

    sr_tt = (
        hists["bkg"]["template"]
        .integrate("region", "sr")
        .integrate("process", "TT")
        .integrate("systematic", "nominal")
    )

    ###
    # Then, shapes in pass/fail
    ###

    sr_ttMass = sr_tt.sum("gentype", "recoil")
    sr_ttMassPass = sr_ttMass.values()[()][
        :, 1
    ]  # get the pass histogram, inclusive in recoil
    sr_ttMassPass = (
        sr_ttMassPass / sr_ttMassPass.sum()
    )  # normalize to the integral to get the shape in pass
    sr_ttMassFail = sr_ttMass.values()[()][
        :, 0
    ]  # get the fail histogram, inclusive in recoil
    sr_ttMassFail = (
        sr_ttMassFail / sr_ttMassFail.sum()
    )  # normalize to the integral to get the shape in fail
    sr_ttShape = {
        "pass": np.array(  # one nuisance per mass shape bin in pass
            [
                rl.IndependentParameter(
                    "sr" + year + "_ttshape_pass_mass%d" % i,
                    b,
                    0,
                    sr_ttMassPass.max() * 2,
                )
                for i, b in enumerate(sr_ttMassPass)
            ]
        ),
        "fail": np.array(  # one nuisance per mass shape bin in fail
            [
                rl.IndependentParameter(
                    "sr" + year + "_ttshape_fail_mass%d" % i,
                    b,
                    0,
                    sr_ttMassFail.max() * 2,
                )
                for i, b in enumerate(sr_ttMassFail)
            ]
        ),
    }

    ###
    # Lastly, let's initialize nuisances per recoil bin
    ###

    sr_ttRecoil = sr_tt.sum("gentype", "fjmass", "ZHbbvsQCD").values()[()][:]
    sr_ttRate = np.array(
        [
            rl.IndependentParameter(
                "sr" + year + "_tt_recoil%d" % i, b, 0, sr_ttRecoil.max() * 2
            )
            for i, b in enumerate(sr_ttRecoil)
        ]
    )

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
    return sr_zjetsShape, sr_zjetsRate, sr_ttShape, sr_ttRate, tt_weight


def computeTFs(hists, year, recoil, category):

    model_id = year + category + "recoil" + str(recoil)

    background = {}
    for r in bkg_hists["template"].identifiers("region"):
        background[str(r)] = bkg_hists["template"].integrate("region", r).sum("gentype")

    def addBtagSyst(process, region, templ):
        btagUp = template(background, process, "btagUp", region, category)[0]
        btagDown = template(background, process, "btagDown", region, category)[0]
        templ.setParamEffect(btag, btagUp, btagDown)

    def addVJetsSyst(process, region, templ):
        ew1Up = template(background, process, "ew1Up", region, category)[0]
        ew1Down = template(background, process, "ew1Down", region, category)[0]
        templ.setParamEffect(ew1, ew1Up, ew1Down)
        ew2GUp = template(background, process, "ew2GUp", region, category)[0]
        ew2GDown = template(background, process, "ew2GDown", region, category)[0]
        templ.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2GUp = templateMass(sr_zjets, "ew2GUp")[0]
    ew2GDown = templateMass(sr_zjets, "ew2GDown")[0]
    sr_zjets.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(sr_zjets, "ew2WUp")[0]
    ew2WDown = templateMass(sr_zjets, "ew2WDown")[0]
    sr_zjets.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(sr_zjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(sr_zjets, "ew2ZDown")[0]
    sr_zjets.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(sr_zjets, "ew3GUp")[0]
    ew3GDown = templateMass(sr_zjets, "ew3GDown")[0]
    sr_zjets.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(sr_zjets, "ew3WUp")[0]
    ew3WDown = templateMass(sr_zjets, "ew3WDown")[0]
    sr_zjets.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(sr_zjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(sr_zjets, "ew3ZDown")[0]
    sr_zjets.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(sr_zjets, "mixUp")[0]
    mixDown = templateMass(sr_zjets, "mixDown")[0]
    sr_zjets.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(sr_zjets, "muFUp")[0]
    muFDown = templateMass(sr_zjets, "muFDown")[0]
    sr_zjets.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(sr_zjets, "muRUp")[0]
    muRDown = templateMass(sr_zjets, "muRDown")[0]
    sr_zjets.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(sr_zjets, "qcd1Up")[0]
    qcd1Down = templateMass(sr_zjets, "qcd1Down")[0]
    sr_zjets.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(sr_zjets, "qcd2Up")[0]
    qcd2Down = templateMass(sr_zjets, "qcd2Down")[0]
    sr_zjets.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(sr_zjets, "qcd3Up")[0]
    qcd3Down = templateMass(sr_zjets, "qcd3Down")[0]
    sr_zjets.setParamEffect(qcd3, qcd3Up, qcd3Down)


    ###
    # Z+jets template
    ###

    sr_zjets = rl.TemplateSample(
        "sr" + model_id + "_zjets",
        rl.Sample.BACKGROUND,
        template(background, "Z+jets", "nominal", "sr", category),
    )
    sr_zjets.setParamEffect(lumi, 1.027)
    sr_zjets.setParamEffect(zjets_norm, 1.4)
    sr_zjets.setParamEffect(trig_met, 1.01)
    sr_zjets.setParamEffect(veto_tau, 1.03)
    sr_zjets.setParamEffect(jec, 1.05)
    sr_zjets.setParamEffect(zhf_fraction, hf_systematic["Z+jets"]["sr"][category])
    addBtagSyst("Z+jets", "sr", sr_zjets)
    addVjetsSyst("Z+jets", "sr", sr_zjets)

    ###
    # W+jets templates
    ###

    sr_wjets = (
        hists["bkg"]["template"]
        .integrate("region", "sr")
        .integrate("process", "W+jets")
        .sum("gentype")
    )
    sr_wjetsMass = rl.TemplateSample(
        "sr" + model_id + "_wjetsMass",
        rl.Sample.BACKGROUND,
        templateMass(sr_wjets, "nominal"),
    )
    sr_wjetsMass.setParamEffect(lumi, 1.027)
    sr_wjetsMass.setParamEffect(wjets_norm, 1.4)
    sr_wjetsMass.setParamEffect(trig_met, 1.01)
    sr_wjetsMass.setParamEffect(veto_tau, 1.03)
    sr_wjetsMass.setParamEffect(jec, 1.05)
    sr_wjetsMass.setParamEffect(whf_fraction, hf_systematic["W+jets"]["sr"][category])
    btagUp = templateMass(sr_wjets, "btagUp")[0]
    btagDown = templateMass(sr_wjets, "btagDown")[0]
    sr_wjetsMass.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateMass(sr_wjets, "ew1Up")[0]
    ew1Down = templateMass(sr_wjets, "ew1Down")[0]
    sr_wjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(sr_wjets, "ew2GUp")[0]
    ew2GDown = templateMass(sr_wjets, "ew2GDown")[0]
    sr_wjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(sr_wjets, "ew2WUp")[0]
    ew2WDown = templateMass(sr_wjets, "ew2WDown")[0]
    sr_wjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(sr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(sr_wjets, "ew2ZDown")[0]
    sr_wjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(sr_wjets, "ew3GUp")[0]
    ew3GDown = templateMass(sr_wjets, "ew3GDown")[0]
    sr_wjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(sr_wjets, "ew3WUp")[0]
    ew3WDown = templateMass(sr_wjets, "ew3WDown")[0]
    sr_wjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(sr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(sr_wjets, "ew3ZDown")[0]
    sr_wjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(sr_wjets, "mixUp")[0]
    mixDown = templateMass(sr_wjets, "mixDown")[0]
    sr_wjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(sr_wjets, "muFUp")[0]
    muFDown = templateMass(sr_wjets, "muFDown")[0]
    sr_wjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(sr_wjets, "muRUp")[0]
    muRDown = templateMass(sr_wjets, "muRDown")[0]
    sr_wjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(sr_wjets, "qcd1Up")[0]
    qcd1Down = templateMass(sr_wjets, "qcd1Down")[0]
    sr_wjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(sr_wjets, "qcd2Up")[0]
    qcd2Down = templateMass(sr_wjets, "qcd2Down")[0]
    sr_wjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(sr_wjets, "qcd3Up")[0]
    qcd3Down = templateMass(sr_wjets, "qcd3Down")[0]
    sr_wjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    sr_wjetsRecoil = rl.TemplateSample(
        "sr" + model_id + "_wjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(sr_wjets, "nominal"),
    )
    sr_wjetsRecoil.setParamEffect(lumi, 1.027)
    sr_wjetsRecoil.setParamEffect(wjets_norm, 1.4)
    sr_wjetsRecoil.setParamEffect(trig_met, 1.01)
    sr_wjetsRecoil.setParamEffect(veto_tau, 1.03)
    sr_wjetsRecoil.setParamEffect(jec, 1.05)
    sr_wjetsRecoil.setParamEffect(whf_fraction, hf_systematic["W+jets"]["sr"][category])
    btagUp = templateRecoil(sr_wjets, "btagUp")[0]
    btagDown = templateRecoil(sr_wjets, "btagDown")[0]
    sr_wjetsRecoil.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateRecoil(sr_wjets, "ew1Up")[0]
    ew1Down = templateRecoil(sr_wjets, "ew1Down")[0]
    sr_wjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(sr_wjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(sr_wjets, "ew2GDown")[0]
    sr_wjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(sr_wjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(sr_wjets, "ew2WDown")[0]
    sr_wjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(sr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(sr_wjets, "ew2ZDown")[0]
    sr_wjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(sr_wjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(sr_wjets, "ew3GDown")[0]
    sr_wjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(sr_wjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(sr_wjets, "ew3WDown")[0]
    sr_wjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(sr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(sr_wjets, "ew3ZDown")[0]
    sr_wjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(sr_wjets, "mixUp")[0]
    mixDown = templateRecoil(sr_wjets, "mixDown")[0]
    sr_wjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(sr_wjets, "muFUp")[0]
    muFDown = templateRecoil(sr_wjets, "muFDown")[0]
    sr_wjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(sr_wjets, "muRUp")[0]
    muRDown = templateRecoil(sr_wjets, "muRDown")[0]
    sr_wjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(sr_wjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(sr_wjets, "qcd1Down")[0]
    sr_wjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(sr_wjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(sr_wjets, "qcd2Down")[0]
    sr_wjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(sr_wjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(sr_wjets, "qcd3Down")[0]
    sr_wjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    wmcr_wjets = (
        hists["bkg"]["template"]
        .integrate("region", "wmcr")
        .integrate("process", "W+jets")
        .sum("gentype")
    )
    wmcr_wjetsMass = rl.TemplateSample(
        "wmcr" + model_id + "_wjetsMass",
        rl.Sample.BACKGROUND,
        templateMass(wmcr_wjets, "nominal"),
    )
    wmcr_wjetsMass.setParamEffect(lumi, 1.027)
    wmcr_wjetsMass.setParamEffect(trig_met, 1.01)
    wmcr_wjetsMass.setParamEffect(veto_tau, 1.03)
    wmcr_wjetsMass.setParamEffect(wjets_norm, 1.4)
    wmcr_wjetsMass.setParamEffect(jec, 1.05)
    wmcr_wjetsMass.setParamEffect(id_mu, 1.02)
    wmcr_wjetsMass.setParamEffect(iso_mu, 1.02)
    wmcr_wjetsMass.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["wmcr"][category]
    )
    btagUp = templateMass(wmcr_wjets, "btagUp")[0]
    btagDown = templateMass(wmcr_wjets, "btagDown")[0]
    wmcr_wjetsMass.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateMass(wmcr_wjets, "ew1Up")[0]
    ew1Down = templateMass(wmcr_wjets, "ew1Down")[0]
    wmcr_wjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(wmcr_wjets, "ew2GUp")[0]
    ew2GDown = templateMass(wmcr_wjets, "ew2GDown")[0]
    wmcr_wjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(wmcr_wjets, "ew2WUp")[0]
    ew2WDown = templateMass(wmcr_wjets, "ew2WDown")[0]
    wmcr_wjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(wmcr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(wmcr_wjets, "ew2ZDown")[0]
    wmcr_wjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(wmcr_wjets, "ew3GUp")[0]
    ew3GDown = templateMass(wmcr_wjets, "ew3GDown")[0]
    wmcr_wjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(wmcr_wjets, "ew3WUp")[0]
    ew3WDown = templateMass(wmcr_wjets, "ew3WDown")[0]
    wmcr_wjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(wmcr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(wmcr_wjets, "ew3ZDown")[0]
    wmcr_wjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(wmcr_wjets, "mixUp")[0]
    mixDown = templateMass(wmcr_wjets, "mixDown")[0]
    wmcr_wjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(wmcr_wjets, "muFUp")[0]
    muFDown = templateMass(wmcr_wjets, "muFDown")[0]
    wmcr_wjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(wmcr_wjets, "muRUp")[0]
    muRDown = templateMass(wmcr_wjets, "muRDown")[0]
    wmcr_wjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(wmcr_wjets, "qcd1Up")[0]
    qcd1Down = templateMass(wmcr_wjets, "qcd1Down")[0]
    wmcr_wjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(wmcr_wjets, "qcd2Up")[0]
    qcd2Down = templateMass(wmcr_wjets, "qcd2Down")[0]
    wmcr_wjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(wmcr_wjets, "qcd3Up")[0]
    qcd3Down = templateMass(wmcr_wjets, "qcd3Down")[0]
    wmcr_wjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    wmcr_wjetsRecoil = rl.TemplateSample(
        "wmcr" + model_id + "_wjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(wmcr_wjets, "nominal"),
    )
    wmcr_wjetsRecoil.setParamEffect(lumi, 1.027)
    wmcr_wjetsRecoil.setParamEffect(trig_met, 1.01)
    wmcr_wjetsRecoil.setParamEffect(veto_tau, 1.03)
    wmcr_wjetsRecoil.setParamEffect(wjets_norm, 1.4)
    wmcr_wjetsRecoil.setParamEffect(jec, 1.05)
    wmcr_wjetsRecoil.setParamEffect(id_mu, 1.02)
    wmcr_wjetsRecoil.setParamEffect(iso_mu, 1.02)
    wmcr_wjetsRecoil.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["wmcr"][category]
    )
    btagUp = templateRecoil(wmcr_wjets, "btagUp")[0]
    btagDown = templateRecoil(wmcr_wjets, "btagDown")[0]
    wmcr_wjetsRecoil.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateRecoil(wmcr_wjets, "ew1Up")[0]
    ew1Down = templateRecoil(wmcr_wjets, "ew1Down")[0]
    wmcr_wjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(wmcr_wjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(wmcr_wjets, "ew2GDown")[0]
    wmcr_wjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(wmcr_wjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(wmcr_wjets, "ew2WDown")[0]
    wmcr_wjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(wmcr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(wmcr_wjets, "ew2ZDown")[0]
    wmcr_wjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(wmcr_wjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(wmcr_wjets, "ew3GDown")[0]
    wmcr_wjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(wmcr_wjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(wmcr_wjets, "ew3WDown")[0]
    wmcr_wjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(wmcr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(wmcr_wjets, "ew3ZDown")[0]
    wmcr_wjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(wmcr_wjets, "mixUp")[0]
    mixDown = templateRecoil(wmcr_wjets, "mixDown")[0]
    wmcr_wjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(wmcr_wjets, "muFUp")[0]
    muFDown = templateRecoil(wmcr_wjets, "muFDown")[0]
    wmcr_wjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(wmcr_wjets, "muRUp")[0]
    muRDown = templateRecoil(wmcr_wjets, "muRDown")[0]
    wmcr_wjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(wmcr_wjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(wmcr_wjets, "qcd1Down")[0]
    wmcr_wjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(wmcr_wjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(wmcr_wjets, "qcd2Down")[0]
    wmcr_wjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(wmcr_wjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(wmcr_wjets, "qcd3Down")[0]
    wmcr_wjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    wecr_wjets = (
        hists["bkg"]["template"]
        .integrate("region", "wecr")
        .integrate("process", "W+jets")
        .sum("gentype")
    )
    wecr_wjetsMass = rl.TemplateSample(
        "wecr" + model_id + "_wjetsMass",
        rl.Sample.BACKGROUND,
        templateMass(wecr_wjets, "nominal"),
    )
    wecr_wjetsMass.setParamEffect(lumi, 1.027)
    wecr_wjetsMass.setParamEffect(trig_e, 1.01)
    wecr_wjetsMass.setParamEffect(veto_tau, 1.03)
    wecr_wjetsMass.setParamEffect(wjets_norm, 1.4)
    wecr_wjetsMass.setParamEffect(jec, 1.05)
    wecr_wjetsMass.setParamEffect(id_e, 1.02)
    wecr_wjetsMass.setParamEffect(reco_e, 1.02)
    wecr_wjetsMass.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["wecr"][category]
    )
    btagUp = templateMass(wecr_wjets, "btagUp")[0]
    btagDown = templateMass(wecr_wjets, "btagDown")[0]
    wecr_wjetsMass.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateMass(wecr_wjets, "ew1Up")[0]
    ew1Down = templateMass(wecr_wjets, "ew1Down")[0]
    wecr_wjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(wecr_wjets, "ew2GUp")[0]
    ew2GDown = templateMass(wecr_wjets, "ew2GDown")[0]
    wecr_wjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(wecr_wjets, "ew2WUp")[0]
    ew2WDown = templateMass(wecr_wjets, "ew2WDown")[0]
    wecr_wjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(wecr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(wecr_wjets, "ew2ZDown")[0]
    wecr_wjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(wecr_wjets, "ew3GUp")[0]
    ew3GDown = templateMass(wecr_wjets, "ew3GDown")[0]
    wecr_wjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(wecr_wjets, "ew3WUp")[0]
    ew3WDown = templateMass(wecr_wjets, "ew3WDown")[0]
    wecr_wjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(wecr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(wecr_wjets, "ew3ZDown")[0]
    wecr_wjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(wecr_wjets, "mixUp")[0]
    mixDown = templateMass(wecr_wjets, "mixDown")[0]
    wecr_wjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(wecr_wjets, "muFUp")[0]
    muFDown = templateMass(wecr_wjets, "muFDown")[0]
    wecr_wjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(wecr_wjets, "muRUp")[0]
    muRDown = templateMass(wecr_wjets, "muRDown")[0]
    wecr_wjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(wecr_wjets, "qcd1Up")[0]
    qcd1Down = templateMass(wecr_wjets, "qcd1Down")[0]
    wecr_wjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(wecr_wjets, "qcd2Up")[0]
    qcd2Down = templateMass(wecr_wjets, "qcd2Down")[0]
    wecr_wjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(wecr_wjets, "qcd3Up")[0]
    qcd3Down = templateMass(wecr_wjets, "qcd3Down")[0]
    wecr_wjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    wecr_wjetsRecoil = rl.TemplateSample(
        "wecr" + model_id + "_wjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(wecr_wjets, "nominal"),
    )
    wecr_wjetsRecoil.setParamEffect(lumi, 1.027)
    wecr_wjetsRecoil.setParamEffect(trig_e, 1.01)
    wecr_wjetsRecoil.setParamEffect(veto_tau, 1.03)
    wecr_wjetsRecoil.setParamEffect(wjets_norm, 1.4)
    wecr_wjetsRecoil.setParamEffect(jec, 1.05)
    wecr_wjetsRecoil.setParamEffect(id_e, 1.02)
    wecr_wjetsRecoil.setParamEffect(reco_e, 1.02)
    wecr_wjetsRecoil.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["wecr"][category]
    )
    btagUp = templateRecoil(wecr_wjets, "btagUp")[0]
    btagDown = templateRecoil(wecr_wjets, "btagDown")[0]
    wecr_wjetsRecoil.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateRecoil(wecr_wjets, "ew1Up")[0]
    ew1Down = templateRecoil(wecr_wjets, "ew1Down")[0]
    wecr_wjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(wecr_wjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(wecr_wjets, "ew2GDown")[0]
    wecr_wjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(wecr_wjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(wecr_wjets, "ew2WDown")[0]
    wecr_wjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(wecr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(wecr_wjets, "ew2ZDown")[0]
    wecr_wjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(wecr_wjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(wecr_wjets, "ew3GDown")[0]
    wecr_wjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(wecr_wjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(wecr_wjets, "ew3WDown")[0]
    wecr_wjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(wecr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(wecr_wjets, "ew3ZDown")[0]
    wecr_wjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(wecr_wjets, "mixUp")[0]
    mixDown = templateRecoil(wecr_wjets, "mixDown")[0]
    wecr_wjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(wecr_wjets, "muFUp")[0]
    muFDown = templateRecoil(wecr_wjets, "muFDown")[0]
    wecr_wjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(wecr_wjets, "muRUp")[0]
    muRDown = templateRecoil(wecr_wjets, "muRDown")[0]
    wecr_wjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(wecr_wjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(wecr_wjets, "qcd1Down")[0]
    wecr_wjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(wecr_wjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(wecr_wjets, "qcd2Down")[0]
    wecr_wjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(wecr_wjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(wecr_wjets, "qcd3Down")[0]
    wecr_wjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    tmcr_wjets = (
        hists["bkg"]["template"]
        .integrate("region", "tmcr")
        .integrate("process", "W+jets")
        .sum("gentype")
    )
    tmcr_wjetsMass = rl.TemplateSample(
        "tmcr" + model_id + "_wjetsMass",
        rl.Sample.BACKGROUND,
        templateRecoil(tmcr_wjets, "nominal"),
    )
    tmcr_wjetsMass.setParamEffect(lumi, 1.027)
    tmcr_wjetsMass.setParamEffect(trig_met, 1.01)
    tmcr_wjetsMass.setParamEffect(veto_tau, 1.03)
    tmcr_wjetsMass.setParamEffect(wjets_norm, 1.4)
    tmcr_wjetsMass.setParamEffect(jec, 1.05)
    tmcr_wjetsMass.setParamEffect(id_mu, 1.02)
    tmcr_wjetsMass.setParamEffect(iso_mu, 1.02)
    tmcr_wjetsMass.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["tmcr"][category]
    )
    btagUp = templateMass(tmcr_wjets, "btagUp")[0]
    btagDown = templateMass(tmcr_wjets, "btagDown")[0]
    tmcr_wjetsMass.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateMass(tmcr_wjets, "ew1Up")[0]
    ew1Down = templateMass(tmcr_wjets, "ew1Down")[0]
    tmcr_wjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(tmcr_wjets, "ew2GUp")[0]
    ew2GDown = templateMass(tmcr_wjets, "ew2GDown")[0]
    tmcr_wjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(tmcr_wjets, "ew2WUp")[0]
    ew2WDown = templateMass(tmcr_wjets, "ew2WDown")[0]
    tmcr_wjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(tmcr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(tmcr_wjets, "ew2ZDown")[0]
    tmcr_wjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(tmcr_wjets, "ew3GUp")[0]
    ew3GDown = templateMass(tmcr_wjets, "ew3GDown")[0]
    tmcr_wjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(tmcr_wjets, "ew3WUp")[0]
    ew3WDown = templateMass(tmcr_wjets, "ew3WDown")[0]
    tmcr_wjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(tmcr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(tmcr_wjets, "ew3ZDown")[0]
    tmcr_wjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(tmcr_wjets, "mixUp")[0]
    mixDown = templateMass(tmcr_wjets, "mixDown")[0]
    tmcr_wjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(tmcr_wjets, "muFUp")[0]
    muFDown = templateMass(tmcr_wjets, "muFDown")[0]
    tmcr_wjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(tmcr_wjets, "muRUp")[0]
    muRDown = templateMass(tmcr_wjets, "muRDown")[0]
    tmcr_wjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(tmcr_wjets, "qcd1Up")[0]
    qcd1Down = templateMass(tmcr_wjets, "qcd1Down")[0]
    tmcr_wjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(tmcr_wjets, "qcd2Up")[0]
    qcd2Down = templateMass(tmcr_wjets, "qcd2Down")[0]
    tmcr_wjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(tmcr_wjets, "qcd3Up")[0]
    qcd3Down = templateMass(tmcr_wjets, "qcd3Down")[0]
    tmcr_wjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    tmcr_wjetsRecoil = rl.TemplateSample(
        "tmcr" + model_id + "_wjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(tmcr_wjets, "nominal"),
    )
    tmcr_wjetsRecoil.setParamEffect(lumi, 1.027)
    tmcr_wjetsRecoil.setParamEffect(trig_met, 1.01)
    tmcr_wjetsRecoil.setParamEffect(veto_tau, 1.03)
    tmcr_wjetsRecoil.setParamEffect(wjets_norm, 1.4)
    tmcr_wjetsRecoil.setParamEffect(jec, 1.05)
    tmcr_wjetsRecoil.setParamEffect(id_mu, 1.02)
    tmcr_wjetsRecoil.setParamEffect(iso_mu, 1.02)
    tmcr_wjetsRecoil.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["tmcr"][category]
    )
    btagUp = templateRecoil(tmcr_wjets, "btagUp")[0]
    btagDown = templateRecoil(tmcr_wjets, "btagDown")[0]
    tmcr_wjetsRecoil.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateRecoil(tmcr_wjets, "ew1Up")[0]
    ew1Down = templateRecoil(tmcr_wjets, "ew1Down")[0]
    tmcr_wjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(tmcr_wjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(tmcr_wjets, "ew2GDown")[0]
    tmcr_wjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(tmcr_wjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(tmcr_wjets, "ew2WDown")[0]
    tmcr_wjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(tmcr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(tmcr_wjets, "ew2ZDown")[0]
    tmcr_wjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(tmcr_wjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(tmcr_wjets, "ew3GDown")[0]
    tmcr_wjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(tmcr_wjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(tmcr_wjets, "ew3WDown")[0]
    tmcr_wjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(tmcr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(tmcr_wjets, "ew3ZDown")[0]
    tmcr_wjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(tmcr_wjets, "mixUp")[0]
    mixDown = templateRecoil(tmcr_wjets, "mixDown")[0]
    tmcr_wjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(tmcr_wjets, "muFUp")[0]
    muFDown = templateRecoil(tmcr_wjets, "muFDown")[0]
    tmcr_wjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(tmcr_wjets, "muRUp")[0]
    muRDown = templateRecoil(tmcr_wjets, "muRDown")[0]
    tmcr_wjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(tmcr_wjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(tmcr_wjets, "qcd1Down")[0]
    tmcr_wjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(tmcr_wjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(tmcr_wjets, "qcd2Down")[0]
    tmcr_wjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(tmcr_wjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(tmcr_wjets, "qcd3Down")[0]
    tmcr_wjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    tecr_wjets = (
        hists["bkg"]["template"]
        .integrate("region", "tecr")
        .integrate("process", "W+jets")
        .sum("gentype")
    )
    tecr_wjetsMass = rl.TemplateSample(
        "tecr" + model_id + "_wjetsMass",
        rl.Sample.BACKGROUND,
        templateMass(tecr_wjets, "nominal"),
    )
    tecr_wjetsMass.setParamEffect(lumi, 1.027)
    tecr_wjetsMass.setParamEffect(trig_e, 1.01)
    tecr_wjetsMass.setParamEffect(veto_tau, 1.03)
    tecr_wjetsMass.setParamEffect(wjets_norm, 1.4)
    tecr_wjetsMass.setParamEffect(jec, 1.05)
    tecr_wjetsMass.setParamEffect(id_e, 1.02)
    tecr_wjetsMass.setParamEffect(reco_e, 1.02)
    tecr_wjetsMass.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["tecr"][category]
    )
    btagUp = templateMass(tecr_wjets, "btagUp")[0]
    btagDown = templateMass(tecr_wjets, "btagDown")[0]
    tecr_wjetsMass.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateMass(tecr_wjets, "ew1Up")[0]
    ew1Down = templateMass(tecr_wjets, "ew1Down")[0]
    tecr_wjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(tecr_wjets, "ew2GUp")[0]
    ew2GDown = templateMass(tecr_wjets, "ew2GDown")[0]
    tecr_wjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(tecr_wjets, "ew2WUp")[0]
    ew2WDown = templateMass(tecr_wjets, "ew2WDown")[0]
    tecr_wjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(tecr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(tecr_wjets, "ew2ZDown")[0]
    tecr_wjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(tecr_wjets, "ew3GUp")[0]
    ew3GDown = templateMass(tecr_wjets, "ew3GDown")[0]
    tecr_wjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(tecr_wjets, "ew3WUp")[0]
    ew3WDown = templateMass(tecr_wjets, "ew3WDown")[0]
    tecr_wjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(tecr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(tecr_wjets, "ew3ZDown")[0]
    tecr_wjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(tecr_wjets, "mixUp")[0]
    mixDown = templateMass(tecr_wjets, "mixDown")[0]
    tecr_wjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(tecr_wjets, "muFUp")[0]
    muFDown = templateMass(tecr_wjets, "muFDown")[0]
    tecr_wjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(tecr_wjets, "muRUp")[0]
    muRDown = templateMass(tecr_wjets, "muRDown")[0]
    tecr_wjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(tecr_wjets, "qcd1Up")[0]
    qcd1Down = templateMass(tecr_wjets, "qcd1Down")[0]
    tecr_wjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(tecr_wjets, "qcd2Up")[0]
    qcd2Down = templateMass(tecr_wjets, "qcd2Down")[0]
    tecr_wjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(tecr_wjets, "qcd3Up")[0]
    qcd3Down = templateMass(tecr_wjets, "qcd3Down")[0]
    tecr_wjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    tecr_wjetsRecoil = rl.TemplateSample(
        "tecr" + model_id + "_wjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(tecr_wjets, "nominal"),
    )
    tecr_wjetsRecoil.setParamEffect(lumi, 1.027)
    tecr_wjetsRecoil.setParamEffect(trig_e, 1.01)
    tecr_wjetsRecoil.setParamEffect(veto_tau, 1.03)
    tecr_wjetsRecoil.setParamEffect(wjets_norm, 1.4)
    tecr_wjetsRecoil.setParamEffect(jec, 1.05)
    tecr_wjetsRecoil.setParamEffect(id_e, 1.02)
    tecr_wjetsRecoil.setParamEffect(reco_e, 1.02)
    tecr_wjetsRecoil.setParamEffect(
        whf_fraction, hf_systematic["W+jets"]["tecr"][category]
    )
    btagUp = templateRecoil(tecr_wjets, "btagUp")[0]
    btagDown = templateRecoil(tecr_wjets, "btagDown")[0]
    tecr_wjetsRecoil.setParamEffect(btag, btagUp, btagDown)
    ew1Up = templateRecoil(tecr_wjets, "ew1Up")[0]
    ew1Down = templateRecoil(tecr_wjets, "ew1Down")[0]
    tecr_wjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(tecr_wjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(tecr_wjets, "ew2GDown")[0]
    tecr_wjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(tecr_wjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(tecr_wjets, "ew2WDown")[0]
    tecr_wjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(tecr_wjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(tecr_wjets, "ew2ZDown")[0]
    tecr_wjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(tecr_wjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(tecr_wjets, "ew3GDown")[0]
    tecr_wjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(tecr_wjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(tecr_wjets, "ew3WDown")[0]
    tecr_wjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(tecr_wjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(tecr_wjets, "ew3ZDown")[0]
    tecr_wjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(tecr_wjets, "mixUp")[0]
    mixDown = templateRecoil(tecr_wjets, "mixDown")[0]
    tecr_wjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(tecr_wjets, "muFUp")[0]
    muFDown = templateRecoil(tecr_wjets, "muFDown")[0]
    tecr_wjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(tecr_wjets, "muRUp")[0]
    muRDown = templateRecoil(tecr_wjets, "muRDown")[0]
    tecr_wjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(tecr_wjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(tecr_wjets, "qcd1Down")[0]
    tecr_wjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(tecr_wjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(tecr_wjets, "qcd2Down")[0]
    tecr_wjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(tecr_wjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(tecr_wjets, "qcd3Down")[0]
    tecr_wjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    ###
    # TT templates
    ###

    sr_tt = (
        hists["bkg"]["template"]
        .integrate("region", "sr")
        .integrate("process", "TT")
        .sum("gentype")
    )
    sr_ttMass = rl.TemplateSample(
        "sr" + model_id + "_ttMass",
        rl.Sample.BACKGROUND,
        templateMass(sr_tt, "nominal"),
    )
    sr_ttMass.setParamEffect(lumi, 1.027)
    sr_ttMass.setParamEffect(tt_norm, 1.2)
    sr_ttMass.setParamEffect(trig_met, 1.01)
    sr_ttMass.setParamEffect(veto_tau, 1.03)
    sr_ttMass.setParamEffect(jec, 1.05)
    btagUp = templateMass(sr_tt, "btagUp")[0]
    btagDown = templateMass(sr_tt, "btagDown")[0]
    sr_ttMass.setParamEffect(btag, btagUp, btagDown)
    sr_ttRecoil = rl.TemplateSample(
        "sr" + model_id + "_ttRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(sr_tt, "nominal"),
    )
    sr_ttRecoil.setParamEffect(lumi, 1.027)
    sr_ttRecoil.setParamEffect(tt_norm, 1.2)
    sr_ttRecoil.setParamEffect(trig_met, 1.01)
    sr_ttRecoil.setParamEffect(veto_tau, 1.03)
    sr_ttRecoil.setParamEffect(jec, 1.05)
    btagUp = templateRecoil(sr_tt, "btagUp")[0]
    btagDown = templateRecoil(sr_tt, "btagDown")[0]
    sr_ttRecoil.setParamEffect(btag, btagUp, btagDown)

    wmcr_tt = (
        hists["bkg"]["template"]
        .integrate("region", "wmcr")
        .integrate("process", "TT")
        .sum("gentype")
    )
    wmcr_ttMass = rl.TemplateSample(
        "wmcr" + model_id + "_ttMass",
        rl.Sample.BACKGROUND,
        templateMass(wmcr_tt, "nominal"),
    )
    wmcr_ttMass.setParamEffect(lumi, 1.027)
    wmcr_ttMass.setParamEffect(trig_met, 1.01)
    wmcr_ttMass.setParamEffect(veto_tau, 1.03)
    wmcr_ttMass.setParamEffect(tt_norm, 1.2)
    wmcr_ttMass.setParamEffect(jec, 1.05)
    wmcr_ttMass.setParamEffect(id_mu, 1.02)
    wmcr_ttMass.setParamEffect(iso_mu, 1.02)
    btagUp = templateMass(wmcr_tt, "btagUp")[0]
    btagDown = templateMass(wmcr_tt, "btagDown")[0]
    wmcr_ttMass.setParamEffect(btag, btagUp, btagDown)
    wmcr_ttRecoil = rl.TemplateSample(
        "wmcr" + model_id + "_ttRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(wmcr_tt, "nominal"),
    )
    wmcr_ttRecoil.setParamEffect(lumi, 1.027)
    wmcr_ttRecoil.setParamEffect(trig_met, 1.01)
    wmcr_ttRecoil.setParamEffect(veto_tau, 1.03)
    wmcr_ttRecoil.setParamEffect(tt_norm, 1.2)
    wmcr_ttRecoil.setParamEffect(jec, 1.05)
    wmcr_ttRecoil.setParamEffect(id_mu, 1.02)
    wmcr_ttRecoil.setParamEffect(iso_mu, 1.02)
    btagUp = templateRecoil(wmcr_tt, "btagUp")[0]
    btagDown = templateRecoil(wmcr_tt, "btagDown")[0]
    wmcr_ttRecoil.setParamEffect(btag, btagUp, btagDown)

    wecr_tt = (
        hists["bkg"]["template"]
        .integrate("region", "wecr")
        .integrate("process", "TT")
        .sum("gentype")
    )
    wecr_ttMass = rl.TemplateSample(
        "wecr" + model_id + "_ttMass",
        rl.Sample.BACKGROUND,
        templateMass(wecr_tt, "nominal"),
    )
    wecr_ttMass.setParamEffect(lumi, 1.027)
    wecr_ttMass.setParamEffect(trig_e, 1.01)
    wecr_ttMass.setParamEffect(veto_tau, 1.03)
    wecr_ttMass.setParamEffect(tt_norm, 1.2)
    wecr_ttMass.setParamEffect(jec, 1.05)
    wecr_ttMass.setParamEffect(id_e, 1.02)
    wecr_ttMass.setParamEffect(reco_e, 1.02)
    btagUp = templateMass(wecr_tt, "btagUp")[0]
    btagDown = templateMass(wecr_tt, "btagDown")[0]
    wecr_ttMass.setParamEffect(btag, btagUp, btagDown)
    wecr_ttRecoil = rl.TemplateSample(
        "wecr" + model_id + "_ttRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(wecr_tt, "nominal"),
    )
    wecr_ttRecoil.setParamEffect(lumi, 1.027)
    wecr_ttRecoil.setParamEffect(trig_e, 1.01)
    wecr_ttRecoil.setParamEffect(veto_tau, 1.03)
    wecr_ttRecoil.setParamEffect(tt_norm, 1.2)
    wecr_ttRecoil.setParamEffect(jec, 1.05)
    wecr_ttRecoil.setParamEffect(id_e, 1.02)
    wecr_ttRecoil.setParamEffect(reco_e, 1.02)
    btagUp = templateRecoil(wecr_tt, "btagUp")[0]
    btagDown = templateRecoil(wecr_tt, "btagDown")[0]
    wecr_ttRecoil.setParamEffect(btag, btagUp, btagDown)

    tmcr_tt = (
        hists["bkg"]["template"]
        .integrate("region", "tmcr")
        .integrate("process", "TT")
        .sum("gentype")
    )
    tmcr_ttMass = rl.TemplateSample(
        "tmcr" + model_id + "_ttMass",
        rl.Sample.BACKGROUND,
        templateMass(tmcr_tt, "nominal"),
    )
    tmcr_ttMass.setParamEffect(lumi, 1.027)
    tmcr_ttMass.setParamEffect(trig_met, 1.01)
    tmcr_ttMass.setParamEffect(veto_tau, 1.03)
    tmcr_ttMass.setParamEffect(tt_norm, 1.2)
    tmcr_ttMass.setParamEffect(jec, 1.05)
    tmcr_ttMass.setParamEffect(id_mu, 1.02)
    tmcr_ttMass.setParamEffect(iso_mu, 1.02)
    btagUp = templateMass(tmcr_tt, "btagUp")[0]
    btagDown = templateMass(tmcr_tt, "btagDown")[0]
    tmcr_ttMass.setParamEffect(btag, btagUp, btagDown)
    tmcr_ttRecoil = rl.TemplateSample(
        "tmcr" + model_id + "_ttRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(tmcr_tt, "nominal"),
    )
    tmcr_ttRecoil.setParamEffect(lumi, 1.027)
    tmcr_ttRecoil.setParamEffect(trig_met, 1.01)
    tmcr_ttRecoil.setParamEffect(veto_tau, 1.03)
    tmcr_ttRecoil.setParamEffect(tt_norm, 1.2)
    tmcr_ttRecoil.setParamEffect(jec, 1.05)
    tmcr_ttRecoil.setParamEffect(id_mu, 1.02)
    tmcr_ttRecoil.setParamEffect(iso_mu, 1.02)
    btagUp = templateRecoil(tmcr_tt, "btagUp")[0]
    btagDown = templateRecoil(tmcr_tt, "btagDown")[0]
    tmcr_ttRecoil.setParamEffect(btag, btagUp, btagDown)

    tecr_tt = (
        hists["bkg"]["template"]
        .integrate("region", "tecr")
        .integrate("process", "TT")
        .sum("gentype")
    )
    tecr_ttMass = rl.TemplateSample(
        "tecr" + model_id + "_ttMass",
        rl.Sample.BACKGROUND,
        templateMass(tecr_tt, "nominal"),
    )
    tecr_ttMass.setParamEffect(lumi, 1.027)
    tecr_ttMass.setParamEffect(trig_e, 1.01)
    tecr_ttMass.setParamEffect(veto_tau, 1.03)
    tecr_ttMass.setParamEffect(tt_norm, 1.2)
    tecr_ttMass.setParamEffect(jec, 1.05)
    tecr_ttMass.setParamEffect(id_e, 1.02)
    tecr_ttMass.setParamEffect(reco_e, 1.02)
    btagUp = templateMass(tecr_tt, "btagUp")[0]
    btagDown = templateMass(tecr_tt, "btagDown")[0]
    tecr_ttMass.setParamEffect(btag, btagUp, btagDown)
    tecr_ttRecoil = rl.TemplateSample(
        "tecr" + model_id + "_ttRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(tecr_tt, "nominal"),
    )
    tecr_ttRecoil.setParamEffect(lumi, 1.027)
    tecr_ttRecoil.setParamEffect(trig_e, 1.01)
    tecr_ttRecoil.setParamEffect(veto_tau, 1.03)
    tecr_ttRecoil.setParamEffect(tt_norm, 1.2)
    tecr_ttRecoil.setParamEffect(jec, 1.05)
    tecr_ttRecoil.setParamEffect(id_e, 1.02)
    tecr_ttRecoil.setParamEffect(reco_e, 1.02)
    btagUp = templateRecoil(tecr_tt, "btagUp")[0]
    btagDown = templateRecoil(tecr_tt, "btagDown")[0]
    tecr_ttRecoil.setParamEffect(btag, btagUp, btagDown)

    ###
    # DY+jets templates
    ###

    zmcr_dyjets = (
        hists["bkg"]["template"]
        .integrate("region", "zmcr")
        .integrate("process", "DY+jets")
        .sum("gentype")
    )
    zmcr_dyjetsMass = rl.TemplateSample(
        "zmcr" + model_id + "_dyjetsMass",
        rl.Sample.BACKGROUND,
        templateMass(zmcr_dyjets, "nominal"),
    )
    zmcr_dyjetsMass.setParamEffect(lumi, 1.027)
    zmcr_dyjetsMass.setParamEffect(trig_met, 1.01)
    zmcr_dyjetsMass.setParamEffect(veto_tau, 1.03)
    zmcr_dyjetsMass.setParamEffect(zjets_norm, 1.4)
    zmcr_dyjetsMass.setParamEffect(jec, 1.05)
    zmcr_dyjetsMass.setParamEffect(id_mu, 1.02)
    zmcr_dyjetsMass.setParamEffect(iso_mu, 1.02)
    zmcr_dyjetsMass.setParamEffect(
        zhf_fraction, hf_systematic["DY+jets"]["zmcr"][category]
    )
    ew1Up = templateMass(zmcr_dyjets, "ew1Up")[0]
    ew1Down = templateMass(zmcr_dyjets, "ew1Down")[0]
    zmcr_dyjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(zmcr_dyjets, "ew2GUp")[0]
    ew2GDown = templateMass(zmcr_dyjets, "ew2GDown")[0]
    zmcr_dyjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(zmcr_dyjets, "ew2WUp")[0]
    ew2WDown = templateMass(zmcr_dyjets, "ew2WDown")[0]
    zmcr_dyjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(zmcr_dyjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(zmcr_dyjets, "ew2ZDown")[0]
    zmcr_dyjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(zmcr_dyjets, "ew3GUp")[0]
    ew3GDown = templateMass(zmcr_dyjets, "ew3GDown")[0]
    zmcr_dyjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(zmcr_dyjets, "ew3WUp")[0]
    ew3WDown = templateMass(zmcr_dyjets, "ew3WDown")[0]
    zmcr_dyjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(zmcr_dyjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(zmcr_dyjets, "ew3ZDown")[0]
    zmcr_dyjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(zmcr_dyjets, "mixUp")[0]
    mixDown = templateMass(zmcr_dyjets, "mixDown")[0]
    zmcr_dyjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(zmcr_dyjets, "muFUp")[0]
    muFDown = templateMass(zmcr_dyjets, "muFDown")[0]
    zmcr_dyjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(zmcr_dyjets, "muRUp")[0]
    muRDown = templateMass(zmcr_dyjets, "muRDown")[0]
    zmcr_dyjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(zmcr_dyjets, "qcd1Up")[0]
    qcd1Down = templateMass(zmcr_dyjets, "qcd1Down")[0]
    zmcr_dyjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(zmcr_dyjets, "qcd2Up")[0]
    qcd2Down = templateMass(zmcr_dyjets, "qcd2Down")[0]
    zmcr_dyjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(zmcr_dyjets, "qcd3Up")[0]
    qcd3Down = templateMass(zmcr_dyjets, "qcd3Down")[0]
    zmcr_dyjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    zmcr_dyjetsRecoil = rl.TemplateSample(
        "zmcr" + model_id + "_dyjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(zmcr_dyjets, "nominal"),
    )
    zmcr_dyjetsRecoil.setParamEffect(lumi, 1.027)
    zmcr_dyjetsRecoil.setParamEffect(trig_met, 1.01)
    zmcr_dyjetsRecoil.setParamEffect(veto_tau, 1.03)
    zmcr_dyjetsRecoil.setParamEffect(zjets_norm, 1.4)
    zmcr_dyjetsRecoil.setParamEffect(jec, 1.05)
    zmcr_dyjetsRecoil.setParamEffect(id_mu, 1.02)
    zmcr_dyjetsRecoil.setParamEffect(iso_mu, 1.02)
    zmcr_dyjetsRecoil.setParamEffect(
        zhf_fraction, hf_systematic["DY+jets"]["zmcr"][category]
    )
    ew1Up = templateRecoil(zmcr_dyjets, "ew1Up")[0]
    ew1Down = templateRecoil(zmcr_dyjets, "ew1Down")[0]
    zmcr_dyjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(zmcr_dyjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(zmcr_dyjets, "ew2GDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(zmcr_dyjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(zmcr_dyjets, "ew2WDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(zmcr_dyjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(zmcr_dyjets, "ew2ZDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(zmcr_dyjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(zmcr_dyjets, "ew3GDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(zmcr_dyjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(zmcr_dyjets, "ew3WDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(zmcr_dyjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(zmcr_dyjets, "ew3ZDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(zmcr_dyjets, "mixUp")[0]
    mixDown = templateRecoil(zmcr_dyjets, "mixDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(zmcr_dyjets, "muFUp")[0]
    muFDown = templateRecoil(zmcr_dyjets, "muFDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(zmcr_dyjets, "muRUp")[0]
    muRDown = templateRecoil(zmcr_dyjets, "muRDown")[0]
    zmcr_dyjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(zmcr_dyjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(zmcr_dyjets, "qcd1Down")[0]
    zmcr_dyjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(zmcr_dyjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(zmcr_dyjets, "qcd2Down")[0]
    zmcr_dyjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(zmcr_dyjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(zmcr_dyjets, "qcd3Down")[0]
    zmcr_dyjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    zecr_dyjets = (
        hists["bkg"]["template"]
        .integrate("region", "zecr")
        .integrate("process", "DY+jets")
        .sum("gentype")
    )
    zecr_dyjetsMass = rl.TemplateSample(
        "zecr" + model_id + "_dyjetsMass",
        rl.Sample.BACKGROUND,
        templateMass(zecr_dyjets, "nominal"),
    )
    zecr_dyjetsMass.setParamEffect(lumi, 1.027)
    zecr_dyjetsMass.setParamEffect(trig_e, 1.01)
    zecr_dyjetsMass.setParamEffect(veto_tau, 1.03)
    zecr_dyjetsMass.setParamEffect(zjets_norm, 1.4)
    zecr_dyjetsMass.setParamEffect(jec, 1.05)
    zecr_dyjetsMass.setParamEffect(id_e, 1.02)
    zecr_dyjetsMass.setParamEffect(reco_e, 1.02)
    zecr_dyjetsMass.setParamEffect(
        zhf_fraction, hf_systematic["DY+jets"]["zecr"][category]
    )
    ew1Up = templateMass(zecr_dyjets, "ew1Up")[0]
    ew1Down = templateMass(zecr_dyjets, "ew1Down")[0]
    zecr_dyjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(zecr_dyjets, "ew2GUp")[0]
    ew2GDown = templateMass(zecr_dyjets, "ew2GDown")[0]
    zecr_dyjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(zecr_dyjets, "ew2WUp")[0]
    ew2WDown = templateMass(zecr_dyjets, "ew2WDown")[0]
    zecr_dyjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(zecr_dyjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(zecr_dyjets, "ew2ZDown")[0]
    zecr_dyjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(zecr_dyjets, "ew3GUp")[0]
    ew3GDown = templateMass(zecr_dyjets, "ew3GDown")[0]
    zecr_dyjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(zecr_dyjets, "ew3WUp")[0]
    ew3WDown = templateMass(zecr_dyjets, "ew3WDown")[0]
    zecr_dyjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(zecr_dyjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(zecr_dyjets, "ew3ZDown")[0]
    zecr_dyjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(zecr_dyjets, "mixUp")[0]
    mixDown = templateMass(zecr_dyjets, "mixDown")[0]
    zecr_dyjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(zecr_dyjets, "muFUp")[0]
    muFDown = templateMass(zecr_dyjets, "muFDown")[0]
    zecr_dyjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(zecr_dyjets, "muRUp")[0]
    muRDown = templateMass(zecr_dyjets, "muRDown")[0]
    zecr_dyjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(zecr_dyjets, "qcd1Up")[0]
    qcd1Down = templateMass(zecr_dyjets, "qcd1Down")[0]
    zecr_dyjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(zecr_dyjets, "qcd2Up")[0]
    qcd2Down = templateMass(zecr_dyjets, "qcd2Down")[0]
    zecr_dyjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(zecr_dyjets, "qcd3Up")[0]
    qcd3Down = templateMass(zecr_dyjets, "qcd3Down")[0]
    zecr_dyjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    zecr_dyjetsRecoil = rl.TemplateSample(
        "zecr" + model_id + "_dyjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(zecr_dyjets, "nominal"),
    )
    zecr_dyjetsRecoil.setParamEffect(lumi, 1.027)
    zecr_dyjetsRecoil.setParamEffect(trig_e, 1.01)
    zecr_dyjetsRecoil.setParamEffect(veto_tau, 1.03)
    zecr_dyjetsRecoil.setParamEffect(zjets_norm, 1.4)
    zecr_dyjetsRecoil.setParamEffect(jec, 1.05)
    zecr_dyjetsRecoil.setParamEffect(id_e, 1.02)
    zecr_dyjetsRecoil.setParamEffect(reco_e, 1.02)
    zecr_dyjetsRecoil.setParamEffect(
        zhf_fraction, hf_systematic["DY+jets"]["zecr"][category]
    )
    ew1Up = templateRecoil(zecr_dyjets, "ew1Up")[0]
    ew1Down = templateRecoil(zecr_dyjets, "ew1Down")[0]
    zecr_dyjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(zecr_dyjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(zecr_dyjets, "ew2GDown")[0]
    zecr_dyjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(zecr_dyjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(zecr_dyjets, "ew2WDown")[0]
    zecr_dyjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(zecr_dyjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(zecr_dyjets, "ew2ZDown")[0]
    zecr_dyjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(zecr_dyjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(zecr_dyjets, "ew3GDown")[0]
    zecr_dyjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(zecr_dyjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(zecr_dyjets, "ew3WDown")[0]
    zecr_dyjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(zecr_dyjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(zecr_dyjets, "ew3ZDown")[0]
    zecr_dyjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(zecr_dyjets, "mixUp")[0]
    mixDown = templateRecoil(zecr_dyjets, "mixDown")[0]
    zecr_dyjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(zecr_dyjets, "muFUp")[0]
    muFDown = templateRecoil(zecr_dyjets, "muFDown")[0]
    zecr_dyjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(zecr_dyjets, "muRUp")[0]
    muRDown = templateRecoil(zecr_dyjets, "muRDown")[0]
    zecr_dyjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(zecr_dyjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(zecr_dyjets, "qcd1Down")[0]
    zecr_dyjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(zecr_dyjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(zecr_dyjets, "qcd2Down")[0]
    zecr_dyjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(zecr_dyjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(zecr_dyjets, "qcd3Down")[0]
    zecr_dyjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    ###
    # G+jets templates
    ###

    gcr_gjets = (
        hists["bkg"]["template"]
        .integrate("region", "gcr")
        .integrate("process", "G+jets")
        .sum("gentype")
    )
    gcr_gjetsMass = rl.TemplateSample(
        "gcr" + model_id + "_gjetsMass",
        rl.Sample.BACKGROUND,
        templateMass(gcr_gjets, "nominal"),
    )
    gcr_gjetsMass.setParamEffect(lumi, 1.027)
    gcr_gjetsMass.setParamEffect(trig_pho, 1.01)
    gcr_gjetsMass.setParamEffect(veto_tau, 1.03)
    gcr_gjetsMass.setParamEffect(gjets_norm, 1.4)
    gcr_gjetsMass.setParamEffect(jec, 1.05)
    gcr_gjetsMass.setParamEffect(id_pho, 1.02)
    gcr_gjetsMass.setParamEffect(ghf_fraction, hf_systematic["G+jets"]["gcr"][category])
    ew1Up = templateMass(gcr_gjets, "ew1Up")[0]
    ew1Down = templateMass(gcr_gjets, "ew1Down")[0]
    gcr_gjetsMass.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateMass(gcr_gjets, "ew2GUp")[0]
    ew2GDown = templateMass(gcr_gjets, "ew2GDown")[0]
    gcr_gjetsMass.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateMass(gcr_gjets, "ew2WUp")[0]
    ew2WDown = templateMass(gcr_gjets, "ew2WDown")[0]
    gcr_gjetsMass.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateMass(gcr_gjets, "ew2ZUp")[0]
    ew2ZDown = templateMass(gcr_gjets, "ew2ZDown")[0]
    gcr_gjetsMass.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateMass(gcr_gjets, "ew3GUp")[0]
    ew3GDown = templateMass(gcr_gjets, "ew3GDown")[0]
    gcr_gjetsMass.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateMass(gcr_gjets, "ew3WUp")[0]
    ew3WDown = templateMass(gcr_gjets, "ew3WDown")[0]
    gcr_gjetsMass.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateMass(gcr_gjets, "ew3ZUp")[0]
    ew3ZDown = templateMass(gcr_gjets, "ew3ZDown")[0]
    gcr_gjetsMass.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateMass(gcr_gjets, "mixUp")[0]
    mixDown = templateMass(gcr_gjets, "mixDown")[0]
    gcr_gjetsMass.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateMass(gcr_gjets, "muFUp")[0]
    muFDown = templateMass(gcr_gjets, "muFDown")[0]
    gcr_gjetsMass.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateMass(gcr_gjets, "muRUp")[0]
    muRDown = templateMass(gcr_gjets, "muRDown")[0]
    gcr_gjetsMass.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateMass(gcr_gjets, "qcd1Up")[0]
    qcd1Down = templateMass(gcr_gjets, "qcd1Down")[0]
    gcr_gjetsMass.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateMass(gcr_gjets, "qcd2Up")[0]
    qcd2Down = templateMass(gcr_gjets, "qcd2Down")[0]
    gcr_gjetsMass.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateMass(gcr_gjets, "qcd3Up")[0]
    qcd3Down = templateMass(gcr_gjets, "qcd3Down")[0]
    gcr_gjetsMass.setParamEffect(qcd3, qcd3Up, qcd3Down)

    gcr_gjetsRecoil = rl.TemplateSample(
        "gcr" + model_id + "_gjetsRecoil",
        rl.Sample.BACKGROUND,
        templateRecoil(gcr_gjets, "nominal"),
    )
    gcr_gjetsRecoil.setParamEffect(lumi, 1.027)
    gcr_gjetsRecoil.setParamEffect(trig_pho, 1.01)
    gcr_gjetsRecoil.setParamEffect(veto_tau, 1.03)
    gcr_gjetsRecoil.setParamEffect(gjets_norm, 1.4)
    gcr_gjetsRecoil.setParamEffect(jec, 1.05)
    gcr_gjetsRecoil.setParamEffect(id_pho, 1.02)
    gcr_gjetsRecoil.setParamEffect(
        ghf_fraction, hf_systematic["G+jets"]["gcr"][category]
    )
    ew1Up = templateRecoil(gcr_gjets, "ew1Up")[0]
    ew1Down = templateRecoil(gcr_gjets, "ew1Down")[0]
    gcr_gjetsRecoil.setParamEffect(ew1, ew1Up, ew1Down)
    ew2GUp = templateRecoil(gcr_gjets, "ew2GUp")[0]
    ew2GDown = templateRecoil(gcr_gjets, "ew2GDown")[0]
    gcr_gjetsRecoil.setParamEffect(ew2G, ew2GUp, ew2GDown)
    ew2WUp = templateRecoil(gcr_gjets, "ew2WUp")[0]
    ew2WDown = templateRecoil(gcr_gjets, "ew2WDown")[0]
    gcr_gjetsRecoil.setParamEffect(ew2W, ew2WUp, ew2WDown)
    ew2ZUp = templateRecoil(gcr_gjets, "ew2ZUp")[0]
    ew2ZDown = templateRecoil(gcr_gjets, "ew2ZDown")[0]
    gcr_gjetsRecoil.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    ew3GUp = templateRecoil(gcr_gjets, "ew3GUp")[0]
    ew3GDown = templateRecoil(gcr_gjets, "ew3GDown")[0]
    gcr_gjetsRecoil.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = templateRecoil(gcr_gjets, "ew3WUp")[0]
    ew3WDown = templateRecoil(gcr_gjets, "ew3WDown")[0]
    gcr_gjetsRecoil.setParamEffect(ew3W, ew3WUp, ew3WDown)
    ew3ZUp = templateRecoil(gcr_gjets, "ew3ZUp")[0]
    ew3ZDown = templateRecoil(gcr_gjets, "ew3ZDown")[0]
    gcr_gjetsRecoil.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    mixUp = templateRecoil(gcr_gjets, "mixUp")[0]
    mixDown = templateRecoil(gcr_gjets, "mixDown")[0]
    gcr_gjetsRecoil.setParamEffect(mix, mixUp, mixDown)
    muFUp = templateRecoil(gcr_gjets, "muFUp")[0]
    muFDown = templateRecoil(gcr_gjets, "muFDown")[0]
    gcr_gjetsRecoil.setParamEffect(muF, muFUp, muFDown)
    muRUp = templateRecoil(gcr_gjets, "muRUp")[0]
    muRDown = templateRecoil(gcr_gjets, "muRDown")[0]
    gcr_gjetsRecoil.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = templateRecoil(gcr_gjets, "qcd1Up")[0]
    qcd1Down = templateRecoil(gcr_gjets, "qcd1Down")[0]
    gcr_gjetsRecoil.setParamEffect(qcd1, qcd1Up, qcd1Down)
    qcd2Up = templateRecoil(gcr_gjets, "qcd2Up")[0]
    qcd2Down = templateRecoil(gcr_gjets, "qcd2Down")[0]
    gcr_gjetsRecoil.setParamEffect(qcd2, qcd2Up, qcd2Down)
    qcd3Up = templateRecoil(gcr_gjets, "qcd3Up")[0]
    qcd3Down = templateRecoil(gcr_gjets, "qcd3Down")[0]
    gcr_gjetsRecoil.setParamEffect(qcd3, qcd3Up, qcd3Down)

    ###
    # Compute TFs
    ###

    sr_wjetsTransferFactor = (
        sr_wjetsMass.getExpectation() * sr_wjetsRecoil.getExpectation()
    ) / (sr_zjetsMass.getExpectation() * sr_zjetsRecoil.getExpectation())
    wmcr_wjetsTransferFactor = (
        wmcr_wjetsMass.getExpectation() * wmcr_wjetsRecoil.getExpectation()
    ) / (sr_wjetsMass.getExpectation() * sr_wjetsRecoil.getExpectation())
    wmcr_ttTransferFactor = (
        wmcr_ttMass.getExpectation() * wmcr_ttRecoil.getExpectation()
    ) / (sr_ttMass.getExpectation() * sr_ttRecoil.getExpectation())
    tmcr_wjetsTransferFactor = (
        tmcr_wjetsMass.getExpectation() * tmcr_wjetsRecoil.getExpectation()
    ) / (sr_wjetsMass.getExpectation() * sr_wjetsRecoil.getExpectation())
    tmcr_ttTransferFactor = (
        tmcr_ttMass.getExpectation() * tmcr_ttRecoil.getExpectation()
    ) / (sr_ttMass.getExpectation() * sr_ttRecoil.getExpectation())
    wecr_wjetsTransferFactor = (
        wecr_wjetsMass.getExpectation() * wecr_wjetsRecoil.getExpectation()
    ) / (sr_wjetsMass.getExpectation() * sr_wjetsRecoil.getExpectation())
    wecr_ttTransferFactor = (
        wecr_ttMass.getExpectation() * wecr_ttRecoil.getExpectation()
    ) / (sr_ttMass.getExpectation() * sr_ttRecoil.getExpectation())
    tecr_wjetsTransferFactor = (
        tecr_wjetsMass.getExpectation() * tecr_wjetsRecoil.getExpectation()
    ) / (sr_wjetsMass.getExpectation() * sr_wjetsRecoil.getExpectation())
    tecr_ttTransferFactor = (
        tecr_ttMass.getExpectation() * tecr_ttRecoil.getExpectation()
    ) / (sr_ttMass.getExpectation() * sr_ttRecoil.getExpectation())
    zmcr_dyjetsTransferFactor = (
        zmcr_dyjetsMass.getExpectation() * zmcr_dyjetsRecoil.getExpectation()
    ) / (sr_zjetsMass.getExpectation() * sr_zjetsRecoil.getExpectation())
    zecr_dyjetsTransferFactor = (
        zecr_dyjetsMass.getExpectation() * zecr_dyjetsRecoil.getExpectation()
    ) / (sr_zjetsMass.getExpectation() * sr_zjetsRecoil.getExpectation())
    gcr_gjetsTransferFactor = (
        gcr_gjetsMass.getExpectation() * gcr_gjetsRecoil.getExpectation()
    ) / (sr_zjetsMass.getExpectation() * sr_zjetsRecoil.getExpectation())

    # sr_wjetsTransferFactor = sr_wjetsRecoil.getExpectation() / sr_zjetsRecoil.getExpectation()
    # wmcr_wjetsTransferFactor = wmcr_wjetsRecoil.getExpectation() / sr_wjetsRecoil.getExpectation()
    # tmcr_wjetsTransferFactor = tmcr_wjetsRecoil.getExpectation() / sr_wjetsRecoil.getExpectation()
    # wecr_wjetsTransferFactor = wecr_wjetsRecoil.getExpectation() / sr_wjetsRecoil.getExpectation()
    # tecr_wjetsTransferFactor = tecr_wjetsRecoil.getExpectation() / sr_wjetsRecoil.getExpectation()
    # zmcr_dyjetsTransferFactor = zmcr_dyjetsRecoil.getExpectation() / sr_zjetsRecoil.getExpectation()
    # zecr_dyjetsTransferFactor = zecr_dyjetsRecoil.getExpectation() / sr_zjetsRecoil.getExpectation()
    # gcr_gjetsTransferFactor = gcr_gjetsRecoil.getExpectation() / sr_zjetsRecoil.getExpectation()

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

    sr.setObservation(template(data, "MET", "data", "sr"))

    ###
    # Z(->nunu)+jets data-driven model
    ###
    sr_zjetsTemplate = template(background, "Z+jets", "nominal", "sr")
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
    # top-antitop data-driven model
    ###

    sr_ttTemplate = template(background, "TT", "nominal", "sr")
    sr_ttObservable = rl.Observable("fjmass", sr_ttTemplate[1])
    sr_tt = rl.ParametericSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields
    )
    sr.addSample(sr_tt)

    ###
    # Other MC-driven processes
    ###

    sr_stTemplate = template(background, "ST", "nominal", "sr")
    sr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, sr_stTemplate)
    sr_st.setParamEffect(lumi, 1.027)
    sr_st.setParamEffect(trig_met, 1.01)
    sr_st.setParamEffect(veto_tau, 1.03)
    sr_st.setParamEffect(st_norm, 1.2)
    sr_st.setParamEffect(jec, 1.05)
    btagUp = template(background, "ST", "btagUp", "sr")[0]
    btagDown = template(background, "ST", "btagDown", "sr")[0]
    sr_st.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_st)

    sr_dyjetsTemplate = template(background, "DY+jets", "nominal", "sr")
    sr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, sr_dyjetsTemplate
    )
    sr_dyjets.setParamEffect(lumi, 1.027)
    sr_dyjets.setParamEffect(trig_met, 1.01)
    sr_dyjets.setParamEffect(veto_tau, 1.03)
    sr_dyjets.setParamEffect(zjets_norm, 1.4)
    sr_dyjets.setParamEffect(jec, 1.05)
    btagUp = template(background, "DY+jets", "btagUp", "sr")[0]
    btagDown = template(background, "DY+jets", "btagDown", "sr")[0]
    btagDown[btagDown <= 0] = 1e-7
    sr_dyjets.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_dyjets)

    sr_vvTemplate = template(background, "VV", "nominal", "sr")
    sr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, sr_vvTemplate)
    sr_vv.setParamEffect(lumi, 1.027)
    sr_vv.setParamEffect(trig_met, 1.01)
    sr_vv.setParamEffect(veto_tau, 1.03)
    sr_vv.setParamEffect(vv_norm, 1.2)
    sr_vv.setParamEffect(jec, 1.05)
    btagUp = template(background, "VV", "btagUp", "sr")[0]
    btagDown = template(background, "VV", "btagDown", "sr")[0]
    btagDown[btagDown <= 0] = 1e-7
    sr_vv.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_vv)

    sr_hbbTemplate = template(background, "Hbb", "nominal", "sr")
    sr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, sr_hbbTemplate)
    sr_hbb.setParamEffect(lumi, 1.027)
    sr_hbb.setParamEffect(trig_met, 1.01)
    sr_hbb.setParamEffect(veto_tau, 1.03)
    sr_hbb.setParamEffect(hbb_norm, 1.2)
    sr_hbb.setParamEffect(jec, 1.05)
    btagUp = template(background, "Hbb", "btagUp", "sr")[0]
    btagDown = template(background, "Hbb", "btagDown", "sr")[0]
    sr_hbb.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_hbb)

    sr_qcdTemplate = template(background, "QCD", "nominal", "sr")
    sr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, sr_qcdTemplate)
    sr_qcd.setParamEffect(lumi, 1.027)
    sr_qcd.setParamEffect(trig_met, 1.01)
    sr_qcd.setParamEffect(veto_tau, 1.03)
    sr_qcd.setParamEffect(qcdsig_norm, 2.0)
    sr_qcd.setParamEffect(jec, 1.05)
    btagUp = template(background, "QCD", "btagUp", "sr")[0]
    btagDown = template(background, "QCD", "btagDown", "sr")[0]
    sr_qcd.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_qcd)

    for s in signal["sr"].identifiers("process"):
        if "Mhs_50" not in str(s):
            continue
        sr_signalTemplate = template(signal, s, "nominal", "sr")
        sr_signal = rl.TemplateSample(
            ch_name + "_" + str(s), rl.Sample.SIGNAL, sr_signalTemplate
        )
        sr_signal.setParamEffect(lumi, 1.027)
        sr_signal.setParamEffect(trig_met, 1.01)
        sr_signal.setParamEffect(veto_tau, 1.03)
        sr_signal.setParamEffect(jec, 1.05)
        btagUp = template(signal, s, "btagUp", "sr")[0]
        btagDown = template(signal, s, "btagDown", "sr")[0]
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

    wmcr.setObservation(template(data, "MET", "data", "wmcr"))

    ###
    # W(->lnu)+jets data-driven model
    ###

    wmcr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, wmcr_wjetsTransferFactor, sr_wjets
    )
    wmcr.addSample(wmcr_wjets)

    ###
    # top-antitop data-driven model
    ###

    wmcr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt
    )
    wmcr.addSample(wmcr_tt)

    ###
    # Other MC-driven processes
    ###

    wmcr_stTemplate = template(background, "ST", "nominal", "wmcr")
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
    btagUp = template(background, "ST", "btagUp", "wmcr")[0]
    btagDown = template(background, "ST", "btagDown", "wmcr")[0]
    wmcr_st.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_st)

    wmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", "wmcr")
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
    btagUp = template(background, "DY+jets", "btagUp", "wmcr")[0]
    btagDown = template(background, "DY+jets", "btagDown", "wmcr")[0]
    wmcr_dyjets.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_dyjets)

    wmcr_vvTemplate = template(background, "VV", "nominal", "wmcr")
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
    btagUp = template(background, "VV", "btagUp", "wmcr")[0]
    btagDown = template(background, "VV", "btagDown", "wmcr")[0]
    wmcr_vv.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_vv)

    wmcr_hbbTemplate = template(background, "Hbb", "nominal", "wmcr")
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
    btagUp = template(background, "Hbb", "btagUp", "wmcr")[0]
    btagDown = template(background, "Hbb", "btagDown", "wmcr")[0]
    wmcr_hbb.setParamEffect(btag, btagUp, btagDown)
    wmcr.addSample(wmcr_hbb)

    wmcr_qcdTemplate = template(background, "QCD", "nominal", "wmcr")
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
    btagUp = template(background, "QCD", "btagUp", "wmcr")[0]
    btagDown = template(background, "QCD", "btagDown", "wmcr")[0]
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

    tmcr.setObservation(template(data, "MET", "data", "tmcr"))

    ###
    # W(->lnu)+jets data-driven model
    ###

    tmcr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, tmcr_wjetsTransferFactor, sr_wjets
    )
    tmcr.addSample(tmcr_wjets)

    ###
    # top-antitop data-driven model
    ###

    tmcr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt
    )
    tmcr.addSample(tmcr_tt)

    ###
    # Other MC-driven processes
    ###

    tmcr_stTemplate = template(background, "ST", "nominal", "tmcr")
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
    btagUp = template(background, "ST", "btagUp", "tmcr")[0]
    btagDown = template(background, "ST", "btagDown", "tmcr")[0]
    tmcr_st.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_st)

    tmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", "tmcr")
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
    btagUp = template(background, "DY+jets", "btagUp", "tmcr")[0]
    btagDown = template(background, "DY+jets", "btagDown", "tmcr")[0]
    tmcr_dyjets.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_dyjets)

    tmcr_vvTemplate = template(background, "VV", "nominal", "tmcr")
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
    btagUp = template(background, "VV", "btagUp", "tmcr")[0]
    btagDown = template(background, "VV", "btagDown", "tmcr")[0]
    tmcr_vv.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_vv)

    tmcr_hbbTemplate = template(background, "Hbb", "nominal", "tmcr")
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
    btagUp = template(background, "Hbb", "btagUp", "tmcr")[0]
    btagDown = template(background, "Hbb", "btagDown", "tmcr")[0]
    tmcr_hbb.setParamEffect(btag, btagUp, btagDown)
    tmcr.addSample(tmcr_hbb)

    tmcr_qcdTemplate = template(background, "QCD", "nominal", "tmcr")
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
    btagUp = template(background, "QCD", "btagUp", "tmcr")[0]
    btagDown = template(background, "QCD", "btagDown", "tmcr")[0]
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
        wecr.setObservation(template(data, "EGamma", "data", "wecr"))
    else:
        wecr.setObservation(template(data, "SingleElectron", "data", "wecr"))

    ###
    # W(->lnu)+jets data-driven model
    ###

    wecr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, wecr_wjetsTransferFactor, sr_wjets
    )
    wecr.addSample(wecr_wjets)

    ###
    # top-antitop data-driven model
    ###

    wecr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt
    )
    wecr.addSample(wecr_tt)

    ###
    # Other MC-driven processes
    ###

    wecr_stTemplate = template(background, "ST", "nominal", "wecr")
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
    btagUp = template(background, "ST", "btagUp", "wecr")[0]
    btagDown = template(background, "ST", "btagDown", "wecr")[0]
    wecr_st.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_st)

    wecr_dyjetsTemplate = template(background, "DY+jets", "nominal", "wecr")
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
    btagUp = template(background, "DY+jets", "btagUp", "wecr")[0]
    btagDown = template(background, "DY+jets", "btagDown", "wecr")[0]
    wecr_dyjets.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_dyjets)

    wecr_vvTemplate = template(background, "VV", "nominal", "wecr")
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
    btagUp = template(background, "VV", "btagUp", "wecr")[0]
    btagDown = template(background, "VV", "btagDown", "wecr")[0]
    wecr_vv.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_vv)

    wecr_hbbTemplate = template(background, "Hbb", "nominal", "wecr")
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
    btagUp = template(background, "Hbb", "btagUp", "wecr")[0]
    btagDown = template(background, "Hbb", "btagDown", "wecr")[0]
    wecr_hbb.setParamEffect(btag, btagUp, btagDown)
    wecr.addSample(wecr_hbb)

    wecr_qcdTemplate = template(background, "QCD", "nominal", "wecr")
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
    btagUp = template(background, "QCD", "btagUp", "wecr")[0]
    btagDown = template(background, "QCD", "btagDown", "wecr")[0]
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
        tecr.setObservation(template(data, "EGamma", "data", "tecr"))
    else:
        tecr.setObservation(template(data, "SingleElectron", "data", "tecr"))

    ###
    # W(->lnu)+jets data-driven model
    ###

    tecr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, tecr_wjetsTransferFactor, sr_wjets
    )
    tecr.addSample(tecr_wjets)

    ###
    # top-antitop data-driven model
    ###

    tecr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt
    )
    tecr.addSample(tecr_tt)

    ###
    # Other MC-driven processes
    ###

    tecr_stTemplate = template(background, "ST", "nominal", "tecr")
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
    btagUp = template(background, "ST", "btagUp", "tecr")[0]
    btagDown = template(background, "ST", "btagDown", "tecr")[0]
    tecr_st.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_st)

    tecr_dyjetsTemplate = template(background, "DY+jets", "nominal", "tecr")
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
    btagUp = template(background, "DY+jets", "btagUp", "tecr")[0]
    btagDown = template(background, "DY+jets", "btagDown", "tecr")[0]
    tecr_dyjets.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_dyjets)

    tecr_vvTemplate = template(background, "VV", "nominal", "tecr")
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
    btagUp = template(background, "VV", "btagUp", "tecr")[0]
    btagDown = template(background, "VV", "btagDown", "tecr")[0]
    tecr_vv.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_vv)

    tecr_hbbTemplate = template(background, "Hbb", "nominal", "tecr")
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
    btagUp = template(background, "Hbb", "btagUp", "tecr")[0]
    btagDown = template(background, "Hbb", "btagDown", "tecr")[0]
    tecr_hbb.setParamEffect(btag, btagUp, btagDown)
    tecr.addSample(tecr_hbb)

    tecr_qcdTemplate = template(background, "QCD", "nominal", "tecr")
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
    btagUp = template(background, "QCD", "btagUp", "tecr")[0]
    btagDown = template(background, "QCD", "btagDown", "tecr")[0]
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

    zmcr.setObservation(template(data, "MET", "data", "zmcr"))

    zmcr_dyjets = rl.TransferFactorSample(
        ch_name + "_dyjets", rl.Sample.BACKGROUND, zmcr_dyjetsTransferFactor, sr_zjets
    )
    zmcr.addSample(zmcr_dyjets)

    ###
    # Other MC-driven processes
    ###

    zmcr_ttTemplate = template(background, "TT", "nominal", "zmcr")
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

    zmcr_stTemplate = template(background, "ST", "nominal", "zmcr")
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

    zmcr_vvTemplate = template(background, "VV", "nominal", "zmcr")
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

    zmcr_hbbTemplate = template(background, "Hbb", "nominal", "zmcr")
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
        zecr.setObservation(template(data, "EGamma", "data", "zecr"))
    else:
        zecr.setObservation(template(data, "SingleElectron", "data", "zecr"))

    zecr_dyjets = rl.TransferFactorSample(
        ch_name + "_dyjets", rl.Sample.BACKGROUND, zecr_dyjetsTransferFactor, sr_zjets
    )
    zecr.addSample(zecr_dyjets)

    ###
    # Other MC-driven processes
    ###

    zecr_ttTemplate = template(background, "TT", "nominal", "zecr")
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

    zecr_stTemplate = template(background, "ST", "nominal", "zecr")
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

    zecr_vvTemplate = template(background, "VV", "nominal", "zecr")
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

    zecr_hbbTemplate = template(background, "Hbb", "nominal", "zecr")
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
        gcr.setObservation(template(data, "EGamma", "data", "gcr"))
    else:
        gcr.setObservation(template(data, "SinglePhoton", "data", "gcr"))

    gcr_gjets = rl.TransferFactorSample(
        ch_name + "_gjets", rl.Sample.BACKGROUND, gcr_gjetsTransferFactor, sr_zjets
    )
    gcr.addSample(gcr_gjets)

    gcr_qcdTemplate = template(background, "QCD", "nominal", "gcr")
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
        sr_ttShape,
        sr_ttRate,
        tt_weight,
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
            sr_ttBinYields = (
                sr_ttShape[category] * sr_ttRate[recoilbin] * tt_weight[category]
            )
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
