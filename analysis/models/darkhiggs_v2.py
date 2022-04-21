from __future__ import print_function, division
from optparse import OptionParser
from collections import defaultdict, OrderedDict
import concurrent.futures
import sys
import os
import rhalphalib as rl
import numpy as np
import pickle
import gzip
import json
from coffea import hist, processor
from coffea.util import load, save
from scipy import stats
import ROOT

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

mass_binning = [
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
recoil_binning_dict = {
    "2018": [250, 310, 370, 470, 590, 3000],
    "2017": [250, 310, 370, 470, 590, 3000],
    "2016": [250, 310, 370, 470, 590, 3000]
}

category_map = {"pass": 1, "fail": 0}

def template(dictionary, process, systematic, recoil, region, category, read_sumw2=False, bkg=False):
    histogram = None
    if bkg:
        histogram = dictionary[region].integrate("process")
    else:
        histogram = dictionary[region].integrate("process", process)
    nominal, sumw2 = histogram.integrate("systematic", "nominal").values(sumw2=True)[()]
    nominal=nominal[recoil, :, category_map[category]]
    sumw2=sumw2[recoil, :, category_map[category]]
    zerobins = nominal <= 0.
    output = nominal
    if "data" not in systematic:
        if "Z+" in str(process):
            output[zerobins] = 1.
            sumw2[zerobins] = 0.
        elif "W+" in str(process) and "fail" in category:
            output[zerobins] = 1.
            sumw2[zerobins] = 0.
        elif "W+" in str(process) and "pass" in category and recoil<4:
            output[zerobins] = 1.
            sumw2[zerobins] = 0.
        elif "TT" in str(process) and "pass" in category  and recoil<4:
            output[zerobins] = 1.
            sumw2[zerobins] = 0.
        else:
            output[zerobins] = 1e-5
            sumw2[zerobins] = 0.
    if "nominal" not in systematic and "data" not in systematic:
        output = histogram.integrate("systematic", systematic).values()[()][recoil, :, category_map[category]]
        output[zerobins] = 1.
        output[~zerobins] /= nominal[~zerobins]
        output[np.isnan(output)] = 1.
    binning = (
        dictionary[region]
        .integrate("process", process)
        .integrate("systematic", systematic)
        .axis("fjmass")
        .edges()
    )
    if read_sumw2:
        return (output, binning, "fjmass", sumw2)
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
    bkg_map["W+HF"] = ("W+HF",)
    bkg_map["W+LF"] = ("W+LF",)
    bkg_map["Z+jets"] = ("Z+*",)
    bkg_map["Z+HF"] = ("Z+HF",)
    bkg_map["Z+LF"] = ("Z+LF",)
    bkg_map["G+jets"] = ("G+*",)
    bkg_map["QCD"] = ("QCD*",)
    data_map["MET"] = ("MET",)
    data_map["SingleElectron"] = ("SingleElectron",)
    data_map["SinglePhoton"] = ("SinglePhoton",)
    data_map["EGamma"] = ("EGamma",)
    for signal in hists['sig']['template'].identifiers('process'):
        if 'mhs' not in str(signal): continue
        sig_map[str(signal)] = (str(signal),)  ## signals

    for key in hists["data"].keys():
        bkg_hists[key] = hists["bkg"][key].group(cats, process, bkg_map)
        signal_hists[key] = hists["sig"][key].group(cats, process, sig_map)
        data_hists[key] = hists["data"][key].group(cats, process, data_map)

    print('initial recoil binning',bkg_hists["template"].axis("recoil").edges())

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

def addBtagSyst(dictionary, recoil, process, region, templ, category):
    btagUp = template(dictionary, process, "btagUp", recoil, region, category)[0]
    btagDown = template(dictionary, process, "btagDown", recoil, region, category)[0]
    templ.setParamEffect(btag, btagUp, btagDown)

def addVJetsSyst(dictionary, recoil, process, region, templ, category):
    def addSyst(dictionary, recoil, process, region, templ, category, syst, string):
        histogram = dictionary[region].integrate("process", process)
        nominal=histogram.integrate("systematic", "nominal").values()[()][recoil, :, category_map[category]]
        up=histogram.integrate("systematic", string+"Up").values()[()][recoil, :, category_map[category]]
        down=histogram.integrate("systematic",string+"Down").values()[()][recoil, :, category_map[category]]
        systUp = np.array( up.sum() / nominal.sum() )
        systUp[np.isnan(systUp)] = 1.
        systUp = systUp.sum()
        templ.setParamEffect(syst, systUp)
    addSyst(dictionary, recoil, process, region, templ, category, ew1, "ew1")
    addSyst(dictionary, recoil, process, region, templ, category, ew2W, "ew2W")
    addSyst(dictionary, recoil, process, region, templ, category, ew2Z, "ew2Z")
    addSyst(dictionary, recoil, process, region, templ, category, ew3W, "ew3W")
    addSyst(dictionary, recoil, process, region, templ, category, ew3Z, "ew3Z")
    addSyst(dictionary, recoil, process, region, templ, category, mix, "mix")
    addSyst(dictionary, recoil, process, region, templ, category, qcd1, "qcd1")
    addSyst(dictionary, recoil, process, region, templ, category, qcd2, "qcd2")
    addSyst(dictionary, recoil, process, region, templ, category, qcd3, "qcd3")

def addCommonRateSyst(region, templ):
    print('Common rate systematics have been applied')
    if region == 'sr':
        templ.setParamEffect(lumi, nlumi)
        templ.setParamEffect(trig_met, ntrig_met)
        templ.setParamEffect(veto_tau, nveto_tau)
        templ.setParamEffect(jec, njec)
    elif region == 'wecr' or region == 'tecr':
        templ.setParamEffect(lumi, nlumi)
        templ.setParamEffect(trig_e, ntrig_e)
        templ.setParamEffect(veto_tau, nveto_tau)
        templ.setParamEffect(jec, njec)
        templ.setParamEffect(id_e, nlepton)
        templ.setParamEffect(reco_e, nlepton)
    elif region == 'wmcr' or region == 'tmcr':
        templ.setParamEffect(lumi, nlumi)
        templ.setParamEffect(trig_met, ntrig_met)
        templ.setParamEffect(veto_tau, nveto_tau)
        templ.setParamEffect(jec, njec)
        templ.setParamEffect(id_mu, nlepton)
        templ.setParamEffect(iso_mu, nlepton)

def otherMCdrivenProcess(dictionary, recoil, region, channel, category, ch_name, param):
    processes = ["ST", "DY+jets", "VV", "Hbb", "QCD"]
    supplementNames = ["_stMC", "_dyjetsMC", "_vvMC", "_hbbMC", "_qcdMC"]
    for i in range(5):
        templ = template(dictionary, processes[i], "nominal", recoil, region, category, read_sumw2=True)
        templSample = rl.TemplateSample(ch_name + supplementNames[i], rl.Sample.BACKGROUND, templ)
        if options.rateSyst:
            addCommonRateSyst(region, templSample)
            if i == 0:
                templSample.setParamEffect(st_norm, nMinor_norm)
            elif i == 1:
                templSample.setParamEffect(zjetsMC_norm, nVjets_norm)
            elif i == 2:
                templSample.setParamEffect(vv_norm, nMinor_norm)
            elif i == 3:
                templSample.setParamEffect(hbb_norm, nMinor_norm)
            elif i == 4:
                if region == 'sr':
                    templSample.setParamEffect(qcdsig_norm, nqcd_norm)
                elif region == 'wecr' or region == 'tecr':
                    templSample.setParamEffect(qcde_norm, nqcd_norm)
                elif region == 'wmcr' or region == 'tmcr':
                    templSample.setParamEffect(qcdmu_norm, nqcd_norm)

        if options.btagVjets:
            if i == 1;
                addBtagSyst(background, recoil, processes[i], region, templSample, category)
                addVJetsSyst(background, recoil, processes[i], region, templSample, category)
            else:
                addBtagSyst(background, recoil, processes[i], region, templSample, category)

        if options.statUncs:
            addBBliteSyst(templSample, param, epsilon=1e-5) ### replace autoMCStats

        channel.addSample(templSample)

def addBBliteSyst(templ, param, epsilon=0):
    for i in range(templ.observable.nbins):
        if templ._nominal[i] <= 0. or templ._sumw2[i] <= 0.:
            continue
        effect_up = np.ones_like(templ._nominal)
        effect_down = np.ones_like(templ._nominal)
        effect_up[i] = (templ._nominal[i] + np.sqrt(templ._sumw2[i]))/templ._nominal[i]
        effect_down[i] = max((templ._nominal[i] - np.sqrt(templ._sumw2[i]))/templ._nominal[i], epsilon)
        templ.setParamEffect(param[i], effect_up, effect_down)

def model(year, recoil, category, s):

    model_id = year + category + "recoil" + str(recoil)
    print(model_id)
    model = rl.Model(str(s) + model_id)

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

    if category == 'pass' and options.sumbkg:
        dataTemplate = template(fake_data, "", "data", recoil, "sr", category, bkg=True)
    else:
        dataTemplate = template(data, "MET", "data", recoil, "sr", category)

    sr.setObservation(dataTemplate)
    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')

    ###
    # Z(->nunu)+jets data-driven model
    ###

    if category == "pass":
        sr_zjetsMC = sr_zjetsMCPass
        sr_zjets = sr_zjetsPass
    else:
        sr_zjetsMC = sr_zjetsMCFail
        sr_zjets = sr_zjetsFail
    sr.addSample(sr_zjets)

    ###
    # W(->lnu)+jets data-driven model
    ###

    if category == "pass":
        sr_wjetsMC = sr_wjetsMCPass
        sr_wjetsTemplate = sr_wjetsMCPassTemplate
        sr_wjets = sr_wjetsPass
        if not (recoil<4):
            sr_wjetsMC = rl.TemplateSample(
                "sr" + model_id + "_wjetsMC",
                rl.Sample.BACKGROUND,
                sr_wjetsTemplate
            )
            if options.rateSyst:
                addCommonRateSyst('sr', sr_wjetsMC)
                sr_wjetsMC.setParamEffect(wjetsMC_norm, nVjets_norm)
            if options.statUncs:
                addBBliteSyst(sr_wjetsMC, param, epsilon=1e-5) ### replace autoMCStats
            if options.btagVjets:
                addBtagSyst(background, recoil, "W+jets", "sr", sr_wjetsMC, category)
                addVJetsSyst(background, recoil, "W+jets", "sr", sr_wjetsMC, category)
            sr_wjets = sr_wjetsMC
    else:
        sr_wjetsTemplate = sr_wjetsMCFailTemplate
        sr_wjetsMC = sr_wjetsMCFail
        sr_wjets = sr_wjetsFail

    sr.addSample(sr_wjets)

    ###
    # top-antitop model
    ###

    sr_ttTemplate = template(background, "TT", "nominal", recoil, "sr", category, read_sumw2=True)
    sr_ttMC = rl.TemplateSample(
        "sr" + model_id + "_ttMC",
        rl.Sample.BACKGROUND,
        sr_ttTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('sr', sr_ttMC)
        sr_ttMC.setParamEffect(tt_norm, nMinor_norm)
    if options.btagVjets:
        addBtagSyst(background, recoil, "TT", "sr", sr_ttMC, category)
    if options.statUncs:
        addBBliteSyst(sr_ttMC, param, epsilon=1e-5) ### replace autoMCStats

    if category == "pass" and recoil<4:
        sigmascale={
            '2016': 1000,
            '2017': 1000,
            '2018': 1000
        }
        sr_ttObservable = rl.Observable("fjmass", sr_ttTemplate[1])
        sr_ttParameters = np.array(
            [
                rl.IndependentParameter(
                    "sr" + year + "_tt_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i,
                    0
                    #0,
                    #-20,
                    #5,
                )
                for i in range(sr_ttObservable.nbins)
            ]
        )
        sr_ttBinYields = sr_ttTemplate[0] * (1 + (0.1/np.maximum(1., np.sqrt(sr_ttTemplate[0]))))**sr_ttParameters
        #sr_ttBinYields = sr_ttTemplate[0] * (1 + (100./np.maximum(1., np.sqrt(sr_ttTemplate[0]))))**sr_ttParameters
        #sr_ttBinYields = sr_ttTemplate[0] * (1 + (sigmascale[year]/np.maximum(1., np.sqrt(sr_ttTemplate[0]))))**sr_ttParameters
        #sr_ttBinYields = np.array([rl.IndependentParameter('tmp', b, 1e-5, sr_ttTemplate[0].max()*2) for b in sr_ttTemplate[0]])

        sr_tt = rl.ParametericSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields
        )
        sr.addSample(sr_tt)
    else:
        sr.addSample(sr_ttMC)

    ###
    # Other MC-driven processes
    ###
    otherMCdrivenProcess(background, recoil, "sr", sr, category, ch_name, param)

    sr_signalTemplate = template(signal, s, "nominal", recoil, "sr", category)
    sr_signal = rl.TemplateSample(
        ch_name + "_" + str(s), rl.Sample.SIGNAL, sr_signalTemplate
    )
    #if options.rateSyst:
    #    addCommonRateSyst('sr', sr_signal)
    #if options.btagVjets:
    #    addBtagSyst(signal, recoil, str(s), "sr", sr_signal, category)
    #    addBBliteSyst(sr_signal, param, epsilon=1e-5) ### replace autoMCStats
    sr.addSample(sr_signal)

    ###
    # End of SR
    ###

    if category=="pass" and not (recoil<4):
        return model

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

    dataTemplate = template(data, "MET", "data", recoil, "wmcr", category)
    wmcr.setObservation(dataTemplate)
    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')

    ###
    # W(->lnu)+jets data-driven model
    ###

    wmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_wjetsMC = rl.TemplateSample(
        "wmcr" + model_id + "_wjetsMC",
        rl.Sample.BACKGROUND,
        wmcr_wjetsTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('wmcr', wmcr_wjetsMC)
        wmcr_wjetsMC.setParamEffect(wjets_norm, nVjets_norm)
    if options.btagVjets:
        addBtagSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)
        addVJetsSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)
    #addBBliteSyst(wmcr_wjetsMC, param, epsilon=1e-5) ### replace autoMCStats
    #wmcr_wjetsMC.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

    #### Transfer Factor
    wmcr_wjetsTransferFactor = wmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
    wmcr_wjets = rl.TransferFactorSample(ch_name + "_wjets", rl.Sample.BACKGROUND, wmcr_wjetsTransferFactor, sr_wjets)
    wmcr.addSample(wmcr_wjets)

    ###
    # top-antitop model
    ###

    wmcr_ttTemplate = template(background, "TT", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_ttMC = rl.TemplateSample(
        "wmcr" + model_id + "_ttMC",
        rl.Sample.BACKGROUND,
        wmcr_ttTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('wmcr', wmcr_ttMC)
        wmcr_ttMC.setParamEffect(tt_norm, nMinor_norm)
    if options.btagVjets:
        addBtagSyst(background, recoil, "TT", "wmcr", wmcr_ttMC, category)
    #addBBliteSyst(wmcr_ttMC, param, epsilon=1e-5) ### replace autoMCStats
    #wmcr_ttMC.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

    if category == "pass":
        #### Transfer Factor
        #wmcr_ttMC.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample
        wmcr_ttTransferFactor = wmcr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        wmcr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt
        )
        wmcr.addSample(wmcr_tt)
    else:
        #addBBliteSyst(wmcr_ttMC, param, epsilon=1e-5) ### replace autoMCStats
        wmcr.addSample(wmcr_ttMC)

    ###
    # Other MC-driven processes
    ###
    otherMCdrivenProcess(background, recoil, "wmcr", wmcr, category, ch_name, param)

    ###
    # End of single muon W control region
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
        dataTemplate = template(data, "EGamma", "data", recoil, "wecr", category)
    else:
        dataTemplate = template(data, "SingleElectron", "data", recoil, "wecr", category)

    wecr.setObservation(dataTemplate)
    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')

    ###
    # W(->lnu)+jets data-driven model
    ###

    wecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_wjetsMC = rl.TemplateSample(
        "wecr" + model_id + "_wjetsMC",
        rl.Sample.BACKGROUND,
        wecr_wjetsTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('wecr', wecr_wjetsMC)
        wecr_wjetsMC.setParamEffect(wjets_norm, nVjets_norm)
    if options.btagVjets:
        addBtagSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)
        addVJetsSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)
    #addBBliteSyst(wecr_wjetsMC, param, epsilon=1e-5) ### replace autoMCStats
    #wecr_wjetsMC.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

    #### Transfer Factor
    wecr_wjetsTransferFactor = wecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
    wecr_wjets = rl.TransferFactorSample(
        ch_name + "_wjets", rl.Sample.BACKGROUND, wecr_wjetsTransferFactor, sr_wjets
    )
    wecr.addSample(wecr_wjets)

    ###
    # top-antitop model
    ###

    wecr_ttTemplate = template(background, "TT", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_ttMC = rl.TemplateSample(
        "wecr" + model_id + "_ttMC",
        rl.Sample.BACKGROUND,
        wecr_ttTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('wecr', wecr_ttMC)
        wecr_ttMC.setParamEffect(tt_norm, nMinor_norm)
    if options.btagVjets:
        addBtagSyst(background, recoil, "TT", "wecr", wecr_ttMC, category)
    #addBBliteSyst(wecr_ttMC, param, epsilon=1e-5) ### replace autoMCStats
    #wecr_ttmc.automcstats(epsilon=1e-5) ### automcstats is used for transferfactorsample

    if category == "pass":
        #### Transfer Factor
        #wecr_ttmc.autoMCStats(epsilon=1e-5) ### autoMCStats is used for transferfactorsample
        wecr_ttTransferFactor = wecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        wecr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt
        )
        wecr.addSample(wecr_tt)
    else:
        #addBBliteSyst(wecr_ttMC, param, epsilon=1e-5) ### replace autoMCStats
        wecr.addSample(wecr_ttMC)

    ###
    # Other MC-driven processes
    ###
    otherMCdrivenProcess(background, recoil, "wecr", wecr, category, ch_name, param)

    ###
    # End of single electron W control region
    ###

    if not (category=="pass"):
        return model

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

    dataTemplate = template(data, "MET", "data", recoil, "tmcr", category)
    tmcr.setObservation(dataTemplate)
    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')

    ###
    # top-antitop model
    ###

    tmcr_ttTemplate = template(background, "TT", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_ttMC = rl.TemplateSample(
        "tmcr" + model_id + "_ttMC",
        rl.Sample.BACKGROUND,
        tmcr_ttTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('tmcr', tmcr_ttMC)
        tmcr_ttMC.setParamEffect(tt_norm, nMinor_norm)
    if options.btagVjets:
        addBtagSyst(background, recoil, "TT", "tmcr", tmcr_ttMC, category)
    #addBBliteSyst(tmcr_ttMC, param, epsilon=1e-5) ### replace autoMCStats
    #tmcr_ttMC.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

    #### Transfer Factor
    tmcr_ttTransferFactor = tmcr_ttMC.getExpectation() / sr_ttMC.getExpectation()
    tmcr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt
    )
    tmcr.addSample(tmcr_tt)

    ###
    # Other MC-driven processes
    ###

    tmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_wjets = rl.TemplateSample(
        ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tmcr_wjetsTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('tmcr', tmcr_wjets)
        tmcr_wjets.setParamEffect(wjetsMC_norm, nVjets_norm)
    if options.btagVjets:
        addBtagSyst(background, recoilbin, "W+jets", "tmcr", tmcr_wjets, category)
        addVJetsSyst(background, recoil, "W+jets", "tmcr", tmcr_wjets, category)
    #addBBliteSyst(tmcr_wjets, param, epsilon=1e-5) ### replace autoMCStats
    tmcr.addSample(tmcr_wjets)

    otherMCdrivenProcess(background, recoil, "tmcr", tmcr, category, ch_name, param)


    ###
    # End of single muon top control region
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
        dataTemplate = template(data, "EGamma", "data", recoil, "tecr", category)
    else:
        dataTemplate = template(data, "SingleElectron", "data", recoil, "tecr", category)

    tecr.setObservation(dataTemplate)
    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')

    ###
    # top-antitop model
    ###

    tecr_ttTemplate = template(background, "TT", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_ttMC = rl.TemplateSample(
        "tecr" + model_id + "_ttMC",
        rl.Sample.BACKGROUND,
        tecr_ttTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('tecr', tecr_ttMC)
        tecr_ttMC.setParamEffect(tt_norm, nMinor_norm)
    if options.btagVjets:
        addBtagSyst(background, recoil, "TT", "tecr", tecr_ttMC, category)
    #addBBliteSyst(tecr_ttMC, param, epsilon=1e-5) ### replace autoMCStats
    #tecr_ttMC.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

    #### Transfer Factor
    tecr_ttTransferFactor = tecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
    tecr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt
    )
    tecr.addSample(tecr_tt)

    ###
    # Other MC-driven processes
    ###

    tecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_wjets = rl.TemplateSample(
        ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tecr_wjetsTemplate
    )
    if options.rateSyst:
        addCommonRateSyst('tecr', tecr_wjets)
        tecr_wjets.setParamEffect(wjetsMC_norm, nVjets_norm)
    if options.btagVjets:
        addBtagSyst(background, recoilbin, "W+jets", "tecr", tecr_wjets, category)
        addVJetsSyst(background, recoil, "W+jets", "tecr", tecr_wjets, category)
    #addBBliteSyst(tecr_wjets, param, epsilon=1e-5) ### replace autoMCStats
    tecr.addSample(tecr_wjets)

    otherMCdrivenProcess(background, recoil, "tecr", tecr, category, ch_name, param)

    ###
    # End of single electron top control region
    ###

    return model


if __name__ == "__main__":
    if not os.path.exists("datacards"):
        os.mkdir("datacards")
    parser = OptionParser()
    parser.add_option("-y", "--year", help="year", dest="year", default="")
    parser.add_option("-b", "--sumbkg", help="replace data to sum of backgrounds", action="store_true", dest="sumbkg")
    parser.add_option("--rateSyst", help="Apply rate systematics to the template", action="store_true", dest="rateSyst")
    parser.add_option("--btagVjets", help="Apply b tagging and V+jets systematics to the template", action="store_true", dest="btagVjets")
    parser.add_option("--statUncs", help="Apply statistical uncertainties to the template", action="store_true", dest="statUncs")
    (options, args) = parser.parse_args()
    year = options.year
    recoil_binning = recoil_binning_dict[year]
    recoilbins = np.array(recoil_binning)
    nrecoil = len(recoilbins) - 1

    ###
    # Extract histograms from input file
    ###

    print("Grouping histograms")
    hists = load("hists/darkhiggs" + year + ".scaled")


    ###
    # Prepare fake data (MUST DO BEFORE REMAPPING!!!!)
    ###

    fake_data_hists = hists['bkg']
    fake_data_hists['template'] = fake_data_hists['template'].rebin('fjmass',hist.Bin('fjmass','Mass', mass_binning))
    fake_data_hists['template'] = fake_data_hists['template'].rebin('recoil',hist.Bin('recoil','Recoil',recoil_binning))

    ###
    # Remapping histograms 
    ###

    hists = remap_histograms(hists)
    data_hists = hists["data"]
    bkg_hists = hists["bkg"]
    signal_hists = hists["sig"]

    ###
    # Preparing histograms for fit
    ##

    data = {}
    fake_data = {}

    #### Use real data in CRs + SR fail while fake data in SR pass
    if options.sumbkg:
        fake_data['sr'] = fake_data_hists["template"].integrate("region", 'sr').sum("gentype")
    for r in data_hists["template"].identifiers("region"):
        data[str(r)] = data_hists["template"].integrate("region", r).sum("gentype")

    background = {}
    for r in bkg_hists["template"].identifiers("region"):
        background[str(r)] = bkg_hists["template"].integrate("region", r).sum("gentype")

    signal = {}
    for r in signal_hists["template"].identifiers("region"):
        signal[str(r)] = signal_hists["template"].integrate("region", r).sum("gentype")

    ### 
    # Fill missing signal mass points for 2018
    ###

    if year == '2018':
        print("Grouping missing signal histograms in 2018")
        tmphists = load("hists/signal" + year + ".scaled")
        tmphists = remap_histograms(tmphists)
        signal_tmphists = tmphists['sig']

        for r in signal_hists["template"].identifiers("region"):
            signal[str(r)] += signal_tmphists["template"].integrate("region", r).sum("gentype")


    ###
    ###
    # Setting up systematics
    ###
    ###
    lumi = rl.NuisanceParameter("lumi" + year, "lnN")
    zjets_norm = rl.NuisanceParameter("zjets_norm", "lnN")
    wjets_norm = rl.NuisanceParameter("wjets_norm", "lnN")
    tt_norm = rl.NuisanceParameter("tt_norm", "lnN")
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
    ew1 = rl.NuisanceParameter("ew1", "lnN")
    ew2W = rl.NuisanceParameter("ew2W", "lnN")
    ew2Z = rl.NuisanceParameter("ew2Z", "lnN")
    ew3W = rl.NuisanceParameter("ew3W", "lnN")
    ew3Z = rl.NuisanceParameter("ew3Z", "lnN")
    mix = rl.NuisanceParameter("mix", "lnN")
    qcd1 = rl.NuisanceParameter("qcd1", "lnN")
    qcd2 = rl.NuisanceParameter("qcd2", "lnN")
    qcd3 = rl.NuisanceParameter("qcd3", "lnN")
    whf_fraction = rl.NuisanceParameter("whf_fraction", "lnN")
    zhf_fraction = rl.NuisanceParameter("zhf_fraction", "lnN")

    ###
    # Set lnN or shape numbers
    ###

    nlumi = 1.027
    ntrig_met = 1.02
    ntrig_e = 1.01
    nveto_tau = 1.03
    njec = 1.05
    nlepton = 1.02 ## id_mu, iso_mu, id_e, reco_e
    nVjets_norm = 1.4 ## wjetsMC_norm, wjets_norm, zjetsMC_norm, zjets_norm, whf_fraction, zhf_fraction
    nMinor_norm = 1.2 ## tt_norm, ttMC_norm, st_norm, vv_norm, hbb_norm
    nqcd_norm = 2.0 ## qcdsig_norm, qcde_norm, qcdmu_norm

    ###
    # Preparing Rhalphabet
    ###

    msdbins = np.array(mass_binning)
    msd = rl.Observable('fjmass', msdbins)
    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(recoilbins[:-1] + 0.3 * np.diff(recoilbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    recoilscaled = (ptpts - 250.) / (3000. - 250.)
    msdscaled = (msdpts - 40.) / (300.0 - 40.)

    tf_dataResidualW = rl.BernsteinPoly("tf_dataResidualW"+year, (1, 1), ['recoil', 'fjmass'], limits=(-100, 100))
    tf_dataResidualW_params = tf_dataResidualW(recoilscaled, msdscaled)
    tf_dataResidualZ = rl.BernsteinPoly("tf_dataResidualZ"+year, (1, 1), ['recoil', 'fjmass'], limits=(-100, 100))
    tf_dataResidualZ_params = tf_dataResidualZ(recoilscaled, msdscaled)

    model_dict = {}
    for recoilbin in range(nrecoil):


        sr_zjetsMCFailTemplate = template(background, "Z+jets", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zjetsMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zjetsMC",
            rl.Sample.BACKGROUND,
            sr_zjetsMCFailTemplate
        )
        if options.rateSyst:
            addCommonRateSyst('sr', sr_zjetsMCFail)
            sr_zjetsMCFail.setParamEffect(zjets_norm, nVjets_norm)
        if options.btagVjets:
            addBtagSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCFail, "fail")
            addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCFail, "fail")
        #addBBliteSyst(sr_zjetsMCFail, param, epsilon=1e-5) ### replace autoMCStats
        #sr_zjetsMCFail.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

        sr_zhfMCFailTemplate = template(background, "Z+HF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zhfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zhfMC",
            rl.Sample.BACKGROUND,
            sr_zhfMCFailTemplate
        )
        if options.rateSyst:
            sr_zhfMCFail.setParamEffect(zhf_fraction, nVjets_norm, )

        sr_zlfMCFailTemplate = template(background, "Z+LF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zlfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zlfMC",
            rl.Sample.BACKGROUND,
            sr_zlfMCFailTemplate
        )
        if options.rateSyst:
            sr_zlfMCFail.setParamEffect(zhf_fraction, 0.95)

        sr_zjetsObservable = rl.Observable("fjmass", sr_zjetsMCFailTemplate[1])
        sr_zjetsParameters = np.array(
            [
                rl.IndependentParameter(
                    "sr" + year + "_zjets_fail_recoil"+str(recoilbin)+"_mass%d" % i,
                    0
                )
                for i in range(sr_zjetsObservable.nbins)
            ]
        )
        sr_zjetsBinYields = sr_zjetsMCFailTemplate[0] * (1 + (10./np.maximum(1., np.sqrt(sr_zjetsMCFailTemplate[0]))))**sr_zjetsParameters
        #sr_zjetsBinYields = np.array([rl.IndependentParameter('tmp', b, 1e-5, sr_zjetsMCFailTemplate[0].max()*2) for b in sr_zjetsMCFailTemplate[0]])

        sr_zjetsFail = rl.ParametericSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin)+ "_zjets",
            rl.Sample.BACKGROUND,
            sr_zjetsObservable,
            sr_zjetsBinYields
        )

        sr_wjetsMCFailTemplate = template(background, "W+jets", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_wjetsMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_wjetsMC",
            rl.Sample.BACKGROUND,
            sr_wjetsMCFailTemplate
        )
        if options.rateSyst:
            addCommonRateSyst('sr', sr_wjetsMCFail)
            sr_wjetsMCFail.setParamEffect(wjets_norm, nVjets_norm)
        if options.btagVjets:
            addBtagSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCFail, "fail")
            addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCFail, "fail")
        #addBBliteSyst(sr_wjetsMCFail, param, epsilon=1e-5) ### replace autoMCStats
        #sr_wjetsMCFail.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

        sr_whfMCFailTemplate = template(background, "W+HF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_whfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_whfMC",
            rl.Sample.BACKGROUND,
            sr_whfMCFailTemplate
        )
        if options.rateSyst:
            sr_whfMCFail.setParamEffect(whf_fraction, nVjets_norm, )

        sr_wlfMCFailTemplate = template(background, "W+LF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_wlfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_wlfMC",
            rl.Sample.BACKGROUND,
            sr_wlfMCFailTemplate
        )
        if options.rateSyst:
            sr_wlfMCFail.setParamEffect(whf_fraction, 0.9)

        sr_wjetsFailTransferFactor = sr_wjetsMCFail.getExpectation() / sr_zjetsMCFail.getExpectation()
        sr_wjetsFail = rl.TransferFactorSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin)+ "_wjets",
            rl.Sample.BACKGROUND,
            sr_wjetsFailTransferFactor,
            sr_zjetsFail
        )

        sr_zjetsMCPassTemplate = template(background, "Z+jets", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_zjetsMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_zjetsMC",
            rl.Sample.BACKGROUND,
            sr_zjetsMCPassTemplate
        )
        if options.rateSyst:
            addCommonRateSyst('sr', sr_zjetsMCPass)
            sr_zjetsMCPass.setParamEffect(zjets_norm, nVjets_norm)
        if options.btagVjets:
            addBtagSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCPass, "pass")
            addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCPass, "pass")
        #addBBliteSyst(sr_zjetsMCPass, param, epsilon=1e-5) ### replace autoMCStats
        #sr_zjetsMCPass.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

        sr_zhfMCPassTemplate = template(background, "Z+HF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_zhfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_zhfMC",
            rl.Sample.BACKGROUND,
            sr_zhfMCPassTemplate
        )
        if options.rateSyst:
            sr_zhfMCPass.setParamEffect(zhf_fraction, nVjets_norm, )

        sr_zlfMCPassTemplate = template(background, "Z+LF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_zlfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_zlfMC",
            rl.Sample.BACKGROUND,
            sr_zlfMCPassTemplate
        )
        if options.rateSyst:
            sr_zlfMCPass.setParamEffect(zhf_fraction, 0.95)

        tf_paramsZdeco = (sr_zlfMCPass.getExpectation()+sr_zhfMCPass.getExpectation()) / (sr_zlfMCFail.getExpectation()+sr_zhfMCFail.getExpectation())
        tf_paramsZ = tf_paramsZdeco #* tf_dataResidualZ_params[recoilbin, :]

        sr_zjetsPass = rl.TransferFactorSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin)+ "_zjets",
            rl.Sample.BACKGROUND,
            tf_paramsZ,
            sr_zjetsFail
        )

        sr_wjetsMCPassTemplate = template(background, "W+jets", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_wjetsMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_wjetsMC",
            rl.Sample.BACKGROUND,
            sr_wjetsMCPassTemplate
        )
        if options.rateSyst:
            addCommonRateSyst('sr', sr_wjetsMCPass)
            sr_wjetsMCPass.setParamEffect(wjets_norm, nVjets_norm)
        if options.btagVjets:
            addBtagSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCPass, "pass")
            addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCPass, "pass")
        #addBBliteSyst(sr_wjetsMCPass, param, epsilon=1e-5) ### replace autoMCStats
        #sr_wjetsMCPass.autoMCStats(epsilon=1e-5) ### autoMCStats is used for TransferFactorSample

        sr_whfMCPassTemplate = template(background, "W+HF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_whfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_whfMC",
            rl.Sample.BACKGROUND,
            sr_whfMCPassTemplate
        )
        if options.rateSyst:
            sr_whfMCPass.setParamEffect(whf_fraction, nVjets_norm)

        sr_wlfMCPassTemplate = template(background, "W+LF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_wlfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_wlfMC",
            rl.Sample.BACKGROUND,
            sr_wlfMCPassTemplate
        )
        if options.rateSyst:
            sr_wlfMCPass.setParamEffect(whf_fraction, 0.9)

        tf_paramsWdeco = (sr_wlfMCPass.getExpectation()+sr_whfMCPass.getExpectation()) / (sr_wlfMCFail.getExpectation()+sr_whfMCFail.getExpectation())
        tf_paramsW = tf_paramsWdeco #* tf_dataResidualW_params[recoilbin, :]

        sr_wjetsPass = rl.TransferFactorSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin)+ "_wjets",
            rl.Sample.BACKGROUND,
            tf_paramsW,
            sr_wjetsFail
        )

        for s in signal["sr"].identifiers("process"):
            print("Signal is:", str(s))
            for category in ["pass", "fail"]:
                qcdpho_norm = rl.NuisanceParameter("qcdpho_norm" + year + category, "lnN")
                qcde_norm = rl.NuisanceParameter("qcde_norm" + year + category, "lnN")
                qcdmu_norm = rl.NuisanceParameter("qcdmu_norm" + year + category, "lnN")
                qcdsig_norm = rl.NuisanceParameter("qcdsig_norm" + year + category, "lnN")
                st_norm = rl.NuisanceParameter("st_norm" + year + category, "lnN")
                ttMC_norm = rl.NuisanceParameter("tt_norm" + year + category, "lnN")
                vv_norm = rl.NuisanceParameter("vv_norm" + year + category, "lnN")
                hbb_norm = rl.NuisanceParameter("hbb_norm" + year + category, "lnN")
                wjetsMC_norm = rl.NuisanceParameter("wjets_norm" + year + category, "lnN")
                zjetsMC_norm = rl.NuisanceParameter("zjets_norm" + year + category, "lnN")

                with open(
                        "data/"
                        + str(s).replace('_','')
                        + "-"
                        + year
                        + "-"
                        + category
                        + "-recoil"
                        + str(recoilbin)
                        + ".model",
                        "wb",
                ) as fout:
                    pickle.dump(model(year, recoilbin, category, s), fout, protocol=2)
