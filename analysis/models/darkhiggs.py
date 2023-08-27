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
import uncertainties as unc  
import uncertainties.unumpy as unumpy 
from coffea import hist, processor
from coffea.util import load, save
from scipy import stats
import ROOT

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

mass_binning = [40., 50., 60., 70., 80., 90., 100., 120., 150., 180., 240., 300.,]
#mass_binning = [40., 50., 60., 70., 80., 90., 100., 120., 160., 200., 300.,] 
recoil_binning = [250., 310., 370., 470., 590., 3000.]
category_map = {"pass": 1, "fail": 0}

def template(dictionary, process, systematic, recoil, region, category, mass, min_value=1e-5, read_sumw2=False):
    histogram = dictionary[region].integrate("process", process)
    nominal, sumw2 = histogram.integrate("systematic", "nominal").values(sumw2=True)[()]
    nominal=nominal[recoil, :, category_map[category]]
    sumw2=sumw2[recoil, :, category_map[category]]
    zerobins = nominal <= 0.
    output = nominal
    if "data" not in systematic:
        output[zerobins] = min_value
        sumw2[zerobins] = 0.
    if "nominal" not in systematic and "data" not in systematic:
        output = histogram.integrate("systematic", systematic).values()[()][recoil, :, category_map[category]]
        output[zerobins] = 1.
        output[~zerobins] /= nominal[~zerobins]
        output[~zerobins] = np.maximum(output[~zerobins], 1e-5)
        output[np.isnan(output)] = 1.
    binning = (
        dictionary[region]
        .integrate("process", process)
        .integrate("systematic", systematic)
        .axis("fjmass")
        .edges()
    )
    if read_sumw2:
        return (output, binning, 'fjmass'+mass, sumw2)
    return (output, binning, 'fjmass'+mass)

def remap_histograms(hists):
    data_hists = {}
    bkg_hists = {}
    signal_hists = {}
    fakedata_map = OrderedDict()
  
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
        
    fakedata_list = []
    for bkg in hists['bkg']['template'].identifiers('process'):
        fakedata_list.append(str(bkg))
    fakedata_map['FakeData'] = (fakedata_list,)
    
    for key in hists["data"].keys():
        bkg_hists[key] = hists["bkg"][key].group(cats, process, bkg_map)
        signal_hists[key] = hists["sig"][key].group(cats, process, sig_map)
        data_hists[key] = hists["data"][key].group(cats, process, data_map)
        data_hists[key] += hists["bkg"][key].group(cats, process, fakedata_map)
    
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

'''
def get_mergedMC_stat_variations(dictionary, recoil, region, category, bkg_list):
    MCbkg = {}
    MCbkg_map = OrderedDict()
    process = hist.Cat("process", "Process", sorting="placement")
    cats = ("process",)
    for bkg in bkg_list:
        MCbkg_map[bkg] = (bkg,)
    MCbkg=dictionary[region].group(cats, process, MCbkg_map)
    merged_obj = MCbkg.integrate("process")
    merged_central, merged_error2 = merged_obj.integrate("systematic", "nominal").values(sumw2=True)[()]
    merged_central=merged_central[recoil, :, category_map[category]]
    merged_error2=merged_error2[recoil, :, category_map[category]]

    return merged_central, merged_error2
'''
def get_mergedMC_stat_variations(dictionary, recoil, region, category, mass, bkg_list):
    templ=template(dictionary, bkg_list[0], "nominal", recoil, region, category, mass, read_sumw2=True)
    merged_central=np.zeros_like(templ[0])
    merged_error2=np.zeros_like(templ[3])
    for bkg in bkg_list:
        templ=template(dictionary, bkg, "nominal", recoil, region, category, mass, read_sumw2=True)
        for i in range(len(templ[0])):
            if templ[0][i] <= 1e-5 or templ[3][i] <= 0.:
                continue
            merged_central[i] += templ[0][i]
            merged_error2[i]  += templ[3][i]
    return merged_central, merged_error2

def addBBliteSyst(templ, param, merged_central, merged_error2, epsilon=1e-5, threshold=0.01):
    for i in range(templ.observable.nbins):
        if merged_central[i] <= 0. or merged_error2[i] <= 0.:
            continue
        if isinstance(templ, rl.TemplateSample) and templ._nominal[i] <= 1e-5:
            continue
        effect_up = np.ones_like(templ._nominal)
        effect_down = np.ones_like(templ._nominal)
        if (np.sqrt(merged_error2[i]) / (merged_central[i] + 1e-12)) < threshold:
            continue
        effect_up[i] = 1.0 + np.sqrt(merged_error2[i])/merged_central[i]
        effect_down[i] = max(epsilon, 1.0 - np.sqrt(merged_error2[i])/merged_central[i])
        #print(templ._name,templ._nominal[i],effect_up[i],effect_down[i])
        templ.setParamEffect(param[i], effect_up, effect_down)

def addMCStatsTFSyst(templ, num, den, epsilon=1e-5, threshold=0.01):
    num=unumpy.uarray(( num[0], np.sqrt(num[3]) ))  
    den=unumpy.uarray(( den[0], np.sqrt(den[3]) ))  
    tf=num/den
    #err=unumpy.std_devs(tf)  
    for i in range(templ.observable.nbins):
        effect_up = np.ones_like(templ._nominal)
        effect_down = np.ones_like(templ._nominal)
        if unumpy.std_devs(tf)[i]/unumpy.nominal_values(tf)[i] < threshold:
            continue
        effect_up[i] = 1.0 + unumpy.std_devs(tf)[i]/unumpy.nominal_values(tf)[i]
        effect_down[i] = max(epsilon, 1.0 - unumpy.std_devs(tf)[i]/unumpy.nominal_values(tf)[i])
        print(templ._name, i, effect_up[i], effect_down[i], tf[i])
        param = rl.NuisanceParameter(templ._name + '_mcstat_bin%i' % i, combinePrior='shape')
        templ.setParamEffect(param, effect_up, effect_down)

def addBtagSyst(dictionary, recoil, process, region, templ, category, mass):
    btagUp = template(dictionary, process, "btagUp", recoil, region, category, mass)[0]
    btagDown = template(dictionary, process, "btagDown", recoil, region, category, mass)[0]
    templ.setParamEffect(btag, btagUp, btagDown)

def addDoubleBtagSyst(dictionary, recoil, process, region, templ, category, mass):
    for syst in dictionary[region].identifiers("systematic"):
        if 'doublebtag' not in str(syst): continue
        if 'Down' in str(syst): continue
        doublebtagUp = template(dictionary, process, str(syst), recoil, region, category, mass)[0]
        doublebtagDown = template(dictionary, process, str(syst).replace('Up','Down'), recoil, region, category, mass)[0]
        templ.setParamEffect(doublebtag[str(syst).replace('Up','')], doublebtagUp, doublebtagDown)

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
    #addSyst(dictionary, recoil, process, region, templ, category, muF, "muF")
    #addSyst(dictionary, recoil, process, region, templ, category, muR, "muR")

def addLumiSyst(templ, year):
    vlumi={
        '2016': 1.01,
        '2017': 1.02,
        '2018': 1.015,
    }
    vlumi_corr={
        '2016': 1.006,
        '2017': 1.009,
        '2018': 1.02,
    }
    vlumi_1718={
        '2017': 1.006,
        '2018': 1.002,
    }
    templ.setParamEffect(lumi, vlumi[year])
    templ.setParamEffect(lumi_corr, vlumi_corr[year])
    if '2016' not in year: templ.setParamEffect(lumi_1718, vlumi_1718[year])

def addPileupSyst(templ):
    templ.setParamEffect(pu, 1.01)

def addMETSyst(templ):
    templ.setParamEffect(met, 1.05)

def addJESSyst(templ):
    templ.setParamEffect(jes, 1.04)

def addPrefiringSyst(templ, year):
    if '2018' not in year: templ.setParamEffect(prefiring, 1.01)

def addSingleTopNormSyst(templ):
    templ.setParamEffect(st_norm, 1.1)

def addTTbarNormSyst(templ):
    templ.setParamEffect(ttMC_norm, 1.1)

def addHbbNormSyst(templ):
    templ.setParamEffect(hbb_norm, 1.2)

def addDibosonNormSyst(templ):
    templ.setParamEffect(vv_norm, 1.2)

def addDrellYanNormSyst(templ):
    templ.setParamEffect(zjetsMC_norm, 1.2)

def addWjetsNormSyst(templ):
    templ.setParamEffect(wjetsMC_norm, 1.4)

def addMETTrigSyst(templ, year):
    if '2016' in year:
        templ.setParamEffect(trig_met, 1.02)
    else:
        templ.setParamEffect(trig_met, 1.01)

def addEleIDSyst(templ, year):
    if '2016' in year:
        templ.setParamEffect(id_e, 1.02)
    else:
        templ.setParamEffect(id_e, 1.03)


    
def model(year, mass, recoil, category):

    model_id = year + category + "mass" + mass+ "recoil" + str(recoil)
    #print(model_id)
    model = rl.Model(model_id)

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

    if category == 'pass' and options.fakedata:
        dataTemplate = template(data, "FakeData", "data", recoil, "sr", category, mass)
    else:
        dataTemplate = template(data, "MET", "data", recoil, "sr", category, mass)
    sr.setObservation(dataTemplate)

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')

    MCbkgList = ["ST", "DY+jets", "VV", "Hbb", "QCD"]
    if isttMC: MCbkgList.append("TT")
    if iswjetsMC: MCbkgList.append("W+jets")
    sr_central, sr_error2 = get_mergedMC_stat_variations(background, recoil, "sr", category, mass, MCbkgList)

    ###
    # Z(->nunu)+jets data-driven model
    ###

    if category == "pass":
        #param = param_pass
        #sr_central = sr_central_pass
        #sr_error2 = sr_error2_pass
        sr_zjets = sr_zjetsPass
        #sr_zjetsBinYields = sr_zjetsPassBinYields
    else:
        #param = param_fail
        #sr_central = sr_central_fail
        #sr_error2 = sr_error2_fail
        sr_zjets = sr_zjetsFail
        #sr_zjetsBinYields = sr_zjetsFailBinYields
    sr.addSample(sr_zjets)

    ###
    # W(->lnu)+jets data-driven model
    ###

    if not iswjetsMC:
        if category == "pass":
            sr_wjets = sr_wjetsPass
            sr_wjetsMC = sr_wjetsMCPass
            sr_wjetsTemplate = sr_wjetsMCPassTemplate
            #sr_wjetsBinYields = sr_wjetsPassBinYields
        else:
            sr_wjets = sr_wjetsFail
            sr_wjetsMC = sr_wjetsMCFail
            sr_wjetsTemplate = sr_wjetsMCFailTemplate
            #sr_wjetsBinYields = sr_wjetsFailBinYields
        sr.addSample(sr_wjets)
        
    ###
    # top-antitop data-driven model
    ###

    if not isttMC:
        sr_ttTemplate = template(background, "TT", "nominal", recoil, "sr", category, mass, min_value=1., read_sumw2=True)
        sr_ttMC = rl.TemplateSample("sr" + model_id + "_ttMC",rl.Sample.BACKGROUND,sr_ttTemplate)
        addMETTrigSyst(sr_ttMC, year)
        addBtagSyst(background, recoil, "TT", "sr", sr_ttMC, category, mass)
        
        sr_ttObservable = rl.Observable("fjmass"+mass, sr_ttTemplate[1])
        sr_ttBinYields = np.array([rl.IndependentParameter(ch_name + "_tt_mu"+str(b), sr_ttTemplate[0][b], 1e-5, sr_ttTemplate[0].max()*2) for b in range(len(sr_ttTemplate[0]))])
        sr_tt = rl.ParametericSample(ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields)
        #addBBliteSyst(sr_tt, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
        sr.addSample(sr_tt)
    
    ###
    # Other MC-driven processes
    ###
    
    if iswjetsMC: 
        sr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "sr", category, mass, read_sumw2=True)
        sr_wjets = rl.TemplateSample( "sr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, sr_wjetsTemplate)
        addLumiSyst(sr_wjets, year)
        addPileupSyst(sr_wjets)
        addPrefiringSyst(sr_wjets, year)
        addMETSyst(sr_wjets)
        addJESSyst(sr_wjets)
        addMETTrigSyst(sr_wjets, year)
        sr_wjets.setParamEffect(veto_tau, nveto_tau)
        addWjetsNormSyst(sr_wjets)
        addBBliteSyst(sr_wjets, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
        addBtagSyst(background, recoil, "W+jets", "sr", sr_wjets, category, mass)
        addVJetsSyst(background, recoil, "W+jets", "sr", sr_wjets, category)
        sr.addSample(sr_wjets)

    if isttMC: 
        sr_ttTemplate = template(background, "TT", "nominal", recoil, "sr", category, mass, read_sumw2=True)
        sr_tt = rl.TemplateSample("sr" + model_id + "_ttMC",rl.Sample.BACKGROUND,sr_ttTemplate)
        addLumiSyst(sr_tt, year)
        addPileupSyst(sr_tt)
        addPrefiringSyst(sr_tt, year)
        addMETSyst(sr_tt)
        addJESSyst(sr_tt)
        addMETTrigSyst(sr_tt, year)
        sr_tt.setParamEffect(veto_tau, nveto_tau)
        addTTbarNormSyst(sr_tt)
        addBtagSyst(background, recoil, "TT", "sr", sr_tt, category, mass)
        addBBliteSyst(sr_tt, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
        sr.addSample(sr_tt)
    
    sr_stTemplate = template(background, "ST", "nominal", recoil, "sr", category, mass, read_sumw2=True)
    sr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, sr_stTemplate)
    addLumiSyst(sr_st, year)
    addPileupSyst(sr_st)
    addPrefiringSyst(sr_st, year)
    addMETSyst(sr_st)
    addJESSyst(sr_st)
    addMETTrigSyst(sr_st, year)
    sr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(sr_st)
    addBBliteSyst(sr_st, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoil, "ST", "sr", sr_st, category, mass)
    sr.addSample(sr_st)

    sr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "sr", category, mass, read_sumw2=True)
    sr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, sr_dyjetsTemplate)
    addLumiSyst(sr_dyjets, year)
    addPileupSyst(sr_dyjets)
    addPrefiringSyst(sr_dyjets, year)
    addMETSyst(sr_dyjets)
    addJESSyst(sr_dyjets)
    addMETTrigSyst(sr_dyjets, year)
    sr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(sr_dyjets)
    addBBliteSyst(sr_dyjets, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoil, "DY+jets", "sr", sr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "sr", sr_dyjets, category)
    sr.addSample(sr_dyjets)

    sr_vvTemplate = template(background, "VV", "nominal", recoil, "sr", category, mass, read_sumw2=True)
    sr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, sr_vvTemplate)
    addLumiSyst(sr_vv, year)
    addPileupSyst(sr_vv)
    addPrefiringSyst(sr_vv, year)
    addMETSyst(sr_vv)
    addJESSyst(sr_vv)
    addMETTrigSyst(sr_vv, year)
    sr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(sr_vv)
    addBBliteSyst(sr_vv, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoil, "VV", "sr", sr_vv, category, mass)
    sr.addSample(sr_vv)

    sr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "sr", category, mass, read_sumw2=True)
    sr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, sr_hbbTemplate)
    addLumiSyst(sr_hbb, year)
    addPileupSyst(sr_hbb)
    addPrefiringSyst(sr_hbb, year)
    addMETSyst(sr_hbb)
    addJESSyst(sr_hbb)
    addMETTrigSyst(sr_hbb, year)
    sr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(sr_hbb)
    addBBliteSyst(sr_hbb, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoil, "Hbb", "sr", sr_hbb, category, mass)
    sr.addSample(sr_hbb)

    sr_qcdTemplate = template(background, "QCD", "nominal", recoil, "sr", category, mass, read_sumw2=True)
    sr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, sr_qcdTemplate)
    addLumiSyst(sr_qcd, year)
    addPileupSyst(sr_qcd)
    addPrefiringSyst(sr_qcd, year)
    addMETSyst(sr_qcd)
    addJESSyst(sr_qcd)
    addMETTrigSyst(sr_qcd, year)
    sr_qcd.setParamEffect(veto_tau, nveto_tau)
    sr_qcd.setParamEffect(qcdsig_norm, nqcd_norm)
    addBBliteSyst(sr_qcd, param, sr_central, sr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoil, "QCD", "sr", sr_qcd, category, mass)
    sr.addSample(sr_qcd)

    if category=="pass": 
        for s in signal["sr"].identifiers("process"):
            sr_signalTemplate = template(signal, s, "nominal", recoil, "sr", category, mass, read_sumw2=True)
            sr_signal = rl.TemplateSample(ch_name + "_" + str(s), rl.Sample.SIGNAL, sr_signalTemplate)
            addLumiSyst(sr_signal, year)
            addPileupSyst(sr_signal)
            addPrefiringSyst(sr_signal, year)
            addMETSyst(sr_signal)
            addJESSyst(sr_signal)
            addMETTrigSyst(sr_signal, year)
            sr_signal.setParamEffect(veto_tau, nveto_tau)
            #addBBliteSyst(sr_signal, param, sr_central, sr_error2, epsilon=1e-5)
            addBtagSyst(signal, recoil, str(s), "sr", sr_signal, category, mass)
            addDoubleBtagSyst(signal, recoil, str(s), "sr", sr_signal, category, mass)
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

    dataTemplate = template(data, "MET", "data", recoil, "wmcr", category, mass)
    wmcr.setObservation(dataTemplate)

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')
    #wmcr_central, wmcr_error2 = get_mergedMC_stat_variations(background, recoil, "wmcr", category, mass, ["W+jets", "TT", "ST", "DY+jets", "VV", "Hbb", "QCD"])

    MCbkgList = ["ST", "DY+jets", "VV", "Hbb", "QCD"]
    if isttMC: MCbkgList.append("TT")
    if iswjetsMC: MCbkgList.append("W+jets")
    wmcr_central, wmcr_error2 = get_mergedMC_stat_variations(background, recoil, "wmcr", category, mass, MCbkgList)
   
    ###
    # W(->lnu)+jets data-driven model
    ###

    if not iswjetsMC:
        wmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wmcr", category, mass, min_value=1., read_sumw2=True)
        wmcr_wjetsMC = rl.TemplateSample("wmcr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wmcr_wjetsTemplate)
        addMETTrigSyst(wmcr_wjetsMC, year)
        wmcr_wjetsMC.setParamEffect(id_mu, nlepton)
        wmcr_wjetsMC.setParamEffect(iso_mu, nlepton)
        addVJetsSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)
        addMCStatsTFSyst(wmcr_wjetsMC, wmcr_wjetsTemplate, sr_wjetsTemplate, epsilon=1e-5)
        
        #### Transfer Factor
        wmcr_wjetsTransferFactor = wmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
        #wmcr_wjets = rl.ParametericSample(ch_name + "_wjets", rl.Sample.BACKGROUND, sr_zjetsObservable, sr_wjetsBinYields*wmcr_wjetsTransferFactor)
        wmcr_wjets = rl.TransferFactorSample(ch_name + "_wjets", rl.Sample.BACKGROUND, wmcr_wjetsTransferFactor, sr_wjets)
        #addBBliteSyst(wmcr_wjets, param, wmcr_central, wmcr_error2, epsilon=1e-5)
        wmcr.addSample(wmcr_wjets)

    ###
    # top-antitop data-driven model
    ###

    if not isttMC:
        wmcr_ttTemplate = template(background, "TT", "nominal", recoil, "wmcr", category, mass, min_value=1., read_sumw2=True)
        wmcr_ttMC = rl.TemplateSample( "wmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wmcr_ttTemplate)
        addMETTrigSyst(wmcr_ttMC, year)
        wmcr_ttMC.setParamEffect(id_mu, nlepton)
        wmcr_ttMC.setParamEffect(iso_mu, nlepton)
        addBtagSyst(background, recoil, "TT", "wmcr", wmcr_ttMC, category, mass)
        addMCStatsTFSyst(wmcr_ttMC, wmcr_ttTemplate, sr_ttTemplate, epsilon=1e-5)
        
        #### Transfer Factor
        wmcr_ttTransferFactor = wmcr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        #wmcr_tt = rl.ParametericSample(ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields*wmcr_ttTransferFactor)
        wmcr_tt = rl.TransferFactorSample(ch_name + "_tt", rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt)
        #addBBliteSyst(wmcr_tt, param, wmcr_central, wmcr_error2, epsilon=1e-5)
        wmcr.addSample(wmcr_tt)
    
    ###
    # Other MC-driven processes
    ###

    if iswjetsMC:
        wmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
        wmcr_wjets = rl.TemplateSample("wmcr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wmcr_wjetsTemplate)
        addLumiSyst(wmcr_wjets, year)
        addPileupSyst(wmcr_wjets)
        addPrefiringSyst(wmcr_wjets, year)
        addJESSyst(wmcr_wjets)
        addMETTrigSyst(wmcr_wjets, year)
        wmcr_wjets.setParamEffect(veto_tau, nveto_tau)
        addWjetsNormSyst(wmcr_wjets)
        wmcr_wjets.setParamEffect(id_mu, nlepton)
        wmcr_wjets.setParamEffect(iso_mu, nlepton)
        addBBliteSyst(wmcr_wjets, param, wmcr_central, wmcr_error2, epsilon=1e-5)
        addBtagSyst(background, recoil, "W+jets", "wmcr", wmcr_wjets, category, mass)
        addVJetsSyst(background, recoil, "W+jets", "wmcr", wmcr_wjets, category)
        wmcr.addSample(wmcr_wjets)

    if isttMC: 
        wmcr_ttTemplate = template(background, "TT", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
        wmcr_tt = rl.TemplateSample( "wmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wmcr_ttTemplate)
        addLumiSyst(wmcr_tt, year)
        addPileupSyst(wmcr_tt)
        addPrefiringSyst(wmcr_tt, year)
        addJESSyst(wmcr_tt)
        addMETTrigSyst(wmcr_tt, year)
        wmcr_tt.setParamEffect(veto_tau, nveto_tau)
        wmcr_tt.setParamEffect(id_mu, nlepton)
        wmcr_tt.setParamEffect(iso_mu, nlepton)
        addTTbarNormSyst(wmcr_tt)
        addBtagSyst(background, recoil, "TT", "wmcr", wmcr_tt, category, mass)
        addBBliteSyst(wmcr_tt, param, wmcr_central, wmcr_error2, epsilon=1e-5) ### replace autoMCStats
        wmcr.addSample(wmcr_tt)
                    
    wmcr_stTemplate = template(background, "ST", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, wmcr_stTemplate)
    addLumiSyst(wmcr_st, year)
    addPileupSyst(wmcr_st)
    addPrefiringSyst(wmcr_st, year)
    addJESSyst(wmcr_st)
    addMETTrigSyst(wmcr_st, year)
    wmcr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(wmcr_st)
    wmcr_st.setParamEffect(id_mu, nlepton)
    wmcr_st.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(wmcr_st, param, wmcr_central, wmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "ST", "wmcr", wmcr_st, category, mass)
    wmcr.addSample(wmcr_st)

    wmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wmcr_dyjetsTemplate)
    addLumiSyst(wmcr_dyjets, year)
    addPileupSyst(wmcr_dyjets)
    addPrefiringSyst(wmcr_dyjets, year)
    addJESSyst(wmcr_dyjets)
    addMETTrigSyst(wmcr_dyjets, year)
    wmcr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(wmcr_dyjets)
    wmcr_dyjets.setParamEffect(id_mu, nlepton)
    wmcr_dyjets.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(wmcr_dyjets, param, wmcr_central, wmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "DY+jets", "wmcr", wmcr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "wmcr", wmcr_dyjets, category)
    wmcr.addSample(wmcr_dyjets)

    wmcr_vvTemplate = template(background, "VV", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, wmcr_vvTemplate)
    addLumiSyst(wmcr_vv, year)
    addPileupSyst(wmcr_vv)
    addPrefiringSyst(wmcr_vv, year)
    addJESSyst(wmcr_vv)
    addMETTrigSyst(wmcr_vv, year)
    wmcr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(wmcr_vv)
    wmcr_vv.setParamEffect(id_mu, nlepton)
    wmcr_vv.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(wmcr_vv, param, wmcr_central, wmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "VV", "wmcr", wmcr_vv, category, mass)
    wmcr.addSample(wmcr_vv)

    wmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, wmcr_hbbTemplate)
    addLumiSyst(wmcr_hbb, year)
    addPileupSyst(wmcr_hbb)
    addPrefiringSyst(wmcr_hbb, year)
    addJESSyst(wmcr_hbb)
    addMETTrigSyst(wmcr_hbb, year)
    wmcr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(wmcr_hbb)
    wmcr_hbb.setParamEffect(id_mu, nlepton)
    wmcr_hbb.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(wmcr_hbb, param, wmcr_central, wmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "Hbb", "wmcr", wmcr_hbb, category, mass)
    wmcr.addSample(wmcr_hbb)

    wmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, wmcr_qcdTemplate)
    addLumiSyst(wmcr_qcd, year)
    addPileupSyst(wmcr_qcd)
    addPrefiringSyst(wmcr_qcd, year)
    addJESSyst(wmcr_qcd)
    addMETTrigSyst(wmcr_qcd, year)
    wmcr_qcd.setParamEffect(veto_tau, nveto_tau)
    wmcr_qcd.setParamEffect(qcdmu_norm, nqcd_norm)
    wmcr_qcd.setParamEffect(id_mu, nlepton)
    wmcr_qcd.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(wmcr_qcd, param, wmcr_central, wmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "QCD", "wmcr", wmcr_qcd, category, mass)
    wmcr.addSample(wmcr_qcd)

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
        dataTemplate = template(data, "EGamma", "data", recoil, "wecr", category, mass)
    else:
        dataTemplate = template(data, "SingleElectron", "data", recoil, "wecr", category, mass)
    wecr.setObservation(dataTemplate)

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')
    #wecr_central, wecr_error2 = get_mergedMC_stat_variations(background, recoil, "wecr", category, mass, ["W+jets", "TT", "ST", "DY+jets", "VV", "Hbb", "QCD"])

    MCbkgList = ["ST", "DY+jets", "VV", "Hbb", "QCD"]
    if isttMC: MCbkgList.append("TT")
    if iswjetsMC: MCbkgList.append("W+jets")
    wecr_central, wecr_error2 = get_mergedMC_stat_variations(background, recoil, "wecr", category, mass, MCbkgList)
    
    ###
    # W(->lnu)+jets data-driven model
    ###

    if not iswjetsMC:
        wecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wecr", category, mass, min_value=1., read_sumw2=True)
        wecr_wjetsMC = rl.TemplateSample("wecr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wecr_wjetsTemplate)
        wecr_wjetsMC.setParamEffect(trig_e, ntrig_e)
        addEleIDSyst(wecr_wjetsMC, year)
        wecr_wjetsMC.setParamEffect(reco_e, nlepton)
        addVJetsSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)
        addMCStatsTFSyst(wecr_wjetsMC, wecr_wjetsTemplate, sr_wjetsTemplate, epsilon=1e-5)
        
        #### Transfer Factor
        wecr_wjetsTransferFactor = wecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
        #wecr_wjets = rl.ParametericSample(ch_name + "_wjets", rl.Sample.BACKGROUND, sr_zjetsObservable, sr_wjetsBinYields*wecr_wjetsTransferFactor)
        wecr_wjets = rl.TransferFactorSample( ch_name + "_wjets", rl.Sample.BACKGROUND, wecr_wjetsTransferFactor, sr_wjets)
        #addBBliteSyst(wecr_wjets, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
        wecr.addSample(wecr_wjets)

    ###
    # top-antitop data-driven model
    ###

    if not isttMC: 
        wecr_ttTemplate = template(background, "TT", "nominal", recoil, "wecr", category, mass, min_value=1., read_sumw2=True)
        wecr_ttMC = rl.TemplateSample("wecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wecr_ttTemplate)
        wecr_ttMC.setParamEffect(trig_e, ntrig_e)
        addEleIDSyst(wecr_ttMC, year)
        wecr_ttMC.setParamEffect(reco_e, nlepton)
        addBtagSyst(background, recoil, "TT", "wecr", wecr_ttMC, category, mass)
        addMCStatsTFSyst(wecr_ttMC, wecr_ttTemplate, sr_ttTemplate, epsilon=1e-5)
        
        #### Transfer Factor
        wecr_ttTransferFactor = wecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        #wecr_tt = rl.ParametericSample(ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields*wecr_ttTransferFactor)        
        wecr_tt = rl.TransferFactorSample( ch_name + "_tt", rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt)
        #addBBliteSyst(wecr_tt, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
        wecr.addSample(wecr_tt)
    
    ###
    # Other MC-driven processes
    ###

    if iswjetsMC:
        wecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
        wecr_wjets = rl.TemplateSample("wecr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wecr_wjetsTemplate)
        addLumiSyst(wecr_wjets, year)
        addPileupSyst(wecr_wjets)
        addPrefiringSyst(wecr_wjets, year)
        addJESSyst(wecr_wjets)
        addEleIDSyst(wecr_wjets, year)
        wecr_wjets.setParamEffect(veto_tau, nveto_tau)
        addWjetsNormSyst(wecr_wjets)
        wecr_wjets.setParamEffect(id_e, nlepton)
        wecr_wjets.setParamEffect(reco_e, nlepton)
        addBBliteSyst(wecr_wjets, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
        addBtagSyst(background, recoil, "W+jets", "wecr", wecr_wjets, category, mass)
        addVJetsSyst(background, recoil, "W+jets", "wecr", wecr_wjets, category)
        wecr.addSample(wecr_wjets)

    if isttMC: 
        wecr_ttTemplate = template(background, "TT", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
        wecr_tt = rl.TemplateSample("wecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wecr_ttTemplate)
        addLumiSyst(wecr_tt, year)
        addPileupSyst(wecr_tt)
        addPrefiringSyst(wecr_tt, year)
        addJESSyst(wecr_tt)
        wecr_tt.setParamEffect(trig_e, ntrig_e)
        wecr_tt.setParamEffect(veto_tau, nveto_tau)
        addEleIDSyst(wecr_tt, year)
        wecr_tt.setParamEffect(reco_e, nlepton)
        addTTbarNormSyst(wecr_tt)
        addBBliteSyst(wecr_tt, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
        addBtagSyst(background, recoil, "TT", "wecr", wecr_tt, category, mass)
        wecr.addSample(wecr_tt)

    wecr_stTemplate = template(background, "ST", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, wecr_stTemplate)
    addLumiSyst(wecr_st, year)
    addPileupSyst(wecr_st)
    addPrefiringSyst(wecr_st, year)
    addJESSyst(wecr_st)
    wecr_st.setParamEffect(trig_e, ntrig_e)
    wecr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(wecr_st)
    addEleIDSyst(wecr_st, year)
    wecr_st.setParamEffect(reco_e, nlepton)
    addBBliteSyst(wecr_st, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "ST", "wecr", wecr_st, category, mass)
    wecr.addSample(wecr_st)

    wecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wecr_dyjetsTemplate)
    addLumiSyst(wecr_dyjets, year)
    addPileupSyst(wecr_dyjets)
    addPrefiringSyst(wecr_dyjets, year)
    addJESSyst(wecr_dyjets)
    wecr_dyjets.setParamEffect(trig_e, ntrig_e)
    wecr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(wecr_dyjets)
    addEleIDSyst(wecr_dyjets, year)
    wecr_dyjets.setParamEffect(reco_e, nlepton)
    addBBliteSyst(wecr_dyjets, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "DY+jets", "wecr", wecr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "wecr", wecr_dyjets, category)
    wecr.addSample(wecr_dyjets)

    wecr_vvTemplate = template(background, "VV", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, wecr_vvTemplate)
    addLumiSyst(wecr_vv, year)
    addPileupSyst(wecr_vv)
    addPrefiringSyst(wecr_vv, year)
    addJESSyst(wecr_vv)
    wecr_vv.setParamEffect(trig_e, ntrig_e)
    wecr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(wecr_vv)
    addEleIDSyst(wecr_vv, year)
    wecr_vv.setParamEffect(reco_e, nlepton)
    addBBliteSyst(wecr_vv, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "VV", "wecr", wecr_vv, category, mass)
    wecr.addSample(wecr_vv)

    wecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, wecr_hbbTemplate)
    addLumiSyst(wecr_hbb, year)
    addPileupSyst(wecr_hbb)
    addPrefiringSyst(wecr_hbb, year)
    addJESSyst(wecr_hbb)
    wecr_hbb.setParamEffect(trig_e, ntrig_e)
    wecr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(wecr_hbb)
    addEleIDSyst(wecr_hbb, year)
    wecr_hbb.setParamEffect(reco_e, nlepton)
    addBBliteSyst(wecr_hbb, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "Hbb", "wecr", wecr_hbb, category, mass)
    wecr.addSample(wecr_hbb)

    wecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, wecr_qcdTemplate)
    addLumiSyst(wecr_qcd, year)
    addPileupSyst(wecr_qcd)
    addPrefiringSyst(wecr_qcd, year)
    addJESSyst(wecr_qcd)
    wecr_qcd.setParamEffect(trig_e, ntrig_e)
    wecr_qcd.setParamEffect(veto_tau, nveto_tau)
    wecr_qcd.setParamEffect(qcde_norm, nqcd_norm)
    addEleIDSyst(wecr_qcd, year)
    wecr_qcd.setParamEffect(reco_e, nlepton)
    addBBliteSyst(wecr_qcd, param, wecr_central, wecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "QCD", "wecr", wecr_qcd, category, mass)
    wecr.addSample(wecr_qcd)

    ###
    # End of single electron W control region
    ###

    if category=="fail": return model

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

    dataTemplate = template(data, "MET", "data", recoil, "tmcr", category, mass)
    tmcr.setObservation(dataTemplate)

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')
    #tmcr_central, tmcr_error2 = get_mergedMC_stat_variations(background, recoil, "tmcr", category, mass, ["TT", "ST", "DY+jets", "VV", "Hbb", "W+jets", "QCD"])

    MCbkgList = ["ST", "DY+jets", "VV", "Hbb", "W+jets", "QCD"]
    if isttMC: MCbkgList.append("TT")
    tmcr_central, tmcr_error2 = get_mergedMC_stat_variations(background, recoil, "tmcr", category, mass, MCbkgList)


    ###
    # top-antitop data-driven model
    ###

    if not isttMC:
        tmcr_ttTemplate = template(background, "TT", "nominal", recoil, "tmcr", category, mass, min_value=1., read_sumw2=True)
        tmcr_ttMC = rl.TemplateSample("tmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tmcr_ttTemplate)
        addMETTrigSyst(tmcr_ttMC, year)
        tmcr_ttMC.setParamEffect(id_mu, nlepton)
        tmcr_ttMC.setParamEffect(iso_mu, nlepton)
        addBtagSyst(background, recoil, "TT", "tmcr", tmcr_ttMC, category, mass)
        addMCStatsTFSyst(tmcr_ttMC, tmcr_ttTemplate, sr_ttTemplate, epsilon=1e-5)
            
        #### Transfer Factor
        tmcr_ttTransferFactor = tmcr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        #tmcr_tt = rl.ParametericSample(ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields*tmcr_ttTransferFactor)        
        tmcr_tt = rl.TransferFactorSample(ch_name + "_tt", rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt)
        #addBBliteSyst(tmcr_tt, param, tmcr_central, tmcr_error2, epsilon=1e-5)
        tmcr.addSample(tmcr_tt)

    ###
    # Other MC-driven processes
    ###

    if isttMC:
        tmcr_ttTemplate = template(background, "TT", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
        tmcr_tt = rl.TemplateSample("tmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tmcr_ttTemplate)
        addLumiSyst(tmcr_tt, year)
        addPileupSyst(tmcr_tt)
        addPrefiringSyst(tmcr_tt, year)
        addJESSyst(tmcr_tt)
        addMETTrigSyst(tmcr_tt, year)
        tmcr_tt.setParamEffect(veto_tau, nveto_tau)
        addTTbarNormSyst(tmcr_tt)
        tmcr_tt.setParamEffect(id_mu, nlepton)
        tmcr_tt.setParamEffect(iso_mu, nlepton)
        addBBliteSyst(tmcr_tt, param, tmcr_central, tmcr_error2, epsilon=1e-5)
        addBtagSyst(background, recoil, "TT", "tmcr", tmcr_tt, category, mass)
        tmcr.addSample(tmcr_tt)

    tmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_wjets = rl.TemplateSample(ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tmcr_wjetsTemplate)
    addLumiSyst(tmcr_wjets, year)
    addPileupSyst(tmcr_wjets)
    addPrefiringSyst(tmcr_wjets, year)
    addJESSyst(tmcr_wjets)
    addMETTrigSyst(tmcr_wjets, year)
    tmcr_wjets.setParamEffect(veto_tau, nveto_tau)
    addWjetsNormSyst(tmcr_wjets)
    tmcr_wjets.setParamEffect(id_mu, nlepton)
    tmcr_wjets.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(tmcr_wjets, param, tmcr_central, tmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "W+jets", "tmcr", tmcr_wjets, category, mass)
    addVJetsSyst(background, recoil, "W+jets", "tmcr", tmcr_wjets, category)
    tmcr.addSample(tmcr_wjets)

    tmcr_stTemplate = template(background, "ST", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, tmcr_stTemplate)
    addLumiSyst(tmcr_st, year)
    addPileupSyst(tmcr_st)
    addPrefiringSyst(tmcr_st, year)
    addJESSyst(tmcr_st)
    addMETTrigSyst(tmcr_st, year)
    tmcr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(tmcr_st)
    tmcr_st.setParamEffect(id_mu, nlepton)
    tmcr_st.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(tmcr_st, param, tmcr_central, tmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "ST", "tmcr", tmcr_st, category, mass)
    tmcr.addSample(tmcr_st)

    tmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tmcr_dyjetsTemplate)
    addLumiSyst(tmcr_dyjets, year)
    addPileupSyst(tmcr_dyjets)
    addPrefiringSyst(tmcr_dyjets, year)
    addJESSyst(tmcr_dyjets)
    addMETTrigSyst(tmcr_dyjets, year)
    tmcr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(tmcr_dyjets)
    tmcr_dyjets.setParamEffect(id_mu, nlepton)
    tmcr_dyjets.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(tmcr_dyjets, param, tmcr_central, tmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "DY+jets", "tmcr", tmcr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "tmcr", tmcr_dyjets, category)
    tmcr.addSample(tmcr_dyjets)

    tmcr_vvTemplate = template(background, "VV", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, tmcr_vvTemplate)
    addLumiSyst(tmcr_vv, year)
    addPileupSyst(tmcr_vv)
    addPrefiringSyst(tmcr_vv, year)
    addJESSyst(tmcr_vv)
    addMETTrigSyst(tmcr_vv, year)
    tmcr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(tmcr_vv)
    tmcr_vv.setParamEffect(id_mu, nlepton)
    tmcr_vv.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(tmcr_vv, param, tmcr_central, tmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "VV", "tmcr", tmcr_vv, category, mass)
    tmcr.addSample(tmcr_vv)

    tmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, tmcr_hbbTemplate)
    addLumiSyst(tmcr_hbb, year)
    addPileupSyst(tmcr_hbb)
    addPrefiringSyst(tmcr_hbb, year)
    addJESSyst(tmcr_hbb)
    addMETTrigSyst(tmcr_hbb, year)
    tmcr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(tmcr_hbb)
    tmcr_hbb.setParamEffect(id_mu, nlepton)
    tmcr_hbb.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(tmcr_hbb, param, tmcr_central, tmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "Hbb", "tmcr", tmcr_hbb, category, mass)
    tmcr.addSample(tmcr_hbb)

    tmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, tmcr_qcdTemplate)
    addLumiSyst(tmcr_qcd, year)
    addPileupSyst(tmcr_qcd)
    addPrefiringSyst(tmcr_qcd, year)
    addJESSyst(tmcr_qcd)
    addMETTrigSyst(tmcr_qcd, year)
    tmcr_qcd.setParamEffect(veto_tau, nveto_tau)
    tmcr_qcd.setParamEffect(qcdmu_norm, nqcd_norm)
    tmcr_qcd.setParamEffect(id_mu, nlepton)
    tmcr_qcd.setParamEffect(iso_mu, nlepton)
    addBBliteSyst(tmcr_qcd, param, tmcr_central, tmcr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "QCD", "tmcr", tmcr_qcd, category, mass)
    tmcr.addSample(tmcr_qcd)

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
        dataTemplate = template(data, "EGamma", "data", recoil, "tecr", category, mass)
    else:
        dataTemplate = template(data, "SingleElectron", "data", recoil, "tecr", category, mass)
    tecr.setObservation(dataTemplate)

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')
    #tecr_central, tecr_error2 = get_mergedMC_stat_variations(background, recoil, "tecr", category, mass, ["TT", "ST", "DY+jets", "VV", "Hbb", "W+jets", "QCD"])

    MCbkgList = ["ST", "DY+jets", "VV", "Hbb", "W+jets", "QCD"]
    if isttMC: MCbkgList.append("TT")
    tecr_central, tecr_error2 = get_mergedMC_stat_variations(background, recoil, "tecr", category, mass, MCbkgList)
    
    ###
    # top-antitop data-driven model
    ###

    if not isttMC:
        tecr_ttTemplate = template(background, "TT", "nominal", recoil, "tecr", category, mass, min_value=1., read_sumw2=True)
        tecr_ttMC = rl.TemplateSample("tecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tecr_ttTemplate)
        tecr_ttMC.setParamEffect(trig_e, ntrig_e)
        addEleIDSyst(tecr_ttMC, year)
        tecr_ttMC.setParamEffect(reco_e, nlepton)
        addBtagSyst(background, recoil, "TT", "tecr", tecr_ttMC, category, mass)
        addMCStatsTFSyst(tecr_ttMC, tecr_ttTemplate, sr_ttTemplate, epsilon=1e-5)
        
        #### Transfer Factor
        tecr_ttTransferFactor = tecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        #tecr_tt = rl.ParametericSample(ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields*tecr_ttTransferFactor)
        tecr_tt = rl.TransferFactorSample(ch_name + "_tt", rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt)
        #addBBliteSyst(tecr_tt, param, tecr_central, tecr_error2, epsilon=1e-5)
        tecr.addSample(tecr_tt)

    ###
    # Other MC-driven processes
    ###

    if isttMC:
        tecr_ttTemplate = template(background, "TT", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
        tecr_tt = rl.TemplateSample("tecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tecr_ttTemplate)
        addLumiSyst(tecr_tt, year)
        addPileupSyst(tecr_tt)
        addPrefiringSyst(tecr_tt, year)
        addJESSyst(tecr_tt)
        tecr_tt.setParamEffect(trig_e, ntrig_e)
        tecr_tt.setParamEffect(veto_tau, nveto_tau)
        addTTbarNormSyst(tecr_tt)
        addEleIDSyst(tecr_tt, year)
        tecr_tt.setParamEffect(reco_e, nlepton)
        addBBliteSyst(tecr_tt, param, tecr_central, tecr_error2, epsilon=1e-5)
        addBtagSyst(background, recoil, "TT", "tecr", tecr_tt, category, mass)
        tecr.addSample(tecr_tt)

    tecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_wjets = rl.TemplateSample(ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tecr_wjetsTemplate)
    addLumiSyst(tecr_wjets, year)
    addPileupSyst(tecr_wjets)
    addPrefiringSyst(tecr_wjets, year)
    addJESSyst(tecr_wjets)
    tecr_wjets.setParamEffect(trig_e, ntrig_e)
    tecr_wjets.setParamEffect(veto_tau, nveto_tau)
    addWjetsNormSyst(tecr_wjets)
    addEleIDSyst(tecr_wjets, year)
    tecr_wjets.setParamEffect(reco_e, nlepton)
    addBBliteSyst(tecr_wjets, param, tecr_central, tecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "W+jets", "tecr", tecr_wjets, category, mass)
    addVJetsSyst(background, recoil, "W+jets", "tecr", tecr_wjets, category)
    tecr.addSample(tecr_wjets)

    tecr_stTemplate = template(background, "ST", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, tecr_stTemplate)
    addLumiSyst(tecr_st, year)
    addPileupSyst(tecr_st)
    addPrefiringSyst(tecr_st, year)
    addJESSyst(tecr_st)
    tecr_st.setParamEffect(trig_e, ntrig_e)
    tecr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(tecr_st)
    addEleIDSyst(tecr_st, year)
    tecr_st.setParamEffect(reco_e, nlepton)
    addBBliteSyst(tecr_st, param, tecr_central, tecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "ST", "tecr", tecr_st, category, mass)
    tecr.addSample(tecr_st)

    tecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tecr_dyjetsTemplate)
    addLumiSyst(tecr_dyjets, year)
    addPileupSyst(tecr_dyjets)
    addPrefiringSyst(tecr_dyjets, year)
    addJESSyst(tecr_dyjets)
    tecr_dyjets.setParamEffect(trig_e, ntrig_e)
    tecr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(tecr_dyjets)
    addEleIDSyst(tecr_dyjets, year)
    tecr_dyjets.setParamEffect(reco_e, nlepton)
    addBBliteSyst(tecr_dyjets, param, tecr_central, tecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "DY+jets", "tecr", tecr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "tecr", tecr_dyjets, category)
    tecr.addSample(tecr_dyjets)

    tecr_vvTemplate = template(background, "VV", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, tecr_vvTemplate)
    addLumiSyst(tecr_vv, year)
    addPileupSyst(tecr_vv)
    addPrefiringSyst(tecr_vv, year)
    addJESSyst(tecr_vv)
    tecr_vv.setParamEffect(trig_e, ntrig_e)
    tecr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(tecr_vv)
    addEleIDSyst(tecr_vv, year)
    tecr_vv.setParamEffect(reco_e, nlepton)
    addBBliteSyst(tecr_vv, param, tecr_central, tecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "VV", "tecr", tecr_vv, category, mass)
    tecr.addSample(tecr_vv)

    tecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, tecr_hbbTemplate)
    addLumiSyst(tecr_hbb, year)
    addPileupSyst(tecr_hbb)
    addPrefiringSyst(tecr_hbb, year)
    addJESSyst(tecr_hbb)
    tecr_hbb.setParamEffect(trig_e, ntrig_e)
    tecr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(tecr_hbb)
    addEleIDSyst(tecr_hbb, year)
    tecr_hbb.setParamEffect(reco_e, nlepton)
    addBBliteSyst(tecr_hbb, param, tecr_central, tecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "Hbb", "tecr", tecr_hbb, category, mass)
    tecr.addSample(tecr_hbb)

    tecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, tecr_qcdTemplate)
    addLumiSyst(tecr_qcd, year)
    addPileupSyst(tecr_qcd)
    addPrefiringSyst(tecr_qcd, year)
    addJESSyst(tecr_qcd)
    tecr_qcd.setParamEffect(trig_e, ntrig_e)
    tecr_qcd.setParamEffect(veto_tau, nveto_tau)
    tecr_qcd.setParamEffect(qcde_norm, nqcd_norm)
    addEleIDSyst(tecr_qcd, year)
    tecr_qcd.setParamEffect(reco_e, nlepton)
    addBBliteSyst(tecr_qcd, param, tecr_central, tecr_error2, epsilon=1e-5) ### replace autoMCStats
    addBtagSyst(background, recoilbin, "QCD", "tecr", tecr_qcd, category, mass)
    tecr.addSample(tecr_qcd)

    ###
    # End of single electron top control region
    ###

    return model


if __name__ == "__main__":
    if not os.path.exists("datacards"):
        os.mkdir("datacards")
    parser = OptionParser()
    parser.add_option("-y", "--year", help="year", dest="year", default="2018")
    parser.add_option("-m", "--mass", help="mass", dest="mass", default="40to300")
    parser.add_option("-f", "--fakedata", help="replace data to sum of backgrounds", action="store_true", dest="fakedata")
    (options, args) = parser.parse_args()
    year = options.year
    mass = options.mass

    #####
    ###
    # Preparing Rhalphabeth
    ###
    #####

    ###
    # Extract histograms from input file and remap
    ###

    hists = load("hists/darkhiggs" + year + ".scaled")
    hists = remap_histograms(hists)

    ###
    # Preparing histograms for Rhalphabeth
    ###

    background = {}
    for r in hists["bkg"]["template"].identifiers("region"):
        background[str(r)] = hists["bkg"]["template"].integrate("region", r)
    
    ###
    # Establishing 2D binning
    ###
    
    recoilbins = np.array(recoil_binning)
    nrecoil = len(recoilbins) - 1
    msdbins = np.array(mass_binning)
    msd = rl.Observable('fjmass', msdbins)
    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(recoilbins[:-1] + 0.3 * np.diff(recoilbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    recoilscaled = (ptpts - recoil_binning[0]) / (recoil_binning[-1] - recoil_binning[0])
    msdscaled = (msdpts - mass_binning[0]) / (mass_binning[-1] - mass_binning[0])

    ###
    # Calculating average pass-to-fail ratio
    ###
    
    def efficiency(pass_templ, fail_templ, qcdmodel):
        qcdpass, qcdfail = 0., 0.
        for recoilbin in range(nrecoil):
            failCh = rl.Channel("recoilbin%d%s" % (recoilbin, 'fail'))
            passCh = rl.Channel("recoilbin%d%s" % (recoilbin, 'pass'))
            qcdmodel.addChannel(failCh)
            qcdmodel.addChannel(passCh)
            failCh.setObservation(fail_templ[recoilbin])
            passCh.setObservation(pass_templ[recoilbin])
            qcdfail += failCh.getObservation().sum()
            qcdpass += passCh.getObservation().sum()
            '''
            try:
                qcdfail = np.append(qcdfail, [failCh.getObservation()], axis = 0)
                qcdpass = np.append(qcdpass, [passCh.getObservation()], axis = 0)
            except:
                qcdfail = np.array([failCh.getObservation()])
                qcdpass = np.array([passCh.getObservation()])
            ''' 

        return qcdpass / qcdfail

    ###
    # Creating first Bernstein polynomial that represents the MC pass-to-fail ratio
    # Incorporates the dependence on mass/recoil (residual tagger correlation, HF fraction)
    # Includes statistical uncertainty by fitting the by-by-bin MC pass-to-fail ratio
    ###

    zjetspass_templ = []
    zjetsfail_templ = []
    for recoilbin in range(nrecoil):
        zjetspass_templ.append(template(background, "Z+jets", "nominal", recoilbin, "sr", "pass", ""))
        zjetsfail_templ.append(template(background, "Z+jets", "nominal", recoilbin, "sr", "fail", ""))

    zjetsmodel = rl.Model("zjetsmodel")
    zjetseff = efficiency(zjetspass_templ, zjetsfail_templ, zjetsmodel)
    tf_MCtemplZ = rl.BernsteinPoly("tf_MCtemplZ"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_MCtemplZ_params = zjetseff * tf_MCtemplZ(recoilscaled, msdscaled)

    wjetspass_templ = []
    wjetsfail_templ = []
    for recoilbin in range(nrecoil):
        wjetspass_templ.append(template(background, "W+jets", "nominal", recoilbin, "sr", "pass", ""))
        wjetsfail_templ.append(template(background, "W+jets", "nominal", recoilbin, "sr", "fail", ""))

    wjetsmodel = rl.Model("wjetsmodel")
    wjetseff = efficiency(wjetspass_templ, wjetsfail_templ, wjetsmodel)
    tf_MCtemplW = rl.BernsteinPoly("tf_MCtemplW"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_MCtemplW_params = wjetseff * tf_MCtemplW(recoilscaled, msdscaled)
    
    ###
    # Prepare model for the MC ratio fit
    ##

    def rhalphabeth(pass_templ, fail_templ, qcdmodel, tf_MCtempl_params):

        for recoilbin in range(nrecoil):
            failCh = qcdmodel['recoilbin%dfail' % recoilbin]
            passCh = qcdmodel['recoilbin%dpass' % recoilbin]
            failObs = failCh.getObservation()
            qcdparams = np.array([rl.IndependentParameter('qcdparam_ptbin%d_msdbin%d' % (recoilbin, i), 0) for i in range(msd.nbins)])
            sigmascale = 10.
            scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**qcdparams
            fail_qcd = rl.ParametericSample('recoilbin'+str(recoilbin)+'fail_'+qcdmodel.name, rl.Sample.BACKGROUND, msd, scaledparams)
            failCh.addSample(fail_qcd)
            pass_qcd = rl.TransferFactorSample('recoilbin'+str(recoilbin)+'pass_'+qcdmodel.name, rl.Sample.BACKGROUND, tf_MCtempl_params[recoilbin, :], fail_qcd)
            passCh.addSample(pass_qcd)

        return qcdmodel

    zjetsmodel = rhalphabeth(zjetspass_templ, zjetsfail_templ, zjetsmodel, tf_MCtemplZ_params)
    wjetsmodel = rhalphabeth(wjetspass_templ, wjetsfail_templ, wjetsmodel, tf_MCtemplW_params)
    
    ###
    # Perform the fit to the bin-by-bin MC ratio
    ###

    def fit(model):
        qcdfit_ws = ROOT.RooWorkspace('qcdfit_ws')
        simpdf, obs = model.renderRoofit(qcdfit_ws)
        qcdfit = simpdf.fitTo(obs,
                            ROOT.RooFit.Extended(True),
                            ROOT.RooFit.SumW2Error(True),
                            ROOT.RooFit.Strategy(2),
                            ROOT.RooFit.Save(),
                            ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                            ROOT.RooFit.PrintLevel(-1),
                            )
        qcdfit_ws.add(qcdfit)
        #if "pytest" not in sys.modules:
        #    qcdfit_ws.writeToFile(os.path.join(str(tmpdir), 'testModel_qcdfit.root'))
        if qcdfit.status() != 0:
            raise RuntimeError('Could not fit qcd')

        return qcdfit

    zjetsfit = fit(zjetsmodel)
    wjetsfit = fit(wjetsmodel)

    ###
    # Use the post-fit values of the Bernstein polynomial coefficients
    ###
    
    def shape(fit, tf_MCtempl):
        param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', fit, param_names)
        tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
        tf_MCtempl_params_final = tf_MCtempl(recoilscaled, msdscaled)

        return tf_MCtempl_params_final

    tf_MCtemplW_params_final = shape(wjetsfit, tf_MCtemplW)
    tf_MCtemplZ_params_final = shape(zjetsfit, tf_MCtemplZ)
    
    ###
    # Create Bernstein polynomials that represent the correction to the MC ratio
    ###
    
    tf_dataResidualW = rl.BernsteinPoly("tf_dataResidualW"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_dataResidualW_params = tf_dataResidualW(recoilscaled, msdscaled)
    tf_dataResidualZ = rl.BernsteinPoly("tf_dataResidualZ"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_dataResidualZ_params = tf_dataResidualZ(recoilscaled, msdscaled)

    #####
    ###
    # End of Rhalphabeth preparation
    ###
    #####
    
    ###
    ###
    # Prepare histograms for the fit
    ###
    ###
    
    ###
    # Split mass range
    ###
    
    
    if '40to' in mass:
        cut = mass.split('40to')[1]
        index = mass_binning.index(int(cut))
        mass_binning = mass_binning[:(index+1)]
        nmass = len(mass_binning) - 1
        tf_MCtemplZ_params_final = tf_MCtemplZ_params_final[:, :nmass]
        tf_dataResidualZ_params = tf_dataResidualZ_params[:, :nmass]
        tf_MCtemplW_params_final = tf_MCtemplW_params_final[:, :nmass]
        tf_dataResidualW_params = tf_dataResidualW_params[:, :nmass]
    if 'to300' in mass:
        nmass = len(mass_binning) - 1
        cut = mass.split('to300')[0]
        index = mass_binning.index(int(cut))
        mass_binning = mass_binning[index:]
        nmass = nmass - (len(mass_binning) - 1)
        tf_MCtemplZ_params_final = tf_MCtemplZ_params_final[:, nmass:]
        tf_dataResidualZ_params = tf_dataResidualZ_params[:, nmass:]
        tf_MCtemplW_params_final = tf_MCtemplW_params_final[:, nmass:]
        tf_dataResidualW_params = tf_dataResidualW_params[:, nmass:]
        
    ###
    # Reload and remap histograms 
    ###

    hists = load("hists/darkhiggs" + year + ".scaled")
    hists = remap_histograms(hists)

    ###
    # Manipulate histograms to be fed to the model
    ###

    signal_hists = hists["sig"]
    signal = {}
    for r in signal_hists["template"].identifiers("region"):
        signal[str(r)] = signal_hists["template"].integrate("region", r)
        
    bkg_hists = hists["bkg"]
    background = {}
    for r in bkg_hists["template"].identifiers("region"):
        background[str(r)] = bkg_hists["template"].integrate("region", r)

    data_hists = hists["data"]
    data = {}
    for r in data_hists["template"].identifiers("region"):
        data[str(r)] = data_hists["template"].integrate("region", r)

    ###
    ###
    # Set up other systematics
    ###
    ###
    
    lumi = rl.NuisanceParameter("lumi" + year, "lnN")
    lumi_corr = rl.NuisanceParameter("lumi_corr", "lnN")
    lumi_1718 = rl.NuisanceParameter("lumi_1718", "lnN")
    pu = rl.NuisanceParameter("pu" + year, "lnN")
    prefiring = rl.NuisanceParameter("prefiring" + year, "lnN")
    id_e = rl.NuisanceParameter("id_e" + year, "lnN")
    id_mu = rl.NuisanceParameter("id_mu" + year, "lnN")
    id_pho = rl.NuisanceParameter("id_pho" + year, "lnN")
    reco_e = rl.NuisanceParameter("reco_e" + year, "lnN")
    iso_mu = rl.NuisanceParameter("iso_mu" + year, "lnN")
    trig_e = rl.NuisanceParameter("trig_e" + year, "lnN")
    trig_met = rl.NuisanceParameter("trig_met" + year, "lnN")
    trig_pho = rl.NuisanceParameter("trig_pho" + year, "lnN")
    veto_tau = rl.NuisanceParameter("veto_tau" + year, "lnN")
    jes = rl.NuisanceParameter("jes" + year, "lnN")
    met = rl.NuisanceParameter("met" + year, "lnN")
    btag = rl.NuisanceParameter("btag" + year, "shape")  # AK4 btag
    ew1 = rl.NuisanceParameter("ew1", "lnN")
    #ew2G = rl.NuisanceParameter("ew2G", "lnN")
    ew2W = rl.NuisanceParameter("ew2W", "lnN")
    ew2Z = rl.NuisanceParameter("ew2Z", "lnN")
    #ew3G = rl.NuisanceParameter("ew3G", "lnN")
    ew3W = rl.NuisanceParameter("ew3W", "lnN")
    ew3Z = rl.NuisanceParameter("ew3Z", "lnN")
    mix = rl.NuisanceParameter("mix", "lnN")
    #muF = rl.NuisanceParameter("muF", "lnN")
    #muR = rl.NuisanceParameter("muR", "lnN")
    qcd1 = rl.NuisanceParameter("qcd1", "lnN")
    qcd2 = rl.NuisanceParameter("qcd2", "lnN")
    qcd3 = rl.NuisanceParameter("qcd3", "lnN")
    doublebtag={}
    for syst in signal['sr'].identifiers('systematic'):
        if 'doublebtag' not in str(syst): continue
        if 'Up' not in str(syst): continue
        doublebtag[str(syst).replace('Up','')]=rl.NuisanceParameter(str(syst).replace('Up','') + year, "shape")
        
    ###
    # Set lnN numbers
    ###

    ntrig_e = 1.01
    nveto_tau = 1.03
    nlepton = 1.01 ## id_mu, iso_mu, reco_e
    nqcd_norm = 2.0 ## qcdsig_norm, qcde_norm, qcdmu_norm

    ###
    ###
    # End of systematics setup
    ###
    ###

    model_dict = {}
    for recoilbin in range(nrecoil):

        
        #####
        ###
        # Z+jets "fail"
        ###  
        #####
        
        sr_zjetsMCFailTemplate = template(background, "Z+jets", "nominal", recoilbin, "sr", "fail", mass, min_value=1., read_sumw2=True)

        '''
        ch_name_pass = "sr" + year + "pass" + "mass" + mass+ "recoil" + str(recoilbin)
        ch_name_fail = "sr" + year + "fail" + "mass" + mass+ "recoil" + str(recoilbin)
        
        nbins = len(sr_zjetsMCFailTemplate[1]) - 1
        param_pass = [None for _ in range(nbins)]
        param_fail = [None for _ in range(nbins)]
        for i in range(nbins):
            param_pass[i] = rl.NuisanceParameter(ch_name_pass + '_mcstat_bin%i' % i, combinePrior='shape')
            param_fail[i] = rl.NuisanceParameter(ch_name_fail + '_mcstat_bin%i' % i, combinePrior='shape')
    
        sr_central_pass, sr_error2_pass = get_mergedMC_stat_variations(background, recoilbin, "sr", "pass", mass, ["Z+jets", "W+jets", "TT", "ST", "DY+jets", "VV", "Hbb", "QCD"])
        sr_central_fail, sr_error2_fail = get_mergedMC_stat_variations(background, recoilbin, "sr", "fail", mass, ["Z+jets", "W+jets", "TT", "ST", "DY+jets", "VV", "Hbb", "QCD"])
        '''

        sr_zjetsMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoilbin) + "_zjetsMC",
            rl.Sample.BACKGROUND,
            sr_zjetsMCFailTemplate
        )
        addMETTrigSyst(sr_zjetsMCFail, year)
        addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCFail, "fail")

        sr_zjetsObservable = rl.Observable("fjmass"+mass, sr_zjetsMCFailTemplate[1])
        sr_zjetsFailBinYields = np.array([rl.IndependentParameter("sr" + year + "fail" + "mass" + mass + "recoil" + str(recoilbin) + "_zjets_mu"+str(b), sr_zjetsMCFailTemplate[0][b], 1e-5, sr_zjetsMCFailTemplate[0].max()*2) for b in range(len(sr_zjetsMCFailTemplate[0]))])

        sr_zjetsFail = rl.ParametericSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoilbin) + "_zjets",
            rl.Sample.BACKGROUND,
            sr_zjetsObservable,
            sr_zjetsFailBinYields
        )
        #addBBliteSyst(sr_zjetsFail, param_fail, sr_central_fail, sr_error2_fail, epsilon=1e-5) ### replace autoMCStats

        #####
        ###
        # W+jets "fail"
        ###  
        #####
      
        sr_wjetsMCFailTemplate = template(background, "W+jets", "nominal", recoilbin, "sr", "fail", mass, min_value=1., read_sumw2=True)
        sr_wjetsMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoilbin) + "_wjetsMC",
            rl.Sample.BACKGROUND,
            sr_wjetsMCFailTemplate
        )
        addMETTrigSyst(sr_wjetsMCFail, year)
        addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCFail, "fail")
        addMCStatsTFSyst(sr_wjetsMCFail, sr_wjetsMCFailTemplate, sr_zjetsMCFailTemplate, epsilon=1e-5)

        sr_wjetsFailTransferFactor = sr_wjetsMCFail.getExpectation() / sr_zjetsMCFail.getExpectation()
        '''
        sr_wjetsFailBinYields = sr_zjetsFailBinYields*sr_wjetsFailTransferFactor
        sr_wjetsFail = rl.ParametericSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoilbin) + "_wjets",
            rl.Sample.BACKGROUND,
            sr_zjetsObservable,
            sr_wjetsFailBinYields
        )
        '''
        sr_wjetsFail = rl.TransferFactorSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoilbin) + "_wjets",
            rl.Sample.BACKGROUND,
            sr_wjetsFailTransferFactor,
            sr_zjetsFail
        )
        #addBBliteSyst(sr_wjetsFail, param_fail, sr_central_fail, sr_error2_fail, epsilon=1e-5) ### replace autoMCStats

        #####
        ###
        # Z+jets "pass"
        ###  
        ##### 

        sr_zjetsMCPassTemplate = template(background, "Z+jets", "nominal", recoilbin, "sr", "pass", mass, min_value=1., read_sumw2=True)
        sr_zjetsMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoilbin) + "_zjetsMC",
            rl.Sample.BACKGROUND,
            sr_zjetsMCPassTemplate
        )
        addMETTrigSyst(sr_zjetsMCPass, year)
        addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCPass, "pass")
        addMCStatsTFSyst(sr_zjetsMCPass, sr_zjetsMCPassTemplate, sr_zjetsMCFailTemplate, epsilon=1e-5)

        tf_MCtemplZ = sr_zjetsMCPass.getExpectation() / sr_zjetsMCFail.getExpectation()
        #tf_paramsZ = zjetseff *tf_MCtemplZ_params_final[recoilbin, :] * tf_dataResidualZ_params[recoilbin, :]
        tf_paramsZ = tf_MCtemplZ * tf_dataResidualZ_params[recoilbin, :]

        '''        
        sr_zjetsPassBinYields=sr_zjetsFailBinYields*tf_paramsZ
        sr_zjetsPass = rl.ParametericSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoilbin) + "_zjets",
            rl.Sample.BACKGROUND,
            sr_zjetsObservable,
            sr_zjetsPassBinYields
        )
        '''
        sr_zjetsPass = rl.TransferFactorSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoilbin) + "_zjets",
            rl.Sample.BACKGROUND,
            tf_paramsZ,
            sr_zjetsFail
        )
        #addBBliteSyst(sr_zjetsPass, param_pass, sr_central_pass, sr_error2_pass, epsilon=1e-5) ### replace autoMCStats
        
        #####
        ###
        # W+jets "pass"
        ###
        #####

        sr_wjetsMCPassTemplate = template(background, "W+jets", "nominal", recoilbin, "sr", "pass", mass, min_value=1., read_sumw2=True)
        sr_wjetsMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoilbin) + "_wjetsMC",
            rl.Sample.BACKGROUND,
            sr_wjetsMCPassTemplate
        )
        addMETTrigSyst(sr_wjetsMCPass, year)
        addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCPass, "pass")
        addMCStatsTFSyst(sr_wjetsMCPass, sr_wjetsMCPassTemplate, sr_wjetsMCFailTemplate, epsilon=1e-5)
        
        tf_MCtemplW = sr_wjetsMCPass.getExpectation() / sr_wjetsMCFail.getExpectation()
        #tf_paramsW = wjetseff * tf_MCtemplW_params_final[recoilbin, :] * tf_dataResidualW_params[recoilbin, :]
        tf_paramsW = tf_MCtemplW * tf_dataResidualW_params[recoilbin, :]
    
        '''
        sr_wjetsPassBinYields = sr_wjetsFailBinYields*tf_paramsW
        sr_wjetsPass = rl.ParametericSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoilbin) + "_wjets",
            rl.Sample.BACKGROUND,
            sr_zjetsObservable,
            sr_wjetsPassBinYields
        )
        '''
        sr_wjetsPass = rl.TransferFactorSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoilbin) + "_wjets",
            rl.Sample.BACKGROUND,
            tf_paramsW,
            sr_wjetsFail
        )
        #addBBliteSyst(sr_wjetsPass, param_pass, sr_central_pass, sr_error2_pass, epsilon=1e-5) ### replace autoMCStats

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
            
            isttMC = ('40to' in mass and not 'to300' in mass) | (category=='fail') | (recoilbin==4)
            iswjetsMC = (recoilbin==4) & (category=='pass')

            with open(
                    "data/models/"
                    + "darkhiggs"
                    + "-"
                    + year
                    + "-"
                    + category
                    + "-mass"
                    + mass
                    + "-recoil"
                    + str(recoilbin)
                    + ".model",
                    "wb",
            ) as fout:
                pickle.dump(model(year, mass, recoilbin, category), fout, protocol=2)
