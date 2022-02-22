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

with open("data/hf_systematic.json") as fin:
    hf_systematic = json.load(fin)

def simple_error_propagation(pw, tw, pw2, tw2, debug=False):
    """Compute errors based on the propagation of uncertainty
    Parameters
    ----------
    pw : numpy.ndarray
        Numerator, or number of (weighted) successes, vectorized
    tw : numpy.ndarray
        Denominator or number of (weighted) trials, vectorized
    pw2 : numpy.ndarray
        Numerator sum of weights squared, vectorized
    tw2 : numpy.ndarray
        Denominator sum of weights squared, vectorized
    """
    dx = np.sqrt(pw2)
    dy = np.sqrt(tw2)
    ratio = tw / pw
    dz = ratio * np.sqrt((dx/pw)**2 + (dy/tw)**2)
    if debug:
        print('================ Propation of uncertainty ================')
        print('ratio:', ratio)
        print('dz:', dz)
        print('================ Propation of uncertainty ================ \n')

    return dz

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

    #print(string, systUp)
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
        systDown = np.array( down.sum() / nominal.sum() )
        systDown[np.isnan(systDown)] = 1.
        systDown = systDown.sum()
        templ.setParamEffect(syst, systUp, systDown)
    addSyst(dictionary, recoil, process, region, templ, category, ew1, "ew1")
    addSyst(dictionary, recoil, process, region, templ, category, ew2W, "ew2W")
    addSyst(dictionary, recoil, process, region, templ, category, ew2Z, "ew2Z")
    addSyst(dictionary, recoil, process, region, templ, category, ew3W, "ew3W")
    addSyst(dictionary, recoil, process, region, templ, category, ew3Z, "ew3Z")
    addSyst(dictionary, recoil, process, region, templ, category, mix, "mix")
    addSyst(dictionary, recoil, process, region, templ, category, qcd1, "qcd1")
    addSyst(dictionary, recoil, process, region, templ, category, qcd2, "qcd2")
    addSyst(dictionary, recoil, process, region, templ, category, qcd3, "qcd3")

def rhalphabeth2D(process, tf_dataResidual_params, ord1, ord2):

    process_map = {
        "W+jets": 'W',
        "Z+jets": 'Z'
        }

    # Build qcd MC pass+fail model and fit to polynomial
    qcdmodel = rl.Model("qcdmodel")
    qcdpass, qcdfail = 0., 0.
    for recoilbin in range(nrecoil):
        failCh = rl.Channel(process_map[process]+"recoil%d%s" % (recoilbin, 'fail'))
        passCh = rl.Channel(process_map[process]+"recoil%d%s" % (recoilbin, 'pass'))
        qcdmodel.addChannel(failCh)
        qcdmodel.addChannel(passCh)
        # mock template
        ptnorm = 1
        #add templates
        failTempl = template(background, process, "nominal", recoilbin, "sr", "fail")
        passTempl = template(background, process, "nominal", recoilbin, "sr", "pass")
        failCh.setObservation(failTempl)
        passCh.setObservation(passTempl)
        qcdfail += failCh.getObservation().sum()
        qcdpass += passCh.getObservation().sum()

    qcdeff = qcdpass / qcdfail
    tf_MCtempl = rl.BernsteinPoly("tf_MCtempl"+process_map[process], (ord1, ord2), ['recoil', 'fjmass'], limits=(0, 10))
    tf_MCtempl_params = qcdeff * tf_MCtempl(recoilscaled, msdscaled)
    for recoilbin in range(nrecoil):
        failCh = qcdmodel[process_map[process]+'recoil%dfail' % recoilbin]
        passCh = qcdmodel[process_map[process]+'recoil%dpass' % recoilbin]
        failObs = failCh.getObservation()
        qcdparams = np.array([rl.IndependentParameter(process_map[process]+'param_recoilbin%d_msdbin%d' % (recoilbin, i), 0) for i in range(msd.nbins)])
        sigmascale = 10.
        scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**qcdparams
        fail_qcd = rl.ParametericSample(process_map[process]+'recoil%dfail_qcd' % recoilbin, rl.Sample.BACKGROUND, msd, scaledparams)
        failCh.addSample(fail_qcd)
        print(tf_MCtempl_params[recoilbin, :])
        pass_qcd = rl.TransferFactorSample(process_map[process]+'recoil%dpass_qcd' % recoilbin, rl.Sample.BACKGROUND, tf_MCtempl_params[recoilbin, :], fail_qcd)
        passCh.addSample(pass_qcd)

        #failCh.mask = validbins[recoilbin]
        #passCh.mask = validbins[recoilbin]

    qcdfit_ws = ROOT.RooWorkspace('qcdfit_ws')
    simpdf, obs = qcdmodel.renderRoofit(qcdfit_ws)
    qcdfit = simpdf.fitTo(obs,
                          ROOT.RooFit.Extended(True),
                          ROOT.RooFit.SumW2Error(True),
                          ROOT.RooFit.Strategy(2),
                          ROOT.RooFit.Save(),
                          ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                          ROOT.RooFit.PrintLevel(-1),
                          )
    qcdfit_ws.add(qcdfit)
    if "pytest" not in sys.modules:
         #qcdfit_ws.writeToFile(os.path.join(str(tmpdir), 'testModel_qcdfit.root'))
        qcdfit_ws.writeToFile(os.path.join(str("models"), process_map[process]+"testModel_qcdfit.root"))
    if qcdfit.status() != 0:
        raise RuntimeError('Could not fit qcd')

    param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
    decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco'+process_map[process], qcdfit, param_names)
    tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
    tf_MCtempl_params_final = tf_MCtempl(recoilscaled, msdscaled)
    tf_params = qcdeff * tf_MCtempl_params_final * tf_dataResidual_params
    return tf_params

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
        sr.setObservation(template(fake_data, "", "data", recoil, "sr", category, bkg=True))
    else:
        sr.setObservation(template(data, "MET", "data", recoil, "sr", category))

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
            #sr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "sr", category, read_sumw2=True)
            sr_wjetsMC = rl.TemplateSample(
                "sr" + model_id + "_wjetsMC",
                rl.Sample.BACKGROUND,
                sr_wjetsTemplate
            )
            sr_wjetsMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
            sr_wjetsMC.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
            sr_wjetsMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
            sr_wjetsMC.setParamEffect(wjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
            sr_wjetsMC.setParamEffect(jec, njec, np.reciprocal(njec))
            #sr_wjetsMC.autoMCStats(shape=True, name="sr"+year+category+"_wjetsMC", epsilon=1e-5)
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
    sr_ttMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    sr_ttMC.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    sr_ttMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    sr_ttMC.setParamEffect(jec, njec, np.reciprocal(njec))
    #sr_ttMC.autoMCStats(shape=True, name="sr"+year+category+"_ttMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "TT", "sr", sr_ttMC, category)

    if category == "pass" and recoil<4:
        sr_ttMC.setParamEffect(tt_norm, nMinor_norm, np.reciprocal(nMinor_norm))
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
                    0,
                    -1*sigmascale[year],
                    sigmascale[year],
                )
                for i in range(sr_ttObservable.nbins)
            ]
        )
        sr_ttBinYields = sr_ttTemplate[0] * (1 + (sigmascale[year]/np.maximum(1., np.sqrt(sr_ttTemplate[0]))))**sr_ttParameters
        #print(sr_ttBinYields)

        sr_tt = rl.ParametericSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields
        )
        sr.addSample(sr_tt)
    else:
        sr_ttMC.setParamEffect(ttMC_norm, nMinor_norm, np.reciprocal(nMinor_norm))
        sr.addSample(sr_ttMC)

    ###
    # Other MC-driven processes
    ###

    sr_stTemplate = template(background, "ST", "nominal", recoil, "sr", category, read_sumw2=True)
    sr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, sr_stTemplate)
    sr_st.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    sr_st.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    sr_st.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    sr_st.setParamEffect(st_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    sr_st.setParamEffect(jec, njec, np.reciprocal(njec))
    #sr_st.autoMCStats(shape=True, name="sr"+year+category+"_stMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "ST", "sr", sr_st, category)
    sr.addSample(sr_st)

    sr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "sr", category, read_sumw2=True)
    sr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, sr_dyjetsTemplate
    )
    sr_dyjets.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    sr_dyjets.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    sr_dyjets.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    sr_dyjets.setParamEffect(zjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    sr_dyjets.setParamEffect(jec, njec, np.reciprocal(njec))
    #sr_dyjets.autoMCStats(shape=True, name="sr"+year+category+"_dyjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "DY+jets", "sr", sr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "sr", sr_dyjets, category)
    sr.addSample(sr_dyjets)

    sr_vvTemplate = template(background, "VV", "nominal", recoil, "sr", category, read_sumw2=True)
    sr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, sr_vvTemplate)
    sr_vv.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    sr_vv.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    sr_vv.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    sr_vv.setParamEffect(vv_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    sr_vv.setParamEffect(jec, njec, np.reciprocal(njec))
    #sr_vv.autoMCStats(shape=True, name="sr"+year+category+"_vvMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "VV", "sr", sr_vv, category)
    sr.addSample(sr_vv)

    sr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "sr", category, read_sumw2=True)
    sr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, sr_hbbTemplate)
    sr_hbb.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    sr_hbb.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    sr_hbb.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    sr_hbb.setParamEffect(hbb_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    sr_hbb.setParamEffect(jec, njec, np.reciprocal(njec))
    #sr_hbb.autoMCStats(shape=True, name="sr"+year+category+"_hbbMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "Hbb", "sr", sr_hbb, category)
    sr.addSample(sr_hbb)

    sr_qcdTemplate = template(background, "QCD", "nominal", recoil, "sr", category, read_sumw2=True)
    sr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, sr_qcdTemplate)
    sr_qcd.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    sr_qcd.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    sr_qcd.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    sr_qcd.setParamEffect(qcdsig_norm, nqcd_norm, np.reciprocal(nqcd_norm))
    sr_qcd.setParamEffect(jec, njec, np.reciprocal(njec))
    #sr_qcd.autoMCStats(shape=True, name="sr"+year+category+"_qcdMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "QCD", "sr", sr_qcd, category)
    sr.addSample(sr_qcd)

    sr_signalTemplate = template(signal, s, "nominal", recoil, "sr", category)
    sr_signal = rl.TemplateSample(
        ch_name + "_" + str(s), rl.Sample.SIGNAL, sr_signalTemplate
    )
    sr_signal.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    sr_signal.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    sr_signal.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    sr_signal.setParamEffect(jec, njec, np.reciprocal(njec))
    addBtagSyst(signal, recoil, str(s), "sr", sr_signal, category)
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

    wmcr.setObservation(template(data, "MET", "data", recoil, "wmcr", category))

    ###
    # W(->lnu)+jets data-driven model
    ###

    wmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_wjetsMC = rl.TemplateSample(
        "wmcr" + model_id + "_wjetsMC",
        rl.Sample.BACKGROUND,
        wmcr_wjetsTemplate
    )
    wmcr_wjetsMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wmcr_wjetsMC.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    wmcr_wjetsMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wmcr_wjetsMC.setParamEffect(wjets_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    wmcr_wjetsMC.setParamEffect(jec, njec, np.reciprocal(njec))
    wmcr_wjetsMC.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    wmcr_wjetsMC.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #wmcr_wjetsMC.autoMCStats(shape=True, name="wmcr"+year+category+"_wjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)
    addVJetsSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)

    #wmcr_wjetsTFstatParameters =  np.array([rl.NuisanceParameter("wmcr_"+year+"_wjetsTFstat_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i, "shape") for i in range(wmcr_wjetsTemplate[0].size)])
    wmcr_wjetsTransferFactor = wmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
    #nominal =  wmcr_wjetsTemplate[0] / sr_wjetsTemplate[0]
    #dz = simple_error_propagation(sr_wjetsTemplate[0], wmcr_wjetsTemplate[0], sr_wjetsTemplate[3], wmcr_wjetsTemplate[3])
    #print('wmcr wjets TF', dz/nominal)
    #wmcr_wjetsTransferFactor = wmcr_wjetsTransferFactor * ( 1. + (dz/nominal)*wmcr_wjetsTFstatParameters )
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
    wmcr_ttMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wmcr_ttMC.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    wmcr_ttMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wmcr_ttMC.setParamEffect(jec, njec, np.reciprocal(njec))
    wmcr_ttMC.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    wmcr_ttMC.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #wmcr_ttMC.autoMCStats(shape=True, name="wmcr"+year+category+"_ttMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "TT", "wmcr", wmcr_ttMC, category)

    if category == "pass":
        wmcr_ttMC.setParamEffect(tt_norm, nMinor_norm, np.reciprocal(nMinor_norm))
        #wmcr_ttTFstatParameters =  np.array([rl.NuisanceParameter("wmcr_"+year+"_ttTFstat_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i, "shape") for i in range(wmcr_ttTemplate[0].size)])
        wmcr_ttTransferFactor = wmcr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        #nominal =  wmcr_ttTemplate[0] / sr_wjetsTemplate[0]
        #dz = simple_error_propagation(sr_ttTemplate[0], wmcr_ttTemplate[0], sr_ttTemplate[3], wmcr_ttTemplate[3])
        #print('wmcr tt TF', dz/nominal)
        #wmcr_ttTransferFactor = wmcr_ttTransferFactor * ( 1. + (dz/nominal)*wmcr_ttTFstatParameters )
        wmcr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt
        )
        wmcr.addSample(wmcr_tt)
    else:
        wmcr_ttMC.setParamEffect(ttMC_norm, nMinor_norm, np.reciprocal(nMinor_norm))
        wmcr.addSample(wmcr_ttMC)

    ###
    # Other MC-driven processes
    ###

    wmcr_stTemplate = template(background, "ST", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, wmcr_stTemplate
        )
    wmcr_st.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wmcr_st.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    wmcr_st.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wmcr_st.setParamEffect(st_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    wmcr_st.setParamEffect(jec, njec, np.reciprocal(njec))
    wmcr_st.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    wmcr_st.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #wmcr_st.autoMCStats(shape=True, name="wmcr"+year+category+"_stMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "ST", "wmcr", wmcr_st, category)
    wmcr.addSample(wmcr_st)

    wmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wmcr_dyjetsTemplate
    )
    wmcr_dyjets.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wmcr_dyjets.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    wmcr_dyjets.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wmcr_dyjets.setParamEffect(zjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    wmcr_dyjets.setParamEffect(jec, njec, np.reciprocal(njec))
    wmcr_dyjets.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    wmcr_dyjets.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #wmcr_dyjets.autoMCStats(shape=True, name="wmcr"+year+category+"_dyjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "DY+jets", "wmcr", wmcr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "wmcr", wmcr_dyjets, category)
    wmcr.addSample(wmcr_dyjets)

    wmcr_vvTemplate = template(background, "VV", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, wmcr_vvTemplate
    )
    wmcr_vv.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wmcr_vv.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    wmcr_vv.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wmcr_vv.setParamEffect(vv_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    wmcr_vv.setParamEffect(jec, njec, np.reciprocal(njec))
    wmcr_vv.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    wmcr_vv.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #wmcr_vv.autoMCStats(shape=True, name="wmcr"+year+category+"_vvMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "VV", "wmcr", wmcr_vv, category)
    wmcr.addSample(wmcr_vv)

    wmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, wmcr_hbbTemplate
    )
    wmcr_hbb.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wmcr_hbb.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    wmcr_hbb.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wmcr_hbb.setParamEffect(hbb_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    wmcr_hbb.setParamEffect(jec, njec, np.reciprocal(njec))
    wmcr_hbb.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    wmcr_hbb.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #wmcr_hbb.autoMCStats(shape=True, name="wmcr"+year+category+"_hbbMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "Hbb", "wmcr", wmcr_hbb, category)
    wmcr.addSample(wmcr_hbb)

    wmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, wmcr_qcdTemplate
    )
    wmcr_qcd.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wmcr_qcd.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    wmcr_qcd.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wmcr_qcd.setParamEffect(qcdmu_norm, nqcd_norm, np.reciprocal(nqcd_norm))
    wmcr_qcd.setParamEffect(jec, njec, np.reciprocal(njec))
    wmcr_qcd.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    wmcr_qcd.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #wmcr_qcd.autoMCStats(shape=True, name="wmcr"+year+category+"_qcdMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "QCD", "wmcr", wmcr_qcd, category)
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
        wecr.setObservation(template(data, "EGamma", "data", recoil, "wecr", category))
    else:
        wecr.setObservation(template(data, "SingleElectron", "data", recoil, "wecr", category))

    ###
    # W(->lnu)+jets data-driven model
    ###

    wecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_wjetsMC = rl.TemplateSample(
        "wecr" + model_id + "_wjetsMC",
        rl.Sample.BACKGROUND,
        wecr_wjetsTemplate
    )
    wecr_wjetsMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wecr_wjetsMC.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    wecr_wjetsMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wecr_wjetsMC.setParamEffect(wjets_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    wecr_wjetsMC.setParamEffect(jec, njec, np.reciprocal(njec))
    wecr_wjetsMC.setParamEffect(id_e, nlepton), np.reciprocal(nlepton)
    wecr_wjetsMC.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #wecr_wjetsMC.autoMCStats(shape=True, name="wecr"+year+category+"_wjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)
    addVJetsSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)

    #wecr_wjetsTFstatParameters =  np.array([rl.NuisanceParameter("wecr_"+year+"_wjetsTFstat_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i, "shape") for i in range(wecr_wjetsTemplate[0].size)])
    wecr_wjetsTransferFactor = wecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
    #nominal =  wecr_wjetsTemplate[0] / sr_wjetsTemplate[0]
    #dz = simple_error_propagation(sr_wjetsTemplate[0], wecr_wjetsTemplate[0], sr_wjetsTemplate[3], wecr_wjetsTemplate[3])
    #print('wecr wjets TF', dz/nominal)
    #wecr_wjetsTransferFactor = wecr_wjetsTransferFactor * ( 1. + (dz/nominal)*wecr_wjetsTFstatParameters )
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
    wecr_ttMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wecr_ttMC.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    wecr_ttMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wecr_ttMC.setParamEffect(jec, njec, np.reciprocal(njec))
    wecr_ttMC.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    wecr_ttMC.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #wecr_ttMC.autoMCStats(shape=True, name="wecr"+year+category+"_ttMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "TT", "wecr", wecr_ttMC, category)

    if category == "pass":
        wecr_ttMC.setParamEffect(tt_norm, nMinor_norm, np.reciprocal(nMinor_norm))
        #wecr_ttTFstatParameters =  np.array([rl.NuisanceParameter("wecr_"+year+"_ttTFstat_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i, "shape") for i in range(wecr_ttTemplate[0].size)])
        wecr_ttTransferFactor = wecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        #nominal =  wecr_ttTemplate[0] / sr_ttTemplate[0]
        #dz = simple_error_propagation(sr_ttTemplate[0], wecr_ttTemplate[0], sr_ttTemplate[3], wecr_ttTemplate[3])
        #print('wecr tt TF', dz/nominal)
        #wecr_ttTransferFactor = wecr_wjetsTransferFactor * ( 1. + (dz/nominal)*wecr_ttTFstatParameters )
        wecr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt
        )
        wecr.addSample(wecr_tt)
    else:
        wecr_ttMC.setParamEffect(ttMC_norm, nMinor_norm, np.reciprocal(nMinor_norm))
        wecr.addSample(wecr_ttMC)

    ###
    # Other MC-driven processes
    ###

    wecr_stTemplate = template(background, "ST", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, wecr_stTemplate
    )
    wecr_st.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wecr_st.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    wecr_st.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wecr_st.setParamEffect(st_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    wecr_st.setParamEffect(jec, njec, np.reciprocal(njec))
    wecr_st.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    wecr_st.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #wecr_st.autoMCStats(shape=True, name="wecr"+year+category+"_stMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "ST", "wecr", wecr_st, category)
    wecr.addSample(wecr_st)

    wecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wecr_dyjetsTemplate
    )
    wecr_dyjets.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wecr_dyjets.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    wecr_dyjets.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wecr_dyjets.setParamEffect(zjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    wecr_dyjets.setParamEffect(jec, njec, np.reciprocal(njec))
    wecr_dyjets.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    wecr_dyjets.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #wecr_dyjets.autoMCStats(shape=True, name="wecr"+year+category+"_dyjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "DY+jets", "wecr", wecr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "wecr", wecr_dyjets, category)
    wecr.addSample(wecr_dyjets)

    wecr_vvTemplate = template(background, "VV", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, wecr_vvTemplate
    )
    wecr_vv.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wecr_vv.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    wecr_vv.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wecr_vv.setParamEffect(vv_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    wecr_vv.setParamEffect(jec, njec, np.reciprocal(njec))
    wecr_vv.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    wecr_vv.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #wecr_vv.autoMCStats(shape=True, name="wecr"+year+category+"_vvMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "VV", "wecr", wecr_vv, category)
    wecr.addSample(wecr_vv)

    wecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, wecr_hbbTemplate
    )
    wecr_hbb.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wecr_hbb.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    wecr_hbb.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wecr_hbb.setParamEffect(hbb_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    wecr_hbb.setParamEffect(jec, njec, np.reciprocal(njec))
    wecr_hbb.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    wecr_hbb.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #wecr_hbb.autoMCStats(shape=True, name="wecr"+year+category+"_hbbMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "Hbb", "wecr", wecr_hbb, category)
    wecr.addSample(wecr_hbb)

    wecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wecr", category, read_sumw2=True)
    wecr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, wecr_qcdTemplate
    )
    wecr_qcd.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    wecr_qcd.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    wecr_qcd.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    wecr_qcd.setParamEffect(qcde_norm, nqcd_norm)
    wecr_qcd.setParamEffect(jec, njec, np.reciprocal(njec))
    wecr_qcd.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    wecr_qcd.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #wecr_qcd.autoMCStats(shape=True, name="wecr"+year+category+"_qcdMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "QCD", "wecr", wecr_qcd, category)
    wecr.addSample(wecr_qcd)

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

    tmcr.setObservation(template(data, "MET", "data", recoil, "tmcr", category))

    ###
    # top-antitop model
    ###

    tmcr_ttTemplate = template(background, "TT", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_ttMC = rl.TemplateSample(
        "tmcr" + model_id + "_ttMC",
        rl.Sample.BACKGROUND,
        tmcr_ttTemplate
    )
    tmcr_ttMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tmcr_ttMC.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    tmcr_ttMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tmcr_ttMC.setParamEffect(tt_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tmcr_ttMC.setParamEffect(jec, njec, np.reciprocal(njec))
    tmcr_ttMC.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    tmcr_ttMC.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #tmcr_ttMC.autoMCStats(shape=True, name="tmcr"+year+category+"_ttMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "TT", "tmcr", tmcr_ttMC, category)

    #tmcr_ttTFstatParameters =  np.array([rl.NuisanceParameter("tmcr_"+year+"_ttTFstat_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i, "shape") for i in range(tmcr_ttTemplate[0].size)])
    tmcr_ttTransferFactor = tmcr_ttMC.getExpectation() / sr_ttMC.getExpectation()
    #nominal =  tmcr_ttTemplate[0] / sr_ttTemplate[0]
    #dz = simple_error_propagation(sr_ttTemplate[0], tmcr_ttTemplate[0], sr_ttTemplate[3], tmcr_ttTemplate[3])
    #print('tmcr tt TF', dz/nominal)
    #tmcr_ttTransferFactor = tmcr_ttTransferFactor * ( 1. + (dz/nominal)*tmcr_ttTFstatParameters )
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
    tmcr_wjets.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tmcr_wjets.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    tmcr_wjets.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tmcr_wjets.setParamEffect(wjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    tmcr_wjets.setParamEffect(jec, njec, np.reciprocal(njec))
    tmcr_wjets.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    tmcr_wjets.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #tmcr_wjets.autoMCStats(shape=True, name="tmcr"+year+category+"_wjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "W+jets", "tmcr", tmcr_wjets, category)
    addVJetsSyst(background, recoil, "W+jets", "tmcr", tmcr_wjets, category)
    tmcr.addSample(tmcr_wjets)

    tmcr_stTemplate = template(background, "ST", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, tmcr_stTemplate
    )
    tmcr_st.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tmcr_st.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    tmcr_st.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tmcr_st.setParamEffect(st_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tmcr_st.setParamEffect(jec, njec, np.reciprocal(njec))
    tmcr_st.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    tmcr_st.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #tmcr_st.autoMCStats(shape=True, name="tmcr"+year+category+"_stMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "ST", "tmcr", tmcr_st, category)
    tmcr.addSample(tmcr_st)

    tmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tmcr_dyjetsTemplate
    )
    tmcr_dyjets.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tmcr_dyjets.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    tmcr_dyjets.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tmcr_dyjets.setParamEffect(zjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    tmcr_dyjets.setParamEffect(jec, njec, np.reciprocal(njec))
    tmcr_dyjets.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    tmcr_dyjets.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #tmcr_dyjets.autoMCStats(shape=True, name="tmcr"+year+category+"_dyjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "DY+jets", "tmcr", tmcr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "tmcr", tmcr_dyjets, category)
    tmcr.addSample(tmcr_dyjets)

    tmcr_vvTemplate = template(background, "VV", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, tmcr_vvTemplate
    )
    tmcr_vv.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tmcr_vv.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    tmcr_vv.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tmcr_vv.setParamEffect(vv_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tmcr_vv.setParamEffect(jec, njec, np.reciprocal(njec))
    tmcr_vv.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    tmcr_vv.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #tmcr_vv.autoMCStats(shape=True, name="tmcr"+year+category+"_vvMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "VV", "tmcr", tmcr_vv, category)
    tmcr.addSample(tmcr_vv)

    tmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, tmcr_hbbTemplate
    )
    tmcr_hbb.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tmcr_hbb.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    tmcr_hbb.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tmcr_hbb.setParamEffect(hbb_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tmcr_hbb.setParamEffect(jec, njec, np.reciprocal(njec))
    tmcr_hbb.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    tmcr_hbb.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #tmcr_hbb.autoMCStats(shape=True, name="tmcr"+year+category+"_hbbMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "Hbb", "tmcr", tmcr_hbb, category)
    tmcr.addSample(tmcr_hbb)

    tmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tmcr", category, read_sumw2=True)
    tmcr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, tmcr_qcdTemplate
    )
    tmcr_qcd.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tmcr_qcd.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
    tmcr_qcd.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tmcr_qcd.setParamEffect(qcdmu_norm, nqcd_norm)
    tmcr_qcd.setParamEffect(jec, njec, np.reciprocal(njec))
    tmcr_qcd.setParamEffect(id_mu, nlepton, np.reciprocal(nlepton))
    tmcr_qcd.setParamEffect(iso_mu, nlepton, np.reciprocal(nlepton))
    #tmcr_qcd.autoMCStats(shape=True, name="tmcr"+year+category+"_qcdMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "QCD", "tmcr", tmcr_qcd, category)
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
        tecr.setObservation(template(data, "EGamma", "data", recoil, "tecr", category))
    else:
        tecr.setObservation(template(data, "SingleElectron", "data", recoil, "tecr", category))

    ###
    # top-antitop model
    ###

    tecr_ttTemplate = template(background, "TT", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_ttMC = rl.TemplateSample(
        "tecr" + model_id + "_ttMC",
        rl.Sample.BACKGROUND,
        tecr_ttTemplate
    )
    tecr_ttMC.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tecr_ttMC.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    tecr_ttMC.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tecr_ttMC.setParamEffect(tt_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tecr_ttMC.setParamEffect(jec, njec, np.reciprocal(njec))
    tecr_ttMC.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    tecr_ttMC.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #tecr_ttMC.autoMCStats(shape=True, name="tecr"+year+category+"_ttMC", epsilon=1e-5)
    addBtagSyst(background, recoil, "TT", "tecr", tecr_ttMC, category)

    #tecr_ttTFstatParameters =  np.array([rl.NuisanceParameter("tecr_"+year+"_ttTFstat_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i, "shape") for i in range(tecr_ttTemplate[0].size)])
    tecr_ttTransferFactor = tecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
    #nominal =  tecr_ttTemplate[0] / sr_ttTemplate[0]
    #dz = simple_error_propagation(sr_ttTemplate[0], tecr_ttTemplate[0], sr_ttTemplate[3], tecr_ttTemplate[3])
    #print('tecr tt TF', dz/nominal)
    #tecr_ttTransferFactor = tecr_ttTransferFactor * ( 1. + (dz/nominal)*tecr_ttTFstatParameters )
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
    tecr_wjets.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tecr_wjets.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    tecr_wjets.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tecr_wjets.setParamEffect(wjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    tecr_wjets.setParamEffect(jec, njec, np.reciprocal(njec))
    tecr_wjets.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    tecr_wjets.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #tecr_wjets.autoMCStats(shape=True, name="tecr"+year+category+"_wjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "W+jets", "tecr", tecr_wjets, category)
    addVJetsSyst(background, recoil, "W+jets", "tecr", tecr_wjets, category)
    tecr.addSample(tecr_wjets)

    tecr_stTemplate = template(background, "ST", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, tecr_stTemplate
    )
    tecr_st.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tecr_st.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    tecr_st.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tecr_st.setParamEffect(st_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tecr_st.setParamEffect(jec, njec, np.reciprocal(njec))
    tecr_st.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    tecr_st.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #tecr_st.autoMCStats(shape=True, name="tecr"+year+category+"_stMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "ST", "tecr", tecr_st, category)
    tecr.addSample(tecr_st)

    tecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tecr_dyjetsTemplate
    )
    tecr_dyjets.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tecr_dyjets.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    tecr_dyjets.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tecr_dyjets.setParamEffect(zjetsMC_norm, nVjets_norm, np.reciprocal(nVjets_norm))
    tecr_dyjets.setParamEffect(jec, njec, np.reciprocal(njec))
    tecr_dyjets.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    tecr_dyjets.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #tecr_dyjets.autoMCStats(shape=True, name="tecr"+year+category+"_dyjetsMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "DY+jets", "tecr", tecr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "tecr", tecr_dyjets, category)
    tecr.addSample(tecr_dyjets)

    tecr_vvTemplate = template(background, "VV", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, tecr_vvTemplate
    )
    tecr_vv.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tecr_vv.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    tecr_vv.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tecr_vv.setParamEffect(vv_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tecr_vv.setParamEffect(jec, njec, np.reciprocal(njec))
    tecr_vv.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    tecr_vv.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #tecr_vv.autoMCStats(shape=True, name="tecr"+year+category+"_vvMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "VV", "tecr", tecr_vv, category)
    tecr.addSample(tecr_vv)

    tecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, tecr_hbbTemplate
    )
    tecr_hbb.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tecr_hbb.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    tecr_hbb.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tecr_hbb.setParamEffect(hbb_norm, nMinor_norm, np.reciprocal(nMinor_norm))
    tecr_hbb.setParamEffect(jec, njec, np.reciprocal(njec))
    tecr_hbb.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    tecr_hbb.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #tecr_hbb.autoMCStats(shape=True, name="tecr"+year+category+"_hbbMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "Hbb", "tecr", tecr_hbb, category)
    tecr.addSample(tecr_hbb)

    tecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tecr", category, read_sumw2=True)
    tecr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, tecr_qcdTemplate
    )
    tecr_qcd.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
    tecr_qcd.setParamEffect(trig_e, ntrig_e, np.reciprocal(ntrig_e))
    tecr_qcd.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
    tecr_qcd.setParamEffect(qcde_norm, nqcd_norm)
    tecr_qcd.setParamEffect(jec, njec, np.reciprocal(njec))
    tecr_qcd.setParamEffect(id_e, nlepton, np.reciprocal(nlepton))
    tecr_qcd.setParamEffect(reco_e, nlepton, np.reciprocal(nlepton))
    #tecr_qcd.autoMCStats(shape=True, name="tecr"+year+category+"_qcdMC", epsilon=1e-5)
    addBtagSyst(background, recoilbin, "QCD", "tecr", tecr_qcd, category)
    tecr.addSample(tecr_qcd)

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

    #if options.sumbkg:
    #    for r in data_hists["template"].identifiers("region"):
    #        if str(r) == 'sr':
    #            data[str(r)] = fake_data_hists["template"].integrate("region", r).sum("gentype")
    #        else:
    #            data[str(r)] = data_hists["template"].integrate("region", r).sum("gentype")
    #    #### Use real data in CRs while fake data in SR
    #else:
    #    for r in data_hists["template"].identifiers("region"):
    #        data[str(r)] = data_hists["template"].integrate("region", r).sum("gentype")

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
    #btag = rl.NuisanceParameter("btag" + year, "shapeN")  # AK4 btag
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
    whf_fraction = rl.NuisanceParameter("whf_fraction", "lnN")
    zhf_fraction = rl.NuisanceParameter("zhf_fraction", "lnN")
    #ghf_fraction = rl.NuisanceParameter("ghf_fraction", "shapeN")

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
    #print(recoilbins)
    #print(ptpts)
    #print(msdbins)
    #print(msdpts)
    recoilscaled = (ptpts - 250.) / (3000. - 250.)
    #    msdpts = np.sqrt(msdpts) * np.sqrt(msdpts)
    msdscaled = (msdpts - 40.) / (300.0 - 40.)
    #print(recoilscaled)
    #print(msdscaled)

    tf_dataResidualW = rl.BernsteinPoly("tf_dataResidualW"+year, (1, 1), ['recoil', 'fjmass'], limits=(-100, 100))
    tf_dataResidualW_params = tf_dataResidualW(recoilscaled, msdscaled)
    tf_dataResidualZ = rl.BernsteinPoly("tf_dataResidualZ"+year, (1, 1), ['recoil', 'fjmass'], limits=(-100, 100))
    tf_dataResidualZ_params = tf_dataResidualZ(recoilscaled, msdscaled)
    #tf_paramsZ = rhalphabeth2D("Z+jets", tf_dataResidual_params, 3, 3)
    #tf_paramsW = rhalphabeth2D("W+jets", tf_dataResidual_params, 3, 2)

    model_dict = {}
    for recoilbin in range(nrecoil):

        sr_zjetsMCFailTemplate = template(background, "Z+jets", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zjetsMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zjetsMC",
            rl.Sample.BACKGROUND,
            sr_zjetsMCFailTemplate
        )
        sr_zjetsMCFail.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
        sr_zjetsMCFail.setParamEffect(zjets_norm, nVjets_norm, np.reciprocal(nVjets_norm))
        sr_zjetsMCFail.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
        sr_zjetsMCFail.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
        sr_zjetsMCFail.setParamEffect(jec, njec, np.reciprocal(njec))
        #sr_zjetsMCFail.autoMCStats(shape=True, name="sr"+year+"fail_zjetsMC", epsilon=1e-5)
        addBtagSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCFail, "fail")
        addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCFail, "fail")

        sr_zhfMCFailTemplate = template(background, "Z+HF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zhfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zhfMC",
            rl.Sample.BACKGROUND,
            sr_zhfMCFailTemplate
        )
        sr_zhfMCFail.setParamEffect(zhf_fraction, nVjets_norm, np.reciprocal(nVjets_norm))
        #sr_zhfMCFail.autoMCStats(shape=True, name="sr"+year+"fail_zhfMC", epsilon=1e-5)

        sr_zlfMCFailTemplate = template(background, "Z+LF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zlfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zlfMC",
            rl.Sample.BACKGROUND,
            sr_zlfMCFailTemplate
        )
        sr_zlfMCFail.setParamEffect(zhf_fraction, 0.95, np.reciprocal(0.95))
        #sr_zlfMCFail.autoMCStats(shape=True, name="sr"+year+"fail_zlfMC", epsilon=1e-5)

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
        sr_wjetsMCFail.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
        sr_wjetsMCFail.setParamEffect(wjets_norm, nVjets_norm, np.reciprocal(nVjets_norm))
        sr_wjetsMCFail.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
        sr_wjetsMCFail.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
        sr_wjetsMCFail.setParamEffect(jec, njec, np.reciprocal(njec))
        #sr_wjetsMCFail.autoMCStats(shape=True, name="sr"+year+"fail_wjetsMC", epsilon=1e-5)
        addBtagSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCFail, "fail")
        addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCFail, "fail")

        sr_whfMCFailTemplate = template(background, "W+HF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_whfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_whfMC",
            rl.Sample.BACKGROUND,
            sr_whfMCFailTemplate
        )
        sr_whfMCFail.setParamEffect(whf_fraction, nVjets_norm, np.reciprocal(nVjets_norm))
        #sr_whfMCFail.autoMCStats(shape=True, name="sr"+year+"fail_whfMC", epsilon=1e-5)

        sr_wlfMCFailTemplate = template(background, "W+LF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_wlfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_wlfMC",
            rl.Sample.BACKGROUND,
            sr_wlfMCFailTemplate
        )
        sr_wlfMCFail.setParamEffect(whf_fraction, 0.9, np.reciprocal(0.9))
        #sr_wlfMCFail.autoMCStats(shape=True, name="sr"+year+"fail_wlfMC", epsilon=1e-5)

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
        sr_zjetsMCPass.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
        sr_zjetsMCPass.setParamEffect(zjets_norm, nVjets_norm, np.reciprocal(nVjets_norm))
        sr_zjetsMCPass.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
        sr_zjetsMCPass.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
        sr_zjetsMCPass.setParamEffect(jec, njec, np.reciprocal(njec))
        #sr_zjetsMCPass.autoMCStats(shape=True, name="sr"+year+"pass_zjetsMC", epsilon=1e-5)
        addBtagSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCPass, "pass")
        addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCPass, "pass")

        sr_zhfMCPassTemplate = template(background, "Z+HF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_zhfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_zhfMC",
            rl.Sample.BACKGROUND,
            sr_zhfMCPassTemplate
        )
        sr_zhfMCPass.setParamEffect(zhf_fraction, nVjets_norm, np.reciprocal(nVjets_norm))
        #sr_zhfMCPass.autoMCStats(shape=True, name="sr"+year+"pass_zhfMC", epsilon=1e-5)

        sr_zlfMCPassTemplate = template(background, "Z+LF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_zlfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_zlfMC",
            rl.Sample.BACKGROUND,
            sr_zlfMCPassTemplate
        )
        sr_zlfMCPass.setParamEffect(zhf_fraction, 0.95, np.reciprocal(0.95))
        #sr_zlfMCPass.autoMCStats(shape=True, name="sr"+year+"pass_zlfMC", epsilon=1e-5)

        tf_paramsZdeco = (sr_zlfMCPass.getExpectation()+sr_zhfMCPass.getExpectation()) / (sr_zlfMCFail.getExpectation()+sr_zhfMCFail.getExpectation())
        tf_paramsZ = tf_paramsZdeco * tf_dataResidualZ_params[recoilbin, :]

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
        sr_wjetsMCPass.setParamEffect(lumi, nlumi, np.reciprocal(nlumi))
        sr_wjetsMCPass.setParamEffect(wjets_norm, nVjets_norm, np.reciprocal(nVjets_norm))
        sr_wjetsMCPass.setParamEffect(trig_met, ntrig_met, np.reciprocal(ntrig_met))
        sr_wjetsMCPass.setParamEffect(veto_tau, nveto_tau, np.reciprocal(nveto_tau))
        sr_wjetsMCPass.setParamEffect(jec, njec, np.reciprocal(njec))
        #sr_wjetsMCPass.autoMCStats(shape=True, name="sr"+year+"pass_wjetsMC", epsilon=1e-5)
        addBtagSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCPass, "pass")
        addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCPass, "pass")

        sr_whfMCPassTemplate = template(background, "W+HF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_whfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_whfMC",
            rl.Sample.BACKGROUND,
            sr_whfMCPassTemplate
        )
        sr_whfMCPass.setParamEffect(whf_fraction, nVjets_norm, np.reciprocal(nVjets_norm))
        #sr_whfMCPass.autoMCStats(shape=True, name="sr"+year+"pass_whfMC", epsilon=1e-5)

        sr_wlfMCPassTemplate = template(background, "W+LF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_wlfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_wlfMC",
            rl.Sample.BACKGROUND,
            sr_wlfMCPassTemplate
        )
        sr_wlfMCPass.setParamEffect(whf_fraction, 0.9, np.reciprocal(0.9))
        #sr_wlfMCPass.autoMCStats(shape=True, name="sr"+year+"pass_wlfMC", epsilon=1e-5)

        tf_paramsWdeco = (sr_wlfMCPass.getExpectation()+sr_whfMCPass.getExpectation()) / (sr_wlfMCFail.getExpectation()+sr_whfMCFail.getExpectation())
        tf_paramsW = tf_paramsWdeco * tf_dataResidualW_params[recoilbin, :]

        sr_wjetsPass = rl.TransferFactorSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin)+ "_wjets",
            rl.Sample.BACKGROUND,
            tf_paramsW,
            sr_wjetsFail
        )

        for s in signal["sr"].identifiers("process"):
            #if "Mz500_mhs90_Mdm250" not in str(s):
            #    continue
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
