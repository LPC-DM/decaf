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
    #0,
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
#recoil_binning = [250, 310, 370, 470, 590, 3000]
recoil_binning_dict = {
    "2018": [250, 310, 370, 470, 3000],
    "2017": [250, 310, 370, 470, 3000],
    "2016": [250, 310, 370, 470, 3000]
}

category_map = {"pass": 1, "fail": 0}

with open("data/hf_systematic.json") as fin:
    hf_systematic = json.load(fin)

def template(dictionary, process, systematic, recoil, region, category, read_sumw2=False):
    histogram = dictionary[region].integrate("process", process)
    nominal, sumw2 = histogram.integrate("systematic", "nominal").values(sumw2=True)[()]
    nominal=nominal[recoil, :, category_map[category]]
    sumw2=sumw2[recoil, :, category_map[category]]
    zerobins = nominal <= 0.
    output = nominal
    if "data" not in systematic:
        output[zerobins] = 1.
        sumw2[zerobins] = 0.
    if "nominal" not in systematic and "data" not in systematic:
        #print('Normalizing',systematic,'histogram of',process,'in region',region)
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

def addSyst(dictionary, recoil, process, region, templ, category, syst, string):
    histogram = dictionary[region].integrate("process", process)
    nominal=histogram.integrate("systematic", "nominal").values()[()][recoil, :, category_map[category]]
    up=histogram.integrate("systematic", string+"Up").values()[()][recoil, :, category_map[category]]
    down=histogram.integrate("systematic",string+"Down").values()[()][recoil, :, category_map[category]]
    systUp = np.array( up.sum() / nominal.sum() )
    systUp[np.isnan(systUp)] = 1.
    systUp = systUp.sum()
    templ.setParamEffect(syst, systUp)
    #print(string, systUp)
def addBtagSyst(dictionary, recoil, process, region, templ, category):
    #btagUp = template(dictionary, process, "btagUp", recoil, region, category)[0]
    #btagDown = template(dictionary, process, "btagDown", recoil, region, category)[0]
    #btagDown = np.array( down.sum() / nominal.sum() )
    #btagDown[np.isnan(btagDown)] = 1.
    #btagDown = btagDown.sum()
    #print(process, btagUp, btagDown)
    #templ.setParamEffect(btag, btagUp, btagDown)
    addSyst(dictionary, recoil, process, region, templ, category, btag, "btag")

def addVJetsSyst(dictionary, recoil, process, region, templ, category):
    addSyst(dictionary, recoil, process, region, templ, category, ew1, "ew1")
    addSyst(dictionary, recoil, process, region, templ, category, ew2W, "ew2W")
    addSyst(dictionary, recoil, process, region, templ, category, ew2Z, "ew2Z")
    addSyst(dictionary, recoil, process, region, templ, category, ew3W, "ew3W")
    addSyst(dictionary, recoil, process, region, templ, category, ew3Z, "ew3Z")
    addSyst(dictionary, recoil, process, region, templ, category, mix, "mix")
    addSyst(dictionary, recoil, process, region, templ, category, qcd1, "qcd1")
    addSyst(dictionary, recoil, process, region, templ, category, qcd2, "qcd2")
    addSyst(dictionary, recoil, process, region, templ, category, qcd3, "qcd3")
'''
def addVJetsSyst(dictionary, recoil, process, region, templ, category):
    ew1Up = template(dictionary, process, "ew1Up", recoil, region, category)[0]
    ew1Down = template(dictionary, process, "ew1Down", recoil, region, category)[0]
    templ.setParamEffect(ew1, ew1Up, ew1Down)
    print("ew1",ew1Up)
    ew2GUp = template(dictionary, process, "ew2GUp", recoil, region, category)[0]
    ew2GDown = template(dictionary, process, "ew2GDown", recoil, region, category)[0]
    templ.setParamEffect(ew2G, ew2GUp, ew2GDown)
    
    ew2WUp = template(dictionary, process, "ew2WUp", recoil, region, category)[0]
    ew2WDown = template(dictionary, process, "ew2WDown", recoil, region, category)[0]
    templ.setParamEffect(ew2W, ew2WUp, ew2WDown)
    print("ew2W",ew2WUp)
    ew2ZUp = template(dictionary, process, "ew2ZUp", recoil, region, category)[0]
    ew2ZDown = template(dictionary, process, "ew2ZDown", recoil, region, category)[0]
    templ.setParamEffect(ew2Z, ew2ZUp, ew2ZDown)
    print("ew2Z",ew2ZUp)
    ew3GUp = template(dictionary, process, "ew3GUp", recoil, region, category)[0]
    ew3GDown = template(dictionary, process, "ew3GDown", recoil, region, category)[0]
    templ.setParamEffect(ew3G, ew3GUp, ew3GDown)
    ew3WUp = template(dictionary, process, "ew3WUp", recoil, region, category)[0]
    ew3WDown = template(dictionary, process, "ew3WDown", recoil, region, category)[0]
    templ.setParamEffect(ew3W, ew3WUp, ew3WDown)
    print("ew3W",ew3WUp)
    ew3ZUp = template(dictionary, process, "ew3ZUp", recoil, region, category)[0]
    ew3ZDown = template(dictionary, process, "ew3ZDown", recoil, region, category)[0]
    templ.setParamEffect(ew3Z, ew3ZUp, ew3ZDown)
    print("ew3Z",ew3ZUp)
    mixUp = template(dictionary, process, "mixUp", recoil, region, category)[0]
    mixDown = template(dictionary, process, "mixDown", recoil, region, category)[0]
    templ.setParamEffect(mix, mixUp, mixDown)
    print("mix",mixUp)
    muFUp = template(dictionary, process, "muFUp", recoil, region, category)[0]
    muFDown = template(dictionary, process, "muFDown", recoil, region, category)[0]
    #templ.setParamEffect(muF, muFUp, muFDown)
    muRUp = template(dictionary, process, "muRUp", recoil, region, category)[0]
    muRDown = template(dictionary, process, "muRDown", recoil, region, category)[0]
    #templ.setParamEffect(muR, muRUp, muRDown)
    qcd1Up = template(dictionary, process, "qcd1Up", recoil, region, category)[0]
    qcd1Down = template(dictionary, process, "qcd1Down", recoil, region, category)[0]
    templ.setParamEffect(qcd1, qcd1Up, qcd1Down)
    print("qcd1",qcd1Up)
    qcd2Up = template(dictionary, process, "qcd2Up", recoil, region, category)[0]
    qcd2Down = template(dictionary, process, "qcd2Down", recoil, region, category)[0]
    templ.setParamEffect(qcd2, qcd2Up, qcd2Down)
    print("qcd2",qcd2Up)
    qcd3Up = template(dictionary, process, "qcd3Up", recoil, region, category)[0]
    qcd3Down = template(dictionary, process, "qcd3Down", recoil, region, category)[0]
    templ.setParamEffect(qcd3, qcd3Up, qcd3Down)
    print("qcd3",qcd3Up)
'''
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

def model(year, recoil, category):

    model_id = year + category + "recoil" + str(recoil)
    print(model_id)
    model = rl.Model("darkhiggs" + model_id)

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
        sr_wjetsTemplate = sr_wjetsMCPassTemplate
        sr_wjetsMC = sr_wjetsMCPass
        sr_wjets = sr_wjetsPass
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
    sr_ttMC.setParamEffect(lumi, 1.027)
    sr_ttMC.setParamEffect(tt_norm, 1.2)
    sr_ttMC.setParamEffect(trig_met, 1.02)
    sr_ttMC.setParamEffect(veto_tau, 1.03)
    sr_ttMC.setParamEffect(jec, 1.05)
    addBtagSyst(background, recoil, "TT", "sr", sr_ttMC, category)

    if category == "pass":
        sr_ttMC.autoMCStats()
        sr_ttObservable = rl.Observable("fjmass", sr_ttTemplate[1])
        sr_ttParameters = np.array(
            [
                rl.IndependentParameter(                                                                                                                                     
                    "sr" + year + "_tt_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i,
                    0,
                    #-100.,
                    #100.,
                )
                for i in range(sr_ttObservable.nbins)
            ]
        )
        sr_ttBinYields = sr_ttTemplate[0] * (1 + (10./np.maximum(1., np.sqrt(sr_ttTemplate[0]))))**sr_ttParameters

        '''
        sr_ttBinYields = np.array(  # one nuisance per mass shape bin in pass                                              
            [
                rl.IndependentParameter(
                    "sr" + year + "_tt_" + category + "_recoil"+str(recoilbin)+"_mass%d" % i,
                    b,
                    0,
                    sr_ttTemplate[0].max() * 10,
                )
                for i, b in enumerate(sr_ttTemplate[0])
            ]
        )
        '''
        sr_tt = rl.ParametericSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields
        )
        sr.addSample(sr_tt)
    else:
        sr.addSample(sr_ttMC)

    ###
    # Other MC-driven processes
    ###

    sr_stTemplate = template(background, "ST", "nominal", recoil, "sr", category)
    sr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, sr_stTemplate)
    sr_st.setParamEffect(lumi, 1.027)
    sr_st.setParamEffect(trig_met, 1.02)
    sr_st.setParamEffect(veto_tau, 1.03)
    sr_st.setParamEffect(st_norm, 1.2)
    sr_st.setParamEffect(jec, 1.05)
    addBtagSyst(background, recoil, "ST", "sr", sr_st, category)
    sr.addSample(sr_st)

    sr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "sr", category)
    sr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, sr_dyjetsTemplate
    )
    sr_dyjets.setParamEffect(lumi, 1.027)
    sr_dyjets.setParamEffect(trig_met, 1.02)
    sr_dyjets.setParamEffect(veto_tau, 1.03)
    sr_dyjets.setParamEffect(zjets_norm, 1.4)
    sr_dyjets.setParamEffect(jec, 1.05)
    addBtagSyst(background, recoil, "DY+jets", "sr", sr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "sr", sr_dyjets, category)
    sr.addSample(sr_dyjets)

    sr_vvTemplate = template(background, "VV", "nominal", recoil, "sr", category)
    sr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, sr_vvTemplate)
    sr_vv.setParamEffect(lumi, 1.027)
    sr_vv.setParamEffect(trig_met, 1.02)
    sr_vv.setParamEffect(veto_tau, 1.03)
    sr_vv.setParamEffect(vv_norm, 1.2)
    sr_vv.setParamEffect(jec, 1.05)
    addBtagSyst(background, recoil, "VV", "sr", sr_vv, category)
    sr.addSample(sr_vv)

    sr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "sr", category)
    sr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, sr_hbbTemplate)
    sr_hbb.setParamEffect(lumi, 1.027)
    sr_hbb.setParamEffect(trig_met, 1.02)
    sr_hbb.setParamEffect(veto_tau, 1.03)
    sr_hbb.setParamEffect(hbb_norm, 1.2)
    sr_hbb.setParamEffect(jec, 1.05)
    addBtagSyst(background, recoil, "Hbb", "sr", sr_hbb, category)
    sr.addSample(sr_hbb)

    sr_qcdTemplate = template(background, "QCD", "nominal", recoil, "sr", category)
    sr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, sr_qcdTemplate)
    sr_qcd.setParamEffect(lumi, 1.027)
    sr_qcd.setParamEffect(trig_met, 1.02)
    sr_qcd.setParamEffect(veto_tau, 1.03)
    sr_qcd.setParamEffect(qcdsig_norm, 2.0)
    sr_qcd.setParamEffect(jec, 1.05)
    addBtagSyst(background, recoil, "QCD", "sr", sr_qcd, category)
    sr.addSample(sr_qcd)

    for s in signal["sr"].identifiers("process"):
        if "Mz500_mhs90_Mdm250" not in str(s):
            continue
        print("Signal is:", str(s))
        sr_signalTemplate = template(signal, s, "nominal", recoil, "sr", category)
        sr_signal = rl.TemplateSample(
            ch_name + "_" + str(s), rl.Sample.SIGNAL, sr_signalTemplate
        )
        sr_signal.setParamEffect(lumi, 1.027)
        sr_signal.setParamEffect(trig_met, 1.02)
        sr_signal.setParamEffect(veto_tau, 1.03)
        sr_signal.setParamEffect(jec, 1.05)
        addBtagSyst(signal, recoil, str(s), "sr", sr_signal, category)
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

    wmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wmcr", category, read_sumw2=True)
    wmcr_wjetsMC = rl.TemplateSample(
        "wmcr" + model_id + "_wjetsMC",
        rl.Sample.BACKGROUND,
        wmcr_wjetsTemplate
    )
    wmcr_wjetsMC.setParamEffect(lumi, 1.027)
    wmcr_wjetsMC.setParamEffect(trig_met, 1.02)
    wmcr_wjetsMC.setParamEffect(veto_tau, 1.03)
    wmcr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    wmcr_wjetsMC.setParamEffect(jec, 1.05)
    wmcr_wjetsMC.setParamEffect(id_mu, 1.02)
    wmcr_wjetsMC.setParamEffect(iso_mu, 1.02)
    wmcr_wjetsMC.autoMCStats()
    addBtagSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)
    addVJetsSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)

    wmcr_wjetsTransferFactor = wmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
    #wmcr_wjetsTransferFactor = wmcr_wjetsTemplate[0] / sr_wjetsTemplate[0]
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
    wmcr_ttMC.setParamEffect(lumi, 1.027)
    wmcr_ttMC.setParamEffect(trig_met, 1.02)
    wmcr_ttMC.setParamEffect(veto_tau, 1.03)
    wmcr_ttMC.setParamEffect(tt_norm, 1.2)
    wmcr_ttMC.setParamEffect(jec, 1.05)
    wmcr_ttMC.setParamEffect(id_mu, 1.02)
    wmcr_ttMC.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoil, "TT", "wmcr", wmcr_ttMC, category)
    
    if category == "pass":
        wmcr_ttMC.autoMCStats()
        wmcr_ttTransferFactor = wmcr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        wmcr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt
        )
        wmcr.addSample(wmcr_tt)
    else:
        wmcr.addSample(wmcr_ttMC)
        
    ###
    # Other MC-driven processes
    ###

    wmcr_stTemplate = template(background, "ST", "nominal", recoil, "wmcr", category)
    wmcr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, wmcr_stTemplate
        )
    wmcr_st.setParamEffect(lumi, 1.027)
    wmcr_st.setParamEffect(trig_met, 1.02)
    wmcr_st.setParamEffect(veto_tau, 1.03)
    wmcr_st.setParamEffect(st_norm, 1.2)
    wmcr_st.setParamEffect(jec, 1.05)
    wmcr_st.setParamEffect(id_mu, 1.02)
    wmcr_st.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "ST", "wmcr", wmcr_st, category)
    wmcr.addSample(wmcr_st)

    wmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wmcr", category)
    wmcr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wmcr_dyjetsTemplate
    )
    wmcr_dyjets.setParamEffect(lumi, 1.027)
    wmcr_dyjets.setParamEffect(trig_met, 1.02)
    wmcr_dyjets.setParamEffect(veto_tau, 1.03)
    wmcr_dyjets.setParamEffect(zjets_norm, 1.4)
    wmcr_dyjets.setParamEffect(jec, 1.05)
    wmcr_dyjets.setParamEffect(id_mu, 1.02)
    wmcr_dyjets.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "DY+jets", "wmcr", wmcr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "wmcr", wmcr_dyjets, category)
    wmcr.addSample(wmcr_dyjets)

    wmcr_vvTemplate = template(background, "VV", "nominal", recoil, "wmcr", category)
    wmcr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, wmcr_vvTemplate
    )
    wmcr_vv.setParamEffect(lumi, 1.027)
    wmcr_vv.setParamEffect(trig_met, 1.02)
    wmcr_vv.setParamEffect(veto_tau, 1.03)
    wmcr_vv.setParamEffect(vv_norm, 1.2)
    wmcr_vv.setParamEffect(jec, 1.05)
    wmcr_vv.setParamEffect(id_mu, 1.02)
    wmcr_vv.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "VV", "wmcr", wmcr_vv, category)
    wmcr.addSample(wmcr_vv)

    wmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wmcr", category)
    wmcr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, wmcr_hbbTemplate
    )
    wmcr_hbb.setParamEffect(lumi, 1.027)
    wmcr_hbb.setParamEffect(trig_met, 1.02)
    wmcr_hbb.setParamEffect(veto_tau, 1.03)
    wmcr_hbb.setParamEffect(hbb_norm, 1.2)
    wmcr_hbb.setParamEffect(jec, 1.05)
    wmcr_hbb.setParamEffect(id_mu, 1.02)
    wmcr_hbb.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "Hbb", "wmcr", wmcr_hbb, category)
    wmcr.addSample(wmcr_hbb)

    wmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wmcr", category)
    wmcr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, wmcr_qcdTemplate
    )
    wmcr_qcd.setParamEffect(lumi, 1.027)
    wmcr_qcd.setParamEffect(trig_met, 1.02)
    wmcr_qcd.setParamEffect(veto_tau, 1.03)
    wmcr_qcd.setParamEffect(qcdmu_norm, 2.0)
    wmcr_qcd.setParamEffect(jec, 1.05)
    wmcr_qcd.setParamEffect(id_mu, 1.02)
    wmcr_qcd.setParamEffect(iso_mu, 1.02)
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
    wecr_wjetsMC.setParamEffect(lumi, 1.027)
    wecr_wjetsMC.setParamEffect(trig_e, 1.01)
    wecr_wjetsMC.setParamEffect(veto_tau, 1.03)
    wecr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    wecr_wjetsMC.setParamEffect(jec, 1.05)
    wecr_wjetsMC.setParamEffect(id_e, 1.02)
    wecr_wjetsMC.setParamEffect(reco_e, 1.02)
    wecr_wjetsMC.autoMCStats()
    addBtagSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)
    addVJetsSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)

    wecr_wjetsTransferFactor = wecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation()
    #wecr_wjetsTransferFactor = wecr_wjetsTemplate[0] / sr_wjetsTemplate[0]
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
    wecr_ttMC.setParamEffect(lumi, 1.027)
    wecr_ttMC.setParamEffect(trig_e, 1.01)
    wecr_ttMC.setParamEffect(veto_tau, 1.03)
    wecr_ttMC.setParamEffect(tt_norm, 1.2)
    wecr_ttMC.setParamEffect(jec, 1.05)
    wecr_ttMC.setParamEffect(id_e, 1.02)
    wecr_ttMC.setParamEffect(reco_e, 1.02)
    addBtagSyst(background, recoil, "TT", "wecr", wecr_ttMC, category)

    if category == "pass":
        wecr_ttMC.autoMCStats()
        wecr_ttTransferFactor = wecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
        wecr_tt = rl.TransferFactorSample(
            ch_name + "_tt", rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt
        )
        wecr.addSample(wecr_tt)
    else:
        wecr.addSample(wecr_ttMC)

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
    addBtagSyst(background, recoilbin, "ST", "wecr", wecr_st, category)
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
    addBtagSyst(background, recoilbin, "DY+jets", "wecr", wecr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "wecr", wecr_dyjets, category)
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
    addBtagSyst(background, recoilbin, "VV", "wecr", wecr_vv, category)
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
    addBtagSyst(background, recoilbin, "Hbb", "wecr", wecr_hbb, category)
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
    addBtagSyst(background, recoilbin, "QCD", "wecr", wecr_qcd, category)
    wecr.addSample(wecr_qcd)

    ###
    # End of single electron W control region
    ###

    if not (category == "pass"):
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
    tmcr_ttMC.setParamEffect(lumi, 1.027)
    tmcr_ttMC.setParamEffect(trig_met, 1.02)
    tmcr_ttMC.setParamEffect(veto_tau, 1.03)
    tmcr_ttMC.setParamEffect(tt_norm, 1.2)
    tmcr_ttMC.setParamEffect(jec, 1.05)
    tmcr_ttMC.setParamEffect(id_mu, 1.02)
    tmcr_ttMC.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoil, "TT", "tmcr", tmcr_ttMC, category)
    tmcr_ttMC.autoMCStats()
    tmcr_ttTransferFactor = tmcr_ttMC.getExpectation() / sr_ttMC.getExpectation() 
    tmcr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt
    )
    tmcr.addSample(tmcr_tt)    
    
    ###
    # Other MC-driven processes
    ###
    
    tmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tmcr", category)
    tmcr_wjets = rl.TemplateSample(
        ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tmcr_wjetsTemplate
    )
    tmcr_wjets.setParamEffect(lumi, 1.027)
    tmcr_wjets.setParamEffect(trig_met, 1.02)
    tmcr_wjets.setParamEffect(veto_tau, 1.03)
    tmcr_wjets.setParamEffect(zjets_norm, 1.4)
    tmcr_wjets.setParamEffect(jec, 1.05)
    tmcr_wjets.setParamEffect(id_mu, 1.02)
    tmcr_wjets.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "W+jets", "tmcr", tmcr_wjets, category)
    addVJetsSyst(background, recoil, "W+jets", "tmcr", tmcr_wjets, category)
    tmcr.addSample(tmcr_wjets)
    
    tmcr_stTemplate = template(background, "ST", "nominal", recoil, "tmcr", category)
    tmcr_st = rl.TemplateSample(
        ch_name + "_stMC", rl.Sample.BACKGROUND, tmcr_stTemplate
    )
    tmcr_st.setParamEffect(lumi, 1.027)
    tmcr_st.setParamEffect(trig_met, 1.02)
    tmcr_st.setParamEffect(veto_tau, 1.03)
    tmcr_st.setParamEffect(st_norm, 1.2)
    tmcr_st.setParamEffect(jec, 1.05)
    tmcr_st.setParamEffect(id_mu, 1.02)
    tmcr_st.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "ST", "tmcr", tmcr_st, category)
    tmcr.addSample(tmcr_st)

    tmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tmcr", category)
    tmcr_dyjets = rl.TemplateSample(
        ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tmcr_dyjetsTemplate
    )
    tmcr_dyjets.setParamEffect(lumi, 1.027)
    tmcr_dyjets.setParamEffect(trig_met, 1.02)
    tmcr_dyjets.setParamEffect(veto_tau, 1.03)
    tmcr_dyjets.setParamEffect(zjets_norm, 1.4)
    tmcr_dyjets.setParamEffect(jec, 1.05)
    tmcr_dyjets.setParamEffect(id_mu, 1.02)
    tmcr_dyjets.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "DY+jets", "tmcr", tmcr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "tmcr", tmcr_dyjets, category)
    tmcr.addSample(tmcr_dyjets)

    tmcr_vvTemplate = template(background, "VV", "nominal", recoil, "tmcr", category)
    tmcr_vv = rl.TemplateSample(
        ch_name + "_vvMC", rl.Sample.BACKGROUND, tmcr_vvTemplate
    )
    tmcr_vv.setParamEffect(lumi, 1.027)
    tmcr_vv.setParamEffect(trig_met, 1.02)
    tmcr_vv.setParamEffect(veto_tau, 1.03)
    tmcr_vv.setParamEffect(vv_norm, 1.2)
    tmcr_vv.setParamEffect(jec, 1.05)
    tmcr_vv.setParamEffect(id_mu, 1.02)
    tmcr_vv.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "VV", "tmcr", tmcr_vv, category)
    tmcr.addSample(tmcr_vv)
    
    tmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tmcr", category)
    tmcr_hbb = rl.TemplateSample(
        ch_name + "_hbbMC", rl.Sample.BACKGROUND, tmcr_hbbTemplate
    )
    tmcr_hbb.setParamEffect(lumi, 1.027)
    tmcr_hbb.setParamEffect(trig_met, 1.02)
    tmcr_hbb.setParamEffect(veto_tau, 1.03)
    tmcr_hbb.setParamEffect(hbb_norm, 1.2)
    tmcr_hbb.setParamEffect(jec, 1.05)
    tmcr_hbb.setParamEffect(id_mu, 1.02)
    tmcr_hbb.setParamEffect(iso_mu, 1.02)
    addBtagSyst(background, recoilbin, "Hbb", "tmcr", tmcr_hbb, category)
    tmcr.addSample(tmcr_hbb)

    tmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tmcr", category)
    tmcr_qcd = rl.TemplateSample(
        ch_name + "_qcdMC", rl.Sample.BACKGROUND, tmcr_qcdTemplate
    )
    tmcr_qcd.setParamEffect(lumi, 1.027)
    tmcr_qcd.setParamEffect(trig_met, 1.02)
    tmcr_qcd.setParamEffect(veto_tau, 1.03)
    tmcr_qcd.setParamEffect(qcdmu_norm, 2.0)
    tmcr_qcd.setParamEffect(jec, 1.05)
    tmcr_qcd.setParamEffect(id_mu, 1.02)
    tmcr_qcd.setParamEffect(iso_mu, 1.02)
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
    tecr_ttMC.setParamEffect(lumi, 1.027)
    tecr_ttMC.setParamEffect(trig_e, 1.01)
    tecr_ttMC.setParamEffect(veto_tau, 1.03)
    tecr_ttMC.setParamEffect(tt_norm, 1.2)
    tecr_ttMC.setParamEffect(jec, 1.05)
    tecr_ttMC.setParamEffect(id_e, 1.02)
    tecr_ttMC.setParamEffect(reco_e, 1.02)
    addBtagSyst(background, recoil, "TT", "tecr", tecr_ttMC, category)
    tecr_ttMC.autoMCStats()
    tecr_ttTransferFactor = tecr_ttMC.getExpectation() / sr_ttMC.getExpectation()
    tecr_tt = rl.TransferFactorSample(
        ch_name + "_tt", rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt
    )
    tecr.addSample(tecr_tt)
    
    ###
    # Other MC-driven processes
    ###

    tecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tecr", category)
    tecr_wjets = rl.TemplateSample(
        ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tecr_wjetsTemplate
    )
    tecr_wjets.setParamEffect(lumi, 1.027)
    tecr_wjets.setParamEffect(trig_e, 1.01)
    tecr_wjets.setParamEffect(veto_tau, 1.03)
    tecr_wjets.setParamEffect(zjets_norm, 1.4)
    tecr_wjets.setParamEffect(jec, 1.05)
    tecr_wjets.setParamEffect(id_e, 1.02)
    tecr_wjets.setParamEffect(reco_e, 1.02)
    addBtagSyst(background, recoilbin, "W+jets", "tecr", tecr_wjets, category)
    addVJetsSyst(background, recoil, "W+jets", "tecr", tecr_wjets, category)
    tecr.addSample(tecr_wjets)

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
    addBtagSyst(background, recoilbin, "ST", "tecr", tecr_st, category)
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
    addBtagSyst(background, recoilbin, "DY+jets", "tecr", tecr_dyjets, category)
    addVJetsSyst(background, recoil, "DY+jets", "tecr", tecr_dyjets, category)
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
    addBtagSyst(background, recoilbin, "VV", "tecr", tecr_vv, category)
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
    addBtagSyst(background, recoilbin, "Hbb", "tecr", tecr_hbb, category)
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
    hists = remap_histograms(hists)
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
    btag = rl.NuisanceParameter("btag" + year, "lnN")  # AK4 btag
    ew1 = rl.NuisanceParameter("ew1", "lnN")
    #ew2G = rl.NuisanceParameter("ew2G", "lnN")
    ew2W = rl.NuisanceParameter("ew2W", "lnN")
    ew2Z = rl.NuisanceParameter("ew2Z", "lnN")
    #ew3G = rl.NuisanceParameter("ew3G", "lnN")
    ew3W = rl.NuisanceParameter("ew3W", "lnN")
    ew3Z = rl.NuisanceParameter("ew3Z", "lnN")
    mix = rl.NuisanceParameter("mix", "lnN")
    muF = rl.NuisanceParameter("muF", "lnN")
    muR = rl.NuisanceParameter("muR", "lnN")
    qcd1 = rl.NuisanceParameter("qcd1", "lnN")
    qcd2 = rl.NuisanceParameter("qcd2", "lnN")
    qcd3 = rl.NuisanceParameter("qcd3", "lnN")
    whf_fraction = rl.NuisanceParameter("whf_fraction", "lnN")
    zhf_fraction = rl.NuisanceParameter("zhf_fraction", "lnN")
    #ghf_fraction = rl.NuisanceParameter("ghf_fraction", "shapeN")

    ###
    # Preparing Rhalphabet
    ###

    msdbins = np.array(mass_binning)
    msd = rl.Observable('fjmass', msdbins)
    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(recoilbins[:-1] + 0.3 * np.diff(recoilbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    print(recoilbins)
    print(ptpts)
    print(msdbins)
    print(msdpts)
    recoilscaled = (ptpts - 250.) / (3000. - 250.)
    #    msdpts = np.sqrt(msdpts) * np.sqrt(msdpts)
    msdscaled = (msdpts - 40.) / (300.0 - 40.)
    print(recoilscaled)
    print(msdscaled)
    
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
        sr_zjetsMCFail.setParamEffect(lumi, 1.027)
        sr_zjetsMCFail.setParamEffect(zjets_norm, 1.4)
        sr_zjetsMCFail.setParamEffect(trig_met, 1.02)
        sr_zjetsMCFail.setParamEffect(veto_tau, 1.03)
        sr_zjetsMCFail.setParamEffect(jec, 1.05)
        sr_zjetsMCFail.autoMCStats()
        addBtagSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCFail, "fail")
        addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCFail, "fail")

        sr_zhfMCFailTemplate = template(background, "Z+HF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zhfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zhfMC",
            rl.Sample.BACKGROUND,
            sr_zhfMCFailTemplate
        )
        sr_zhfMCFail.setParamEffect(zhf_fraction, 1.4)
        sr_zhfMCFail.autoMCStats()

        sr_zlfMCFailTemplate = template(background, "Z+LF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_zlfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_zlfMC",
            rl.Sample.BACKGROUND,
            sr_zlfMCFailTemplate
        )
        sr_zlfMCFail.autoMCStats()

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
        sr_wjetsMCFail.setParamEffect(lumi, 1.027)
        sr_wjetsMCFail.setParamEffect(wjets_norm, 1.4)
        sr_wjetsMCFail.setParamEffect(trig_met, 1.02)
        sr_wjetsMCFail.setParamEffect(veto_tau, 1.03)
        sr_wjetsMCFail.setParamEffect(jec, 1.05)
        sr_wjetsMCFail.autoMCStats()
        addBtagSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCFail, "fail")
        addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCFail, "fail")

        sr_whfMCFailTemplate = template(background, "W+HF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_whfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_whfMC",
            rl.Sample.BACKGROUND,
            sr_whfMCFailTemplate
        )
        sr_whfMCFail.setParamEffect(whf_fraction, 1.4)
        sr_whfMCFail.autoMCStats()

        sr_wlfMCFailTemplate = template(background, "W+LF", "nominal", recoilbin, "sr", "fail", read_sumw2=True)
        sr_wlfMCFail = rl.TemplateSample(
            "sr" + year + "fail" + "recoil" + str(recoilbin) + "_wlfMC",
            rl.Sample.BACKGROUND,
            sr_wlfMCFailTemplate
        )
        sr_wlfMCFail.autoMCStats()

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
        sr_zjetsMCPass.setParamEffect(lumi, 1.027)
        sr_zjetsMCPass.setParamEffect(zjets_norm, 1.4)
        sr_zjetsMCPass.setParamEffect(trig_met, 1.02)
        sr_zjetsMCPass.setParamEffect(veto_tau, 1.03)
        sr_zjetsMCPass.setParamEffect(jec, 1.05)
        sr_zjetsMCPass.autoMCStats()
        addBtagSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCPass, "pass")
        addVJetsSyst(background, recoilbin, "Z+jets", "sr", sr_zjetsMCPass, "pass")

        sr_zhfMCPassTemplate = template(background, "Z+HF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_zhfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_zhfMC",
            rl.Sample.BACKGROUND,
            sr_zhfMCPassTemplate
        )
        sr_zhfMCPass.setParamEffect(zhf_fraction, 1.4)
        sr_zhfMCPass.autoMCStats()

        sr_zlfMCPassTemplate = template(background, "Z+LF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_zlfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_zlfMC",
            rl.Sample.BACKGROUND,
            sr_zlfMCPassTemplate
        )
        sr_zlfMCPass.autoMCStats()

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
        sr_wjetsMCPass.setParamEffect(lumi, 1.027)
        sr_wjetsMCPass.setParamEffect(wjets_norm, 1.4)
        sr_wjetsMCPass.setParamEffect(trig_met, 1.02)
        sr_wjetsMCPass.setParamEffect(veto_tau, 1.03)
        sr_wjetsMCPass.setParamEffect(jec, 1.05)
        sr_wjetsMCPass.autoMCStats()
        addBtagSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCPass, "pass")
        addVJetsSyst(background, recoilbin, "W+jets", "sr", sr_wjetsMCPass, "pass")

        sr_whfMCPassTemplate = template(background, "W+HF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_whfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_whfMC",
            rl.Sample.BACKGROUND,
            sr_whfMCPassTemplate
        )
        sr_whfMCPass.setParamEffect(whf_fraction, 1.4)
        sr_whfMCPass.autoMCStats()

        sr_wlfMCPassTemplate = template(background, "W+LF", "nominal", recoilbin, "sr", "pass", read_sumw2=True)
        sr_wlfMCPass = rl.TemplateSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin) + "_wlfMC",
            rl.Sample.BACKGROUND,
            sr_wlfMCPassTemplate
        )
        sr_wlfMCPass.autoMCStats()

        tf_paramsWdeco = (sr_wlfMCPass.getExpectation()+sr_whfMCPass.getExpectation()) / (sr_wlfMCFail.getExpectation()+sr_whfMCFail.getExpectation())
        tf_paramsW = tf_paramsWdeco * tf_dataResidualW_params[recoilbin, :]
    
        sr_wjetsPass = rl.TransferFactorSample(
            "sr" + year + "pass" + "recoil" + str(recoilbin)+ "_wjets",
            rl.Sample.BACKGROUND,
            tf_paramsW,
            sr_wjetsFail
        )

        for category in ["pass", "fail"]:
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
