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

mass_binning = [40, 50, 60, 70, 80, 90, 100, 120, 150, 180, 240, 300,]

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

if __name__ == "__main__":
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

    ###
    # Remapping histograms 
    ###

    hists = remap_histograms(hists)
    bkg_hists = hists["bkg"]

    background = {}
    for r in bkg_hists["template"].identifiers("region"):
        background[str(r)] = bkg_hists["template"].integrate("region", r).sum("gentype")

    ###
    # Preparing Rhalphabet
    ###

    msdbins = np.array(mass_binning)
    msd = rl.Observable('fjmass', msdbins)
    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(recoilbins[:-1] + 0.3 * np.diff(recoilbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    recoilscaled = (ptpts - 250.) / (3000. - 250.)
    msdscaled = (msdpts - 40.) / (300.0 - 40.)

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

        return qcdpass / qcdfail

    zjetspass_templ = []
    zjetsfail_templ = []
    for recoilbin in range(nrecoil):
        zjetspass_templ.append(template(background, "Z+jets", "nominal", recoilbin, "sr", "pass"))
        zjetsfail_templ.append(template(background, "Z+jets", "nominal", recoilbin, "sr", "fail"))

    zjetsmodel = rl.Model("zjetsmodel")
    zjetseff = efficiency(zjetspass_templ, zjetsfail_templ, zjetsmodel)
    tf_MCtemplZ = rl.BernsteinPoly("tf_MCtemplZ", (1, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_MCtemplZ_params = zjetseff * tf_MCtemplZ(recoilscaled, msdscaled)

    wjetspass_templ = []
    wjetsfail_templ = []
    for recoilbin in range(nrecoil):
        wjetspass_templ.append(template(background, "W+jets", "nominal", recoilbin, "sr", "pass"))
        wjetsfail_templ.append(template(background, "W+jets", "nominal", recoilbin, "sr", "fail"))

    wjetsmodel = rl.Model("wjetsmodel")
    wjetseff = efficiency(wjetspass_templ, wjetsfail_templ, wjetsmodel)
    tf_MCtemplW = rl.BernsteinPoly("tf_MCtemplW", (1, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_MCtemplW_params = wjetseff * tf_MCtemplW(recoilscaled, msdscaled)

    def model(pass_templ, fail_templ, qcdmodel, tf_MCtempl_params):

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

    zjetsmodel = model(zjetspass_templ, zjetsfail_templ, zjetsmodel, tf_MCtemplZ_params)
    wjetsmodel = model(wjetspass_templ, wjetsfail_templ, wjetsmodel, tf_MCtemplW_params)

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

    def shape(fit, tf_MCtempl):
        param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', fit, param_names)
        tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
        tf_MCtempl_params_final = tf_MCtempl(recoilscaled, msdscaled)

        return tf_MCtempl_params_final

    tf_MCtemplW_params_final = shape(wjetsfit, tf_MCtemplW)
    tf_dataResidualW = rl.BernsteinPoly("tf_dataResidualW"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_dataResidualW_params = tf_dataResidualW(recoilscaled, msdscaled)

    tf_MCtemplZ_params_final = shape(zjetsfit, tf_MCtemplZ)
    tf_dataResidualZ = rl.BernsteinPoly("tf_dataResidualZ"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_dataResidualZ_params = tf_dataResidualZ(recoilscaled, msdscaled)

    TFs = {
      'tf_MCtemplW_params_final' = tf_MCtemplW_params_final,
      'tf_dataResidualW_params' = tf_dataResidualW_params,
      'tf_MCtemplZ_params_final' = tf_MCtemplZ_params_final,
      'tf_dataResidualZ_params' = tf_dataResidualZ_params
    }
    save(TFs, 'data/darkhiggs_rhalphabethTFs.coffea')
      
