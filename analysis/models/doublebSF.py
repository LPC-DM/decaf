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

gentype_map = ['bb', 'b', 'cc', 'c', 'other']

category_map = {
        "pass": 1,
        "fail": 0
        }

bbtagger_eff = {
        "2016": 0.938,
        "2017": 0.947,
        "2018": 0.945
        }

### category: pass/fail flag
def template(dictionary, process, gentype, category, read_sumw2=False):
    histogram = dictionary[gentype].integrate("process", process)
    jp, sumw2 = histogram.values(sumw2=True)[()]
    jp = jp[:, category_map[category]]
    sumw2 = sumw2[:, category_map[category]]
    zerobins = jp <= 0.
    output = jp

    output[zerobins] = 1.
    sumw2[zerobins] = 0.

    binning = (
        dictionary[gentype].integrate("process", process).axis("btagJP").edges()
    )

    if read_sumw2:
        return (output, binning, "btagJP", sumw2)
    return (output, binning, "btagJP")

### s: process in the signal region
def model(year, category):

    model_id = "btagJP" + year + category
    model = rl.Model(model_id)

    ###
    ###
    # Signal region
    ###
    ###

    ch_name = model_id
    sr = rl.Channel(ch_name)
    model.addChannel(sr)

    ###
    # Add data distribution to the channel
    ###

    sr.setObservation(template(data, "BTagMu", "bb", category))

    ###
    # QCD sig process
    ###

    #sr_genbb_Template = template(signal, "QCD", "bb", category)
    #sr_genbb = rl.TemplateSample(ch_name + "_genbb", rl.Sample.SIGNAL, sr_genbb_Template)
    #sr_genbb.setParamEffect(lumi, 1.027)
    #sr_genbb.setParamEffect(pu, 1.05)
    #sr_genbb.setParamEffect(jes, 1.02)
    #sr_genbb.setParamEffect(frac_bb, 1.5)
    #sr.addSample(sr_genbb)

    sr_genbb_Template = template(signal, "QCD", "bb", category)
    sr_genbb_Observable = rl.Observable("btagJP", sr_genbb_Template[1])
    sr_genbb_BinYields = sr_genbb_Template[0] * weight[category]
    sr_genbb = rl.ParametericSample(ch_name + "_genbb", rl.Sample.SIGNAL, sr_genbb_Observable, sr_genbb_BinYields)
    sr_genbb.setParamEffect(lumi, 1.027)
    sr_genbb.setParamEffect(pu, 1.05)
    sr_genbb.setParamEffect(jes, 1.02)
    sr_genbb.setParamEffect(frac_bb, 1.5)
    sr.addSample(sr_genbb)

    ###
    # QCD bkg processes
    ###

    sr_genb_Template = template(background, "QCD", "b", category)
    sr_genb = rl.TemplateSample(ch_name + "_genb", rl.Sample.BACKGROUND, sr_genb_Template)
    sr_genb.setParamEffect(lumi, 1.027)
    sr_genb.setParamEffect(pu, 1.05)
    sr_genb.setParamEffect(jes, 1.02)
    sr_genb.setParamEffect(frac_b, 1.5)
    sr.addSample(sr_genb)

    sr_genc_Template = template(background, "QCD", "c", category)
    sr_genc = rl.TemplateSample(ch_name + "_genc", rl.Sample.BACKGROUND, sr_genc_Template)
    sr_genc.setParamEffect(lumi, 1.027)
    sr_genc.setParamEffect(pu, 1.05)
    sr_genc.setParamEffect(jes, 1.02)
    sr_genc.setParamEffect(frac_c, 1.5)
    sr.addSample(sr_genc)

    sr_gencc_Template = template(background, "QCD", "cc", category)
    sr_gencc = rl.TemplateSample(ch_name + "_gencc", rl.Sample.BACKGROUND, sr_gencc_Template)
    sr_gencc.setParamEffect(lumi, 1.027)
    sr_gencc.setParamEffect(pu, 1.05)
    sr_gencc.setParamEffect(jes, 1.02)
    sr_gencc.setParamEffect(frac_cc, 1.5)
    sr.addSample(sr_gencc)

    sr_genother_Template = template(background, "QCD", "other", category)
    sr_genother = rl.TemplateSample(ch_name + "_genother", rl.Sample.BACKGROUND, sr_genother_Template)
    sr_genother.setParamEffect(lumi, 1.027)
    sr_genother.setParamEffect(pu, 1.05)
    sr_genother.setParamEffect(jes, 1.02)
    sr_genother.setParamEffect(frac_other, 1.5)
    sr.addSample(sr_genother)

    return model

if __name__ == "__main__":
    if not os.path.exists("datacards"):
        os.mkdir("datacards")
    parser = OptionParser()
    parser.add_option("-y", "--year", help="year", dest="year", default="")
    #parser.add_option("-n", "--name", help="name of category (pass or fail)", dest="name", default="")
    (options, args) = parser.parse_args()
    year = options.year
    #name = options.name

    ###
    # Extract histograms from input file
    ###

    print("Grouping histograms")
    hists = load("hists/doublebSF" + year + ".scaled")
    data_hists = hists["data"]
    bkg_hists = hists["bkg"]

    ###
    # Preparing histograms for fit
    ##

    data = {}
    data['bb'] = data_hists["template"].sum("gentype", overflow='all')

    background = {}
    signal = {}

    for i in range(5):
        if gentype_map[i] == 'bb':
            signal[str(gentype_map[i])] = bkg_hists["template"].integrate("gentype", 0)
        else:
            background[str(gentype_map[i])] = bkg_hists["template"].integrate("gentype", i)

    #print(data)
    #print(signal)
    #print(background)

    ###
    ###
    # Setting up systematics
    ###
    ###

    lumi = rl.NuisanceParameter("lumi" + year, "lnN")
    pu = rl.NuisanceParameter("pu" + year, "lnN")
    jes = rl.NuisanceParameter("jes" + year, "lnN")

    #### fractional systematics (assume 50%)
    #frac_bb = rl.NuisanceParameter("frac_bb" + year, "lnN")
    #frac_b = rl.NuisanceParameter("frac_b" + year + str(name), "lnN")
    #frac_cc = rl.NuisanceParameter("frac_cc" + year + str(name), "lnN")
    #frac_c = rl.NuisanceParameter("frac_c" + year + str(name), "lnN")
    #frac_other = rl.NuisanceParameter("frac_other" + year + str(name), "lnN")

    #### SF weight
    sf = rl.IndependentParameter("sf" + year, 1.0, 0.01, 1.0 / bbtagger_eff[year])
    weight = {
            "pass": sf,
            "fail": (1 - (sf*bbtagger_eff[year])) / (1 - bbtagger_eff[year])
            }

    for category in ["pass", "fail"]:
        #### fractional systematics (assume 50%)
        frac_bb = rl.NuisanceParameter("frac_bb" + year, "lnN")
        frac_b = rl.NuisanceParameter("frac_b" + year + category, "lnN")
        frac_cc = rl.NuisanceParameter("frac_cc" + year + category, "lnN")
        frac_c = rl.NuisanceParameter("frac_c" + year + category, "lnN")
        frac_other = rl.NuisanceParameter("frac_other" + year + category, "lnN")
        with open(
            "data/doublebSF-"
            + year
            + "-"
            + category
            + ".model",
            "wb",
        ) as fout:
            pickle.dump(model(year, category), fout, protocol=2)
