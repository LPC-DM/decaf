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

category_map = {"pass": 1, "fail": 0}
gentype_map = ['bb', 'b', 'cc', 'c', 'other']

### category: pass/fail flag
def template(dictionary, process, gentype, category, read_sumw2=False):
    histogram = dictionary[gentype].integrate("process", process)
    jp, sumw2 = histogram.values(sumw2=True)[()]
    jp = jp[category_map[category]]
    sumw2 = sumw2[category_map[category]]
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
def model(year, btagJP, category, s):

    model_id = year + category + "btagJP"
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

    sr.setObservation(template(data, "BTagMu", "bb", category))

    ###
    # QCD sig process
    ###

    sr_genbb_Template = template(signal, "QCD", "bb", category)
    sr_genbb = rl.TemplateSample(ch_name + "_genbb", rl.Sample.SIGNAL, sr_genbb_Template)
    sr.addSample(sr_genbb)

    ###
    # QCD bkg processes
    ###

    sr_genb_Template = template(background, "QCD", "b", category)
    sr_genb = rl.TemplateSample(ch_name + "_genb", rl.Sample.BACKGROUND, sr_genb_Template)
    sr.addSample(sr_genb)

    sr_genc_Template = template(background, "QCD", "c", category)
    sr_genc = rl.TemplateSample(ch_name + "_genc", rl.Sample.BACKGROUND, sr_genc_Template)
    sr.addSample(sr_genc)

    sr_gencc_Template = template(background, "QCD", "cc", category)
    sr_gencc = rl.TemplateSample(ch_name + "_gencc", rl.Sample.BACKGROUND, sr_gencc_Template)
    sr.addSample(sr_gencc)

    sr_genother_Template = template(background, "QCD", "other", category)
    sr_genother = rl.TemplateSample(ch_name + "_genother", rl.Sample.BACKGROUND, sr_genother_Template)
    sr.addSample(sr_genother)

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

    print(data)
    print(signal)
    print(background)
