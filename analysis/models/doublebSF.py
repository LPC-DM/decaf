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
    for r in data_hists["template"].identifiers("process"):
        data[str(r)] = data_hists["template"].sum("process", "gentype", overflow='all')

    background = {}
    signal = {}

    for i in range(5):
        if gentype_map[i] == 'bb':
            signal[str(gentype_map[i])] = bkg_hists["template"].integrate("gentype", 0).sum("process")
        else:
            background[str(gentype_map[i])] = bkg_hists["template"].integrate("gentype", i).sum("process")
