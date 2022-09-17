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

binning = {
        "2016": [-0.8, -0.4, 0.,  0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.5, 3.2],
        "2017": [-0.8, -0.4, 0.,  0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.5, 3.2],
        "2018": [-0.8, -0.4, 0.,  0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.5, 3.2]
        }

### category: pass/fail flag
def template(dictionary, process, gentype, category, read_sumw2=False):
    histogram = dictionary[gentype].integrate("process", process)
    nominal, sumw2 = histogram.values(sumw2=True)[()]
    nominal = nominal[:, category_map[category]]
    sumw2 = sumw2[:, category_map[category]]
    zerobins = nominal <= 0.
    output = nominal
    if "data" not in gentype:
        output[zerobins] = 1e-5
        sumw2[zerobins] = 0.
    binning = (dictionary[gentype].integrate("process", process).axis("svmass").edges())
    if read_sumw2:
        return (output, binning, "svmass", sumw2)
    return (output, binning, "svmass")

def get_mergedMC_stat_variations(dictionary, category, mc_list):
    templ=template(dictionary, 'QCD', mc_list[0], category, read_sumw2=True)
    merged_central=np.zeros_like(templ[0])
    merged_error2=np.zeros_like(templ[3])
    for mc in mc_list:
        templ=template(dictionary, 'QCD', mc, category, read_sumw2=True)
        for i in range(len(templ[0])):
            if templ[0][i] <= 1e-5 or templ[3][i] <= 0.:
                continue
            merged_central[i] += templ[0][i]
            merged_error2[i]  += templ[3][i]
    return merged_central, merged_error2

def addBBliteSyst(templ, param, merged_central, merged_error2, epsilon=0):
    for i in range(templ.observable.nbins):
        if merged_central[i] <= 0. or merged_error2[i] <= 0.:
            continue
        if templ._nominal[i] <= 1e-5:
            continue
        effect_up = np.ones_like(templ._nominal)
        effect_down = np.ones_like(templ._nominal)
        effect_up[i] = 1.0 + np.sqrt(merged_error2[i])/merged_central[i]
        print('Effect up = ',effect_up[i])
        effect_down[i] = max(epsilon, 1.0 - np.sqrt(merged_error2[i])/merged_central[i])
        print('Effect down = ', effect_down[i])
        print('Central value',templ._nominal[i])
        templ.setParamEffect(param[i], effect_up, effect_down)
        
def model(year, category):

    model_id = "svmass" + year + category
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

    dataTemplate = template(data, "BTagMu", "data", category)
    sr.setObservation(dataTemplate)
    
    ###
    # Other MC processes
    ###

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')
        
    mc_list = ['bb', 'b', 'cc', 'c', 'other']
    total_yields, total_error2 = get_mergedMC_stat_variations(mc, category, mc_list)

    ###
    # QCD sig process
    ###

    ##### Template to use bb stat uncertainties
    sr_genbb_Template = template(mc, "QCD", "bb", category, read_sumw2=True)
    sr_genbb = rl.TemplateSample(ch_name + "_genbb", rl.Sample.SIGNAL, sr_genbb_Template)
    sr_genbb.setParamEffect(lumi, nlumi)
    sr_genbb.setParamEffect(pu, npu)
    sr_genbb.setParamEffect(jes, njes)
    sr_genbb.setParamEffect(frac, nfrac)
    sr_genbb.setParamEffect(sf_weight, weight[category])
    addBBliteSyst(sr_genbb, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genbb)

    ###
    # QCD bkg processes
    ###

    sr_genb_Template = template(mc, "QCD", "b", category, read_sumw2=True)
    sr_genb = rl.TemplateSample(ch_name + "_genb", rl.Sample.BACKGROUND, sr_genb_Template)
    sr_genb.setParamEffect(lumi, nlumi)
    sr_genb.setParamEffect(pu, npu)
    sr_genb.setParamEffect(jes, njes)
    sr_genb.setParamEffect(frac, nfrac)
    addBBliteSyst(sr_genb, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genb)

    sr_genc_Template = template(mc, "QCD", "c", category, read_sumw2=True)
    sr_genc = rl.TemplateSample(ch_name + "_genc", rl.Sample.BACKGROUND, sr_genc_Template)
    sr_genc.setParamEffect(lumi, nlumi)
    sr_genc.setParamEffect(pu, npu)
    sr_genc.setParamEffect(jes, njes)
    sr_genc.setParamEffect(frac, nfrac)
    addBBliteSyst(sr_genc, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genc)

    sr_gencc_Template = template(mc, "QCD", "cc", category, read_sumw2=True)
    sr_gencc = rl.TemplateSample(ch_name + "_gencc", rl.Sample.BACKGROUND, sr_gencc_Template)
    sr_gencc.setParamEffect(lumi, nlumi)
    sr_gencc.setParamEffect(pu, npu)
    sr_gencc.setParamEffect(jes, njes)
    sr_gencc.setParamEffect(frac, nfrac)
    addBBliteSyst(sr_gencc, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_gencc)

    sr_genother_Template = template(mc, "QCD", "other", category, read_sumw2=True)
    sr_genother = rl.TemplateSample(ch_name + "_genother", rl.Sample.BACKGROUND, sr_genother_Template)
    sr_genother.setParamEffect(lumi, nlumi)
    sr_genother.setParamEffect(pu, npu)
    sr_genother.setParamEffect(jes, njes)
    sr_genother.setParamEffect(frac, nfrac)
    addBBliteSyst(sr_genother, param, total_yields, total_error2, epsilon=1e-5)
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
    ###
    # Setting up common systematics
    ###
    ###

    lumi = rl.NuisanceParameter("lumi" + year, "lnN")
    pu = rl.NuisanceParameter("pu" + year, "lnN")
    jes = rl.NuisanceParameter("jes" + year, "lnN")

    #### fractional systematics (assume 50%)
    frac = rl.NuisanceParameter("frac" + year, "lnN")

    ###
    # Set lnN or shape parameters
    ### 

    nlumi = 1.027
    npu = 1.05
    njes = 1.02
    nfrac = 1.5

    #### SF weight (TemplateSample version) ####
    sf = rl.IndependentParameter("sf" + year, 1.0, 0.01, 1.0 / bbtagger_eff[year])
    weight = {
            "pass": rl.DependentParameter("weight", "{0}", sf),
            "fail": rl.DependentParameter("weight", "(1-({0}*%f))/(1-%f)" % (bbtagger_eff[year], bbtagger_eff[year]), sf)
            }
    sf_weight = rl.IndependentParameter("sf_weight" + year, 1.0)

    for category in ["pass", "fail"]:

        ###
        # Extract histograms from input file
        ###

        print("Extracting histograms for", year, category)
        hists = load("hists/doublebsf" + year + ".scaled")
        data_hists = hists["data"]
        bkg_hists = hists["bkg"]

        #### Setting up fractional systematics (assume 50%)
        frac_b = rl.NuisanceParameter("frac_b" + year + category, "lnN")
        frac_cc = rl.NuisanceParameter("frac_cc" + year + category, "lnN")
        frac_c = rl.NuisanceParameter("frac_c" + year + category, "lnN")
        frac_other = rl.NuisanceParameter("frac_other" + year + category, "lnN")


        ###
        # Rebin templates for fit 
        ##
        data_hists["svtemplate"] = data_hists["svtemplate"].rebin("svmass", hist.Bin("svmass", "svmass", binning[year]))
        bkg_hists["svtemplate"] = bkg_hists["svtemplate"].rebin("svmass", hist.Bin("svmass", "svmass", binning[year]))

        ###
        # Preparing histograms for fit
        ##
        data = {}
        data['data'] = data_hists["svtemplate"].sum("gentype", overflow='all')

        mc = {}
        for i in range(5):
            mc[str(gentype_map[i])] = bkg_hists["svtemplate"].integrate("gentype", i)

        with open(
            "data/doublebsf-"
            + year
            + "-"
            + category
            + ".model",
            "wb",
        ) as fout:
            pickle.dump(model(year, category), fout, protocol=2)
