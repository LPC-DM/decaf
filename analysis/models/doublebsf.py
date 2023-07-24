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
pt_binning = {
        "2016": [350, 450, 500, 600, 2500],
        "2017": [350, 450, 500, 600, 2500],
        "2018": [350, 450, 500, 600, 2500]
        }
labels = {
    'QCD-$\mu$ (bb)':'bb',
    'QCD-$\mu$ (cc)':'cc',
    'QCD-$\mu$ (b)':'b',
    'QCD-$\mu$ (c)':'c',
    'QCD-$\mu$ (l)':'other',
}


### category: pass/fail flag
def template(dictionary, process, category, pt, read_sumw2=False):
    histogram = dictionary[process]#.integrate("process", process)
    nominal, sumw2 = histogram.values(sumw2=True)[()]
    nominal = nominal[:, pt, category_map[category]]
    sumw2 = sumw2[:, pt, category_map[category]]
    zerobins = nominal <= 0.
    output = nominal
    if "data" not in process:
        output[zerobins] = 1e-5
        sumw2[zerobins] = 0.
    binning = (dictionary[process].axis("svmass").edges())
    if read_sumw2:
        return (output, binning, "svmass", sumw2)
    return (output, binning, "svmass")

def get_mergedMC_stat_variations(dictionary, category, pt, bkg_list):
    templ=template(dictionary, bkg_list[0], category, pt, read_sumw2=True)
    merged_central=np.zeros_like(templ[0])
    merged_error2=np.zeros_like(templ[3])
    for bkg in bkg_list:
        templ=template(dictionary, bkg, category, pt, read_sumw2=True)
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
        effect_down[i] = max(epsilon, 1.0 - np.sqrt(merged_error2[i])/merged_central[i])
        templ.setParamEffect(param[i], effect_up, effect_down)
            
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

def addPrefiringSyst(templ, year):
    if '2018' not in year: templ.setParamEffect(prefiring, 1.01)
        
def model(year, category, pt):

    model_id = "svmass" + year + category + 'pt' + str(pt)
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

    dataTemplate = template(data, "data", category, pt)
    sr.setObservation(dataTemplate)
    
    ###
    # Other MC processes
    ###

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')
        
    mc_list = ['QCD-$\mu$ (bb)', 'QCD-$\mu$ (b)', 'QCD-$\mu$ (cc)', 'QCD-$\mu$ (c)', 'QCD-$\mu$ (l)']
    total_yields, total_error2 = get_mergedMC_stat_variations(mc, category, pt, mc_list)

    ###
    # QCD sig process
    ###

    ##### Template to use bb stat uncertainties
    sr_genbb_Template = template(mc, 'QCD-$\mu$ (bb)', category, pt, read_sumw2=True)
    sr_genbb = rl.TemplateSample(ch_name + "_genbb", rl.Sample.SIGNAL, sr_genbb_Template)
    addLumiSyst(sr_genbb, year)
    addPileupSyst(sr_genbb)
    addPrefiringSyst(sr_genbb, year)
    sr_genbb.setParamEffect(jes, 1.04)
    sr_genbb.setParamEffect(frac_bb, 1.2)
    sr_genbb.setParamEffect(doublebtag_weight['QCD-$\mu$ (bb)'], weight['QCD-$\mu$ (bb)'][category])
    addBBliteSyst(sr_genbb, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genbb)

    ###
    # QCD bkg processes
    ###

    sr_genb_Template = template(mc, 'QCD-$\mu$ (b)', category, pt, read_sumw2=True)
    sr_genb = rl.TemplateSample(ch_name + "_genb", rl.Sample.BACKGROUND, sr_genb_Template)
    addLumiSyst(sr_genb, year)
    addPileupSyst(sr_genb)
    addPrefiringSyst(sr_genb, year)
    sr_genb.setParamEffect(jes, 1.04)
    sr_genb.setParamEffect(frac_b, 1.2)
    sr_genb.setParamEffect(doublebtag_weight['QCD-$\mu$ (b)'], weight['QCD-$\mu$ (b)'][category])
    addBBliteSyst(sr_genb, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genb)

    sr_gencc_Template = template(mc, 'QCD-$\mu$ (cc)', category, pt, read_sumw2=True)
    sr_gencc = rl.TemplateSample(ch_name + "_gencc", rl.Sample.BACKGROUND, sr_gencc_Template)
    addLumiSyst(sr_gencc, year)
    addPileupSyst(sr_gencc)
    addPrefiringSyst(sr_gencc, year)
    sr_gencc.setParamEffect(jes, 1.04)
    sr_gencc.setParamEffect(frac_cc, 1.2)
    sr_gencc.setParamEffect(doublebtag_weight['QCD-$\mu$ (cc)'], weight['QCD-$\mu$ (cc)'][category])
    addBBliteSyst(sr_gencc, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_gencc)

    sr_genc_Template = template(mc, 'QCD-$\mu$ (c)', category, pt, read_sumw2=True)
    sr_genc = rl.TemplateSample(ch_name + "_genc", rl.Sample.BACKGROUND, sr_genc_Template)
    addLumiSyst(sr_genc, year)
    addPileupSyst(sr_genc)
    addPrefiringSyst(sr_genc, year)
    sr_genc.setParamEffect(jes, 1.04)
    sr_genc.setParamEffect(frac_c, 1.2)
    sr_genc.setParamEffect(doublebtag_weight['QCD-$\mu$ (c)'], weight['QCD-$\mu$ (c)'][category])
    addBBliteSyst(sr_genc, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genc)

    sr_genother_Template = template(mc, 'QCD-$\mu$ (l)', category, pt, read_sumw2=True)
    sr_genother = rl.TemplateSample(ch_name + "_genother", rl.Sample.BACKGROUND, sr_genother_Template)
    addLumiSyst(sr_genother, year)
    addPileupSyst(sr_genother)
    addPrefiringSyst(sr_genother, year)
    sr_genother.setParamEffect(jes, 1.04)
    sr_genother.setParamEffect(frac_other, 1.2)
    sr_genother.setParamEffect(doublebtag_weight['QCD-$\mu$ (l)'], weight['QCD-$\mu$ (l)'][category])
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
    # Extract histograms from input file
    ###

    #print("Extracting histograms for", year, category)
    hists = load("hists/doublebsf" + year + ".scaled")
    data_hists = hists["data"]
    if year == '2018':
        for k in data_hists:
                data_hists[k].scale(2.)
    bkg_hists = hists["bkg"]


    ###
    # Rebin templates for fit 
    ##
    data_hists["template"] = data_hists["template"].rebin("fj1pt", hist.Bin("fj1pt", "fj1pt", pt_binning[year]))
    bkg_hists["template"] = bkg_hists["template"].rebin("fj1pt", hist.Bin("fj1pt", "fj1pt", pt_binning[year]))

    ###
    # Preparing histograms for fit
    ##
    data = {}
    data['data'] = data_hists["template"].integrate('process','BTagMu')

    mc = {}
    for process in bkg_hists["template"].identifiers('process'):
        mc[str(process)] = bkg_hists["template"].integrate('process', process)

    ###
    ###
    # Setting up common systematics
    ###
    ###

    lumi = rl.NuisanceParameter("lumi" + year, "lnN")
    lumi_corr = rl.NuisanceParameter("lumi_corr", "lnN")
    lumi_1718 = rl.NuisanceParameter("lumi_1718", "lnN")
    pu = rl.NuisanceParameter("pu" + year, "lnN")
    prefiring = rl.NuisanceParameter("prefiring" + year, "lnN")
    jes = rl.NuisanceParameter("jes" + year, "lnN")
    
    #### fractional systematics (assume 50%)
    #frac_bb = rl.NuisanceParameter("frac_bb" + year, "lnN")

    ptbins = np.array(pt_binning[year])
    npt = len(ptbins) - 1
    for ptbin in range(npt):
        print(ptbin)

        frac_b = rl.NuisanceParameter("frac_b_pt" + str(ptbin) + year, "lnN")
        frac_c = rl.NuisanceParameter("frac_c_pt" + str(ptbin) + year, "lnN")
        frac_bb = rl.NuisanceParameter("frac_bb_pt" + str(ptbin) + year, "lnN")
        frac_cc = rl.NuisanceParameter("frac_cc_pt" + str(ptbin) + year, "lnN")
        frac_other = rl.NuisanceParameter("frac_other_pt" + str(ptbin) + year, "lnN")

        ###
        # Calculating efficiencies
        ###

        eff={}
        for k in mc:
            num=mc[k].integrate('svmass').values()[()][ptbin,1]
            den=mc[k].integrate('svmass').sum('ZHbbvsQCD').values()[()][ptbin]
            eff[k] = np.nan_to_num(num/den)
            print(k,num,den)


        #### SF weight (TemplateSample version) ####
        sf={}
        weight={}
        doublebtag_weight={}
        for k in mc:
            sf[k] = rl.IndependentParameter("sf"+ labels[k] + year + 'pt' + str(ptbin), 1.0, 0.01, 1.0 / eff[k])
            weight[k] = {
                "pass": rl.DependentParameter("weight"+labels[k], "{0}", sf[k]),
                "fail": rl.DependentParameter("weight"+labels[k], "(1-({0}*%f))/(1-%f)" % (eff[k], eff[k]), sf[k])
            }
            doublebtag_weight[k] = rl.IndependentParameter("doublebtag_weight"+ labels[k] + year + 'pt' + str(ptbin), 1.0)
        
        for category in ["pass", "fail"]:

                with open(
                    "data/models/doublebsf-"
                    + year
                    + "-"
                    + category
                    + "-pt"
                    + str(ptbin)
                    + ".model",
                    "wb",
        ) as fout:
                    pickle.dump(model(year, category, ptbin), fout, protocol=2)
