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

category_map = {
        "pass": 1,
        "fail": 0
        }
pt_binning = {
        "2016": [350, 400, 500, 2500],
        "2017": [350, 400, 500, 2500],
        "2018": [350, 400, 500, 2500]
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
        
    mc_list = ['QCD-$\mu$ (b+bb)', 'QCD-$\mu$ (c+cc)', 'QCD-$\mu$ (l)']
    total_yields, total_error2 = get_mergedMC_stat_variations(mc, category, pt, mc_list)

    ###
    # QCD processes
    ###

    sr_genbb_Template = template(mc, 'QCD-$\mu$ (b+bb)', category, pt, read_sumw2=True)
    sr_genbb = rl.TemplateSample(ch_name + "_genbb", rl.Sample.SIGNAL, sr_genbb_Template)
    addLumiSyst(sr_genbb, year)
    addPileupSyst(sr_genbb)
    addPrefiringSyst(sr_genbb, year)
    #sr_genbb.setParamEffect(frac_bb, 1.2)
    sr_genbb.setParamEffect(jes, 1.04)
    sr_genbb.setParamEffect(doublebtag_weight['bs'], weight['bs'][category])
    addBBliteSyst(sr_genbb, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genbb)

    sr_gencc_Template = template(mc, 'QCD-$\mu$ (c+cc)', category, pt, read_sumw2=True)
    sr_gencc = rl.TemplateSample(ch_name + "_gencc", rl.Sample.SIGNAL, sr_gencc_Template)
    addLumiSyst(sr_gencc, year)
    addPileupSyst(sr_gencc)
    addPrefiringSyst(sr_gencc, year)
    #sr_gencc.setParamEffect(frac_cc, 1.2)
    sr_gencc.setParamEffect(jes, 1.04)
    sr_gencc.setParamEffect(doublebtag_weight['cs'], weight['cs'][category])
    addBBliteSyst(sr_gencc, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_gencc)

    sr_genother_Template = template(mc, 'QCD-$\mu$ (l)', category, pt, read_sumw2=True)
    sr_genother = rl.TemplateSample(ch_name + "_genother", rl.Sample.SIGNAL, sr_genother_Template)
    addLumiSyst(sr_genother, year)
    addPileupSyst(sr_genother)
    addPrefiringSyst(sr_genother, year)
    #sr_genother.setParamEffect(frac_other, 1.2)
    sr_genother.setParamEffect(jes, 1.04)
    sr_genother.setParamEffect(doublebtag_weight['other'], weight['other'][category])
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
    
    process = hist.Cat("process", "Process", sorting='placement')
    cats = ("process",)
    bkg_map = OrderedDict()
    bkg_map['QCD-$\mu$ (b+bb)'] = (['QCD-$\mu$ (bb)','QCD-$\mu$ (b)'],)
    bkg_map['QCD-$\mu$ (c+cc)'] = (['QCD-$\mu$ (cc)','QCD-$\mu$ (c)'],)
    bkg_map['QCD-$\mu$ (l)'] = ('QCD-$\mu$ (l)',)
    bkg_hists["template"] = bkg_hists["template"].group(cats, process, bkg_map)
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
    frac_bb = rl.NuisanceParameter("frac_bb" + year, "lnN")
    frac_cc = rl.NuisanceParameter("frac_cc" + year, "lnN")
    frac_other = rl.NuisanceParameter("frac_other" + year, "lnN")
    
    ptbins = np.array(pt_binning[year])
    npt = len(ptbins) - 1
    for ptbin in range(npt):
        print(ptbin)

        ###
        # Calculating efficiencies
        ###

        eff={}
        
        num=mc['QCD-$\mu$ (b+bb)'].integrate('svmass').values()[()][ptbin,1]
        den=mc['QCD-$\mu$ (b+bb)'].integrate('svmass').sum('ZHbbvsQCD').values()[()][ptbin]
        eff['bs'] = np.nan_to_num(num/den)

        num=mc['QCD-$\mu$ (c+cc)'].integrate('svmass').values()[()][ptbin,1]
        den=mc['QCD-$\mu$ (c+cc)'].integrate('svmass').sum('ZHbbvsQCD').values()[()][ptbin]
        eff['cs'] = np.nan_to_num(num/den)

        num=mc['QCD-$\mu$ (l)'].integrate('svmass').values()[()][ptbin,1]
        den=mc['QCD-$\mu$ (l)'].integrate('svmass').sum('ZHbbvsQCD').values()[()][ptbin]
        eff['other'] = np.nan_to_num(num/den)

        #### SF weight (TemplateSample version) ####
        sf={}
        weight={}
        doublebtag_weight={}
        for k in eff:
            sf[k] = rl.IndependentParameter("sf"+ k + year + 'pt' + str(ptbin), 1.0, 0.01, 1.0 / eff[k])
            weight[k] = {
                "pass": rl.DependentParameter("weight"+k, "{0}", sf[k]),
                "fail": rl.DependentParameter("weight"+k, "(1-({0}*%f))/(1-%f)" % (eff[k], eff[k]), sf[k])
            }
            doublebtag_weight[k] = rl.IndependentParameter("doublebtag_weight"+ k + year + 'pt' + str(ptbin), 1.0)
        
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
