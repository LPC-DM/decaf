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
        "2016": [350, 400, 450, 500, 600, 1250],
        "2017": [350, 400, 450, 500, 600, 1250],
        "2018": [350, 400, 450, 500, 600, 1250]
        }
weight_dict={
    '2016':[0.80982663, 0.85100921, 0.86931336, 0.89259707, 0.7907999 ],
    '2017':[0.86549006, 0.90644119, 0.93265111, 0.9682362,  0.87836352],
    '2018':[0.70589636, 0.9585072,  1.02255474, 1.07915829, 1.00106306]
}

### category: pass/fail flag
def template(dictionary, process, gentype, category, pt, read_sumw2=False):
    histogram = dictionary[gentype].integrate("process", process)
    if "data" not in gentype:
        histogram.scale(weight_dict[year][pt])
    nominal, sumw2 = histogram.values(sumw2=True)[()]
    nominal = nominal[:, pt, category_map[category]]
    sumw2 = sumw2[:, pt, category_map[category]]
    zerobins = nominal <= 0.
    output = nominal
    if "data" not in gentype:
        output[zerobins] = 1e-5
        sumw2[zerobins] = 0.
    binning = (dictionary[gentype].integrate("process", process).axis("svmass").edges())
    if read_sumw2:
        return (output, binning, "svmass", sumw2)
    return (output, binning, "svmass")

def get_mergedMC_stat_variations(dictionary, category, pt, mc_list):
    templ=template(dictionary, 'QCD', mc_list[0], category, pt, read_sumw2=True)
    merged_central=np.zeros_like(templ[0])
    merged_error2=np.zeros_like(templ[3])
    for mc in mc_list:
        templ=template(dictionary, 'QCD', mc, category, pt, read_sumw2=True)
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

    dataTemplate = template(data, "BTagMu", "data", category, pt)
    sr.setObservation(dataTemplate)
    
    ###
    # Other MC processes
    ###

    nbins = len(dataTemplate[1]) - 1
    param = [None for _ in range(nbins)]
    for i in range(nbins):
        param[i] = rl.NuisanceParameter(ch_name + '_mcstat_bin%i' % i, combinePrior='shape')
        
    mc_list = ['bb', 'b', 'cc', 'c', 'other']
    total_yields, total_error2 = get_mergedMC_stat_variations(mc, category, pt, mc_list)

    ###
    # QCD sig process
    ###

    ##### Template to use bb stat uncertainties
    sr_genbb_Template = template(mc, "QCD", "bb", category, pt, read_sumw2=True)
    sr_genbb = rl.TemplateSample(ch_name + "_genbb", rl.Sample.SIGNAL, sr_genbb_Template)
    addLumiSyst(sr_genbb, year)
    addPileupSyst(sr_genbb)
    addPrefiringSyst(sr_genbb, year)
    sr_genbb.setParamEffect(jes, 1.04)
    sr_genbb.setParamEffect(frac_b, 1.2)
    sr_genbb.setParamEffect(sf_weight['bb'], weight['bb'][category])
    addBBliteSyst(sr_genbb, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genbb)

    ###
    # QCD bkg processes
    ###

    sr_genb_Template = template(mc, "QCD", "b", category, pt, read_sumw2=True)
    sr_genb = rl.TemplateSample(ch_name + "_genb", rl.Sample.BACKGROUND, sr_genb_Template)
    addLumiSyst(sr_genb, year)
    addPileupSyst(sr_genb)
    addPrefiringSyst(sr_genb, year)
    sr_genb.setParamEffect(jes, 1.04)
    sr_genb.setParamEffect(frac_b, 1.2)
    sr_genb.setParamEffect(sf_weight['b'], weight['b'][category])
    addBBliteSyst(sr_genb, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genb)

    sr_genc_Template = template(mc, "QCD", "c", category, pt, read_sumw2=True)
    sr_genc = rl.TemplateSample(ch_name + "_genc", rl.Sample.BACKGROUND, sr_genc_Template)
    addLumiSyst(sr_genc, year)
    addPileupSyst(sr_genc)
    addPrefiringSyst(sr_genc, year)
    sr_genc.setParamEffect(jes, 1.04)
    sr_genc.setParamEffect(frac_c, 1.2)
    sr_genc.setParamEffect(sf_weight['c'], weight['c'][category])
    addBBliteSyst(sr_genc, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_genc)

    sr_gencc_Template = template(mc, "QCD", "cc", category, pt, read_sumw2=True)
    sr_gencc = rl.TemplateSample(ch_name + "_gencc", rl.Sample.BACKGROUND, sr_gencc_Template)
    addLumiSyst(sr_gencc, year)
    addPileupSyst(sr_gencc)
    addPrefiringSyst(sr_gencc, year)
    sr_gencc.setParamEffect(jes, 1.04)
    sr_gencc.setParamEffect(frac_c, 1.2)
    sr_gencc.setParamEffect(sf_weight['cc'], weight['cc'][category])
    addBBliteSyst(sr_gencc, param, total_yields, total_error2, epsilon=1e-5)
    sr.addSample(sr_gencc)

    sr_genother_Template = template(mc, "QCD", "other", category, pt, read_sumw2=True)
    sr_genother = rl.TemplateSample(ch_name + "_genother", rl.Sample.BACKGROUND, sr_genother_Template)
    addLumiSyst(sr_genother, year)
    addPileupSyst(sr_genother)
    addPrefiringSyst(sr_genother, year)
    sr_genother.setParamEffect(jes, 1.04)
    sr_genother.setParamEffect(frac_other, 1.2)
    sr_genother.setParamEffect(sf_weight['other'], weight['other'][category])
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
    data['data'] = data_hists["template"].sum("gentype", overflow='all')

    mc = {}
    for i in range(5):
        mc[str(gentype_map[i])] = bkg_hists["template"].integrate("gentype", i)

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

        frac_b = rl.NuisanceParameter("frac_b_pt" + str(ptbin), "lnN")
        frac_c = rl.NuisanceParameter("frac_c_pt" + str(ptbin), "lnN")
        frac_other = rl.NuisanceParameter("frac_other_pt" + str(ptbin), "lnN")

        ###
        # Calculating efficiencies
        ###

        eff={}
        for i in range(5):
            num=mc[str(gentype_map[i])].integrate('svmass').integrate('process').values()[()][ptbin,1]
            den=mc[str(gentype_map[i])].integrate('svmass').integrate('process').sum('ZHbbvsQCD').values()[()][ptbin]
            eff[str(gentype_map[i])] = np.nan_to_num(num/den)
            print(gentype_map[i],eff[str(gentype_map[i])])
            print(num,den)


        #### SF weight (TemplateSample version) ####
        sf={}
        weight={}
        sf_weight={}
        for i in range(5):
            sf[str(gentype_map[i])] = rl.IndependentParameter("sf"+ str(gentype_map[i]) + year + 'pt' + str(ptbin), 1.0, 0.01, 1.0 / eff[str(gentype_map[i])])
            weight[str(gentype_map[i])] = {
                "pass": rl.DependentParameter("weight"+str(gentype_map[i]), "{0}", sf[str(gentype_map[i])]),
                "fail": rl.DependentParameter("weight"+str(gentype_map[i]), "(1-({0}*%f))/(1-%f)" % (eff[str(gentype_map[i])], eff[str(gentype_map[i])]), sf[str(gentype_map[i])])
            }
            sf_weight[str(gentype_map[i])] = rl.IndependentParameter("sf_weight"+ str(gentype_map[i]) + year + 'pt' + str(ptbin), 1.0)
        
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
