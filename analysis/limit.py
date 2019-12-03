from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import sys
import os
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import ROOT
import gzip
from coffea import hist, processor 

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

def expo_sample(norm, scale, obs):
    cdf = scipy.stats.expon.cdf(scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)

def template(hist, name):
    return (hist.values(overflow='all')[()], hist.axis(name).edges(overflow='all'), name)

def darkhiggs_model(tmpdir):

    ###
    #
    # Loading MC histograms
    #
    ###

    mc_hists={}
    pd = []
    year = '2018'
    dirname = 'pods/darkhiggs' + year
    for filename in os.listdir(dirname):
        if 'MET' in filename or 'SingleElectron' in filename or 'SinglePhoton' in filename or 'EGamma' in filename: continue
        if '.pkl.gz' in filename:
            if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])
            with gzip.open(dirname+'/'+filename) as fin:
                hin = pickle.load(fin)
                for k in hin.keys():
                    if k in mc_hists: mc_hists[k]+=hin[k]
                    else: mc_hists[k]=hin[k]
    ##
    # Defining primary datasets (pdataset) to aggregate all the histograms that belong to a single sample
    ##

    pdataset = hist.Cat("pdataset", "pdataset", sorting='placement')
    pdataset_cats = ("dataset",)
    pdataset_map = OrderedDict()
    for pdi in pd:
        pdataset_map[pdi] = (pdi+"*",)
    for key in mc_hists.keys():
        mc_hists[key] = mc_hists[key].group(pdataset_cats, pdataset, pdataset_map)

    ###
    # Rescaling MC histograms using the xsec weight
    ###

    scale={}
    for pdi in mc_hists['sumw'].identifiers('pdataset'):
        scale[pdi]=mc_hists['sumw'].integrate('pdataset', pdi).values(overflow='all')[()][1]
    for key in mc_hists.keys():
        if key=='sumw': continue
        for pdi in mc_hists[key].identifiers('pdataset'):
            mc_hists[key].scale({pdi:1/scale[pdi]},axis='pdataset')

    ###
    # Defining 'process', to aggregate different samples into a single process
    ##

    process = hist.Cat("process", "Process", sorting='placement')
    bkg_cats = ("pdataset",)
    bkg_map = OrderedDict()
    bkg_map["Hbb"] = ("*HToBB*")
    bkg_map["DY"] = ("DYJets*",)
    bkg_map["Diboson"] = ("*_TuneCP5_13TeV-pythia8",)
    bkg_map["ST"] = ("ST*",)
    bkg_map["TT"] = ("TT*",)
    bkg_map["Wjets"] = ("WJets*",)
    bkg_map["ZJets"] = ("ZJetsToNuNu*",)   ## temporarily 
    bkg_hists = {}

    signal_cats = ("pdataset",)
    signal_map = OrderedDict() ### for signal samples
    signal_map["Mhs_50"] = ("*Mhs_50*",)  ## signals
    signal_map["Mhs_70"] = ("*Mhs_70*",)
    signal_map["Mhs_90"] = ("*Mhs_90*",)
    signal_map["MonoJet"] = ("MonoJet*",)  ## signals
    signal_map["MonoW"] = ("MonoW*",)    ## signals
    signal_map["MonoZ"] = ("MonoZ*",)    ## signals
    signal_hists = {}

    ###
    # Storing signal and background histograms
    ###

    for key in mc_hists.keys():
        signal_hists[key] = mc_hists[key].group(signal_cats, process, signal_map)
        bkg_hists[key] = mc_hists[key].group(bkg_cats, process, bkg_map)

    ###
    #
    # Loading data histograms
    #
    ###

    data_hists={}
    for filename in os.listdir(dirname):
        if 'MET' in filename or 'SingleElectron' in filename or 'SinglePhoton' in filename or 'EGamma' in filename:
            if '.pkl.gz' in filename:
                with gzip.open(dirname+'/'+filename) as fin:
                    hin = pickle.load(fin)
                    for k in hin.keys():
                        if k in data_hists: data_hists[k]+=hin[k]
                        else: data_hists[k]=hin[k]

    data_map = OrderedDict()
    data_map["MET"] = ("MET*", )
    data_map["SingleElectron"] = ("EGamma*", )
    data_map["SinglePhoton"] = ("EGamma*", )
    data_cats = ("dataset",)
    for key in data_hists.keys():
        data_hists[key] = data_hists[key].group(data_cats, process, data_map)

    model = rl.Model("DarkHiggs")

    # lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    jec = rl.NuisanceParameter('CMS_jec', 'shape')
    ele_id_eff = rl.NuisanceParameter('CMS_ele_id_eff', 'shape')
    pho_id_eff = rl.NuisanceParameter('CMS_pho_id_eff', 'shape')
    gamma_to_z_ewk = rl.NuisanceParameter('Theory_gamma_z_ewk', 'shape')

    recoilbins = np.linspace(300, 1200, 13)
    recoil = rl.Observable('recoil', recoilbins)

    signalCh = rl.Channel("signalCh")
    model.addChannel(signalCh)

    #    zvvTemplate = expo_sample(1000, 400, recoil)
    zvvHist = signal_hists['recoil'].integrate('jet_selection','baggy').integrate('region','iszeroL').integrate('process', 'Mhs_90')
    zvvTemplate = template(zvvHist, 'recoil')
    zvvJetsMC = rl.TemplateSample('zvvJetsMC', rl.Sample.BACKGROUND, zvvTemplate)
    zvvJetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(zvvHist.axis('recoil').edges(overflow='all'))-1))

    # these parameters are large, should probably log-transform them
    zvvBinYields = np.array([rl.IndependentParameter('tmp', b, 0, zvvTemplate[0].max()*2) for b in zvvTemplate[0]])  # name will be changed by ParametericSample
    zvvJets = rl.ParametericSample('signalCh_zvvJets', rl.Sample.BACKGROUND, recoil, zvvBinYields)
    signalCh.addSample(zvvJets)

    dmTemplate = expo_sample(100, 800, recoil)
    dmSample = rl.TemplateSample('signalCh_someDarkMatter', rl.Sample.SIGNAL, dmTemplate)
    signalCh.addSample(dmSample)

    signalCh.setObservation(expo_sample(1000, 400, recoil))

    zllCh = rl.Channel("zllCh")
    model.addChannel(zllCh)

    zllTemplate = expo_sample(1000*6.6/20, 400, recoil)
    zllJetsMC = rl.TemplateSample('zllJetsMC', rl.Sample.BACKGROUND, zllTemplate)
    zllJetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
    zllJetsMC.setParamEffect(ele_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins), np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

    zllTransferFactor = zllJetsMC.getExpectation() / zvvJetsMC.getExpectation()
    zllJets = rl.TransferFactorSample('zllCh_zllJets', rl.Sample.BACKGROUND, zllTransferFactor, zvvJets)
    zllCh.addSample(zllJets)

    otherbkgTemplate = expo_sample(200, 250, recoil)
    otherbkg = rl.TemplateSample('zllCh_otherbkg', rl.Sample.BACKGROUND, otherbkgTemplate)
    otherbkg.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=recoil.nbins))
    zllCh.addSample(otherbkg)

    zllCh.setObservation(expo_sample(1200, 380, recoil))

    gammaCh = rl.Channel("gammaCh")
    model.addChannel(gammaCh)

    gammaTemplate = expo_sample(2000, 450, recoil)
    gammaJetsMC = rl.TemplateSample('gammaJetsMC', rl.Sample.BACKGROUND, gammaTemplate)
    gammaJetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
    gammaJetsMC.setParamEffect(pho_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

    gammaTransferFactor = gammaJetsMC.getExpectation() / zvvJetsMC.getExpectation()
    gammaJets = rl.TransferFactorSample('gammaCh_gammaJets', rl.Sample.BACKGROUND, gammaTransferFactor, zvvJets)
    gammaJets.setParamEffect(gamma_to_z_ewk, np.linspace(1.01, 1.05, recoil.nbins))
    gammaCh.addSample(gammaJets)

    gammaCh.setObservation(expo_sample(2000, 450, recoil))

    with open(os.path.join(str(tmpdir), 'monojetModel.pkl'), "wb") as fout:
        pickle.dump(model, fout)

    model.renderCombine(os.path.join(str(tmpdir), 'monojetModel'))


if __name__ == '__main__':
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    darkhiggs_model('tmp')
