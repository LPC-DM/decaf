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
from coffea.util import load

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

def expo_sample(norm, scale, obs):
    cdf = scipy.stats.expon.cdf(scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)

def template(hist, name):
    return (hist.values(overflow='all')[()], hist.axis(name).edges(overflow='all'), name)

def darkhiggs_model(tmpdir):

    hists = load('pods/hists_darkhiggs2018.coffea')
    signal_hists = hists['signal']
    bkg_hists    = hists['bkg']
    data_hists   = hists['data']

    model = rl.Model("DarkHiggs")

    # lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    jec = rl.NuisanceParameter('CMS_jec', 'shape')
    ele_id_eff = rl.NuisanceParameter('CMS_ele_id_eff', 'shape')
    pho_id_eff = rl.NuisanceParameter('CMS_pho_id_eff', 'shape')
    gamma_to_z_ewk = rl.NuisanceParameter('Theory_gamma_z_ewk', 'shape')

    signalCh = rl.Channel("signalCh")
    model.addChannel(signalCh)

    #    zvvTemplate = expo_sample(1000, 400, recoil)
    zvvHist = signal_hists['recoil'].integrate('jet_selection','baggy').integrate('region','iszeroL').integrate('process', 'Mhs_90')
    zvvTemplate = template(zvvHist, 'recoil')
    zvvJetsMC = rl.TemplateSample('zvvJetsMC', rl.Sample.BACKGROUND, zvvTemplate)
    zvvJetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(zvvHist.axis('recoil').edges(overflow='all'))-1))

    # these parameters are large, should probably log-transform them
    zvvBinYields = np.array([rl.IndependentParameter('tmp', b, 0, zvvTemplate[0].max()*2) for b in zvvTemplate[0]])  # name will be changed by ParametericSample
    recoil = rl.Observable('recoil', zvvHist.axis('recoil').edges(overflow='all'))
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
