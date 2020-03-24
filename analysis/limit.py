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
    print(hist.values(overflow='all')[()])
    return (hist.values(overflow='all')[()], hist.axis(name).edges(overflow='all'), name)

def darkhiggs_model(tmpdir,mass,category):

    model = rl.Model('darkhiggs_'+mass+'_'+category)

    binning_map = {
        'mass0': {
            'monohs' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 550.0, 640.0, 740.0, 1250.0],
            'monojet' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 1250.0]
        },
        'mass1': {
            'monohs' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 1250.0],
            'monojet' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 1250.0]
        },
        'mass2': {
            'monohs' : [250.0, 280.0, 310.0, 340.0, 430.0, 1250.0],
            'monojet' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 1250.0]
        },
        'mass3': {
            'monohs' : [250.0, 280.0, 310.0, 340.0, 400.0, 430.0, 470.0, 1250.0],
            'monojet' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 640.0, 1250.0]
        },
        'mass4': {
            'monohs' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 1250.0],
            'monojet' : [250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 740.0, 900.0, 1250.0]
        }
    }
    
    process_map = OrderedDict()
    #process_map["Hbb_merged"] = ("Hbb_merged*",)
    #process_map["Hbb_unmerged"] = ("Hbb_unmerged*",)
    process_map["Hbb"] = ("Hbb*",)    
    process_map["DY"] = ("DY*",)
    #process_map["VVbb"] = ("VVbb*",)
    #process_map["VV"] = ("VV",)
    process_map["VV"] = ("VV*",) 
    #process_map["ST_merged"] = ("ST_merged*",)
    #process_map["ST_unmerged"] = ("ST_unmerged*",)
    process_map["ST"] = ("ST*",) 
    #process_map["TT_merged"] = ("TT_merged*",)
    #process_map["TT_unmerged"] = ("TT_unmerged*",)
    process_map["TT"] = ("TT*",)  
    process_map["WJets"] = ("WJets*",)
    process_map["ZJets"] = ("ZJets*",)
    process_map["GJets"] = ("GJets*",)

    wjets_eff = {
        'mass0': 0.092,
        'mass1': 0.054,
        'mass2': 0.044,
        'mass3': 0.045,
        'mass4': 0.042
    }

    zjets_eff = {
        'mass0': 0.092,
        'mass1': 0.054,
        'mass2': 0.050,
        'mass3': 0.054,
        'mass4': 0.054,
    }

    ###
    #Extract histograms from input file
    ###

    hists = load('hists/darkhiggs2018.scaled')
    recoil = {}
    for r in hists['recoil'].identifiers('region'):
        print(r)
        if category not in str(r) or mass not in str(r): continue
        print(r,category,mass)
        recoil[str(r).split("_")[0]]=hists['recoil'].integrate('region',r).rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning_map[mass][category]))
    print(recoil.keys())

    ###
    ###
    # Setting up rate systematics
    ###
    ###

    ###
    # Luminosity
    ###

    lumi = rl.NuisanceParameter('lumi', 'lnN')

    ###
    # MET bin migration
    ###

    met = rl.NuisanceParameter('met', 'lnN')

    ###
    # Cross section of MC-driven processes
    ###

    QCDe_Norm = rl.NuisanceParameter('QCDe_Norm', 'lnN')
    QCDmu_Norm = rl.NuisanceParameter('QCDmu_Norm', 'lnN')
    QCDsig_Norm = rl.NuisanceParameter('QCDsig_Norm', 'lnN')
    stop_Norm = rl.NuisanceParameter('stop_Norm', 'lnN')
    VV_Norm = rl.NuisanceParameter('VV_Norm', 'lnN')
    Hbb_Norm = rl.NuisanceParameter('Hbb_Norm', 'lnN')
    dy_Norm = rl.NuisanceParameter('dy_Norm', 'lnN') #only in signal region

    ###
    # Lepton/photon ID uncertainties 
    ###

    id_e = rl.NuisanceParameter('id_e', 'lnN')
    id_mu = rl.NuisanceParameter('id_mu', 'lnN')
    id_pho = rl.NuisanceParameter('id_pho', 'lnN')

    ###
    # Trigger efficiency
    ###

    trig_e = rl.NuisanceParameter('trig_e', 'lnN')
    trig_met = rl.NuisanceParameter('trig_met', 'lnN')

    ###
    # DeepAk15 signal scale factor and mistag rate for MC-driven processes
    ###

    sf_deepAK15 = rl.NuisanceParameter('sf_deepAK15', 'lnN')
    mistag_deepAK15 = rl.NuisanceParameter('mistag_deepAK15', 'lnN')

    ###
    # Tau veto
    ###

    veto_tau = rl.NuisanceParameter('veto_tau', 'lnN')


    ###
    # Setting tagger efficiency and scale factor for in-situ calculation
    ###

    wjets_sf = rl.IndependentParameter('wjets_sf', 1., 0, 10)
    wjets_weight = (1-(wjets_sf*wjets_eff[mass]))/(1-wjets_eff[mass])
    if 'monohs' in category: wjets_weight = wjets_sf

    zjets_sf = rl.IndependentParameter('zjets_sf', 1., 0, 10)
    zjets_weight = (1-(zjets_sf*zjets_eff[mass]))/(1-zjets_eff[mass])
    if 'monohs' in category: zjets_weight = zjets_sf

    ttbar_eff = 1.
    ttbar_sf = rl.IndependentParameter('ttbar_sf', 1., 0, 10)
    ttbar_weight = (1-(ttbar_sf*ttbar_eff))/(1-ttbar_eff)
    if 'monohs' in category: ttbar_weight = ttbar_sf

    ###
    ###
    # Shape systematics
    ###
    ###

    ###
    # JEC/JER
    ###
    
    jec = rl.NuisanceParameter('jec', 'shape')
    jer = rl.NuisanceParameter('jer', 'shape')
    btag = rl.NuisanceParameter('btag', 'shape') #AK4 btag
    gamma_to_z_ewk = rl.NuisanceParameter('Theory_gamma_z_ewk', 'shape')

    ###
    ###
    # Signal region
    ###
    ###

    sr = rl.Channel(category+'-sr')
    model.addChannel(sr)

    ###
    # Add data distribution to the channel
    ###

    sr.setObservation(template(recoil['sr'].integrate('process', 'MET'), 'recoil'))

    ###
    # Z(->nunu)+jets data-driven model
    ###

    sr_zvvHist = recoil['sr'].integrate('process', 'ZJets')
    sr_zvvTemplate = template(sr_zvvHist, 'recoil')
    sr_zvvMC =  rl.TemplateSample('sr_zvvMC', rl.Sample.BACKGROUND, sr_zvvTemplate)
    #sr_zvvMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_zvvHist.axis('recoil').edges(overflow='all'))-1))
    
    sr_zvvBinYields = np.array([rl.IndependentParameter(category+'-sr_zvv_bin_%d' % i, b, 0, sr_zvvTemplate[0].max()*2) for i,b in enumerate(sr_zvvTemplate[0])]) * zjets_weight
    sr_zvvObservable = rl.Observable('recoil', sr_zvvHist.axis('recoil').edges(overflow='all'))
    sr_zvv = rl.ParametericSample(category+'-sr_zvv', rl.Sample.BACKGROUND, sr_zvvObservable, sr_zvvBinYields)

    sr.addSample(sr_zvv)

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    sr_wjetsHist = recoil['sr'].integrate('process', 'Wjets')
    sr_wjetsTemplate = template(sr_wjetsHist, 'recoil')
    sr_wjetsMC =  rl.TemplateSample('sr_wjetsMC', rl.Sample.BACKGROUND, sr_wjetsTemplate)
    #sr_wjetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_wjetsHist.axis('recoil').edges(overflow='all'))-1))

    sr_wjetsBinYields = np.array([rl.IndependentParameter(category+'-sr_wjets_bin_%d' % i,b,0,sr_wjetsTemplate[0].max()*2) for i,b in enumerate(sr_wjetsTemplate[0])]) * wjets_weight
    sr_wjetsObservable = rl.Observable('recoil', sr_wjetsHist.axis('recoil').edges(overflow='all'))
    sr_wjets = rl.ParametericSample(category+'-sr_wjets', rl.Sample.BACKGROUND, sr_wjetsObservable, sr_wjetsBinYields)
    sr.addSample(sr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    sr_ttbarHist = recoil['sr'].integrate('process', 'TT')
    sr_ttbarTemplate = template(sr_ttbarHist, 'recoil')
    sr_ttbarMC =  rl.TemplateSample('sr_ttbarMC', rl.Sample.BACKGROUND, sr_ttbarTemplate)
    sr_ttbarMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_ttbarHist.axis('recoil').edges(overflow='all'))-1))

    # these parameters are large, should probably log-transform them
    sr_ttbarBinYields = np.array([rl.IndependentParameter(category+'-sr_ttbar_bin_%d' % i,b,0,sr_ttbarTemplate[0].max()*2) for i,b in enumerate(sr_ttbarTemplate[0])]) * ttbar_weight
    sr_ttbarObservable = rl.Observable('recoil', sr_ttbarHist.axis('recoil').edges(overflow='all'))
    sr_ttbar = rl.ParametericSample(category+'-sr_ttbar', rl.Sample.BACKGROUND, sr_ttbarObservable, sr_ttbarBinYields)
    sr.addSample(sr_ttbar)

    ###
    # Other MC-driven processes
    ###

    sr_singletopHist = recoil['sr'].integrate('process', 'ST')
    sr_singletopTemplate = template(sr_singletopHist, 'recoil')
    sr_singletop = rl.TemplateSample(category+'-sr_singletop', rl.Sample.BACKGROUND, sr_singletopTemplate)
    sr_singletop.setParamEffect(lumi, 1.027)
    sr_singletop.setParamEffect(stop_Norm, 1.2)
    sr_singletop.setParamEffect(trig_met, 1.01)
    sr_singletop.setParamEffect(veto_tau, 1.03)
    sr_singletop.setParamEffect(met, 1.05)
    #sr_singletop.setParamEffect(mistag_deepAK15, )
    
    sr.addSample(sr_singletop)

    sr_dyHist = recoil['sr'].integrate('process', 'DY')
    sr_dyTemplate = template(sr_dyHist, 'recoil')
    sr_dy = rl.TemplateSample(category+'-sr_dy', rl.Sample.BACKGROUND, sr_dyTemplate)
    sr_dy.setParamEffect(lumi, 1.027)
    sr_dy.setParamEffect(dy_Norm, 1.2)
    sr_dy.setParamEffect(trig_met, 1.01)
    sr_dy.setParamEffect(veto_tau, 1.03)
    sr_dy.setParamEffect(met, 1.05)
    sr.addSample(sr_dy)

    sr_dibosonHist = recoil['sr'].integrate('process', 'Diboson')
    sr_dibosonTemplate = template(sr_dibosonHist, 'recoil')
    sr_diboson = rl.TemplateSample(category+'-sr_diboson', rl.Sample.BACKGROUND, sr_dibosonTemplate)
    sr_diboson.setParamEffect(lumi, 1.027)
    sr_diboson.setParamEffect(VV_Norm, 1.2)
    sr_diboson.setParamEffect(trig_met, 1.01)
    sr_diboson.setParamEffect(veto_tau, 1.03)
    sr_diboson.setParamEffect(met, 1.05)
    sr.addSample(sr_diboson)

    sr_higgsHist = recoil['sr'].integrate('process', 'Hbb')
    sr_higgsTemplate = template(sr_higgsHist, 'recoil')
    sr_higgs = rl.TemplateSample(category+'-sr_higgs', rl.Sample.BACKGROUND, sr_higgsTemplate)
    sr_higgs.setParamEffect(lumi, 1.027)
    sr_higgs.setParamEffect(Hbb_Norm, 1.2)
    sr_higgs.setParamEffect(trig_met, 1.01)
    sr_higgs.setParamEffect(veto_tau, 1.03)
    sr_higgs.setParamEffect(met, 1.05)
    sr.addSample(sr_higgs)

    for signal in recoil['sr'].identifiers('process'):
        if 'Mono' not in str(signal): continue
        sr_dmHist = recoil['sr'].integrate('process', signal)
        sr_dmTemplate = template(sr_dmHist, 'recoil')
        sr_dm = rl.TemplateSample(category+'-sr_'+str(signal), rl.Sample.SIGNAL, sr_dmTemplate)
        sr_dm.setParamEffect(lumi, 1.027)
        sr_dm.setParamEffect(trig_met, 1.01)
        sr_dm.setParamEffect(veto_tau, 1.03)
        sr_dm.setParamEffect(met, 1.05)
        sr.addSample(sr_dm)

    ###
    # End of SR
    ###

    ###
    ###
    # Single Lepton Control Regions
    ###
    ###

    cr={}

    ttbarHist = {}
    ttbarTemplate = {}
    ttbarMC = {}
    ttbarTransferFactor = {}
    ttbar = {}

    wjetsHist = {}
    wjetsTemplate = {}
    wjetsMC = {}
    wjetsTransferFactor = {}
    wjets = {}

    singletopHist = {}
    singletopTemplate = {}
    singletop = {}

    dyHist = {}
    dyTemplate = {}
    dyMC = {}
    dyTransferFactor = {}
    dy = {}

    dibosonHist = {}
    dibosonTemplate = {}
    diboson = {}

    higgsHist = {}
    higgsTemplate = {}
    higgs = {}

    for p in ['t','w']:
        for l in ['e','m']:
            cr[p+l]=rl.Channel(category+'-'+p+l+'cr')
            model.addChannel(cr[p+l])
            print(p+l)
            print(recoil[p+l+'cr'].identifiers('process'))
            if 'e' in l: cr[p+l].setObservation(template(recoil[p+l+'cr'].integrate('process', 'SingleElectron'), 'recoil'))
            else: cr[p+l].setObservation(template(recoil[p+l+'cr'].integrate('process', 'MET'), 'recoil'))   


            ttbarHist[p+l] = recoil[p+l+'cr'].integrate('process', 'TT')
            ttbarTemplate[p+l] = template(ttbarHist[p+l], 'recoil')
            ttbarMC[p+l] =  rl.TemplateSample(p+l+'cr_ttbarMC', rl.Sample.BACKGROUND, ttbarTemplate[p+l])
            #ttbarMC[p+l].setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
            #ttbarMC[p+l].setParamEffect(ele_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins), np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

            ttbarTransferFactor[p+l] = ttbarMC[p+l].getExpectation() / sr_ttbarMC.getExpectation()
            ttbar[p+l] = rl.TransferFactorSample(category+'-'+p+l+'cr_ttbar', rl.Sample.BACKGROUND, ttbarTransferFactor[p+l], sr_ttbar)
            cr[p+l].addSample(ttbar[p+l])

            wjetsHist[p+l] = recoil[p+l+'cr'].integrate('process', 'Wjets')
            wjetsTemplate[p+l] = template(wjetsHist[p+l], 'recoil')
            wjetsMC[p+l] =  rl.TemplateSample(p+l+'cr_wjetsMC', rl.Sample.BACKGROUND, wjetsTemplate[p+l])
            #wjetsMC[p+l].setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
            #wjetsMC[p+l].setParamEffect(ele_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins), np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

            wjetsTransferFactor[p+l] = wjetsMC[p+l].getExpectation() / sr_wjetsMC.getExpectation()
            wjets[p+l] = rl.TransferFactorSample(category+'-'+p+l+'cr_wjets', rl.Sample.BACKGROUND, wjetsTransferFactor[p+l], sr_wjets)
            cr[p+l].addSample(wjets[p+l])

            singletopHist[p+l] = recoil[p+l+'cr'].integrate('process', 'ST')
            singletopTemplate[p+l] = template(singletopHist[p+l], 'recoil')
            singletop[p+l] = rl.TemplateSample(category+'-'+p+l+'cr_singletop', rl.Sample.BACKGROUND, singletopTemplate[p+l])
            cr[p+l].addSample(singletop[p+l])
            
            dyHist[p+l] = recoil[p+l+'cr'].integrate('process', 'DY')
            dyTemplate[p+l] = template(dyHist[p+l], 'recoil')
            dy[p+l] = rl.TemplateSample(category+'-'+p+l+'cr_dy', rl.Sample.BACKGROUND, dyTemplate[p+l])
            cr[p+l].addSample(dy[p+l])

            dibosonHist[p+l] = recoil[p+l+'cr'].integrate('process', 'Diboson')
            dibosonTemplate[p+l] = template(dibosonHist[p+l], 'recoil')
            diboson[p+l] = rl.TemplateSample(category+'-'+p+l+'cr_diboson', rl.Sample.BACKGROUND, dibosonTemplate[p+l])
            cr[p+l].addSample(diboson[p+l])

            higgsHist[p+l] = recoil[p+l+'cr'].integrate('process', 'Hbb')
            higgsTemplate[p+l] = template(higgsHist[p+l], 'recoil')
            higgs[p+l] = rl.TemplateSample(category+'-'+p+l+'cr_higgs', rl.Sample.BACKGROUND, higgsTemplate[p+l])
            cr[p+l].addSample(higgs[p+l])
    ###
    # End of Single Lepton CR
    ###

    ###
    ###
    # Double Lepton Control Regions
    ###
    ###

    for ll in ['ze','zm']:

        cr[ll] = rl.Channel(category+'-'+ll+'cr')
        model.addChannel(cr[ll])
        print(ll)
        print(recoil[ll+'cr'].identifiers('process'))
        if 'e' in ll: cr[ll].setObservation(template(recoil[ll+'cr'].integrate('process', 'SingleElectron'), 'recoil'))
        else: cr[ll].setObservation(template(recoil[ll+'cr'].integrate('process', 'MET'), 'recoil'))   
        
        dyHist[ll] = recoil[ll+'cr'].integrate('process', 'DY')
        dyTemplate[ll] = template(dyHist[ll], 'recoil')
        dyMC[ll] = rl.TemplateSample(ll+'cr_dy', rl.Sample.BACKGROUND, dyTemplate[ll])
        #zllJetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
        #zllJetsMC.setParamEffect(ele_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins), np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

        dyTransferFactor[ll] = dyMC[ll].getExpectation() / sr_zvvMC.getExpectation()
        dy[ll] = rl.TransferFactorSample(category+'-'+ll+'cr_dy', rl.Sample.BACKGROUND, dyTransferFactor[ll], sr_zvv)
        cr[ll].addSample(dy[ll])

        ttbarHist[ll] = recoil[ll+'cr'].integrate('process', 'TT')
        ttbarTemplate[ll] = template(ttbarHist[ll], 'recoil')
        ttbar[ll] =  rl.TemplateSample(category+'-'+ll+'cr_ttbarMC', rl.Sample.BACKGROUND, ttbarTemplate[ll])
        cr[ll].addSample(ttbar[ll])

        singletopHist[ll] = recoil[ll+'cr'].integrate('process', 'ST')
        singletopTemplate[ll] = template(singletopHist[ll], 'recoil')
        singletop[ll] = rl.TemplateSample(category+'-'+ll+'cr_singletop', rl.Sample.BACKGROUND, singletopTemplate[ll])
        cr[ll].addSample(singletop[ll])
        
        dibosonHist[ll] = recoil[ll+'cr'].integrate('process', 'Diboson')
        dibosonTemplate[ll] = template(dibosonHist[ll], 'recoil')
        diboson[ll] = rl.TemplateSample(category+'-'+ll+'cr_diboson', rl.Sample.BACKGROUND, dibosonTemplate[ll])
        cr[ll].addSample(diboson[ll])

        higgsHist[ll] = recoil[ll+'cr'].integrate('process', 'Hbb')
        higgsTemplate[ll] = template(higgsHist[ll], 'recoil')
        higgs[ll] = rl.TemplateSample(category+'-'+ll+'cr_higgs', rl.Sample.BACKGROUND, higgsTemplate[ll])
        cr[ll].addSample(higgs[ll])

    ###
    # End of Double Lepton CR
    ###

    ###
    ###
    # Single Photon Control Region
    ###
    ###

    gcr = rl.Channel(category+'-gcr')
    model.addChannel(gcr)

    gcr.setObservation(template(recoil['gcr'].integrate('process', 'SinglePhoton'), 'recoil'))

    gcr_gjetsHist = recoil['gcr'].integrate('process', 'Gjets')
    gcr_gjetsTemplate = template(gcr_gjetsHist, 'recoil')
    gcr_gjetsMC = rl.TemplateSample('gjetsMC', rl.Sample.BACKGROUND, gcr_gjetsTemplate)
    #gcr_gjetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
    #gcr_gjetsMC.setParamEffect(pho_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

    gcr_gjetsTransferFactor = gcr_gjetsMC.getExpectation() / sr_zvvMC.getExpectation()
    gcr_gjets = rl.TransferFactorSample(category+'-gcr_gjets', rl.Sample.BACKGROUND, gcr_gjetsTransferFactor, sr_zvv)
    #gammaJets.setParamEffect(gamma_to_z_ewk, np.linspace(1.01, 1.05, recoil.nbins))
    gcr.addSample(gcr_gjets)

    with open(os.path.join(str(tmpdir), 'darkhiggsModel.pkl'), "wb") as fout:
        pickle.dump(model, fout)

    model.renderCombine(os.path.join(str(tmpdir), 'darkhiggsModel/'+mass))


if __name__ == '__main__':
    if not os.path.exists('datacards'):
        os.mkdir('datacards')
    mass='mass0'
    category='monojet'
    darkhiggs_model('datacards',mass,category)
