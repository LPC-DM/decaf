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

def darkhiggs_model(tmpdir,mass,category,year):

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
    
    ###
    #Extract histograms from input file
    ###

    hists = load('hists/darkhiggs'+year+'.scaled')
    
    ###
    # Regrouping histograms
    ###
    
    process = hist.Cat("process", "Process", sorting='placement')
    cats = ("process",)
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
    process_map["MET"]   = ("MET*",)
    process_map["SingleElectron"]   = ("SingleElectron*",)
    process_map["SinglePhoton"]   = ("SinglePhoton*",)
    
    for key in hists.keys():
        hists[key] = hists[key].group(cats, process, process_map)

    ###
    # Preparing histograms for fit
    ##

    recoil = {}
    for r in hists['recoil'].identifiers('region'):
        if category not in str(r) or mass not in str(r): continue
        #print(r,category,mass)
        #print('Before rebin',hists['recoil'].integrate('region',r).values(overflow='all'))
        recoil[str(r).split("_")[0]]=hists['recoil'].integrate('region',r).rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning_map[mass][category]))
        #print('After rebin',recoil[str(r).split("_")[0]].values(overflow='all'))

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

    #met = rl.NuisanceParameter('met', 'lnN')

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
    # Electron reco
    ###

    reco_e = rl.NuisanceParameter('reco_e', 'lnN')

    ###
    # Muon isolation
    ###

    iso_m = rl.NuisanceParameter('reco_e', 'lnN')

    ###
    # Trigger efficiency
    ###

    trig_e = rl.NuisanceParameter('trig_e', 'lnN')
    trig_met = rl.NuisanceParameter('trig_met', 'lnN')

    ###
    # DeepAk15 signal scale factor and mistag rate for MC-driven processes
    ###

    #sf_deepAK15 = rl.NuisanceParameter('sf_deepAK15', 'lnN')
    #mistag_deepAK15 = rl.NuisanceParameter('mistag_deepAK15', 'lnN')

    ###
    # Tau veto
    ###

    veto_tau = rl.NuisanceParameter('veto_tau', 'lnN')

    ###
    # AK b-tagging of iso jet 0-tag efficiencies
    ###

    0tag_eff = {
        'whf': 0.86,
        'wlf': 0.90,
        'zhf': 0.80,
        'zlf': 0.90,
        'ttbqq': 1.,
        'ttqq': 1.,
        'ttother': 1.,
        'stbqq': 1.,
        'stqq':1.,
        'stother': 1.,
        'vvbb': 1,
        'vvqq': 1,
        'vvother': 1,
        'hbb': 1,
        'hother': 1
    }

    ###
    # Defining W/Z/gamma+jets heavy flavor fractions and their corrective k-factors
    ###

    whf_fraction = 0.18
    zhf_fraction = 0.09
    ghf_fraction = 0.12

    whf_k = rl.IndependentParameter('whf_k', 1., 0, 1/whf_fraction)
    zhf_k = rl.IndependentParameter('zhf_k', 1., 0, 1/zhf_fraction)
    ghf_k = rl.IndependentParameter('ghf_k', 1., 0, 1/ghf_fraction)

    ###
    # Taking into account the varying HF fraction to adjust the overall efficiency of ak4 btagging of iso jets
    ###

    whf_0tag_eff = 0.86
    wlf_0tag_eff = 0.90

    wj_0tag_eff = wlf_0tag_eff*(1 - whf_fraction) + whf_0tag_eff*whf_fraction
    wj_0tag_sfxeff = wlf_0tag_eff*(1 - whf_k*whf_fraction) + whf_0tag_eff*whf_k*whf_fraction

    wjets_0tag_weight = wj_0tag_sfxeff / wj_0tag_eff
    wjets_1tag_weight = (1 - wj_0tag_sfxeff) / (1 - wj_0tag_eff)

    zhf_0tag_eff = 0.80
    zlf_0tag_eff = 0.90

    zj_0tag_eff = zlf_0tag_eff*(1 - zhf_fraction) + zhf_0tag_eff*zhf_fraction
    zj_0tag_sfxeff = zlf_0tag_eff*(1 - zhf_k*zhf_fraction) + zhf_0tag_eff*zhf_k*zhf_fraction

    zjets_0tag_weight = zj_0tag_sfxeff / zj_0tag_eff


    ###
    # Setting tagger efficiency and scale factor for in-situ calculation
    ###

    whf_deepak15_eff = 0.1
    wlf_deepak15_eff = 0.04

    whf_deepak15_sf = rl.IndependentParameter('whf_deepak15_sf', 1., 0, 1/whf_deepak15_eff)
    wlf_deepak15_sf = rl.IndependentParameter('wlf_deepak15_sf', 1., 0, 1/wlf_deepak15_eff)

    wj_deepak15_sfxeff = wlf_deepak15_sf*wlf_deepak15_eff*(1-whf_k*whf_fraction) + whf_deepak15_sf*whf_deepak15_eff*whf_k*whf_fraction
    wj_deepak15_eff = wlf_deepak15_eff*(1-whf_fraction) + whf_deepak15_eff*whf_fraction

    wjets_deepak15_weight = (1 - wj_deepak15_sfxeff)/(1 - wj_deepak15_eff)
    if 'monohs' in category: wjets_deepak15_weight = wj_deepak15_sfxeff/wj_deepak15_eff

    zhf_deepak15_eff = 0.04
    zlf_deepak15_eff = 0.05

    zhf_deepak15_sf = rl.IndependentParameter('zhf_deepak15_sf', 1., 0, 1/zhf_deepak15_eff)
    zlf_deepak15_sf = rl.IndependentParameter('zlf_deepak15_sf', 1., 0, 1/zlf_deepak15_eff)

    zj_deepak15_sfxeff = zlf_deepak15_sf*zlf_deepak15_eff*(1-zhf_k*zhf_fraction) + zhf_deepak15_sf*zhf_deepak15_eff*zhf_k*zhf_fraction
    zj_deepak15_eff = zlf_deepak15_eff*(1-zhf_fraction) + zhf_deepak15_eff*zhf_fraction

    zjets_deepak15_weight = (1 - zj_deepak15_sfxeff)/(1 - zj_deepak15_eff)
    if 'monohs' in category: zjets_deepak15_weight = zj_deepak15_sfxeff/zj_deepak15_eff

    ghf_deepak15_eff = 0.03
    glf_deepak15_eff = 0.005

    ghf_deepak15_sf = rl.IndependentParameter('ghf_deepak15_sf', 1., 0, 1/ghf_deepak15_eff)
    glf_deepak15_sf = rl.IndependentParameter('glf_deepak15_sf', 1., 0, 1/glf_deepak15_eff)

    gj_deepak15_sfxeff = glf_deepak15_sf*glf_deepak15_eff*(1-ghf_k*ghf_fraction) + ghf_deepak15_sf*ghf_deepak15_eff*ghf_k*ghf_fraction
    gj_deepak15_eff = glf_deepak15_eff*(1-ghf_fraction) + ghf_deepak15_eff*ghf_fraction

    gjets_deepak15_weight = (1 - gj_deepak15_sfxeff)/(1 - gj_deepak15_eff)
    if 'monohs' in category: gjets_deepak15_weight = gj_deepak15_sfxeff/gj_deepak15_eff

    bqq_eff = 0.6
    qq_eff = 0.3
    bb_eff = 0.9
    other_eff = 0.3

    bqq_sf = rl.IndependentParameter('bqq_sf', 1., 0, 1/bqq_eff)
    qq_sf = rl.IndependentParameter('qq_sf', 1., 0, 1/qq_eff)
    bb_sf = rl.IndependentParameter('qq_sf', 1., 0, 1/bb_eff)
    other_sf = rl.IndependentParameter('other_sf', 1., 0, 1/other_eff)

    tt_bqq_fraction = {
        '0tag': {
            'mass0': 0.04,
            'mass1': 0.06,
            'mass2': 0.11,
            'mass3': 0.19,
            'mass4': 0.6
        },
        '1tag': {
            'mass0': 0.014,
            'mass1': 0.04,
            'mass2': 0.1,
            'mass3': 0.13,
            'mass4': 0.54
        }
    }

    tt_qq_fraction = {
        '0tag': {
            'mass0': 0.04,
            'mass1': 0.06,
            'mass2': 0.11,
            'mass3': 0.19,
            'mass4': 0.6
        },
        '1tag': {
            'mass0': 0.014,
            'mass1': 0.04,
            'mass2': 0.1,
            'mass3': 0.13,
            'mass4': 0.54
        }
    }

    tt_0tag_sfxeff = bqq_sf*bqq_eff*tt_bqq_fraction['0tag'][mass] + qq_sf*qq_eff*tt_qq_fraction['0tag'][mass] + other_sf*other_eff*(1 - tt_bqq_fraction['0tag'][mass] - tt_qq_fraction['0tag'][mass])
    tt_0tag_eff = bqq_eff*tt_bqq_fraction['0tag'][mass] + qq_eff*tt_qq_fraction['0tag'][mass] + other_eff*(1 - tt_bqq_fraction['0tag'][mass] - tt_qq_fraction['0tag'][mass])
    tt_1tag_sfxeff = bqq_sf*bqq_eff*tt_bqq_fraction['1tag'][mass] + qq_sf*qq_eff*tt_qq_fraction['1tag'][mass] + other_sf*other_eff*(1 - tt_bqq_fraction['1tag'][mass] - tt_qq_fraction['1tag'][mass])
    tt_1tag_eff =  bqq_eff*tt_bqq_fraction['1tag'][mass] + qq_eff*tt_qq_fraction['1tag'][mass] + other_eff*(1 - tt_bqq_fraction['1tag'][mass] - tt_qq_fraction['1tag'][mass])

    tt_0tag_weight = (1 - tt_0tag_sfxeff)/(1 - tt_0tag_eff)
    if 'monohs' in category: tt_0tag_weight = tt_0tag_sfxeff / tt_0tag_eff
    tt_1tag_weight = (1 - tt_1tag_sfxeff)/(1 - tt_1tag_eff)
    if 'monohs' in category: tt_1tag_weight = tt_1tag_sfxeff / tt_1tag_eff

    ###
    ###
    # Shape systematics
    ###
    ###

    ###
    # JEC/JER
    ###
    
    #jec = rl.NuisanceParameter('jec', 'shape')
    #jer = rl.NuisanceParameter('jer', 'shape')
    btag = rl.NuisanceParameter('btag', 'shape') #AK4 btag
    gamma_to_z_ewk = rl.NuisanceParameter('Theory_gamma_z_ewk', 'shape')

    ###
    ###
    # Signal region
    ###
    ###

    ch_name = 'sr-'+mass+'-'+category
    sr = rl.Channel(ch_name)
    model.addChannel(sr)

    ###
    # Add data distribution to the channel
    ###

    sr.setObservation(template(recoil['sr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))

    ###
    # Z(->nunu)+jets data-driven model
    ###

    sr_zvvHist = recoil['sr'].integrate('process', 'ZJets').integrate('systematic','nominal')
    sr_zvvTemplate = template(sr_zvvHist, 'recoil')
    sr_zvvMC =  rl.TemplateSample(ch_name+'_zvvMC', rl.Sample.BACKGROUND, sr_zvvTemplate)
    #sr_zvvMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_zvvHist.axis('recoil').edges(overflow='all'))-1))
    
    sr_zvvBinYields = np.array([rl.IndependentParameter(ch_name+'_zvv_bin_%d' % i, b, 0, sr_zvvTemplate[0].max()*2) for i,b in enumerate(sr_zvvTemplate[0])]) 
    sr_zvvBinYields = sr_zvvBinYields * zjets_deepak15_weight * zjets_0tag_weight
    sr_zvvObservable = rl.Observable('recoil', sr_zvvHist.axis('recoil').edges(overflow='all'))
    sr_zvv = rl.ParametericSample(ch_name+'_zvv', rl.Sample.BACKGROUND, sr_zvvObservable, sr_zvvBinYields)

    sr.addSample(sr_zvv)

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    sr_wjetsHist = recoil['sr'].integrate('process', 'WJets').integrate('systematic','nominal')
    sr_wjetsTemplate = template(sr_wjetsHist, 'recoil')
    sr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, sr_wjetsTemplate)
    #sr_wjetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_wjetsHist.axis('recoil').edges(overflow='all'))-1))

    sr_wjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_wjets_bin_%d' % i,b,0,sr_wjetsTemplate[0].max()*2) for i,b in enumerate(sr_wjetsTemplate[0])]) 
    sr_wjetsBinYields = sr_wjetsBinYields * wjets_deepak15_weight * wjets_0tag_weight
    sr_wjetsObservable = rl.Observable('recoil', sr_wjetsHist.axis('recoil').edges(overflow='all'))
    sr_wjets = rl.ParametericSample(ch_name+'_wjets', rl.Sample.BACKGROUND, sr_wjetsObservable, sr_wjetsBinYields)
    sr.addSample(sr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    sr_ttbarHist = recoil['sr'].integrate('process', 'TT').integrate('systematic','nominal')
    sr_ttbarTemplate = template(sr_ttbarHist, 'recoil')
    sr_ttbarMC =  rl.TemplateSample(ch_name+'_ttbarMC', rl.Sample.BACKGROUND, sr_ttbarTemplate)
    #sr_ttbarMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_ttbarHist.axis('recoil').edges(overflow='all'))-1))

    # these parameters are large, should probably log-transform them
    sr_ttbarBinYields = np.array([rl.IndependentParameter(ch_name+'_ttbar_bin_%d' % i,b,0,sr_ttbarTemplate[0].max()*2) for i,b in enumerate(sr_ttbarTemplate[0])]) * tt_0tag_weight
    sr_ttbarObservable = rl.Observable('recoil', sr_ttbarHist.axis('recoil').edges(overflow='all'))
    sr_ttbar = rl.ParametericSample(ch_name+'_ttbar', rl.Sample.BACKGROUND, sr_ttbarObservable, sr_ttbarBinYields)
    sr.addSample(sr_ttbar)

    ###
    # Other MC-driven processes
    ###

    sr_singletopHist = recoil['sr'].integrate('process', 'ST').integrate('systematic','nominal')
    sr_singletopTemplate = template(sr_singletopHist, 'recoil')
    sr_singletop = rl.TemplateSample(ch_name+'_singletop', rl.Sample.BACKGROUND, sr_singletopTemplate)
    sr_singletop.setParamEffect(lumi, 1.027)
    sr_singletop.setParamEffect(stop_Norm, 1.2)
    sr_singletop.setParamEffect(trig_met, 1.01)
    sr_singletop.setParamEffect(veto_tau, 1.03)
    #sr_singletop.setParamEffect(met, 1.05)
    sr.addSample(sr_singletop)

    sr_dyHist = recoil['sr'].integrate('process', 'DY').integrate('systematic','nominal')
    sr_dyTemplate = template(sr_dyHist, 'recoil')
    sr_dy = rl.TemplateSample(ch_name+'_dy', rl.Sample.BACKGROUND, sr_dyTemplate)
    sr_dy.setParamEffect(lumi, 1.027)
    sr_dy.setParamEffect(dy_Norm, 1.2)
    sr_dy.setParamEffect(trig_met, 1.01)
    sr_dy.setParamEffect(veto_tau, 1.03)
    #sr_dy.setParamEffect(met, 1.05)
    sr.addSample(sr_dy)

    sr_dibosonHist = recoil['sr'].integrate('process', 'VV').integrate('systematic','nominal')
    sr_dibosonTemplate = template(sr_dibosonHist, 'recoil')
    sr_diboson = rl.TemplateSample(ch_name+'_diboson', rl.Sample.BACKGROUND, sr_dibosonTemplate)
    sr_diboson.setParamEffect(lumi, 1.027)
    sr_diboson.setParamEffect(VV_Norm, 1.2)
    sr_diboson.setParamEffect(trig_met, 1.01)
    sr_diboson.setParamEffect(veto_tau, 1.03)
    #sr_diboson.setParamEffect(met, 1.05)
    sr.addSample(sr_diboson)

    sr_higgsHist = recoil['sr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    sr_higgsTemplate = template(sr_higgsHist, 'recoil')
    sr_higgs = rl.TemplateSample(ch_name+'_higgs', rl.Sample.BACKGROUND, sr_higgsTemplate)
    sr_higgs.setParamEffect(lumi, 1.027)
    sr_higgs.setParamEffect(Hbb_Norm, 1.2)
    sr_higgs.setParamEffect(trig_met, 1.01)
    sr_higgs.setParamEffect(veto_tau, 1.03)
    #sr_higgs.setParamEffect(met, 1.05)
    sr.addSample(sr_higgs)

    for signal in recoil['sr'].identifiers('process'):
        if 'Mono' not in str(signal): continue
        sr_dmHist = recoil['sr'].integrate('process', signal).integrate('systematic','nominal')
        sr_dmTemplate = template(sr_dmHist, 'recoil')
        sr_dm = rl.TemplateSample(ch_name+'_'+str(signal), rl.Sample.SIGNAL, sr_dmTemplate)
        sr_dm.setParamEffect(lumi, 1.027)
        sr_dm.setParamEffect(trig_met, 1.01)
        sr_dm.setParamEffect(veto_tau, 1.03)
        #sr_dm.setParamEffect(met, 1.05)
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
            ch_name = p+l+'cr-'+mass+'-'+category
            cr[p+l]=rl.Channel(ch_name)
            model.addChannel(cr[p+l])
            if 'e' in l: cr[p+l].setObservation(template(recoil[p+l+'cr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))
            else: cr[p+l].setObservation(template(recoil[p+l+'cr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))   


            ttbarHist[p+l] = recoil[p+l+'cr'].integrate('process', 'TT').integrate('systematic','nominal')
            ttbarTemplate[p+l] = template(ttbarHist[p+l], 'recoil')
            ttbarMC[p+l] =  rl.TemplateSample(ch_name+'_ttbarMC', rl.Sample.BACKGROUND, ttbarTemplate[p+l])
            #ttbarMC[p+l].setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
            #ttbarMC[p+l].setParamEffect(ele_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins), np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

            ttbarTransferFactor[p+l] = ttbarMC[p+l].getExpectation() / sr_ttbarMC.getExpectation()
            ttbar[p+l] = rl.TransferFactorSample(ch_name+'_ttbar', rl.Sample.BACKGROUND, ttbarTransferFactor[p+l], sr_ttbar)
            cr[p+l].addSample(ttbar[p+l])

            wjetsHist[p+l] = recoil[p+l+'cr'].integrate('process', 'WJets').integrate('systematic','nominal')
            wjetsTemplate[p+l] = template(wjetsHist[p+l], 'recoil')
            wjetsMC[p+l] =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, wjetsTemplate[p+l])
            #wjetsMC[p+l].setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
            #wjetsMC[p+l].setParamEffect(ele_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins), np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

            wjetsTransferFactor[p+l] = wjetsMC[p+l].getExpectation() / sr_wjetsMC.getExpectation()
            wjets[p+l] = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, wjetsTransferFactor[p+l], sr_wjets)
            cr[p+l].addSample(wjets[p+l])

            singletopHist[p+l] = recoil[p+l+'cr'].integrate('process', 'ST').integrate('systematic','nominal')
            singletopTemplate[p+l] = template(singletopHist[p+l], 'recoil')
            singletop[p+l] = rl.TemplateSample(ch_name+'_singletop', rl.Sample.BACKGROUND, singletopTemplate[p+l])
            cr[p+l].addSample(singletop[p+l])
            
            dyHist[p+l] = recoil[p+l+'cr'].integrate('process', 'DY').integrate('systematic','nominal')
            dyTemplate[p+l] = template(dyHist[p+l], 'recoil')
            dy[p+l] = rl.TemplateSample(ch_name+'_dy', rl.Sample.BACKGROUND, dyTemplate[p+l])
            cr[p+l].addSample(dy[p+l])

            dibosonHist[p+l] = recoil[p+l+'cr'].integrate('process', 'VV').integrate('systematic','nominal')
            dibosonTemplate[p+l] = template(dibosonHist[p+l], 'recoil')
            diboson[p+l] = rl.TemplateSample(ch_name+'_diboson', rl.Sample.BACKGROUND, dibosonTemplate[p+l])
            cr[p+l].addSample(diboson[p+l])

            higgsHist[p+l] = recoil[p+l+'cr'].integrate('process', 'Hbb').integrate('systematic','nominal')
            higgsTemplate[p+l] = template(higgsHist[p+l], 'recoil')
            higgs[p+l] = rl.TemplateSample(ch_name+'_higgs', rl.Sample.BACKGROUND, higgsTemplate[p+l])
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

        ch_name = ll+'cr-'+mass+'-'+category
        cr[ll] = rl.Channel(ch_name)
        model.addChannel(cr[ll])
        if 'e' in ll: cr[ll].setObservation(template(recoil[ll+'cr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))
        else: cr[ll].setObservation(template(recoil[ll+'cr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))   
        
        dyHist[ll] = recoil[ll+'cr'].integrate('process', 'DY').integrate('systematic','nominal')
        dyTemplate[ll] = template(dyHist[ll], 'recoil')
        dyMC[ll] = rl.TemplateSample(ch_name+'_dyMC', rl.Sample.BACKGROUND, dyTemplate[ll])
        #zllJetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
        #zllJetsMC.setParamEffect(ele_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins), np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

        dyTransferFactor[ll] = dyMC[ll].getExpectation() / sr_zvvMC.getExpectation()
        dy[ll] = rl.TransferFactorSample(ch_name+'_dy', rl.Sample.BACKGROUND, dyTransferFactor[ll], sr_zvv)
        cr[ll].addSample(dy[ll])

        ttbarHist[ll] = recoil[ll+'cr'].integrate('process', 'TT').integrate('systematic','nominal')
        ttbarTemplate[ll] = template(ttbarHist[ll], 'recoil')
        ttbar[ll] =  rl.TemplateSample(ch_name+'_ttbar', rl.Sample.BACKGROUND, ttbarTemplate[ll])
        cr[ll].addSample(ttbar[ll])

        singletopHist[ll] = recoil[ll+'cr'].integrate('process', 'ST').integrate('systematic','nominal')
        singletopTemplate[ll] = template(singletopHist[ll], 'recoil')
        singletop[ll] = rl.TemplateSample(ch_name+'_singletop', rl.Sample.BACKGROUND, singletopTemplate[ll])
        cr[ll].addSample(singletop[ll])
        
        dibosonHist[ll] = recoil[ll+'cr'].integrate('process', 'VV').integrate('systematic','nominal')
        dibosonTemplate[ll] = template(dibosonHist[ll], 'recoil')
        diboson[ll] = rl.TemplateSample(ch_name+'_diboson', rl.Sample.BACKGROUND, dibosonTemplate[ll])
        cr[ll].addSample(diboson[ll])

        higgsHist[ll] = recoil[ll+'cr'].integrate('process', 'Hbb').integrate('systematic','nominal')
        higgsTemplate[ll] = template(higgsHist[ll], 'recoil')
        higgs[ll] = rl.TemplateSample(ch_name+'_higgs', rl.Sample.BACKGROUND, higgsTemplate[ll])
        cr[ll].addSample(higgs[ll])

    ###
    # End of Double Lepton CR
    ###

    ###
    ###
    # Single Photon Control Region
    ###
    ###

    ch_name = 'gcr-'+mass+'-'+category
    gcr = rl.Channel(ch_name)
    model.addChannel(gcr)

    gcr.setObservation(template(recoil['gcr'].integrate('process', 'SinglePhoton').integrate('systematic','nominal'), 'recoil'))

    gcr_gjetsHist = recoil['gcr'].integrate('process', 'GJets').integrate('systematic','nominal')
    gcr_gjetsTemplate = template(gcr_gjetsHist, 'recoil')
    gcr_gjetsMC = rl.TemplateSample(ch_name+'_gjetsMC', rl.Sample.BACKGROUND, gcr_gjetsTemplate)
    #gcr_gjetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.05, size=recoil.nbins))
    #gcr_gjetsMC.setParamEffect(pho_id_eff, np.random.normal(loc=1, scale=0.02, size=recoil.nbins))

    gcr_gjetsTransferFactor = gcr_gjetsMC.getExpectation() / sr_zvvMC.getExpectation()
    gcr_gjets = rl.TransferFactorSample(ch_name+'_gjets', rl.Sample.BACKGROUND, gcr_gjetsTransferFactor, sr_zvv)
    #gammaJets.setParamEffect(gamma_to_z_ewk, np.linspace(1.01, 1.05, recoil.nbins))
    gcr.addSample(gcr_gjets)

    with open(os.path.join(str(tmpdir), 'darkhiggsModel'+year+'.pkl'), "wb") as fout:
        pickle.dump(model, fout)

    model.renderCombine(os.path.join(str(tmpdir), 'darkhiggsModel'+year+'/'+mass))


if __name__ == '__main__':
    if not os.path.exists('datacards'):
        os.mkdir('datacards')
    mass='mass0'
    category='monojet'
    year='2018'
    darkhiggs_model('datacards',mass,category,year)
