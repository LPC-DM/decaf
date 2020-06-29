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

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

def model(year,mass,category):
    
    def template(dictionary, process, systematic, region):
        print('Generating template for',process,'in',region)
        output=dictionary[region].integrate('process', process).integrate('systematic',systematic).sum('fjmass').values()[()][:,0]
        if 'monohs' in category:
            output=dictionary[region].integrate('process', process).integrate('systematic',systematic).values()[()][:,mass,1]
        binning=dictionary[region].integrate('process', process).integrate('systematic',systematic).axis('recoil').edges()
        return (output, binning, 'recoil')

    model_id=year+category
    if mass is not None: model_id=year+'mass'+str(mass)+category
    print(model_id)
    model = rl.Model('darkhiggs'+model_id)

    whf_fraction=0.18
    zhf_fraction=0.10
    ghf_fraction=0.12

    whf_k = rl.IndependentParameter('whf_k', 1., 0.6, 1.4)
    zhf_k = rl.IndependentParameter('zhf_k', 1., 0.6, 1.4)
    ghf_k = rl.IndependentParameter('ghf_k', 1., 0.6, 1.4)


    gentypes = {
        'Hbb': ['xbb','vqq','wcq','b','bb','bc','garbage','other','tbcq','tbqq'],
        'DY+HF': ['b','bb','c','cc','garbage','other'],
        'DY+LF': ['garbage','other'],
        'G+HF': ['b','bb','c','cc','garbage','other'],
        'G+LF': ['garbage','other'],
        'VV': ['xbb','vqq','wcq','b','c','garbage','other','zcc'],
        'ST': ['vqq','wcq','b','bb','bc','c','garbage','other','tbcq','tbqq'],
        'TT': ['vqq','wcq','b','bb','bc','c','garbage','other','tbcq','tbqq'],
        'W+HF': ['b','bb','c','cc','garbage','other'],
        'W+LF': ['garbage','other'],
        'QCD': ['b','bb','c','cc','garbage','other'],
        'Mhs': ['xbb','b','bb','other'],
        'MonoJet': ['c','cc','other'],
        'MonoW': ['vqq','wcq','c','other'],
        'MonoZ': ['vqq','c','cc','other','zcc']
    }

    with open('data/'+year+'_deepak15_pass_eff.json') as fin:
        deepak15_pass_eff = json.load(fin)

    deepak15_pass_sf = {
        'xbb'    : rl.IndependentParameter('xbb_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['xbb']), 
        'vqq'    : rl.IndependentParameter('vqq_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['vqq']),
        'wcq'    : rl.IndependentParameter('wcq_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['wcq']),
        'b'      : rl.IndependentParameter('b_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['b']),
        'bb'     : rl.IndependentParameter('bb_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['bb']), 
        'bc'     : rl.IndependentParameter('bc_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['bc']),
        'c'      : rl.IndependentParameter('c_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['c']),
        'cc'     : rl.IndependentParameter('cc_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['cc']),
        'garbage': 1.,
        'other'  : rl.IndependentParameter('other_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['other']),
        'tbcq'   : rl.IndependentParameter('tbcq_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['tbcq']),
        'tbqq'   : rl.IndependentParameter('tbqq_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['tbqq']),
        'zcc'    : rl.IndependentParameter('zcc_sf_'+year, 1., 0.01, 1/deepak15_pass_eff['zcc'])
    }
    
    with open('data/'+year+'_deepak4_0tag_eff.json') as fin:
        deepak4_0tag_gentype_eff = json.load(fin)

    with open('data/'+year+'_deepak4_0tag_process_eff.json') as fin:
        deepak4_0tag_process_eff = json.load(fin)

    with open('data/'+year+'_abs_signal_fractions.json') as fin:
        abs_signal_fractions = json.load(fin)

    with open('data/'+year+'_mass_signal_fractions_modulation.json') as fin:
        signal_mass_modulation = json.load(fin)

    signal_weight={}
    for process in abs_signal_fractions.keys():
        signal_weight[process]=0.
        for gentype in abs_signal_fractions[process].keys():
            if gentype not in gentypes[process.split('_')[0]]: continue
            print('Extracting',gentype,'fraction for',process)
            fraction=abs_signal_fractions[process][gentype]
            if mass is not None: fraction*=signal_mass_modulation[process][gentype]['mass'+str(mass)]
            if 'monohs' in category:
                signal_weight[process]+=deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*fraction
            elif 'monojet' in category:
                signal_weight[process]+=(1-deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*fraction
                
    with open('data/'+year+'_abs_bkg_fractions_tot.json') as fin:
        abs_fractions = json.load(fin)

    with open('data/'+year+'_mass_bkg_fractions_modulation_tot.json') as fin:
        mass_modulation = json.load(fin)

    deepak15_weight={}
    deepak15_weight['0tag']={}
    deepak15_weight['1tag']={}
    deepak15_weight['notag']={}
    for process in ['Hbb','VV','ST','QCD','TT','DY+HF','DY+LF','W+HF','W+LF']:
        deepak15_weight['0tag'][process]=0.
        deepak15_weight['1tag'][process]=0.
        deepak15_weight['notag'][process]=0.
        for gentype in abs_fractions[process].keys():
            if gentype not in gentypes[process]: continue
            print('Extracting',gentype,'fraction for',process )
            fraction=abs_fractions[process][gentype]
            if mass is not None: fraction*=mass_modulation[process][gentype]['mass'+str(mass)]
            if 'monohs' in category:
                deepak15_weight['0tag'][process] += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*deepak4_0tag_gentype_eff[process][gentype]*fraction
                deepak15_weight['1tag'][process] += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*(1 - deepak4_0tag_gentype_eff[process][gentype])*fraction
                deepak15_weight['notag'][process] += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*fraction
            elif 'monojet' in category:
                deepak15_weight['0tag'][process] += (1-deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*deepak4_0tag_gentype_eff[process][gentype]*fraction
                deepak15_weight['1tag'][process] += (1-deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*(1 - deepak4_0tag_gentype_eff[process][gentype])*fraction
                deepak15_weight['notag'][process] += (1-deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*fraction
        deepak15_weight['0tag'][process]=np.nan_to_num(deepak15_weight['0tag'][process]/deepak4_0tag_process_eff[process])
        deepak15_weight['1tag'][process]=np.nan_to_num(deepak15_weight['1tag'][process]/(1 - deepak4_0tag_process_eff[process]))

    for process in ['G+HF','G+LF']:
        deepak15_weight['notag'][process]=0.
        for gentype in abs_fractions[process].keys():
            if gentype not in gentypes[process]: continue
            print('Extracting',gentype,'fraction for',process )
            fraction=abs_fractions[process][gentype]
            if mass is not None: fraction*=mass_modulation[process][gentype]['mass'+str(mass)]
            if 'monohs' in category:
                deepak15_weight['notag'][process] += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*fraction
            elif 'monojet' in category:
                deepak15_weight['notag'][process] += (1-deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*fraction

    hf_fraction_weight={}
    hf_fraction_weight['0tag']={}
    hf_fraction_weight['1tag']={}
    hf_fraction_weight['notag']={}
    hf_fraction_weight['0tag']['W+jets'] = np.nan_to_num(deepak15_weight['0tag']['W+HF']*(deepak4_0tag_process_eff['W+HF']/deepak4_0tag_process_eff['W+jets'])*whf_k*whf_fraction)
    hf_fraction_weight['0tag']['W+jets'] += np.nan_to_num(deepak15_weight['0tag']['W+LF']*(deepak4_0tag_process_eff['W+LF']/deepak4_0tag_process_eff['W+jets'])*(1 - whf_k*whf_fraction))
    hf_fraction_weight['1tag']['W+jets'] = np.nan_to_num(deepak15_weight['1tag']['W+HF']*((1-deepak4_0tag_process_eff['W+HF'])/(1-deepak4_0tag_process_eff['W+jets']))*whf_k*whf_fraction) 
    hf_fraction_weight['1tag']['W+jets'] += np.nan_to_num(deepak15_weight['1tag']['W+LF']*((1-deepak4_0tag_process_eff['W+LF'])/(1-deepak4_0tag_process_eff['W+jets']))*(1 - whf_k*whf_fraction))
    hf_fraction_weight['notag']['W+jets'] = deepak15_weight['notag']['W+HF']*whf_k*whf_fraction
    hf_fraction_weight['notag']['W+jets'] += deepak15_weight['notag']['W+LF']*(1 - whf_k*whf_fraction)

    hf_fraction_weight['0tag']['DY+jets'] = np.nan_to_num(deepak15_weight['0tag']['DY+HF']*(deepak4_0tag_process_eff['DY+HF']/deepak4_0tag_process_eff['DY+jets'])*zhf_k*zhf_fraction)
    hf_fraction_weight['0tag']['DY+jets'] += np.nan_to_num(deepak15_weight['0tag']['DY+LF']*(deepak4_0tag_process_eff['DY+LF']/deepak4_0tag_process_eff['DY+jets'])*(1 - zhf_k*zhf_fraction))
    hf_fraction_weight['1tag']['DY+jets'] = np.nan_to_num(deepak15_weight['1tag']['DY+HF']*((1-deepak4_0tag_process_eff['DY+HF'])/(1-deepak4_0tag_process_eff['DY+jets']))*zhf_k*zhf_fraction) 
    hf_fraction_weight['1tag']['DY+jets'] += np.nan_to_num(deepak15_weight['1tag']['DY+LF']*((1-deepak4_0tag_process_eff['DY+LF'])/(1-deepak4_0tag_process_eff['DY+jets']))*(1 - zhf_k*zhf_fraction))
    hf_fraction_weight['notag']['DY+jets'] = deepak15_weight['notag']['DY+HF']*zhf_k*zhf_fraction
    hf_fraction_weight['notag']['DY+jets'] += deepak15_weight['notag']['DY+LF']*(1 - zhf_k*zhf_fraction)

    hf_fraction_weight['notag']['G+jets'] = deepak15_weight['notag']['G+HF']*ghf_k*ghf_fraction
    hf_fraction_weight['notag']['G+jets'] += deepak15_weight['notag']['G+LF']*(1 - ghf_k*ghf_fraction)

    data_hists   = hists['data']
    bkg_hists    = hists['bkg']
    signal_hists = hists['sig']

    '''
    binning = {
        'mass0': data_hists['recoil'].integrate('region','sr').integrate('systematic','nominal').integrate('process','MET').sum('gentype').axis('recoil').edges(),
        'mass1': data_hists['recoil'].integrate('region','sr').integrate('systematic','nominal').integrate('process','MET').sum('gentype').axis('recoil').edges(),
        'mass2': data_hists['recoil'].integrate('region','sr').integrate('systematic','nominal').integrate('process','MET').sum('gentype').axis('recoil').edges(),
        'mass3': data_hists['recoil'].integrate('region','sr').integrate('systematic','nominal').integrate('process','MET').sum('gentype').axis('recoil').edges(),
        'mass4': data_hists['recoil'].integrate('region','sr').integrate('systematic','nominal').integrate('process','MET').sum('gentype').axis('recoil').edges(),
        'nomass': data_hists['recoil'].integrate('region','sr').integrate('systematic','nominal').integrate('process','MET').sum('gentype').axis('recoil').edges()
    }    
    '''

    ###
    # Preparing histograms for fit
    ##

    data = {}
    for r in data_hists['template'].identifiers('region'):
        data[str(r)]=data_hists['template'].integrate('region',r).sum('gentype')

    background = {}
    for r in bkg_hists['template'].identifiers('region'):
        background[str(r)]=bkg_hists['template'].integrate('region',r).sum('gentype')

    signal = {}
    for r in bkg_hists['template'].identifiers('region'):
        signal[str(r)]=signal_hists['template'].integrate('region',r).sum('gentype')

    ###
    ###
    # Setting up systematics
    ###
    ###

    lumi = rl.NuisanceParameter('lumi_'+year, 'lnN')
    qcdpho_norm = rl.NuisanceParameter('qcdpho_norm', 'lnN')
    qcde_norm = rl.NuisanceParameter('qcde_norm', 'lnN')
    qcdmu_norm = rl.NuisanceParameter('qcdmu_norm', 'lnN')
    qcdsig_norm = rl.NuisanceParameter('qcdsig_norm', 'lnN')
    st_norm = rl.NuisanceParameter('st_norm', 'lnN')
    tt_norm = rl.NuisanceParameter('tt_norm', 'lnN')
    vv_norm = rl.NuisanceParameter('vv_norm', 'lnN')
    hbb_norm = rl.NuisanceParameter('hbb_norm', 'lnN')
    dyjets_norm = rl.NuisanceParameter('dyjets_norm', 'lnN') 
    zjets_norm = rl.NuisanceParameter('zjets_norm', 'lnN')
    wjets_norm = rl.NuisanceParameter('wjets_norm', 'lnN')
    gjets_norm = rl.NuisanceParameter('gjets_norm', 'lnN')
    id_e = rl.NuisanceParameter('id_e_'+year, 'lnN')
    id_mu = rl.NuisanceParameter('id_mu_'+year, 'lnN')
    id_pho = rl.NuisanceParameter('id_pho_'+year, 'lnN')
    reco_e = rl.NuisanceParameter('reco_e_'+year, 'lnN')
    iso_mu = rl.NuisanceParameter('iso_mu_'+year, 'lnN')
    trig_e = rl.NuisanceParameter('trig_e_'+year, 'lnN')
    trig_met = rl.NuisanceParameter('trig_met_'+year, 'lnN')
    trig_pho = rl.NuisanceParameter('trig_pho_'+year, 'lnN')
    veto_tau = rl.NuisanceParameter('veto_tau_'+year, 'lnN')
    jec = rl.NuisanceParameter('jec_'+year, 'lnN')
    btag = rl.NuisanceParameter('btag_'+year, 'shape') #AK4 btag

    ###
    ###
    # Signal region
    ###
    ###

    ch_name = 'sr'+model_id
    sr = rl.Channel(ch_name)
    model.addChannel(sr)

    ###
    # Add data distribution to the channel
    ###

    sr.setObservation(template(data,'MET','nominal','sr'))

    ###
    # Z(->nunu)+jets data-driven model
    ###

    sr_zjetsTemplate = template(background,'Z+jets','nominal','sr')
    sr_zjetsMC =  rl.TemplateSample('sr'+year+'_zjetsMC', rl.Sample.BACKGROUND, sr_zjetsTemplate)
    sr_zjetsMC.setParamEffect(lumi, 1.027)
    sr_zjetsMC.setParamEffect(trig_met, 1.01)
    sr_zjetsMC.setParamEffect(veto_tau, 1.03)
    sr_zjetsMC.setParamEffect(zjets_norm, 1.4)
    sr_zjetsMC.setParamEffect(jec, 1.05)
    btagUp=template(background,'Z+jets','btagUp','sr')[0]
    btagDown=template(background,'Z+jets','btagDown','sr')[0]
    sr_zjetsMC.setParamEffect(btag, btagUp, btagDown)
    sr_zjetsBinYields = np.array([rl.IndependentParameter('sr'+year+'_zjets_bin_%d' % i, b, 0, sr_zjetsTemplate[0].max()*2) for i,b in enumerate(sr_zjetsTemplate[0])]) 
    sr_zjetsObservable = rl.Observable('recoil', sr_zjetsTemplate[1])
    sr_zjets = rl.ParametericSample('sr'+year+'_zjets', rl.Sample.BACKGROUND, sr_zjetsObservable, sr_zjetsBinYields)
    sr_zjetsBinYields = sr_zjetsBinYields * hf_fraction_weight['0tag']['DY+jets']
    sr_zjetsXweight = rl.ParametericSample(ch_name+'_zjets', rl.Sample.BACKGROUND, sr_zjetsObservable, sr_zjetsBinYields)
    sr.addSample(sr_zjetsXweight)

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    sr_wjetsTemplate = template(background,'W+jets','nominal','sr')
    sr_wjetsMC =  rl.TemplateSample('sr'+year+'_wjetsMC', rl.Sample.BACKGROUND, sr_wjetsTemplate)
    sr_wjetsMC.setParamEffect(lumi, 1.027)
    sr_wjetsMC.setParamEffect(trig_met, 1.01)
    sr_wjetsMC.setParamEffect(veto_tau, 1.03)
    sr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    sr_wjetsMC.setParamEffect(jec, 1.05)
    btagUp=template(background,'W+jets','btagUp','sr')[0]
    btagDown=template(background,'W+jets','btagDown','sr')[0]
    sr_wjetsMC.setParamEffect(btag, btagUp, btagDown)
    sr_wjetsBinYields = np.array([rl.IndependentParameter('sr'+year+'_wjets_bin_%d' % i,b,0,sr_wjetsTemplate[0].max()*2) for i,b in enumerate(sr_wjetsTemplate[0])]) 
    sr_wjetsObservable = rl.Observable('recoil', sr_wjetsTemplate[1])
    sr_wjets = rl.ParametericSample('sr'+year+'_wjets', rl.Sample.BACKGROUND, sr_wjetsObservable, sr_wjetsBinYields)
    sr_wjetsBinYields = sr_wjetsBinYields * hf_fraction_weight['0tag']['W+jets']
    sr_wjetsXweight = rl.ParametericSample(ch_name+'_wjets', rl.Sample.BACKGROUND, sr_wjetsObservable, sr_wjetsBinYields)
    sr.addSample(sr_wjetsXweight)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    sr_ttTemplate = template(background,'TT','nominal','sr')
    sr_ttMC =  rl.TemplateSample('sr'+year+'_ttMC', rl.Sample.BACKGROUND, sr_ttTemplate)
    sr_ttMC.setParamEffect(lumi, 1.027)
    sr_ttMC.setParamEffect(trig_met, 1.01)
    sr_ttMC.setParamEffect(veto_tau, 1.03)
    sr_ttMC.setParamEffect(tt_norm, 1.4)
    sr_ttMC.setParamEffect(jec, 1.05)
    btagUp=template(background,'TT','btagUp','sr')[0]
    btagDown=template(background,'TT','btagDown','sr')[0]
    sr_ttMC.setParamEffect(btag, btagUp, btagDown)
    sr_ttBinYields = np.array([rl.IndependentParameter('sr'+year+'_tt_bin_%d' % i,b,0,sr_ttTemplate[0].max()*2) for i,b in enumerate(sr_ttTemplate[0])])
    sr_ttObservable = rl.Observable('recoil', sr_ttTemplate[1])
    sr_tt = rl.ParametericSample('sr'+year+'_tt', rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields)
    sr_ttBinYields = sr_ttBinYields * deepak15_weight['0tag']['TT']
    sr_ttXweight = rl.ParametericSample(ch_name+'_tt', rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields)
    sr.addSample(sr_ttXweight)

    ###
    # Other MC-driven processes
    ###

    sr_stTemplate = template(background,'ST','nominal','sr')
    sr_st=rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, sr_stTemplate)
    sr_st.setParamEffect(lumi, 1.027)
    sr_st.setParamEffect(trig_met, 1.01)
    sr_st.setParamEffect(veto_tau, 1.03)
    sr_st.setParamEffect(st_norm, 1.2)
    sr_st.setParamEffect(jec, 1.05)
    btagUp=template(background,'ST','btagUp','sr')[0]
    btagDown=template(background,'ST','btagDown','sr')[0]
    sr_st.setParamEffect(btag, btagUp, btagDown)
    sr_st.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['ST'])
    sr.addSample(sr_st)

    sr_dyjetsTemplate = template(background,'DY+jets','nominal','sr')
    sr_dyjets=rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, sr_dyjetsTemplate)
    sr_dyjets.setParamEffect(lumi, 1.027)
    sr_dyjets.setParamEffect(trig_met, 1.01)
    sr_dyjets.setParamEffect(veto_tau, 1.03)
    sr_dyjets.setParamEffect(dyjets_norm, 1.4)
    sr_dyjets.setParamEffect(jec, 1.05)
    btagUp=template(background,'DY+jets','btagUp','sr')[0]
    btagDown=template(background,'DY+jets','btagDown','sr')[0]
    sr_dyjets.setParamEffect(btag, btagUp, btagDown)
    sr_dyjets.setParamEffect(deepak15_pass_sf['other'], hf_fraction_weight['0tag']['DY+jets'])
    sr.addSample(sr_dyjets)

    sr_vvTemplate = template(background,'VV','nominal','sr')
    sr_vv=rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, sr_vvTemplate)
    sr_vv.setParamEffect(lumi, 1.027)
    sr_vv.setParamEffect(trig_met, 1.01)
    sr_vv.setParamEffect(veto_tau, 1.03)
    sr_vv.setParamEffect(vv_norm, 1.2)
    sr_vv.setParamEffect(jec, 1.05)
    btagUp=template(background,'VV','btagUp','sr')[0]
    btagDown=template(background,'VV','btagDown','sr')[0]
    sr_vv.setParamEffect(btag, btagUp, btagDown)
    sr_vv.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['VV'])
    sr.addSample(sr_vv)

    sr_hbbTemplate = template(background,'Hbb','nominal','sr')
    sr_hbb=rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, sr_hbbTemplate)
    sr_hbb.setParamEffect(lumi, 1.027)
    sr_hbb.setParamEffect(trig_met, 1.01)
    sr_hbb.setParamEffect(veto_tau, 1.03)
    sr_hbb.setParamEffect(hbb_norm, 1.2)
    sr_hbb.setParamEffect(jec, 1.05)
    btagUp=template(background,'Hbb','btagUp','sr')[0]
    btagDown=template(background,'Hbb','btagDown','sr')[0]
    sr_hbb.setParamEffect(btag, btagUp, btagDown)
    sr_hbb.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['Hbb'])
    sr.addSample(sr_hbb)

    sr_qcdTemplate = template(background,'QCD','nominal','sr')
    sr_qcd=rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, sr_qcdTemplate)
    sr_qcd.setParamEffect(lumi, 1.027)
    sr_qcd.setParamEffect(trig_met, 1.01)
    sr_qcd.setParamEffect(veto_tau, 1.03)
    sr_qcd.setParamEffect(qcdsig_norm, 2.0)
    sr_qcd.setParamEffect(jec, 1.05)
    btagUp=template(background,'QCD','btagUp','sr')[0]
    btagDown=template(background,'QCD','btagDown','sr')[0]
    sr_qcd.setParamEffect(btag, btagUp, btagDown)
    sr_qcd.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['QCD'])
    sr.addSample(sr_qcd)

    for s in signal['sr'].identifiers('process'):
        sr_signalTemplate = template(signal,s,'nominal','sr')
        sr_signal=rl.TemplateSample(ch_name+'_'+str(s), rl.Sample.SIGNAL, sr_signalTemplate)
        sr_signal.setParamEffect(lumi, 1.027)
        sr_signal.setParamEffect(trig_met, 1.01)
        sr_signal.setParamEffect(veto_tau, 1.03)
        sr_signal.setParamEffect(jec, 1.05)
        btagUp=template(signal, s,'btagUp','sr')[0]
        btagDown=template(signal, s,'btagDown','sr')[0]
        sr_signal.setParamEffect(btag, btagUp, btagDown)
        sr_signal.setParamEffect(deepak15_pass_sf['other'], signal_weight[str(s)])
        sr.addSample(sr_signal)


    ###
    # End of SR
    ###

    

    ###
    ###
    # Single muon W control region
    ###
    ###

    ch_name = 'wmcr'+model_id
    wmcr = rl.Channel(ch_name)
    model.addChannel(wmcr)

    ###
    # Add data distribution to the channel
    ###

    wmcr.setObservation(template(data,'MET','nominal','wmcr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    wmcr_wjetsTemplate = template(background,'W+jets','nominal','wmcr')
    wmcr_wjetsMC =  rl.TemplateSample('wmcr'+year+'_wjetsMC', rl.Sample.BACKGROUND, wmcr_wjetsTemplate)
    wmcr_wjetsMC.setParamEffect(lumi, 1.027)
    wmcr_wjetsMC.setParamEffect(trig_met, 1.01)
    wmcr_wjetsMC.setParamEffect(veto_tau, 1.03)
    wmcr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    wmcr_wjetsMC.setParamEffect(jec, 1.05)
    wmcr_wjetsMC.setParamEffect(id_mu, 1.02)
    wmcr_wjetsMC.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'W+jets','btagUp','wmcr')[0]
    btagDown=template(background,'W+jets','btagDown','wmcr')[0]
    wmcr_wjetsMC.setParamEffect(btag, btagUp, btagDown)
    wmcr_wjetsTransferFactor = wmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['0tag']['W+jets']
    wmcr_wjetsXweight = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, wmcr_wjetsTransferFactor, sr_wjets)
    wmcr.addSample(wmcr_wjetsXweight)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    wmcr_ttTemplate = template(background,'TT','nominal','wmcr')
    wmcr_ttMC =  rl.TemplateSample('wmcr'+year+'_ttMC', rl.Sample.BACKGROUND, wmcr_ttTemplate)
    wmcr_ttMC.setParamEffect(lumi, 1.027)
    wmcr_ttMC.setParamEffect(trig_met, 1.01)
    wmcr_ttMC.setParamEffect(veto_tau, 1.03)
    wmcr_ttMC.setParamEffect(tt_norm, 1.4)
    wmcr_ttMC.setParamEffect(jec, 1.05)
    wmcr_ttMC.setParamEffect(id_mu, 1.02)
    wmcr_ttMC.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'TT','btagUp','wmcr')[0]
    btagDown=template(background,'TT','btagDown','wmcr')[0]
    wmcr_ttMC.setParamEffect(btag, btagUp, btagDown)
    wmcr_ttTransferFactor = wmcr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['0tag']['TT']
    wmcr_ttXweight = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt)
    wmcr.addSample(wmcr_ttXweight)

    ###
    # Other MC-driven processes
    ###

    wmcr_stTemplate = template(background,'ST','nominal','wmcr')
    wmcr_st=rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, wmcr_stTemplate)
    wmcr_st.setParamEffect(lumi, 1.027)
    wmcr_st.setParamEffect(trig_met, 1.01)
    wmcr_st.setParamEffect(veto_tau, 1.03)
    wmcr_st.setParamEffect(st_norm, 1.2)
    wmcr_st.setParamEffect(jec, 1.05)
    wmcr_st.setParamEffect(id_mu, 1.02)
    wmcr_st.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'ST','btagUp','wmcr')[0]
    btagDown=template(background,'ST','btagDown','wmcr')[0]
    wmcr_st.setParamEffect(btag, btagUp, btagDown)
    wmcr_st.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['ST'])
    wmcr.addSample(wmcr_st)

    wmcr_dyjetsTemplate = template(background,'DY+jets','nominal','wmcr')
    wmcr_dyjets=rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, wmcr_dyjetsTemplate)
    wmcr_dyjets.setParamEffect(lumi, 1.027)
    wmcr_dyjets.setParamEffect(trig_met, 1.01)
    wmcr_dyjets.setParamEffect(veto_tau, 1.03)
    wmcr_dyjets.setParamEffect(dyjets_norm, 1.4)
    wmcr_dyjets.setParamEffect(jec, 1.05)
    wmcr_dyjets.setParamEffect(id_mu, 1.02)
    wmcr_dyjets.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'DY+jets','btagUp','wmcr')[0]
    btagDown=template(background,'DY+jets','btagDown','wmcr')[0]
    wmcr_dyjets.setParamEffect(btag, btagUp, btagDown)
    wmcr_dyjets.setParamEffect(deepak15_pass_sf['other'], hf_fraction_weight['0tag']['DY+jets'])
    wmcr.addSample(wmcr_dyjets)

    wmcr_vvTemplate = template(background,'VV','nominal','wmcr')
    wmcr_vv=rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, wmcr_vvTemplate)
    wmcr_vv.setParamEffect(lumi, 1.027)
    wmcr_vv.setParamEffect(trig_met, 1.01)
    wmcr_vv.setParamEffect(veto_tau, 1.03)
    wmcr_vv.setParamEffect(vv_norm, 1.2)
    wmcr_vv.setParamEffect(jec, 1.05)
    wmcr_vv.setParamEffect(id_mu, 1.02)
    wmcr_vv.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'VV','btagUp','wmcr')[0]
    btagDown=template(background,'VV','btagDown','wmcr')[0]
    wmcr_vv.setParamEffect(btag, btagUp, btagDown)
    wmcr_vv.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['VV'])
    wmcr.addSample(wmcr_vv)

    wmcr_hbbTemplate = template(background,'Hbb','nominal','wmcr')
    wmcr_hbb=rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, wmcr_hbbTemplate)
    wmcr_hbb.setParamEffect(lumi, 1.027)
    wmcr_hbb.setParamEffect(trig_met, 1.01)
    wmcr_hbb.setParamEffect(veto_tau, 1.03)
    wmcr_hbb.setParamEffect(hbb_norm, 1.2)
    wmcr_hbb.setParamEffect(jec, 1.05)
    wmcr_hbb.setParamEffect(id_mu, 1.02)
    wmcr_hbb.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'Hbb','btagUp','wmcr')[0]
    btagDown=template(background,'Hbb','btagDown','wmcr')[0]
    wmcr_hbb.setParamEffect(btag, btagUp, btagDown)
    wmcr_hbb.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['Hbb'])
    wmcr.addSample(wmcr_hbb)

    wmcr_qcdTemplate = template(background,'QCD','nominal','wmcr')
    wmcr_qcd=rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, wmcr_qcdTemplate)
    wmcr_qcd.setParamEffect(lumi, 1.027)
    wmcr_qcd.setParamEffect(trig_met, 1.01)
    wmcr_qcd.setParamEffect(veto_tau, 1.03)
    wmcr_qcd.setParamEffect(qcdmu_norm, 2.0)
    wmcr_qcd.setParamEffect(jec, 1.05)
    wmcr_qcd.setParamEffect(id_mu, 1.02)
    wmcr_qcd.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'QCD','btagUp','wmcr')[0]
    btagDown=template(background,'QCD','btagDown','wmcr')[0]
    wmcr_qcd.setParamEffect(btag, btagUp, btagDown)
    wmcr_qcd.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['QCD'])
    wmcr.addSample(wmcr_qcd)

    ###
    # End of single muon W control region
    ###

    ###
    ###
    # Single muon top control region
    ###
    ###

    ch_name = 'tmcr'+model_id
    tmcr = rl.Channel(ch_name)
    model.addChannel(tmcr)

    ###
    # Add data distribution to the channel
    ###

    tmcr.setObservation(template(data,'MET','nominal','tmcr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    tmcr_wjetsTemplate = template(background,'W+jets','nominal','tmcr')
    tmcr_wjetsMC =  rl.TemplateSample('tmcr'+year+'_wjetsMC', rl.Sample.BACKGROUND, tmcr_wjetsTemplate)
    tmcr_wjetsMC.setParamEffect(lumi, 1.027)
    tmcr_wjetsMC.setParamEffect(trig_met, 1.01)
    tmcr_wjetsMC.setParamEffect(veto_tau, 1.03)
    tmcr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    tmcr_wjetsMC.setParamEffect(jec, 1.05)
    tmcr_wjetsMC.setParamEffect(id_mu, 1.02)
    tmcr_wjetsMC.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'W+jets','btagUp','tmcr')[0]
    btagDown=template(background,'W+jets','btagDown','tmcr')[0]
    tmcr_wjetsMC.setParamEffect(btag, btagUp, btagDown)
    tmcr_wjetsTransferFactor = tmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['1tag']['W+jets']
    tmcr_wjetsXweight = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, tmcr_wjetsTransferFactor, sr_wjets)
    tmcr.addSample(tmcr_wjetsXweight)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    tmcr_ttTemplate = template(background,'TT','nominal','tmcr')
    tmcr_ttMC =  rl.TemplateSample('tmcr'+year+'_ttMC', rl.Sample.BACKGROUND, tmcr_ttTemplate)
    tmcr_ttMC.setParamEffect(lumi, 1.027)
    tmcr_ttMC.setParamEffect(trig_met, 1.01)
    tmcr_ttMC.setParamEffect(veto_tau, 1.03)
    tmcr_ttMC.setParamEffect(tt_norm, 1.4)
    tmcr_ttMC.setParamEffect(jec, 1.05)
    tmcr_ttMC.setParamEffect(id_mu, 1.02)
    tmcr_ttMC.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'TT','btagUp','tmcr')[0]
    btagDown=template(background,'TT','btagDown','tmcr')[0]
    tmcr_ttMC.setParamEffect(btag, btagUp, btagDown)
    tmcr_ttTransferFactor = tmcr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['1tag']['TT']
    tmcr_ttXweight = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt)
    tmcr.addSample(tmcr_ttXweight)

    ###
    # Other MC-driven processes
    ###

    tmcr_stTemplate = template(background,'ST','nominal','tmcr')
    tmcr_st=rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, tmcr_stTemplate)
    tmcr_st.setParamEffect(lumi, 1.027)
    tmcr_st.setParamEffect(trig_met, 1.01)
    tmcr_st.setParamEffect(veto_tau, 1.03)
    tmcr_st.setParamEffect(st_norm, 1.2)
    tmcr_st.setParamEffect(jec, 1.05)
    tmcr_st.setParamEffect(id_mu, 1.02)
    tmcr_st.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'ST','btagUp','tmcr')[0]
    btagDown=template(background,'ST','btagDown','tmcr')[0]
    tmcr_st.setParamEffect(btag, btagUp, btagDown)
    tmcr_st.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['ST'])
    tmcr.addSample(tmcr_st)

    tmcr_dyjetsTemplate = template(background,'DY+jets','nominal','tmcr')
    tmcr_dyjets=rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, tmcr_dyjetsTemplate)
    tmcr_dyjets.setParamEffect(lumi, 1.027)
    tmcr_dyjets.setParamEffect(trig_met, 1.01)
    tmcr_dyjets.setParamEffect(veto_tau, 1.03)
    tmcr_dyjets.setParamEffect(dyjets_norm, 1.4)
    tmcr_dyjets.setParamEffect(jec, 1.05)
    tmcr_dyjets.setParamEffect(id_mu, 1.02)
    tmcr_dyjets.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'DY+jets','btagUp','tmcr')[0]
    btagDown=template(background,'DY+jets','btagDown','tmcr')[0]
    tmcr_dyjets.setParamEffect(btag, btagUp, btagDown)
    tmcr_dyjets.setParamEffect(deepak15_pass_sf['other'], hf_fraction_weight['1tag']['DY+jets'])
    tmcr.addSample(tmcr_dyjets)

    tmcr_vvTemplate = template(background,'VV','nominal','tmcr')
    tmcr_vv=rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, tmcr_vvTemplate)
    tmcr_vv.setParamEffect(lumi, 1.027)
    tmcr_vv.setParamEffect(trig_met, 1.01)
    tmcr_vv.setParamEffect(veto_tau, 1.03)
    tmcr_vv.setParamEffect(vv_norm, 1.2)
    tmcr_vv.setParamEffect(jec, 1.05)
    tmcr_vv.setParamEffect(id_mu, 1.02)
    tmcr_vv.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'VV','btagUp','tmcr')[0]
    btagDown=template(background,'VV','btagDown','tmcr')[0]
    tmcr_vv.setParamEffect(btag, btagUp, btagDown)
    tmcr_vv.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['VV'])
    tmcr.addSample(tmcr_vv)

    tmcr_hbbTemplate = template(background,'Hbb','nominal','tmcr')
    tmcr_hbb=rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, tmcr_hbbTemplate)
    tmcr_hbb.setParamEffect(lumi, 1.027)
    tmcr_hbb.setParamEffect(trig_met, 1.01)
    tmcr_hbb.setParamEffect(veto_tau, 1.03)
    tmcr_hbb.setParamEffect(hbb_norm, 1.2)
    tmcr_hbb.setParamEffect(jec, 1.05)
    tmcr_hbb.setParamEffect(id_mu, 1.02)
    tmcr_hbb.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'Hbb','btagUp','tmcr')[0]
    btagDown=template(background,'Hbb','btagDown','tmcr')[0]
    tmcr_hbb.setParamEffect(btag, btagUp, btagDown)
    tmcr_hbb.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['Hbb'])
    tmcr.addSample(tmcr_hbb)

    tmcr_qcdTemplate = template(background,'QCD','nominal','tmcr')
    tmcr_qcd=rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, tmcr_qcdTemplate)
    tmcr_qcd.setParamEffect(lumi, 1.027)
    tmcr_qcd.setParamEffect(trig_met, 1.01)
    tmcr_qcd.setParamEffect(veto_tau, 1.03)
    tmcr_qcd.setParamEffect(qcdmu_norm, 2.0)
    tmcr_qcd.setParamEffect(jec, 1.05)
    tmcr_qcd.setParamEffect(id_mu, 1.02)
    tmcr_qcd.setParamEffect(iso_mu, 1.02)
    btagUp=template(background,'QCD','btagUp','tmcr')[0]
    btagDown=template(background,'QCD','btagDown','tmcr')[0]
    tmcr_qcd.setParamEffect(btag, btagUp, btagDown)
    tmcr_qcd.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['QCD'])
    tmcr.addSample(tmcr_qcd)

    ###
    # End of single muon top control region
    ###

    ###
    ###
    # Single electron W control region
    ###
    ###

    ch_name = 'wecr'+model_id
    wecr = rl.Channel(ch_name)
    model.addChannel(wecr)

    ###
    # Add data distribution to the channel
    ###

    wecr.setObservation(template(data,'SingleElectron','nominal','wecr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    wecr_wjetsTemplate = template(background,'W+jets','nominal','wecr')
    wecr_wjetsMC =  rl.TemplateSample('wecr'+year+'_wjetsMC', rl.Sample.BACKGROUND, wecr_wjetsTemplate)
    wecr_wjetsMC.setParamEffect(lumi, 1.027)
    wecr_wjetsMC.setParamEffect(trig_e, 1.01)
    wecr_wjetsMC.setParamEffect(veto_tau, 1.03)
    wecr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    wecr_wjetsMC.setParamEffect(jec, 1.05)
    wecr_wjetsMC.setParamEffect(id_e, 1.02)
    wecr_wjetsMC.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'W+jets','btagUp','wecr')[0]
    btagDown=template(background,'W+jets','btagDown','wecr')[0]
    wecr_wjetsMC.setParamEffect(btag, btagUp, btagDown)
    wecr_wjetsTransferFactor = wecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['0tag']['W+jets']
    wecr_wjetsXweight = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, wecr_wjetsTransferFactor, sr_wjets)
    wecr.addSample(wecr_wjetsXweight)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    #wecr_ttHist = background['wecr'].integrate('process', 'TT').integrate('systematic','nominal')
    wecr_ttTemplate = template(background,'TT','nominal','wecr')
    wecr_ttMC =  rl.TemplateSample('wecr'+year+'_ttMC', rl.Sample.BACKGROUND, wecr_ttTemplate)
    wecr_ttMC.setParamEffect(lumi, 1.027)
    wecr_ttMC.setParamEffect(trig_e, 1.01)
    wecr_ttMC.setParamEffect(veto_tau, 1.03)
    wecr_ttMC.setParamEffect(tt_norm, 1.4)
    wecr_ttMC.setParamEffect(jec, 1.05)
    wecr_ttMC.setParamEffect(id_e, 1.02)
    wecr_ttMC.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'TT','btagUp','wecr')[0]
    btagDown=template(background,'TT','btagDown','wecr')[0]
    wecr_ttMC.setParamEffect(btag, btagUp, btagDown)
    wecr_ttTransferFactor = wecr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['0tag']['TT']
    wecr_ttXweight = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt)
    wecr.addSample(wecr_ttXweight)

    ###
    # Other MC-driven processes
    ###

    wecr_stTemplate = template(background,'ST','nominal','wecr')
    wecr_st=rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, wecr_stTemplate)
    wecr_st.setParamEffect(lumi, 1.027)
    wecr_st.setParamEffect(trig_e, 1.01)
    wecr_st.setParamEffect(veto_tau, 1.03)
    wecr_st.setParamEffect(st_norm, 1.2)
    wecr_st.setParamEffect(jec, 1.05)
    wecr_st.setParamEffect(id_e, 1.02)
    wecr_st.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'ST','btagUp','wecr')[0]
    btagDown=template(background,'ST','btagDown','wecr')[0]
    wecr_st.setParamEffect(btag, btagUp, btagDown)
    wecr_st.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['ST'])
    wecr.addSample(wecr_st)

    wecr_dyjetsTemplate = template(background,'DY+jets','nominal','wecr')
    wecr_dyjets=rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, wecr_dyjetsTemplate)
    wecr_dyjets.setParamEffect(lumi, 1.027)
    wecr_dyjets.setParamEffect(trig_e, 1.01)
    wecr_dyjets.setParamEffect(veto_tau, 1.03)
    wecr_dyjets.setParamEffect(dyjets_norm, 1.4)
    wecr_dyjets.setParamEffect(jec, 1.05)
    wecr_dyjets.setParamEffect(id_e, 1.02)
    wecr_dyjets.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'DY+jets','btagUp','wecr')[0]
    btagDown=template(background,'DY+jets','btagDown','wecr')[0]
    wecr_dyjets.setParamEffect(btag, btagUp, btagDown)
    wecr_dyjets.setParamEffect(deepak15_pass_sf['other'], hf_fraction_weight['0tag']['DY+jets'])
    wecr.addSample(wecr_dyjets)

    wecr_vvTemplate = template(background,'VV','nominal','wecr')
    wecr_vv=rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, wecr_vvTemplate)
    wecr_vv.setParamEffect(lumi, 1.027)
    wecr_vv.setParamEffect(trig_e, 1.01)
    wecr_vv.setParamEffect(veto_tau, 1.03)
    wecr_vv.setParamEffect(vv_norm, 1.2)
    wecr_vv.setParamEffect(jec, 1.05)
    wecr_vv.setParamEffect(id_e, 1.02)
    wecr_vv.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'VV','btagUp','wecr')[0]
    btagDown=template(background,'VV','btagDown','wecr')[0]
    wecr_vv.setParamEffect(btag, btagUp, btagDown)
    wecr_vv.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['VV'])
    wecr.addSample(wecr_vv)

    wecr_hbbTemplate = template(background,'Hbb','nominal','wecr')
    wecr_hbb=rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, wecr_hbbTemplate)
    wecr_hbb.setParamEffect(lumi, 1.027)
    wecr_hbb.setParamEffect(trig_e, 1.01)
    wecr_hbb.setParamEffect(veto_tau, 1.03)
    wecr_hbb.setParamEffect(hbb_norm, 1.2)
    wecr_hbb.setParamEffect(jec, 1.05)
    wecr_hbb.setParamEffect(id_e, 1.02)
    wecr_hbb.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'Hbb','btagUp','wecr')[0]
    btagDown=template(background,'Hbb','btagDown','wecr')[0]
    wecr_hbb.setParamEffect(btag, btagUp, btagDown)
    wecr_hbb.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['Hbb'])
    wecr.addSample(wecr_hbb)

    wecr_qcdTemplate = template(background,'QCD','nominal','wecr')
    wecr_qcd=rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, wecr_qcdTemplate)
    wecr_qcd.setParamEffect(lumi, 1.027)
    wecr_qcd.setParamEffect(trig_e, 1.01)
    wecr_qcd.setParamEffect(veto_tau, 1.03)
    wecr_qcd.setParamEffect(qcdmu_norm, 2.0)
    wecr_qcd.setParamEffect(jec, 1.05)
    wecr_qcd.setParamEffect(id_e, 1.02)
    wecr_qcd.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'QCD','btagUp','wecr')[0]
    btagDown=template(background,'QCD','btagDown','wecr')[0]
    wecr_qcd.setParamEffect(btag, btagUp, btagDown)
    wecr_qcd.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['0tag']['QCD'])
    wecr.addSample(wecr_qcd)

    ###
    # End of single electron W control region
    ###

    ###
    ###
    # Single electron top control region
    ###
    ###

    ch_name = 'tecr'+model_id
    tecr = rl.Channel(ch_name)
    model.addChannel(tecr)

    ###
    # Add data distribution to the channel
    ###

    tecr.setObservation(template(data,'SingleElectron','nominal','tecr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    tecr_wjetsTemplate = template(background,'W+jets','nominal','tecr')
    tecr_wjetsMC =  rl.TemplateSample('tecr'+year+'_wjetsMC', rl.Sample.BACKGROUND, tecr_wjetsTemplate)
    tecr_wjetsMC.setParamEffect(lumi, 1.027)
    tecr_wjetsMC.setParamEffect(trig_e, 1.01)
    tecr_wjetsMC.setParamEffect(veto_tau, 1.03)
    tecr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    tecr_wjetsMC.setParamEffect(jec, 1.05)
    tecr_wjetsMC.setParamEffect(id_e, 1.02)
    tecr_wjetsMC.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'W+jets','btagUp','tecr')[0]
    btagDown=template(background,'W+jets','btagDown','tecr')[0]
    tecr_wjetsMC.setParamEffect(btag, btagUp, btagDown)
    tecr_wjetsTransferFactor = tecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['1tag']['W+jets']
    tecr_wjetsXweight = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, tecr_wjetsTransferFactor, sr_wjets)
    tecr.addSample(tecr_wjetsXweight)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    tecr_ttTemplate = template(background,'TT','nominal','tecr')
    tecr_ttMC =  rl.TemplateSample('tecr'+year+'_ttMC', rl.Sample.BACKGROUND, tecr_ttTemplate)
    tecr_ttMC.setParamEffect(lumi, 1.027)
    tecr_ttMC.setParamEffect(trig_e, 1.01)
    tecr_ttMC.setParamEffect(veto_tau, 1.03)
    tecr_ttMC.setParamEffect(tt_norm, 1.4)
    tecr_ttMC.setParamEffect(jec, 1.05)
    tecr_ttMC.setParamEffect(id_e, 1.02)
    tecr_ttMC.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'TT','btagUp','tecr')[0]
    btagDown=template(background,'TT','btagDown','tecr')[0]
    tecr_ttMC.setParamEffect(btag, btagUp, btagDown)
    tecr_ttTransferFactor = tecr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['1tag']['TT']
    tecr_ttXweight = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt)
    tecr.addSample(tecr_ttXweight)

    ###
    # Other MC-driven processes
    ###

    tecr_stTemplate = template(background,'ST','nominal','tecr')
    tecr_st=rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, tecr_stTemplate)
    tecr_st.setParamEffect(lumi, 1.027)
    tecr_st.setParamEffect(trig_e, 1.01)
    tecr_st.setParamEffect(veto_tau, 1.03)
    tecr_st.setParamEffect(st_norm, 1.2)
    tecr_st.setParamEffect(jec, 1.05)
    tecr_st.setParamEffect(id_e, 1.02)
    tecr_st.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'ST','btagUp','tecr')[0]
    btagDown=template(background,'ST','btagDown','tecr')[0]
    tecr_st.setParamEffect(btag, btagUp, btagDown)
    tecr_st.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['ST'])
    tecr.addSample(tecr_st)

    tecr_dyjetsTemplate = template(background,'DY+jets','nominal','tecr')
    tecr_dyjets=rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, tecr_dyjetsTemplate)
    tecr_dyjets.setParamEffect(lumi, 1.027)
    tecr_dyjets.setParamEffect(trig_e, 1.01)
    tecr_dyjets.setParamEffect(veto_tau, 1.03)
    tecr_dyjets.setParamEffect(dyjets_norm, 1.4)
    tecr_dyjets.setParamEffect(jec, 1.05)
    tecr_dyjets.setParamEffect(id_e, 1.02)
    tecr_dyjets.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'DY+jets','btagUp','tecr')[0]
    btagDown=template(background,'DY+jets','btagDown','tecr')[0]
    tecr_dyjets.setParamEffect(btag, btagUp, btagDown)
    tecr_dyjets.setParamEffect(deepak15_pass_sf['other'], hf_fraction_weight['1tag']['DY+jets'])
    tecr.addSample(tecr_dyjets)

    tecr_vvTemplate = template(background,'VV','nominal','tecr')
    tecr_vv=rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, tecr_vvTemplate)
    tecr_vv.setParamEffect(lumi, 1.027)
    tecr_vv.setParamEffect(trig_e, 1.01)
    tecr_vv.setParamEffect(veto_tau, 1.03)
    tecr_vv.setParamEffect(vv_norm, 1.2)
    tecr_vv.setParamEffect(jec, 1.05)
    tecr_vv.setParamEffect(id_e, 1.02)
    tecr_vv.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'VV','btagUp','tecr')[0]
    btagDown=template(background,'VV','btagDown','tecr')[0]
    tecr_vv.setParamEffect(btag, btagUp, btagDown)
    tecr_vv.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['VV'])
    tecr.addSample(tecr_vv)

    tecr_hbbTemplate = template(background,'Hbb','nominal','tecr')
    tecr_hbb=rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, tecr_hbbTemplate)
    tecr_hbb.setParamEffect(lumi, 1.027)
    tecr_hbb.setParamEffect(trig_e, 1.01)
    tecr_hbb.setParamEffect(veto_tau, 1.03)
    tecr_hbb.setParamEffect(hbb_norm, 1.2)
    tecr_hbb.setParamEffect(jec, 1.05)
    tecr_hbb.setParamEffect(id_e, 1.02)
    tecr_hbb.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'Hbb','btagUp','tecr')[0]
    btagDown=template(background,'Hbb','btagDown','tecr')[0]
    tecr_hbb.setParamEffect(btag, btagUp, btagDown)
    tecr_hbb.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['Hbb'])
    tecr.addSample(tecr_hbb)

    tecr_qcdTemplate = template(background,'QCD','nominal','tecr')
    tecr_qcd=rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, tecr_qcdTemplate)
    tecr_qcd.setParamEffect(lumi, 1.027)
    tecr_qcd.setParamEffect(trig_e, 1.01)
    tecr_qcd.setParamEffect(veto_tau, 1.03)
    tecr_qcd.setParamEffect(qcdmu_norm, 2.0)
    tecr_qcd.setParamEffect(jec, 1.05)
    tecr_qcd.setParamEffect(id_e, 1.02)
    tecr_qcd.setParamEffect(reco_e, 1.02)
    btagUp=template(background,'QCD','btagUp','tecr')[0]
    btagDown=template(background,'QCD','btagDown','tecr')[0]
    tecr_qcd.setParamEffect(btag, btagUp, btagDown)
    tecr_qcd.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['1tag']['QCD'])
    tecr.addSample(tecr_qcd)

    ###
    # End of single electron top control region
    ###

    ###
    ###
    # Double muon control region
    ###
    ###

    ch_name = 'zmcr'+model_id
    zmcr = rl.Channel(ch_name)
    model.addChannel(zmcr)

    ###
    # Add data distribution to the channel
    ###

    zmcr.setObservation(template(data,'MET','nominal','zmcr'))

    zmcr_dyjetsTemplate = template(background,'DY+jets','nominal','zmcr')
    zmcr_dyjetsMC =  rl.TemplateSample('zmcr'+year+'_dyjetsMC', rl.Sample.BACKGROUND, zmcr_dyjetsTemplate)
    zmcr_dyjetsMC.setParamEffect(lumi, 1.027)
    zmcr_dyjetsMC.setParamEffect(trig_met, 1.01)
    zmcr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    zmcr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    zmcr_dyjetsMC.setParamEffect(jec, 1.05)
    zmcr_dyjetsMC.setParamEffect(id_mu, 1.02)
    zmcr_dyjetsMC.setParamEffect(iso_mu, 1.02)
    zmcr_dyjetsTransferFactor = zmcr_dyjetsMC.getExpectation() / sr_zjetsMC.getExpectation() * hf_fraction_weight['notag']['DY+jets']
    zmcr_dyjetsXweight = rl.TransferFactorSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, zmcr_dyjetsTransferFactor, sr_zjets)
    zmcr.addSample(zmcr_dyjetsXweight)

    ###
    # Other MC-driven processes
    ###

    zmcr_ttTemplate = template(background,'TT','nominal','zmcr')
    zmcr_tt=rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, zmcr_ttTemplate)
    zmcr_tt.setParamEffect(lumi, 1.027)
    zmcr_tt.setParamEffect(trig_met, 1.01)
    zmcr_tt.setParamEffect(veto_tau, 1.03)
    zmcr_tt.setParamEffect(tt_norm, 1.4)
    zmcr_tt.setParamEffect(jec, 1.05)
    zmcr_tt.setParamEffect(id_mu, 1.02)
    zmcr_tt.setParamEffect(iso_mu, 1.02)
    zmcr_tt.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['TT'])
    zmcr.addSample(zmcr_tt)

    zmcr_stTemplate = template(background,'ST','nominal','zmcr')
    zmcr_st=rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, zmcr_stTemplate)
    zmcr_st.setParamEffect(lumi, 1.027)
    zmcr_st.setParamEffect(trig_met, 1.01)
    zmcr_st.setParamEffect(veto_tau, 1.03)
    zmcr_st.setParamEffect(st_norm, 1.2)
    zmcr_st.setParamEffect(jec, 1.05)
    zmcr_st.setParamEffect(id_mu, 1.02)
    zmcr_st.setParamEffect(iso_mu, 1.02)
    zmcr_st.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['ST'])
    zmcr.addSample(zmcr_st)

    zmcr_vvTemplate = template(background,'VV','nominal','zmcr')
    zmcr_vv=rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, zmcr_vvTemplate)
    zmcr_vv.setParamEffect(lumi, 1.027)
    zmcr_vv.setParamEffect(trig_met, 1.01)
    zmcr_vv.setParamEffect(veto_tau, 1.03)
    zmcr_vv.setParamEffect(vv_norm, 1.2)
    zmcr_vv.setParamEffect(jec, 1.05)
    zmcr_vv.setParamEffect(id_mu, 1.02)
    zmcr_vv.setParamEffect(iso_mu, 1.02)
    zmcr_vv.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['VV'])
    zmcr.addSample(zmcr_vv)

    zmcr_hbbTemplate = template(background,'Hbb','nominal','zmcr')
    zmcr_hbb=rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, zmcr_hbbTemplate)
    zmcr_hbb.setParamEffect(lumi, 1.027)
    zmcr_hbb.setParamEffect(trig_met, 1.01)
    zmcr_hbb.setParamEffect(veto_tau, 1.03)
    zmcr_hbb.setParamEffect(hbb_norm, 1.2)
    zmcr_hbb.setParamEffect(jec, 1.05)
    zmcr_hbb.setParamEffect(id_mu, 1.02)
    zmcr_hbb.setParamEffect(iso_mu, 1.02)
    zmcr_hbb.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['Hbb'])
    zmcr.addSample(zmcr_hbb)

    ###
    # End of double muon control region
    ###

    ###
    ###
    # Double electron control region
    ###
    ###

    ch_name = 'zecr'+model_id
    zecr = rl.Channel(ch_name)
    model.addChannel(zecr)

    ###
    # Add data distribution to the channel
    ###

    zecr.setObservation(template(data,'SingleElectron','nominal','zecr'))

    zecr_dyjetsTemplate = template(background,'DY+jets','nominal','zecr')
    zecr_dyjetsMC =  rl.TemplateSample('zecr'+year+'_dyjetsMC', rl.Sample.BACKGROUND, zecr_dyjetsTemplate)
    zecr_dyjetsMC.setParamEffect(lumi, 1.027)
    zecr_dyjetsMC.setParamEffect(trig_e, 1.01)
    zecr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    zecr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    zecr_dyjetsMC.setParamEffect(jec, 1.05)
    zecr_dyjetsMC.setParamEffect(id_e, 1.02)
    zecr_dyjetsMC.setParamEffect(reco_e, 1.02)
    zecr_dyjetsTransferFactor = zecr_dyjetsMC.getExpectation() / sr_zjetsMC.getExpectation() * hf_fraction_weight['notag']['DY+jets']
    zecr_dyjetsXweight = rl.TransferFactorSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, zecr_dyjetsTransferFactor, sr_zjets)
    zecr.addSample(zecr_dyjetsXweight)

    ###
    # Other MC-driven processes
    ###

    zecr_ttTemplate = template(background,'TT','nominal','zecr')
    zecr_tt=rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, zecr_ttTemplate)
    zecr_tt.setParamEffect(lumi, 1.027)
    zecr_tt.setParamEffect(trig_e, 1.01)
    zecr_tt.setParamEffect(veto_tau, 1.03)
    zecr_tt.setParamEffect(tt_norm, 1.4)
    zecr_tt.setParamEffect(jec, 1.05)
    zecr_tt.setParamEffect(id_e, 1.02)
    zecr_tt.setParamEffect(reco_e, 1.02)
    zecr_tt.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['TT'])
    zecr.addSample(zecr_tt)

    zecr_stTemplate = template(background,'ST','nominal','zecr')
    zecr_st=rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, zecr_stTemplate)
    zecr_st.setParamEffect(lumi, 1.027)
    zecr_st.setParamEffect(trig_e, 1.01)
    zecr_st.setParamEffect(veto_tau, 1.03)
    zecr_st.setParamEffect(st_norm, 1.2)
    zecr_st.setParamEffect(jec, 1.05)
    zecr_st.setParamEffect(id_e, 1.02)
    zecr_st.setParamEffect(reco_e, 1.02)
    zecr_st.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['ST'])
    zecr.addSample(zecr_st)

    zecr_vvTemplate = template(background,'VV','nominal','zecr')
    zecr_vv=rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, zecr_vvTemplate)
    zecr_vv.setParamEffect(lumi, 1.027)
    zecr_vv.setParamEffect(trig_e, 1.01)
    zecr_vv.setParamEffect(veto_tau, 1.03)
    zecr_vv.setParamEffect(vv_norm, 1.2)
    zecr_vv.setParamEffect(jec, 1.05)
    zecr_vv.setParamEffect(id_e, 1.02)
    zecr_vv.setParamEffect(reco_e, 1.02)
    zecr_vv.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['VV'])
    zecr.addSample(zecr_vv)

    zecr_hbbTemplate = template(background,'Hbb','nominal','zecr')
    zecr_hbb=rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, zecr_hbbTemplate)
    zecr_hbb.setParamEffect(lumi, 1.027)
    zecr_hbb.setParamEffect(trig_e, 1.01)
    zecr_hbb.setParamEffect(veto_tau, 1.03)
    zecr_hbb.setParamEffect(hbb_norm, 1.2)
    zecr_hbb.setParamEffect(jec, 1.05)
    zecr_hbb.setParamEffect(id_e, 1.02)
    zecr_hbb.setParamEffect(reco_e, 1.02)
    zecr_hbb.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['Hbb'])
    zecr.addSample(zecr_hbb)

    ###
    # End of double electron control region
    ###

    ###
    ###
    # Single photon control region
    ###
    ###

    ch_name = 'gcr'+model_id
    gcr = rl.Channel(ch_name)
    model.addChannel(gcr)

    ###
    # Add data distribution to the channel
    ###

    gcr.setObservation(template(data,'SinglePhoton','nominal','gcr'))

    gcr_gjetsTemplate = template(background,'G+jets','nominal','gcr')
    gcr_gjetsMC =  rl.TemplateSample('gcr'+year+'_gjetsMC', rl.Sample.BACKGROUND, gcr_gjetsTemplate)
    gcr_gjetsMC.setParamEffect(lumi, 1.027)
    gcr_gjetsMC.setParamEffect(trig_pho, 1.01)
    gcr_gjetsMC.setParamEffect(veto_tau, 1.03)
    gcr_gjetsMC.setParamEffect(gjets_norm, 1.4)
    gcr_gjetsMC.setParamEffect(jec, 1.05)
    gcr_gjetsMC.setParamEffect(id_pho, 1.02)
    gcr_gjetsTransferFactor = gcr_gjetsMC.getExpectation() / sr_zjetsMC.getExpectation() * hf_fraction_weight['notag']['G+jets']
    gcr_gjetsXweight = rl.TransferFactorSample(ch_name+'_gjets', rl.Sample.BACKGROUND, gcr_gjetsTransferFactor, sr_zjets)
    gcr.addSample(gcr_gjetsXweight)

    gcr_qcdTemplate = template(background,'QCD','nominal','gcr')
    gcr_qcd=rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, gcr_qcdTemplate)
    gcr_qcd.setParamEffect(lumi, 1.027)
    gcr_qcd.setParamEffect(trig_pho, 1.01)
    gcr_qcd.setParamEffect(veto_tau, 1.03)
    gcr_qcd.setParamEffect(qcdpho_norm, 2.0)
    gcr_qcd.setParamEffect(jec, 1.05)
    gcr_qcd.setParamEffect(id_pho, 1.02)
    gcr_qcd.setParamEffect(deepak15_pass_sf['other'], deepak15_weight['notag']['QCD'])
    gcr.addSample(gcr_qcd)

    return model

if __name__ == '__main__':
    if not os.path.exists('datacards'):
        os.mkdir('datacards')
    parser = OptionParser()
    parser.add_option('-g', '--category', help='category', dest='category', default='')
    parser.add_option('-y', '--year', help='year', dest='year', default='')
    (options, args) = parser.parse_args()
        
    ###
    #Extract histograms from input file
    ###

    print('Grouping histograms')
    hists = load('hists/darkhiggs'+options.year+'.scaled')
    
    data_hists   = {}
    bkg_hists    = {}
    signal_hists = {}

    process = hist.Cat("process", "Process", sorting='placement')
    cats = ("process",)
    sig_map = OrderedDict()
    bkg_map = OrderedDict()
    data_map = OrderedDict()
    bkg_map["Hbb"] = ("Hbb*",)
    bkg_map["DY+jets"] = ("DY+*",)
    bkg_map["VV"] = (["WW","WZ","ZZ"],)
    bkg_map["ST"] = ("ST*",)
    bkg_map["TT"] = ("TT*",)
    bkg_map["W+jets"] = ("W+*",)
    bkg_map["Z+jets"] = ("Z+*",)
    bkg_map["G+jets"] = ("G+*",)
    bkg_map["QCD"] = ("QCD*",)
    sig_map["Mhs_50"] = ("*Mhs_50*",)  ## signals
    sig_map["Mhs_70"] = ("*Mhs_70*",)
    sig_map["Mhs_90"] = ("*Mhs_90*",)
    sig_map["MonoJet"] = ("MonoJet*",)  ## signals
    sig_map["MonoW"] = ("MonoW*",)    ## signals
    sig_map["MonoZ"] = ("MonoZ*",)    ## signals
    data_map["MET"] = ("MET*", )
    data_map["SingleElectron"] = ("EGamma*", )
    data_map["SinglePhoton"] = ("EGamma*", )
    
    for key in hists['data'].keys():
        bkg_hists[key] = hists['bkg'][key].group(cats, process, bkg_map)
        signal_hists[key] = hists['sig'][key].group(cats, process, sig_map)
        data_hists[key] = hists['data'][key].group(cats, process, data_map)

    hists={
        'bkg': bkg_hists,
        'sig': signal_hists,
        'data': data_hists
    }

    model_dict={}
    for category in ['monojet','monohs']:
        if options.category and options.category not in category: continue
        if category=='monojet':
            with open('data/darkhiggs'+options.year+'-'+category+'.model', "wb") as fout:
                pickle.dump(model(options.year,None,category), fout, protocol=2)
        elif category=='monohs':
            for mass in [0,1,2,3,4]:
                with open('data/darkhiggs'+options.year+'-'+category+'-mass'+str(mass)+'.model', "wb") as fout:
                    pickle.dump(model(options.year,mass,category), fout, protocol=2)

