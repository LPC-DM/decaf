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
import json
from coffea import hist, processor 
from coffea.util import load, save

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

def expo_sample(norm, scale, obs):
    cdf = scipy.stats.expon.cdf(scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)

def template(hist, name):
    #return (hist.values(overflow='all')[()], hist.axis(name).edges(overflow='all'), name)
    return (hist.values()[()], hist.axis(name).edges(), name)

def darkhiggs_model(tmpdir,mass,category,year,grouping):
    
    model = rl.Model('darkhiggs_'+mass+'_'+category)

    whf_fraction=0.18
    zhf_fraction=0.10
    ghf_fraction=0.12

    whf_k = rl.IndependentParameter('whf_k', 1., 0, 1/whf_fraction)
    zhf_k = rl.IndependentParameter('zhf_k', 1., 0, 1/zhf_fraction)
    ghf_k = rl.IndependentParameter('ghf_k', 1., 0, 1/ghf_fraction)

    deepak15_pass_eff = {
        'xbb'    : 0.88, 
        'vqq'    : 0.06, 
        'wcq'    : 0.21, 
        'b'      : 0.59, 
        'bb'     : 0.78, 
        'bc'     : 0.73, 
        'c'      : 0.16, 
        'cc'     : 0.25, 
        'garbage': 0.0, 
        'other'  : 0.06, 
        'tbcq'   : 0.68, 
        'tbqq'   : 0.55, 
        'zcc'    : 0.40
    }
    deepak15_pass_sf = {
        'xbb'    : rl.IndependentParameter('xbb_sf', 1., 0, 1/deepak15_pass_eff['xbb']), 
        'vqq'    : rl.IndependentParameter('vqq_sf', 1., 0, 1/deepak15_pass_eff['vqq']),
        'wcq'    : rl.IndependentParameter('wcq_sf', 1., 0, 1/deepak15_pass_eff['wcq']),
        'b'      : rl.IndependentParameter('b_sf', 1., 0, 1/deepak15_pass_eff['b']),
        'bb'     : rl.IndependentParameter('bb_sf', 1., 0, 1/deepak15_pass_eff['bb']), 
        'bc'     : rl.IndependentParameter('bc_sf', 1., 0, 1/deepak15_pass_eff['bc']),
        'c'      : rl.IndependentParameter('c_sf', 1., 0, 1/deepak15_pass_eff['c']),
        'cc'     : rl.IndependentParameter('cc_sf', 1., 0, 1/deepak15_pass_eff['cc']),
        'garbage': 1.,
        'other'  : rl.IndependentParameter('other_sf', 1., 0, 1/deepak15_pass_eff['other']),
        'tbcq'   : rl.IndependentParameter('tbcq_sf', 1., 0, 1/deepak15_pass_eff['tbcq']),
        'tbqq'   : rl.IndependentParameter('tbqq_sf', 1., 0, 1/deepak15_pass_eff['tbqq']),
        'zcc'    : rl.IndependentParameter('zcc_sf', 1., 0, 1/deepak15_pass_eff['zcc'])
    }
    
    deepak4_0tag_component_eff = {
        'Hbb': {
            'xbb': 0.8697046312331598, 
            'vqq': 0.23124434902057903, 
            'wcq': 0.22047080192919324, 
            'b': 0.6820164771169962, 
            'bb': 0.46661580939496744, 
            'bc': 0.38229026101737873, 
            'c': 0.38635373928707456, 
            'cc': 0.46294460945495375, 
            'garbage': 0.67443302378121, 
            'other': 0.4463215680021725, 
            'tbcq': 0.29126844822716513, 
            'tbqq': 0.30888680936982366, 
            'zcc': 0.0
        }, 
        'Z+HF': {
            'xbb': 0.0, 
            'vqq': 0.0, 
            'wcq': 0.0, 
            'b': 0.7979449881691513, 
            'bb': 0.8961869119285696, 
            'bc': 0.4910183479912289, 
            'c': 0.8616372499954866, 
            'cc': 0.8949353436835303, 
            'garbage': 0.7666641925645035, 
            'other': 0.6522792754261665, 
            'tbcq': 0.0, 
            'tbqq': 0.0, 
            'zcc': 0.0
        }, 
        'Z+LF': {
            'xbb': 0.0, 
            'vqq': 0.0, 
            'wcq': 0.0, 
            'b': 0.0, 
            'bb': 0.0, 
            'bc': 0.0, 
            'c': 0.0, 
            'cc': 0.0, 
            'garbage': 0.887874059807142, 
            'other': 0.89619464074004, 
            'tbcq': 0.0, 
            'tbqq': 0.0, 
            'zcc': 0.0
        }, 
        'G+HF': {
            'xbb': 0.0, 
            'vqq': 0.0, 
            'wcq': 0.0, 
            'b': 0.7979449881691513, 
            'bb': 0.8961869119285696, 
            'bc': 0.4910183479912289, 
            'c': 0.8616372499954866, 
            'cc': 0.8949353436835303, 
            'garbage': 0.7666641925645035, 
            'other': 0.6522792754261665, 
            'tbcq': 0.0, 
            'tbqq': 0.0, 
            'zcc': 0.0
        }, 
        'G+LF': {
            'xbb': 0.0, 
            'vqq': 0.0, 
            'wcq': 0.0, 
            'b': 0.0, 
            'bb': 0.0, 
            'bc': 0.0, 
            'c': 0.0, 
            'cc': 0.0, 
            'garbage': 0.887874059807142, 
            'other': 0.89619464074004, 
            'tbcq': 0.0, 
            'tbqq': 0.0, 
            'zcc': 0.0
        }, 
        'VV': {
            'xbb': 0.9419428942394055, 
            'vqq': 0.9280086922176735, 
            'wcq': 0.9336817821532177, 
            'b': 0.595296958248878, 
            'bb': 0.0, 
            'bc': 1.0, 
            'c': 0.8906290333285863, 
            'cc': 1.0, 
            'garbage': 0.8695515160067608, 
            'other': 0.845611807805552, 
            'tbcq': 0.0, 
            'tbqq': 0.0, 
            'zcc': 0.9531669398387308
        }, 
        'ST': {
            'xbb': 0.0, 
            'vqq': 0.5717995037496204, 
            'wcq': 0.5691996356188856, 
            'b': 0.6224105767039395, 
            'bb': 0.8363685746865235, 
            'bc': 0.710128136756667, 
            'c': 0.4489667441833009, 
            'cc': 0.44850369816000596, 
            'garbage': 0.6088520617501884, 
            'other': 0.4324488197731728, 
            'tbcq': 0.8051769837186212, 
            'tbqq': 0.8161799711888433, 
            'zcc': 0.0
        }, 
        'TT': {
            'xbb': 0.0, 
            'vqq': 0.4459557368488155, 
            'wcq': 0.450844924955211, 
            'b': 0.4917130631084224, 
            'bb': 0.8314067750944469, 
            'bc': 0.5031580464547177, 
            'c': 0.3037376518772207, 
            'cc': 0.3042761840537053, 
            'garbage': 0.41293445024465525, 
            'other': 0.28015111794952063, 
            'tbcq': 0.5920703770856897, 
            'tbqq': 0.5885231523899408, 
            'zcc': 0.0
        }, 
        'W+HF': {
            'xbb': 0.0, 
            'vqq': 0.0, 
            'wcq': 0.0, 
            'b': 0.6830977711065136, 
            'bb': 0.9015729268035118, 
            'bc': 0.7316564932962332, 
            'c': 0.9062705207112773, 
            'cc': 0.9011852284020715, 
            'garbage': 0.831577774188298, 
            'other': 0.821589706002926, 
            'tbcq': 0.0, 
            'tbqq': 0.0, 
            'zcc': 0.0
        }, 
        'W+LF': {
            'xbb': 0.0, 
            'vqq': 0.0, 
            'wcq': 0.0, 
            'b': 0.0, 
            'bb': 0.0, 
            'bc': 0.0, 
            'c': 0.0, 
            'cc': 0.0, 
            'garbage': 0.8799004964320957, 
            'other': 0.906565500257813, 
            'tbcq': 0.0, 
            'tbqq': 0.0,
            'zcc': 0.0
        },
        'QCD': {
            'xbb': 0.0, 
            'vqq': 0.0, 
            'wcq': 0.0, 
            'b': 0.8052218395379344, 
            'bb': 0.9676263397285139, 
            'bc': 1.0, 
            'c': 0.8882376783888627, 
            'cc': 0.7150587923851417, 
            'garbage': 0.8563298154080022, 
            'other': 0.838726390208842, 
            'tbcq': 0.0, 
            'tbqq': 0.0, 
            'zcc': 0.0
        }
    }

    deepak4_0tag_process_eff = {
        'Hbb': 0.740183629513201, 
        'Z+HF': 0.7834895739514862, 
        'Z+LF': 0.8920487242439581, 
        'G+HF': 0.7834895739514862,
        'G+LF': 0.8920487242439581,
        'VV': 0.8899188652047255, 
        'ST': 0.6497776622766339, 
        'TT': 0.5140433389517574, 
        'W+HF': 0.8587729399436966, 
        'W+LF': 0.9030101162281751, 
        'QCD': 0.8425783767140694,
        'W+jets':0.8951594995622084,
        'Z+jets':0.8820319256963761,
        'G+jets':0.8820319256963761
    }

    with open('data/signal_fractions.json') as fin:
        signal_fractions = json.load(fin)
        
    signal_weight={}
    for process in signal_fractions.keys():
        for component in signal_fractions[process].keys():
            print('Extracting',component,'fraction for',process)
            
            try:
                weight
            except:
                weight = deepak15_pass_sf[component]*deepak15_pass_eff[component]*np.array(signal_fractions[process][component][mass])
            else:
                weight += deepak15_pass_sf[component]*deepak15_pass_eff[component]*np.array(signal_fractions[process][component][mass])

            if 'monojet' in category:
                try:
                    weight
                except:
                    weight = (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*np.array(signal_fractions[process][component][mass])
                else:
                    weight += (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*np.array(signal_fractions[process][component][mass])

        signal_weight[process]=weight
                
    with open('data/fractions.json') as fin:
        fractions = json.load(fin)

    deepak15_weight={}
    deepak15_weight['0tag']={}
    deepak15_weight['1tag']={}
    deepak15_weight['notag']={}
    for process in ['Hbb','VV','ST','QCD','TT','Z+HF','Z+LF','W+HF','W+LF','G+HF','G+LF']:
        for component in fractions[process].keys():
            print('Extracting',component,'fraction for',process )

            try:
                weight_0tag
            except:
                weight_0tag = deepak15_pass_sf[component]*deepak15_pass_eff[component]*deepak4_0tag_component_eff[process][component]*np.array(fractions[process][component][mass])
            else:
                weight_0tag += deepak15_pass_sf[component]*deepak15_pass_eff[component]*deepak4_0tag_component_eff[process][component]*np.array(fractions[process][component][mass])

            try:
                weight_1tag
            except:
                weight_1tag = deepak15_pass_sf[component]*deepak15_pass_eff[component]*(1 - deepak4_0tag_component_eff[process][component])*np.array(fractions[process][component][mass])
            else:
                weight_1tag += deepak15_pass_sf[component]*deepak15_pass_eff[component]*(1 - deepak4_0tag_component_eff[process][component])*np.array(fractions[process][component][mass])

            try:
                weight_notag
            except:
                weight_notag = deepak15_pass_sf[component]*deepak15_pass_eff[component]*np.array(fractions[process][component][mass])
            else:
                weight_notag += deepak15_pass_sf[component]*deepak15_pass_eff[component]*np.array(fractions[process][component][mass])
                
            if 'monojet' in category: 
                
                try:
                    weight_0tag
                except:
                    weight_0tag = (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*deepak4_0tag_component_eff[process][component]*np.array(fractions[process][component][mass])
                else:
                    weight_0tag += (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*deepak4_0tag_component_eff[process][component]*np.array(fractions[process][component][mass])

                try:
                    weight_1tag
                except:
                    weight_1tag = (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*(1 - deepak4_0tag_component_eff[process][component])*np.array(fractions[process][component][mass])
                else:
                    weight_1tag += (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*(1 - deepak4_0tag_component_eff[process][component])*np.array(fractions[process][component][mass])

                try:
                    weight_notag
                except:
                    weight_notag = (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*np.array(fractions[process][component][mass])
                else:
                    weight_notag += (1 - deepak15_pass_sf[component]*deepak15_pass_eff[component])*np.array(fractions[process][component][mass])

        deepak15_weight['0tag'][process]=weight_0tag/deepak4_0tag_process_eff[process]
        deepak15_weight['1tag'][process]=weight_1tag/(1 - deepak4_0tag_process_eff[process])
        deepak15_weight['notag'][process]=weight_notag

    hf_fraction_weight={}
    hf_fraction_weight['0tag']={}
    hf_fraction_weight['1tag']={}
    hf_fraction_weight['notag']={}

    hf_fraction_weight['0tag']['W+jets'] = deepak15_weight['0tag']['W+HF']*(deepak4_0tag_process_eff['W+HF']/deepak4_0tag_process_eff['W+jets'])*whf_k*whf_fraction 
    hf_fraction_weight['0tag']['W+jets'] += deepak15_weight['0tag']['W+LF']*(deepak4_0tag_process_eff['W+LF']/deepak4_0tag_process_eff['W+jets'])*(1 - whf_k*whf_fraction)
    hf_fraction_weight['1tag']['W+jets'] = deepak15_weight['1tag']['W+HF']*((1-deepak4_0tag_process_eff['W+HF'])/(1-deepak4_0tag_process_eff['W+jets']))*whf_k*whf_fraction 
    hf_fraction_weight['1tag']['W+jets'] += deepak15_weight['1tag']['W+LF']*((1-deepak4_0tag_process_eff['W+HF'])/(1-deepak4_0tag_process_eff['W+jets']))*(1 - whf_k*whf_fraction)
    hf_fraction_weight['notag']['W+jets'] = deepak15_weight['notag']['W+HF']*whf_k*whf_fraction
    hf_fraction_weight['notag']['W+jets'] += deepak15_weight['notag']['W+LF']*(1 - whf_k*whf_fraction)

    hf_fraction_weight['0tag']['Z+jets'] = deepak15_weight['0tag']['Z+HF']*(deepak4_0tag_process_eff['Z+HF']/deepak4_0tag_process_eff['Z+jets'])*zhf_k*zhf_fraction 
    hf_fraction_weight['0tag']['Z+jets'] += deepak15_weight['0tag']['Z+LF']*(deepak4_0tag_process_eff['Z+LF']/deepak4_0tag_process_eff['Z+jets'])*(1 - zhf_k*zhf_fraction)
    hf_fraction_weight['1tag']['Z+jets'] = deepak15_weight['1tag']['Z+HF']*((1-deepak4_0tag_process_eff['Z+HF'])/(1-deepak4_0tag_process_eff['Z+jets']))*zhf_k*zhf_fraction 
    hf_fraction_weight['1tag']['Z+jets'] += deepak15_weight['1tag']['Z+LF']*((1-deepak4_0tag_process_eff['Z+HF'])/(1-deepak4_0tag_process_eff['Z+jets']))*(1 - zhf_k*zhf_fraction)
    hf_fraction_weight['notag']['Z+jets'] = deepak15_weight['notag']['Z+HF']*zhf_k*zhf_fraction
    hf_fraction_weight['notag']['Z+jets'] += deepak15_weight['notag']['Z+LF']*(1 - zhf_k*zhf_fraction)

    hf_fraction_weight['0tag']['G+jets'] = deepak15_weight['0tag']['G+HF']*(deepak4_0tag_process_eff['G+HF']/deepak4_0tag_process_eff['G+jets'])*ghf_k*ghf_fraction 
    hf_fraction_weight['0tag']['G+jets'] += deepak15_weight['0tag']['G+LF']*(deepak4_0tag_process_eff['G+LF']/deepak4_0tag_process_eff['G+jets'])*(1 - ghf_k*ghf_fraction)
    hf_fraction_weight['1tag']['G+jets'] = deepak15_weight['1tag']['G+HF']*((1-deepak4_0tag_process_eff['G+HF'])/(1-deepak4_0tag_process_eff['G+jets']))*ghf_k*ghf_fraction 
    hf_fraction_weight['1tag']['G+jets'] += deepak15_weight['1tag']['G+LF']*((1-deepak4_0tag_process_eff['G+HF'])/(1-deepak4_0tag_process_eff['G+jets']))*(1 - ghf_k*ghf_fraction)
    hf_fraction_weight['notag']['G+jets'] = deepak15_weight['notag']['G+HF']*ghf_k*ghf_fraction
    hf_fraction_weight['notag']['G+jets'] += deepak15_weight['notag']['G+LF']*(1 - ghf_k*ghf_fraction)
    
        
    ###
    #Extract histograms from input file
    ###

    if grouping:

        print('Grouping histograms')
        hists = load('hists/darkhiggs'+year+'.plot')

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
        bkg_map["VV"] = ("VV*",)
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
        data_map["SingleElectron"] = ("SingleElectron*", )
        data_map["SinglePhoton"] = ("SinglePhoton*", )

        for key in hists['data'].keys():
            bkg_hists[key] = hists['bkg'][key].group(cats, process, bkg_map)
            signal_hists[key] = hists['sig'][key].group(cats, process, sig_map)
            data_hists[key] = hists['data'][key].group(cats, process, data_map)

        hists={
            'bkg': bkg_hists,
            'sig': signal_hists,
            'data': data_hists
        }
        save(hists,'hists/darkhiggs'+year+'.limit')

    hists = load('hists/darkhiggs'+year+'.limit')
    data_hists   = hists['data']
    bkg_hists    = hists['bkg']
    signal_hists = hists['sig']

    binning = {
        'mass0': data_hists['recoil'].integrate('gentype','data').integrate('systematic','nominal').integrate('region','sr_mass0').integrate('process','MET').axis('recoil').edges(),
        'mass1': data_hists['recoil'].integrate('gentype','data').integrate('systematic','nominal').integrate('region','sr_mass1').integrate('process','MET').axis('recoil').edges(),
        'mass2': data_hists['recoil'].integrate('gentype','data').integrate('systematic','nominal').integrate('region','sr_mass2').integrate('process','MET').axis('recoil').edges(),
        'mass3': data_hists['recoil'].integrate('gentype','data').integrate('systematic','nominal').integrate('region','sr_mass3').integrate('process','MET').axis('recoil').edges(),
        'mass4': data_hists['recoil'].integrate('gentype','data').integrate('systematic','nominal').integrate('region','sr_mass4').integrate('process','MET').axis('recoil').edges()
    }    

    ###
    # Preparing histograms for fit
    ##

    data = {}
    for r in data_hists['recoil'].identifiers('region'):
        if category not in str(r): continue
        if mass not in str(r): continue
        print(r,category,mass)
        data[str(r).split("_")[0]]=data_hists['recoil'].integrate('region',r).integrate('gentype').rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning[mass]))

    background = {}
    for r in data.keys():
        background[r]=bkg_hists['recoil'].integrate('region',r).integrate('gentype').rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning[mass]))

    signal = {}
    for r in data.keys():
        signal[r]=signal_hists['recoil'].integrate('region',r).integrate('gentype').rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning[mass]))

    ###
    ###
    # Setting up systematics
    ###
    ###

    lumi = rl.NuisanceParameter('lumi', 'lnN')

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

    id_e = rl.NuisanceParameter('id_e', 'lnN')
    id_mu = rl.NuisanceParameter('id_mu', 'lnN')
    id_pho = rl.NuisanceParameter('id_pho', 'lnN')

    reco_e = rl.NuisanceParameter('reco_e', 'lnN')

    iso_mu = rl.NuisanceParameter('iso_mu', 'lnN')

    trig_e = rl.NuisanceParameter('trig_e', 'lnN')
    trig_met = rl.NuisanceParameter('trig_met', 'lnN')
    trig_pho = rl.NuisanceParameter('trig_pho', 'lnN')

    veto_tau = rl.NuisanceParameter('veto_tau', 'lnN')

    jec = rl.NuisanceParameter('jec', 'lnN')

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

    sr.setObservation(template(data['sr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))

    ###
    # Z(->nunu)+jets data-driven model
    ###

    sr_zjetsHist = background['sr'].integrate('process', 'Z+jets').integrate('systematic','nominal')
    sr_zjetsTemplate = template(sr_zjetsHist, 'recoil')
    sr_zjetsMC =  rl.TemplateSample(ch_name+'_zjetsMC', rl.Sample.BACKGROUND, sr_zjetsTemplate)
    sr_zjetsMC.setParamEffect(lumi, 1.027)
    sr_zjetsMC.setParamEffect(trig_met, 1.01)
    sr_zjetsMC.setParamEffect(veto_tau, 1.03)
    sr_zjetsMC.setParamEffect(zjets_norm, 1.4)
    sr_zjetsMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'Z+jets').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'Z+jets').integrate('systematic','btagDown')
    sr_zjetsMC.setParamEffect(btag, btagUp, btagDown)
    
    sr_zjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_zjets_bin_%d' % i, b, 0, sr_zjetsTemplate[0].max()*2) for i,b in enumerate(sr_zjetsTemplate[0])]) 
    sr_zjetsBinYields = sr_zjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    sr_zjetsObservable = rl.Observable('recoil', sr_zjetsHist.axis('recoil').edges())
    sr_zjets = rl.ParametericSample(ch_name+'_zjets', rl.Sample.BACKGROUND, sr_zjetsObservable, sr_zjetsBinYields)

    sr.addSample(sr_zjets)

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    sr_wjetsHist = background['sr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    sr_wjetsTemplate = template(sr_wjetsHist, 'recoil')
    sr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, sr_wjetsTemplate)
    sr_wjetsMC.setParamEffect(lumi, 1.027)
    sr_wjetsMC.setParamEffect(trig_met, 1.01)
    sr_wjetsMC.setParamEffect(veto_tau, 1.03)
    sr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    sr_wjetsMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'W+jets').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'W+jets').integrate('systematic','btagDown')
    sr_wjetsMC.setParamEffect(btag, btagUp, btagDown)

    sr_wjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_wjets_bin_%d' % i,b,0,sr_wjetsTemplate[0].max()*2) for i,b in enumerate(sr_wjetsTemplate[0])]) 
    sr_wjetsBinYields = sr_wjetsBinYields * hf_fraction_weight['0tag']['W+jets']
    sr_wjetsObservable = rl.Observable('recoil', sr_wjetsHist.axis('recoil').edges())
    sr_wjets = rl.ParametericSample(ch_name+'_wjets', rl.Sample.BACKGROUND, sr_wjetsObservable, sr_wjetsBinYields)

    sr.addSample(sr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    sr_ttHist = background['sr'].integrate('process', 'TT').integrate('systematic','nominal')
    sr_ttTemplate = template(sr_ttHist, 'recoil')
    sr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, sr_ttTemplate)
    sr_ttMC.setParamEffect(lumi, 1.027)
    sr_ttMC.setParamEffect(trig_met, 1.01)
    sr_ttMC.setParamEffect(veto_tau, 1.03)
    sr_ttMC.setParamEffect(tt_norm, 1.4)
    sr_ttMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'TT').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'TT').integrate('systematic','btagDown')
    sr_ttMC.setParamEffect(btag, btagUp, btagDown)

    sr_ttBinYields = np.array([rl.IndependentParameter(ch_name+'_tt_bin_%d' % i,b,0,sr_ttTemplate[0].max()*2) for i,b in enumerate(sr_ttTemplate[0])])
    sr_ttBinYields = sr_ttBinYields * deepak15_weight['0tag']['TT']
    sr_ttObservable = rl.Observable('recoil', sr_ttHist.axis('recoil').edges())
    sr_tt = rl.ParametericSample(ch_name+'_tt', rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields)

    sr.addSample(sr_tt)

    ###
    # Other MC-driven processes
    ###

    sr_stHist = background['sr'].integrate('process', 'ST').integrate('systematic','nominal')
    sr_stTemplate = template(sr_stHist, 'recoil')
    sr_stMC =  rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, sr_stTemplate)
    sr_stMC.setParamEffect(lumi, 1.027)
    sr_stMC.setParamEffect(trig_met, 1.01)
    sr_stMC.setParamEffect(veto_tau, 1.03)
    sr_stMC.setParamEffect(st_norm, 1.2)
    sr_stMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'ST').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'ST').integrate('systematic','btagDown')
    sr_stMC.setParamEffect(btag, btagUp, btagDown)

    sr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,sr_stTemplate[0].max()*2) for i,b in enumerate(sr_stTemplate[0])])
    sr_stBinYields = sr_stBinYields * deepak15_weight['0tag']['ST']
    sr_stObservable = rl.Observable('recoil', sr_stHist.axis('recoil').edges())
    sr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, sr_stObservable, sr_stBinYields)

    sr.addSample(sr_st)

    sr_dyjetsHist = background['sr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    sr_dyjetsTemplate = template(sr_dyjetsHist, 'recoil')
    sr_dyjetsMC =  rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, sr_dyjetsTemplate)
    sr_dyjetsMC.setParamEffect(lumi, 1.027)
    sr_dyjetsMC.setParamEffect(trig_met, 1.01)
    sr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    sr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    sr_dyjetsMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'DY+jets').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'DY+jets').integrate('systematic','btagDown')
    sr_dyjetsMC.setParamEffect(btag, btagUp, btagDown)

    sr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,sr_dyjetsTemplate[0].max()*2) for i,b in enumerate(sr_dyjetsTemplate[0])])
    sr_dyjetsBinYields = sr_dyjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    sr_dyjetsObservable = rl.Observable('recoil', sr_dyjetsHist.axis('recoil').edges())
    sr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, sr_dyjetsObservable, sr_dyjetsBinYields)

    sr.addSample(sr_dyjets)

    sr_vvHist = background['sr'].integrate('process', 'VV').integrate('systematic','nominal')
    sr_vvTemplate = template(sr_vvHist, 'recoil')
    sr_vvMC =  rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, sr_vvTemplate)
    sr_vvMC.setParamEffect(lumi, 1.027)
    sr_vvMC.setParamEffect(trig_met, 1.01)
    sr_vvMC.setParamEffect(veto_tau, 1.03)
    sr_vvMC.setParamEffect(vv_norm, 1.2)
    sr_vvMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'VV').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'VV').integrate('systematic','btagDown')
    sr_vvMC.setParamEffect(btag, btagUp, btagDown)

    sr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,sr_vvTemplate[0].max()*2) for i,b in enumerate(sr_vvTemplate[0])])
    sr_vvBinYields = sr_vvBinYields * deepak15_weight['0tag']['VV']
    sr_vvObservable = rl.Observable('recoil', sr_vvHist.axis('recoil').edges())
    sr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, sr_vvObservable, sr_vvBinYields)

    sr.addSample(sr_vv)

    sr_hbbHist = background['sr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    sr_hbbTemplate = template(sr_hbbHist, 'recoil')
    sr_hbbMC =  rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, sr_hbbTemplate)
    sr_hbbMC.setParamEffect(lumi, 1.027)
    sr_hbbMC.setParamEffect(trig_met, 1.01)
    sr_hbbMC.setParamEffect(veto_tau, 1.03)
    sr_hbbMC.setParamEffect(hbb_norm, 1.2)
    sr_hbbMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'Hbb').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'Hbb').integrate('systematic','btagDown')
    sr_hbbMC.setParamEffect(btag, btagUp, btagDown)

    sr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,sr_hbbTemplate[0].max()*2) for i,b in enumerate(sr_hbbTemplate[0])])
    sr_hbbBinYields = sr_hbbBinYields * deepak15_weight['0tag']['Hbb']
    sr_hbbObservable = rl.Observable('recoil', sr_hbbHist.axis('recoil').edges())
    sr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, sr_hbbObservable, sr_hbbBinYields)

    sr.addSample(sr_hbb)

    sr_qcdHist = background['sr'].integrate('process', 'QCD').integrate('systematic','nominal')
    sr_qcdTemplate = template(sr_qcdHist, 'recoil')
    sr_qcdMC =  rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, sr_qcdTemplate)
    sr_qcdMC.setParamEffect(lumi, 1.027)
    sr_qcdMC.setParamEffect(trig_met, 1.01)
    sr_qcdMC.setParamEffect(veto_tau, 1.03)
    sr_qcdMC.setParamEffect(qcdsig_norm, 2.0)
    sr_qcdMC.setParamEffect(jec, 1.05)
    btagUp=background['sr'].integrate('process', 'QCD').integrate('systematic','btagUp')
    btagDown=background['sr'].integrate('process', 'QCD').integrate('systematic','btagDown')
    sr_qcdMC.setParamEffect(btag, btagUp, btagDown)

    sr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,sr_qcdTemplate[0].max()*2) for i,b in enumerate(sr_qcdTemplate[0])])
    sr_qcdBinYields = sr_qcdBinYields * deepak15_weight['0tag']['QCD']
    sr_qcdObservable = rl.Observable('recoil', sr_qcdHist.axis('recoil').edges())
    sr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, sr_qcdObservable, sr_qcdBinYields)

    sr.addSample(sr_qcd)

    for s in signal['sr'].identifiers('process'):
        print(s)
        sr_signalHist = signal['sr'].integrate('process', s).integrate('systematic','nominal')
        sr_signalTemplate = template(sr_signalHist, 'recoil')
        sr_signalMC =  rl.TemplateSample(ch_name+'_signalMC', rl.Sample.BACKGROUND, sr_signalTemplate)
        sr_signalMC.setParamEffect(lumi, 1.027)
        sr_signalMC.setParamEffect(trig_met, 1.01)
        sr_signalMC.setParamEffect(veto_tau, 1.03)
        sr_signalMC.setParamEffect(jec, 1.05)
        btagUp=signal['sr'].integrate('process', s).integrate('systematic','btagUp')
        btagDown=signal['sr'].integrate('process', s).integrate('systematic','btagDown')
        sr_signalMC.setParamEffect(btag, btagUp, btagDown)

        sr_signalBinYields = np.array([rl.IndependentParameter(ch_name+'_signal_bin_%d' % i,b,0,sr_signalTemplate[0].max()*2) for i,b in enumerate(sr_signalTemplate[0])])
        sr_signalBinYields = sr_signalBinYields * signal_weight[str(s)]
        sr_signalObservable = rl.Observable('recoil', sr_signalHist.axis('recoil').edges())
        sr_signal = rl.ParametericSample(ch_name+'_signal', rl.Sample.BACKGROUND, sr_signalObservable, sr_signalBinYields)

    ###
    # End of SR
    ###

    ###
    ###
    # Single muon W control region
    ###
    ###

    ch_name = 'wmcr-'+mass+'-'+category
    wmcr = rl.Channel(ch_name)
    model.addChannel(wmcr)

    ###
    # Add data distribution to the channel
    ###

    wmcr.setObservation(template(data['wmcr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    wmcr_wjetsHist = background['wmcr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    wmcr_wjetsTemplate = template(wmcr_wjetsHist, 'recoil')
    wmcr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, wmcr_wjetsTemplate)
    wmcr_wjetsMC.setParamEffect(lumi, 1.027)
    wmcr_wjetsMC.setParamEffect(trig_met, 1.01)
    wmcr_wjetsMC.setParamEffect(veto_tau, 1.03)
    wmcr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    wmcr_wjetsMC.setParamEffect(jec, 1.05)
    wmcr_wjetsMC.setParamEffect(id_mu, 1.02)
    wmcr_wjetsMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['wmcr'].integrate('process', 'W+jets').integrate('systematic','btagUp')
    btagDown=background['wmcr'].integrate('process', 'W+jets').integrate('systematic','btagDown')
    wmcr_wjetsMC.setParamEffect(btag, btagUp, btagDown)

    wmcr_wjetsTransferFactor = wmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['0tag']['W+jets']
    wmcr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, wmcr_wjetsTransferFactor, sr_wjets)

    wmcr.addSample(wmcr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    wmcr_ttHist = background['wmcr'].integrate('process', 'TT').integrate('systematic','nominal')
    wmcr_ttTemplate = template(wmcr_ttHist, 'recoil')
    wmcr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, wmcr_ttTemplate)
    wmcr_ttMC.setParamEffect(lumi, 1.027)
    wmcr_ttMC.setParamEffect(trig_met, 1.01)
    wmcr_ttMC.setParamEffect(veto_tau, 1.03)
    wmcr_ttMC.setParamEffect(tt_norm, 1.4)
    wmcr_ttMC.setParamEffect(jec, 1.05)
    wmcr_ttMC.setParamEffect(id_mu, 1.02)
    wmcr_ttMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['wmcr'].integrate('process', 'TT').integrate('systematic','btagUp')
    btagDown=background['wmcr'].integrate('process', 'TT').integrate('systematic','btagDown')
    wmcr_ttMC.setParamEffect(btag, btagUp, btagDown)

    wmcr_ttTransferFactor = wmcr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['0tag']['TT']
    wmcr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt)

    wmcr.addSample(wmcr_tt)

    ###
    # Other MC-driven processes
    ###

    wmcr_stHist = background['wmcr'].integrate('process', 'ST').integrate('systematic','nominal')
    wmcr_stTemplate = template(wmcr_stHist, 'recoil')
    wmcr_stMC =  rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, wmcr_stTemplate)
    wmcr_stMC.setParamEffect(lumi, 1.027)
    wmcr_stMC.setParamEffect(trig_met, 1.01)
    wmcr_stMC.setParamEffect(veto_tau, 1.03)
    wmcr_stMC.setParamEffect(st_norm, 1.2)
    wmcr_stMC.setParamEffect(jec, 1.05)
    wmcr_stMC.setParamEffect(id_mu, 1.02)
    wmcr_stMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['wmcr'].integrate('process', 'ST').integrate('systematic','btagUp')
    btagDown=background['wmcr'].integrate('process', 'ST').integrate('systematic','btagDown')
    wmcr_stMC.setParamEffect(btag, btagUp, btagDown)

    wmcr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,wmcr_stTemplate[0].max()*2) for i,b in enumerate(wmcr_stTemplate[0])])
    wmcr_stBinYields = wmcr_stBinYields * deepak15_weight['0tag']['ST']
    wmcr_stObservable = rl.Observable('recoil', wmcr_stHist.axis('recoil').edges())
    wmcr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, wmcr_stObservable, wmcr_stBinYields)

    wmcr.addSample(wmcr_st)

    wmcr_dyjetsHist = background['wmcr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    wmcr_dyjetsTemplate = template(wmcr_dyjetsHist, 'recoil')
    wmcr_dyjetsMC =  rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, wmcr_dyjetsTemplate)
    wmcr_dyjetsMC.setParamEffect(lumi, 1.027)
    wmcr_dyjetsMC.setParamEffect(trig_met, 1.01)
    wmcr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    wmcr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    wmcr_dyjetsMC.setParamEffect(jec, 1.05)
    wmcr_dyjetsMC.setParamEffect(id_mu, 1.02)
    wmcr_dyjetsMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['wmcr'].integrate('process', 'DY+jets').integrate('systematic','btagUp')
    btagDown=background['wmcr'].integrate('process', 'DY+jets').integrate('systematic','btagDown')
    wmcr_dyjetsMC.setParamEffect(btag, btagUp, btagDown)

    wmcr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,wmcr_dyjetsTemplate[0].max()*2) for i,b in enumerate(wmcr_dyjetsTemplate[0])])
    wmcr_dyjetsBinYields = wmcr_dyjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    wmcr_dyjetsObservable = rl.Observable('recoil', wmcr_dyjetsHist.axis('recoil').edges())
    wmcr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, wmcr_dyjetsObservable, wmcr_dyjetsBinYields)

    wmcr.addSample(wmcr_dyjets)

    wmcr_vvHist = background['wmcr'].integrate('process', 'VV').integrate('systematic','nominal')
    wmcr_vvTemplate = template(wmcr_vvHist, 'recoil')
    wmcr_vvMC =  rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, wmcr_vvTemplate)
    wmcr_vvMC.setParamEffect(lumi, 1.027)
    wmcr_vvMC.setParamEffect(trig_met, 1.01)
    wmcr_vvMC.setParamEffect(veto_tau, 1.03)
    wmcr_vvMC.setParamEffect(vv_norm, 1.2)
    wmcr_vvMC.setParamEffect(jec, 1.05)
    wmcr_vvMC.setParamEffect(id_mu, 1.02)
    wmcr_vvMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['wmcr'].integrate('process', 'VV').integrate('systematic','btagUp')
    btagDown=background['wmcr'].integrate('process', 'VV').integrate('systematic','btagDown')
    wmcr_vvMC.setParamEffect(btag, btagUp, btagDown)

    wmcr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,wmcr_vvTemplate[0].max()*2) for i,b in enumerate(wmcr_vvTemplate[0])])
    wmcr_vvBinYields = wmcr_vvBinYields * deepak15_weight['0tag']['VV']
    wmcr_vvObservable = rl.Observable('recoil', wmcr_vvHist.axis('recoil').edges())
    wmcr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, wmcr_vvObservable, wmcr_vvBinYields)

    wmcr.addSample(wmcr_vv)

    wmcr_hbbHist = background['wmcr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    wmcr_hbbTemplate = template(wmcr_hbbHist, 'recoil')
    wmcr_hbbMC =  rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, wmcr_hbbTemplate)
    wmcr_hbbMC.setParamEffect(lumi, 1.027)
    wmcr_hbbMC.setParamEffect(trig_met, 1.01)
    wmcr_hbbMC.setParamEffect(veto_tau, 1.03)
    wmcr_hbbMC.setParamEffect(hbb_norm, 1.2)
    wmcr_hbbMC.setParamEffect(jec, 1.05)
    wmcr_hbbMC.setParamEffect(id_mu, 1.02)
    wmcr_hbbMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['wmcr'].integrate('process', 'Hbb').integrate('systematic','btagUp')
    btagDown=background['wmcr'].integrate('process', 'Hbb').integrate('systematic','btagDown')
    wmcr_hbbMC.setParamEffect(btag, btagUp, btagDown)

    wmcr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,wmcr_hbbTemplate[0].max()*2) for i,b in enumerate(wmcr_hbbTemplate[0])])
    wmcr_hbbBinYields = wmcr_hbbBinYields * deepak15_weight['0tag']['Hbb']
    wmcr_hbbObservable = rl.Observable('recoil', wmcr_hbbHist.axis('recoil').edges())
    wmcr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, wmcr_hbbObservable, wmcr_hbbBinYields)

    wmcr.addSample(wmcr_hbb)

    wmcr_qcdHist = background['wmcr'].integrate('process', 'QCD').integrate('systematic','nominal')
    wmcr_qcdTemplate = template(wmcr_qcdHist, 'recoil')
    wmcr_qcdMC =  rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, wmcr_qcdTemplate)
    wmcr_qcdMC.setParamEffect(lumi, 1.027)
    wmcr_qcdMC.setParamEffect(trig_met, 1.01)
    wmcr_qcdMC.setParamEffect(veto_tau, 1.03)
    wmcr_qcdMC.setParamEffect(qcdmu_norm, 2.0)
    wmcr_qcdMC.setParamEffect(jec, 1.05)
    wmcr_qcdMC.setParamEffect(id_mu, 1.02)
    wmcr_qcdMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['wmcr'].integrate('process', 'QCD').integrate('systematic','btagUp')
    btagDown=background['wmcr'].integrate('process', 'QCD').integrate('systematic','btagDown')
    wmcr_qcdMC.setParamEffect(btag, btagUp, btagDown)

    wmcr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,wmcr_qcdTemplate[0].max()*2) for i,b in enumerate(wmcr_qcdTemplate[0])])
    wmcr_qcdBinYields = wmcr_qcdBinYields * deepak15_weight['0tag']['QCD']
    wmcr_qcdObservable = rl.Observable('recoil', wmcr_qcdHist.axis('recoil').edges())
    wmcr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, wmcr_qcdObservable, wmcr_qcdBinYields)

    wmcr.addSample(wmcr_qcd)

    ###
    # End of single muon W control region
    ###

    ###
    ###
    # Single muon top control region
    ###
    ###

    ch_name = 'tmcr-'+mass+'-'+category
    tmcr = rl.Channel(ch_name)
    model.addChannel(tmcr)

    ###
    # Add data distribution to the channel
    ###

    tmcr.setObservation(template(data['tmcr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    tmcr_wjetsHist = background['tmcr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    tmcr_wjetsTemplate = template(tmcr_wjetsHist, 'recoil')
    tmcr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, tmcr_wjetsTemplate)
    tmcr_wjetsMC.setParamEffect(lumi, 1.027)
    tmcr_wjetsMC.setParamEffect(trig_met, 1.01)
    tmcr_wjetsMC.setParamEffect(veto_tau, 1.03)
    tmcr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    tmcr_wjetsMC.setParamEffect(jec, 1.05)
    tmcr_wjetsMC.setParamEffect(id_mu, 1.02)
    tmcr_wjetsMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['tmcr'].integrate('process', 'W+jets').integrate('systematic','btagUp')
    btagDown=background['tmcr'].integrate('process', 'W+jets').integrate('systematic','btagDown')
    tmcr_wjetsMC.setParamEffect(btag, btagUp, btagDown)

    tmcr_wjetsTransferFactor = tmcr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['1tag']['W+jets']
    tmcr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, tmcr_wjetsTransferFactor, sr_wjets)

    tmcr.addSample(tmcr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    tmcr_ttHist = background['tmcr'].integrate('process', 'TT').integrate('systematic','nominal')
    tmcr_ttTemplate = template(tmcr_ttHist, 'recoil')
    tmcr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, tmcr_ttTemplate)
    tmcr_ttMC.setParamEffect(lumi, 1.027)
    tmcr_ttMC.setParamEffect(trig_met, 1.01)
    tmcr_ttMC.setParamEffect(veto_tau, 1.03)
    tmcr_ttMC.setParamEffect(tt_norm, 1.4)
    tmcr_ttMC.setParamEffect(jec, 1.05)
    tmcr_ttMC.setParamEffect(id_mu, 1.02)
    tmcr_ttMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['tmcr'].integrate('process', 'TT').integrate('systematic','btagUp')
    btagDown=background['tmcr'].integrate('process', 'TT').integrate('systematic','btagDown')
    tmcr_ttMC.setParamEffect(btag, btagUp, btagDown)

    tmcr_ttTransferFactor = tmcr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['1tag']['TT']
    tmcr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt)

    tmcr.addSample(tmcr_tt)

    ###
    # Other MC-driven processes
    ###

    tmcr_stHist = background['tmcr'].integrate('process', 'ST').integrate('systematic','nominal')
    tmcr_stTemplate = template(tmcr_stHist, 'recoil')
    tmcr_stMC =  rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, tmcr_stTemplate)
    tmcr_stMC.setParamEffect(lumi, 1.027)
    tmcr_stMC.setParamEffect(trig_met, 1.01)
    tmcr_stMC.setParamEffect(veto_tau, 1.03)
    tmcr_stMC.setParamEffect(st_norm, 1.2)
    tmcr_stMC.setParamEffect(jec, 1.05)
    tmcr_stMC.setParamEffect(id_mu, 1.02)
    tmcr_stMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['tmcr'].integrate('process', 'ST').integrate('systematic','btagUp')
    btagDown=background['tmcr'].integrate('process', 'ST').integrate('systematic','btagDown')
    tmcr_stMC.setParamEffect(btag, btagUp, btagDown)

    tmcr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,tmcr_stTemplate[0].max()*2) for i,b in enumerate(tmcr_stTemplate[0])])
    tmcr_stBinYields = tmcr_stBinYields * deepak15_weight['1tag']['ST']
    tmcr_stObservable = rl.Observable('recoil', tmcr_stHist.axis('recoil').edges())
    tmcr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, tmcr_stObservable, tmcr_stBinYields)

    tmcr.addSample(tmcr_st)

    tmcr_dyjetsHist = background['tmcr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    tmcr_dyjetsTemplate = template(tmcr_dyjetsHist, 'recoil')
    tmcr_dyjetsMC =  rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, tmcr_dyjetsTemplate)
    tmcr_dyjetsMC.setParamEffect(lumi, 1.027)
    tmcr_dyjetsMC.setParamEffect(trig_met, 1.01)
    tmcr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    tmcr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    tmcr_dyjetsMC.setParamEffect(jec, 1.05)
    tmcr_dyjetsMC.setParamEffect(id_mu, 1.02)
    tmcr_dyjetsMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['tmcr'].integrate('process', 'DY+jets').integrate('systematic','btagUp')
    btagDown=background['tmcr'].integrate('process', 'DY+jets').integrate('systematic','btagDown')
    tmcr_dyjetsMC.setParamEffect(btag, btagUp, btagDown)

    tmcr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,tmcr_dyjetsTemplate[0].max()*2) for i,b in enumerate(tmcr_dyjetsTemplate[0])])
    tmcr_dyjetsBinYields = tmcr_dyjetsBinYields * hf_fraction_weight['1tag']['Z+jets']
    tmcr_dyjetsObservable = rl.Observable('recoil', tmcr_dyjetsHist.axis('recoil').edges())
    tmcr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, tmcr_dyjetsObservable, tmcr_dyjetsBinYields)

    tmcr.addSample(tmcr_dyjets)

    tmcr_vvHist = background['tmcr'].integrate('process', 'VV').integrate('systematic','nominal')
    tmcr_vvTemplate = template(tmcr_vvHist, 'recoil')
    tmcr_vvMC =  rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, tmcr_vvTemplate)
    tmcr_vvMC.setParamEffect(lumi, 1.027)
    tmcr_vvMC.setParamEffect(trig_met, 1.01)
    tmcr_vvMC.setParamEffect(veto_tau, 1.03)
    tmcr_vvMC.setParamEffect(vv_norm, 1.2)
    tmcr_vvMC.setParamEffect(jec, 1.05)
    tmcr_vvMC.setParamEffect(id_mu, 1.02)
    tmcr_vvMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['tmcr'].integrate('process', 'VV').integrate('systematic','btagUp')
    btagDown=background['tmcr'].integrate('process', 'VV').integrate('systematic','btagDown')
    tmcr_vvMC.setParamEffect(btag, btagUp, btagDown)

    tmcr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,tmcr_vvTemplate[0].max()*2) for i,b in enumerate(tmcr_vvTemplate[0])])
    tmcr_vvBinYields = tmcr_vvBinYields * deepak15_weight['1tag']['VV']
    tmcr_vvObservable = rl.Observable('recoil', tmcr_vvHist.axis('recoil').edges())
    tmcr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, tmcr_vvObservable, tmcr_vvBinYields)

    tmcr.addSample(tmcr_vv)

    tmcr_hbbHist = background['tmcr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    tmcr_hbbTemplate = template(tmcr_hbbHist, 'recoil')
    tmcr_hbbMC =  rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, tmcr_hbbTemplate)
    tmcr_hbbMC.setParamEffect(lumi, 1.027)
    tmcr_hbbMC.setParamEffect(trig_met, 1.01)
    tmcr_hbbMC.setParamEffect(veto_tau, 1.03)
    tmcr_hbbMC.setParamEffect(hbb_norm, 1.2)
    tmcr_hbbMC.setParamEffect(jec, 1.05)
    tmcr_hbbMC.setParamEffect(id_mu, 1.02)
    tmcr_hbbMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['tmcr'].integrate('process', 'Hbb').integrate('systematic','btagUp')
    btagDown=background['tmcr'].integrate('process', 'Hbb').integrate('systematic','btagDown')
    tmcr_hbbMC.setParamEffect(btag, btagUp, btagDown)

    tmcr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,tmcr_hbbTemplate[0].max()*2) for i,b in enumerate(tmcr_hbbTemplate[0])])
    tmcr_hbbBinYields = tmcr_hbbBinYields * deepak15_weight['1tag']['Hbb']
    tmcr_hbbObservable = rl.Observable('recoil', tmcr_hbbHist.axis('recoil').edges())
    tmcr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, tmcr_hbbObservable, tmcr_hbbBinYields)

    tmcr.addSample(tmcr_hbb)

    tmcr_qcdHist = background['tmcr'].integrate('process', 'QCD').integrate('systematic','nominal')
    tmcr_qcdTemplate = template(tmcr_qcdHist, 'recoil')
    tmcr_qcdMC =  rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, tmcr_qcdTemplate)
    tmcr_qcdMC.setParamEffect(lumi, 1.027)
    tmcr_qcdMC.setParamEffect(trig_met, 1.01)
    tmcr_qcdMC.setParamEffect(veto_tau, 1.03)
    tmcr_qcdMC.setParamEffect(qcdmu_norm, 2.0)
    tmcr_qcdMC.setParamEffect(jec, 1.05)
    tmcr_qcdMC.setParamEffect(id_mu, 1.02)
    tmcr_qcdMC.setParamEffect(iso_mu, 1.02)
    btagUp=background['tmcr'].integrate('process', 'QCD').integrate('systematic','btagUp')
    btagDown=background['tmcr'].integrate('process', 'QCD').integrate('systematic','btagDown')
    tmcr_qcdMC.setParamEffect(btag, btagUp, btagDown)

    tmcr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,tmcr_qcdTemplate[0].max()*2) for i,b in enumerate(tmcr_qcdTemplate[0])])
    tmcr_qcdBinYields = tmcr_qcdBinYields * deepak15_weight['1tag']['QCD']
    tmcr_qcdObservable = rl.Observable('recoil', tmcr_qcdHist.axis('recoil').edges())
    tmcr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, tmcr_qcdObservable, tmcr_qcdBinYields)

    tmcr.addSample(tmcr_qcd)

    ###
    # End of single muon top control region
    ###

    ###
    ###
    # Single electron W control region
    ###
    ###

    ch_name = 'wecr-'+mass+'-'+category
    wecr = rl.Channel(ch_name)
    model.addChannel(wecr)

    ###
    # Add data distribution to the channel
    ###

    wecr.setObservation(template(data['wecr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    wecr_wjetsHist = background['wecr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    wecr_wjetsTemplate = template(wecr_wjetsHist, 'recoil')
    wecr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, wecr_wjetsTemplate)
    wecr_wjetsMC.setParamEffect(lumi, 1.027)
    wecr_wjetsMC.setParamEffect(trig_e, 1.01)
    wecr_wjetsMC.setParamEffect(veto_tau, 1.03)
    wecr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    wecr_wjetsMC.setParamEffect(jec, 1.05)
    wecr_wjetsMC.setParamEffect(id_e, 1.02)
    wecr_wjetsMC.setParamEffect(reco_e, 1.02)
    btagUp=background['wecr'].integrate('process', 'W+jets').integrate('systematic','btagUp')
    btagDown=background['wecr'].integrate('process', 'W+jets').integrate('systematic','btagDown')
    wecr_wjetsMC.setParamEffect(btag, btagUp, btagDown)

    wecr_wjetsTransferFactor = wecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['0tag']['W+jets']
    wecr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, wecr_wjetsTransferFactor, sr_wjets)

    wecr.addSample(wecr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    wecr_ttHist = background['wecr'].integrate('process', 'TT').integrate('systematic','nominal')
    wecr_ttTemplate = template(wecr_ttHist, 'recoil')
    wecr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, wecr_ttTemplate)
    wecr_ttMC.setParamEffect(lumi, 1.027)
    wecr_ttMC.setParamEffect(trig_e, 1.01)
    wecr_ttMC.setParamEffect(veto_tau, 1.03)
    wecr_ttMC.setParamEffect(tt_norm, 1.4)
    wecr_ttMC.setParamEffect(jec, 1.05)
    wecr_ttMC.setParamEffect(id_e, 1.02)
    wecr_ttMC.setParamEffect(reco_e, 1.02)
    btagUp=background['wecr'].integrate('process', 'TT').integrate('systematic','btagUp')
    btagDown=background['wecr'].integrate('process', 'TT').integrate('systematic','btagDown')
    wecr_ttMC.setParamEffect(btag, btagUp, btagDown)

    wecr_ttTransferFactor = wecr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['0tag']['TT']
    wecr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt)

    wecr.addSample(wecr_tt)

    ###
    # Other MC-driven processes
    ###

    wecr_stHist = background['wecr'].integrate('process', 'ST').integrate('systematic','nominal')
    wecr_stTemplate = template(wecr_stHist, 'recoil')
    wecr_stMC =  rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, wecr_stTemplate)
    wecr_stMC.setParamEffect(lumi, 1.027)
    wecr_stMC.setParamEffect(trig_e, 1.01)
    wecr_stMC.setParamEffect(veto_tau, 1.03)
    wecr_stMC.setParamEffect(st_norm, 1.2)
    wecr_stMC.setParamEffect(jec, 1.05)
    wecr_stMC.setParamEffect(id_e, 1.02)
    wecr_stMC.setParamEffect(reco_e, 1.02)
    btagUp=background['wecr'].integrate('process', 'ST').integrate('systematic','btagUp')
    btagDown=background['wecr'].integrate('process', 'ST').integrate('systematic','btagDown')
    wecr_stMC.setParamEffect(btag, btagUp, btagDown)

    wecr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,wecr_stTemplate[0].max()*2) for i,b in enumerate(wecr_stTemplate[0])])
    wecr_stBinYields = wecr_stBinYields * deepak15_weight['0tag']['ST']
    wecr_stObservable = rl.Observable('recoil', wecr_stHist.axis('recoil').edges())
    wecr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, wecr_stObservable, wecr_stBinYields)

    wecr.addSample(wecr_st)

    wecr_dyjetsHist = background['wecr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    wecr_dyjetsTemplate = template(wecr_dyjetsHist, 'recoil')
    wecr_dyjetsMC =  rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, wecr_dyjetsTemplate)
    wecr_dyjetsMC.setParamEffect(lumi, 1.027)
    wecr_dyjetsMC.setParamEffect(trig_e, 1.01)
    wecr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    wecr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    wecr_dyjetsMC.setParamEffect(jec, 1.05)
    wecr_dyjetsMC.setParamEffect(id_e, 1.02)
    wecr_dyjetsMC.setParamEffect(reco_e, 1.02)
    btagUp=background['wecr'].integrate('process', 'DY+jets').integrate('systematic','btagUp')
    btagDown=background['wecr'].integrate('process', 'DY+jets').integrate('systematic','btagDown')
    wecr_dyjetsMC.setParamEffect(btag, btagUp, btagDown)

    wecr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,wecr_dyjetsTemplate[0].max()*2) for i,b in enumerate(wecr_dyjetsTemplate[0])])
    wecr_dyjetsBinYields = wecr_dyjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    wecr_dyjetsObservable = rl.Observable('recoil', wecr_dyjetsHist.axis('recoil').edges())
    wecr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, wecr_dyjetsObservable, wecr_dyjetsBinYields)

    wecr.addSample(wecr_dyjets)

    wecr_vvHist = background['wecr'].integrate('process', 'VV').integrate('systematic','nominal')
    wecr_vvTemplate = template(wecr_vvHist, 'recoil')
    wecr_vvMC =  rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, wecr_vvTemplate)
    wecr_vvMC.setParamEffect(lumi, 1.027)
    wecr_vvMC.setParamEffect(trig_e, 1.01)
    wecr_vvMC.setParamEffect(veto_tau, 1.03)
    wecr_vvMC.setParamEffect(vv_norm, 1.2)
    wecr_vvMC.setParamEffect(jec, 1.05)
    wecr_vvMC.setParamEffect(id_e, 1.02)
    wecr_vvMC.setParamEffect(reco_e, 1.02)
    btagUp=background['wecr'].integrate('process', 'VV').integrate('systematic','btagUp')
    btagDown=background['wecr'].integrate('process', 'VV').integrate('systematic','btagDown')
    wecr_vvMC.setParamEffect(btag, btagUp, btagDown)

    wecr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,wecr_vvTemplate[0].max()*2) for i,b in enumerate(wecr_vvTemplate[0])])
    wecr_vvBinYields = wecr_vvBinYields * deepak15_weight['0tag']['VV']
    wecr_vvObservable = rl.Observable('recoil', wecr_vvHist.axis('recoil').edges())
    wecr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, wecr_vvObservable, wecr_vvBinYields)

    wecr.addSample(wecr_vv)

    wecr_hbbHist = background['wecr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    wecr_hbbTemplate = template(wecr_hbbHist, 'recoil')
    wecr_hbbMC =  rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, wecr_hbbTemplate)
    wecr_hbbMC.setParamEffect(lumi, 1.027)
    wecr_hbbMC.setParamEffect(trig_e, 1.01)
    wecr_hbbMC.setParamEffect(veto_tau, 1.03)
    wecr_hbbMC.setParamEffect(hbb_norm, 1.2)
    wecr_hbbMC.setParamEffect(jec, 1.05)
    wecr_hbbMC.setParamEffect(id_e, 1.02)
    wecr_hbbMC.setParamEffect(reco_e, 1.02)
    btagUp=background['wecr'].integrate('process', 'Hbb').integrate('systematic','btagUp')
    btagDown=background['wecr'].integrate('process', 'Hbb').integrate('systematic','btagDown')
    wecr_hbbMC.setParamEffect(btag, btagUp, btagDown)

    wecr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,wecr_hbbTemplate[0].max()*2) for i,b in enumerate(wecr_hbbTemplate[0])])
    wecr_hbbBinYields = wecr_hbbBinYields * deepak15_weight['0tag']['Hbb']
    wecr_hbbObservable = rl.Observable('recoil', wecr_hbbHist.axis('recoil').edges())
    wecr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, wecr_hbbObservable, wecr_hbbBinYields)

    wecr.addSample(wecr_hbb)

    wecr_qcdHist = background['wecr'].integrate('process', 'QCD').integrate('systematic','nominal')
    wecr_qcdTemplate = template(wecr_qcdHist, 'recoil')
    wecr_qcdMC =  rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, wecr_qcdTemplate)
    wecr_qcdMC.setParamEffect(lumi, 1.027)
    wecr_qcdMC.setParamEffect(trig_e, 1.01)
    wecr_qcdMC.setParamEffect(veto_tau, 1.03)
    wecr_qcdMC.setParamEffect(qcdmu_norm, 2.0)
    wecr_qcdMC.setParamEffect(jec, 1.05)
    wecr_qcdMC.setParamEffect(id_e, 1.02)
    wecr_qcdMC.setParamEffect(reco_e, 1.02)
    btagUp=background['wecr'].integrate('process', 'QCD').integrate('systematic','btagUp')
    btagDown=background['wecr'].integrate('process', 'QCD').integrate('systematic','btagDown')
    wecr_qcdMC.setParamEffect(btag, btagUp, btagDown)

    wecr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,wecr_qcdTemplate[0].max()*2) for i,b in enumerate(wecr_qcdTemplate[0])])
    wecr_qcdBinYields = wecr_qcdBinYields * deepak15_weight['0tag']['QCD']
    wecr_qcdObservable = rl.Observable('recoil', wecr_qcdHist.axis('recoil').edges())
    wecr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, wecr_qcdObservable, wecr_qcdBinYields)

    wecr.addSample(wecr_qcd)

    ###
    # End of single electron W control region
    ###

    ###
    ###
    # Single electron top control region
    ###
    ###

    ch_name = 'tecr-'+mass+'-'+category
    tecr = rl.Channel(ch_name)
    model.addChannel(tecr)

    ###
    # Add data distribution to the channel
    ###

    tecr.setObservation(template(data['tecr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    tecr_wjetsHist = background['tecr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    tecr_wjetsTemplate = template(tecr_wjetsHist, 'recoil')
    tecr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, tecr_wjetsTemplate)
    tecr_wjetsMC.setParamEffect(lumi, 1.027)
    tecr_wjetsMC.setParamEffect(trig_e, 1.01)
    tecr_wjetsMC.setParamEffect(veto_tau, 1.03)
    tecr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    tecr_wjetsMC.setParamEffect(jec, 1.05)
    tecr_wjetsMC.setParamEffect(id_e, 1.02)
    tecr_wjetsMC.setParamEffect(reco_e, 1.02)
    btagUp=background['tecr'].integrate('process', 'W+jets').integrate('systematic','btagUp')
    btagDown=background['tecr'].integrate('process', 'W+jets').integrate('systematic','btagDown')
    tecr_wjetsMC.setParamEffect(btag, btagUp, btagDown)

    tecr_wjetsTransferFactor = tecr_wjetsMC.getExpectation() / sr_wjetsMC.getExpectation() * hf_fraction_weight['1tag']['W+jets']
    tecr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, tecr_wjetsTransferFactor, sr_wjets)

    tecr.addSample(tecr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    tecr_ttHist = background['tecr'].integrate('process', 'TT').integrate('systematic','nominal')
    tecr_ttTemplate = template(tecr_ttHist, 'recoil')
    tecr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, tecr_ttTemplate)
    tecr_ttMC.setParamEffect(lumi, 1.027)
    tecr_ttMC.setParamEffect(trig_e, 1.01)
    tecr_ttMC.setParamEffect(veto_tau, 1.03)
    tecr_ttMC.setParamEffect(tt_norm, 1.4)
    tecr_ttMC.setParamEffect(jec, 1.05)
    tecr_ttMC.setParamEffect(id_e, 1.02)
    tecr_ttMC.setParamEffect(reco_e, 1.02)
    btagUp=background['tecr'].integrate('process', 'TT').integrate('systematic','btagUp')
    btagDown=background['tecr'].integrate('process', 'TT').integrate('systematic','btagDown')
    tecr_ttMC.setParamEffect(btag, btagUp, btagDown)

    tecr_ttTransferFactor = tecr_ttMC.getExpectation() / sr_ttMC.getExpectation() * deepak15_weight['1tag']['TT']
    tecr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt)

    tecr.addSample(tecr_tt)

    ###
    # Other MC-driven processes
    ###

    tecr_stHist = background['tecr'].integrate('process', 'ST').integrate('systematic','nominal')
    tecr_stTemplate = template(tecr_stHist, 'recoil')
    tecr_stMC =  rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, tecr_stTemplate)
    tecr_stMC.setParamEffect(lumi, 1.027)
    tecr_stMC.setParamEffect(trig_e, 1.01)
    tecr_stMC.setParamEffect(veto_tau, 1.03)
    tecr_stMC.setParamEffect(st_norm, 1.2)
    tecr_stMC.setParamEffect(jec, 1.05)
    tecr_stMC.setParamEffect(id_e, 1.02)
    tecr_stMC.setParamEffect(reco_e, 1.02)
    btagUp=background['tecr'].integrate('process', 'ST').integrate('systematic','btagUp')
    btagDown=background['tecr'].integrate('process', 'ST').integrate('systematic','btagDown')
    tecr_stMC.setParamEffect(btag, btagUp, btagDown)

    tecr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,tecr_stTemplate[0].max()*2) for i,b in enumerate(tecr_stTemplate[0])])
    tecr_stBinYields = tecr_stBinYields * deepak15_weight['1tag']['ST']
    tecr_stObservable = rl.Observable('recoil', tecr_stHist.axis('recoil').edges())
    tecr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, tecr_stObservable, tecr_stBinYields)

    tecr.addSample(tecr_st)

    tecr_dyjetsHist = background['tecr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    tecr_dyjetsTemplate = template(tecr_dyjetsHist, 'recoil')
    tecr_dyjetsMC =  rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, tecr_dyjetsTemplate)
    tecr_dyjetsMC.setParamEffect(lumi, 1.027)
    tecr_dyjetsMC.setParamEffect(trig_e, 1.01)
    tecr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    tecr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    tecr_dyjetsMC.setParamEffect(jec, 1.05)
    tecr_dyjetsMC.setParamEffect(id_e, 1.02)
    tecr_dyjetsMC.setParamEffect(reco_e, 1.02)
    btagUp=background['tecr'].integrate('process', 'DY+jets').integrate('systematic','btagUp')
    btagDown=background['tecr'].integrate('process', 'DY+jets').integrate('systematic','btagDown')
    tecr_dyjetsMC.setParamEffect(btag, btagUp, btagDown)

    tecr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,tecr_dyjetsTemplate[0].max()*2) for i,b in enumerate(tecr_dyjetsTemplate[0])])
    tecr_dyjetsBinYields = tecr_dyjetsBinYields * hf_fraction_weight['1tag']['Z+jets']
    tecr_dyjetsObservable = rl.Observable('recoil', tecr_dyjetsHist.axis('recoil').edges())
    tecr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, tecr_dyjetsObservable, tecr_dyjetsBinYields)

    tecr.addSample(tecr_dyjets)

    tecr_vvHist = background['tecr'].integrate('process', 'VV').integrate('systematic','nominal')
    tecr_vvTemplate = template(tecr_vvHist, 'recoil')
    tecr_vvMC =  rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, tecr_vvTemplate)
    tecr_vvMC.setParamEffect(lumi, 1.027)
    tecr_vvMC.setParamEffect(trig_e, 1.01)
    tecr_vvMC.setParamEffect(veto_tau, 1.03)
    tecr_vvMC.setParamEffect(vv_norm, 1.2)
    tecr_vvMC.setParamEffect(jec, 1.05)
    tecr_vvMC.setParamEffect(id_e, 1.02)
    tecr_vvMC.setParamEffect(reco_e, 1.02)
    btagUp=background['tecr'].integrate('process', 'VV').integrate('systematic','btagUp')
    btagDown=background['tecr'].integrate('process', 'VV').integrate('systematic','btagDown')
    tecr_vvMC.setParamEffect(btag, btagUp, btagDown)

    tecr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,tecr_vvTemplate[0].max()*2) for i,b in enumerate(tecr_vvTemplate[0])])
    tecr_vvBinYields = tecr_vvBinYields * deepak15_weight['1tag']['VV']
    tecr_vvObservable = rl.Observable('recoil', tecr_vvHist.axis('recoil').edges())
    tecr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, tecr_vvObservable, tecr_vvBinYields)

    tecr.addSample(tecr_vv)

    tecr_hbbHist = background['tecr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    tecr_hbbTemplate = template(tecr_hbbHist, 'recoil')
    tecr_hbbMC =  rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, tecr_hbbTemplate)
    tecr_hbbMC.setParamEffect(lumi, 1.027)
    tecr_hbbMC.setParamEffect(trig_e, 1.01)
    tecr_hbbMC.setParamEffect(veto_tau, 1.03)
    tecr_hbbMC.setParamEffect(hbb_norm, 1.2)
    tecr_hbbMC.setParamEffect(jec, 1.05)
    tecr_hbbMC.setParamEffect(id_e, 1.02)
    tecr_hbbMC.setParamEffect(reco_e, 1.02)
    btagUp=background['tecr'].integrate('process', 'Hbb').integrate('systematic','btagUp')
    btagDown=background['tecr'].integrate('process', 'Hbb').integrate('systematic','btagDown')
    tecr_hbbMC.setParamEffect(btag, btagUp, btagDown)

    tecr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,tecr_hbbTemplate[0].max()*2) for i,b in enumerate(tecr_hbbTemplate[0])])
    tecr_hbbBinYields = tecr_hbbBinYields * deepak15_weight['1tag']['Hbb']
    tecr_hbbObservable = rl.Observable('recoil', tecr_hbbHist.axis('recoil').edges())
    tecr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, tecr_hbbObservable, tecr_hbbBinYields)

    tecr.addSample(tecr_hbb)

    tecr_qcdHist = background['tecr'].integrate('process', 'QCD').integrate('systematic','nominal')
    tecr_qcdTemplate = template(tecr_qcdHist, 'recoil')
    tecr_qcdMC =  rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, tecr_qcdTemplate)
    tecr_qcdMC.setParamEffect(lumi, 1.027)
    tecr_qcdMC.setParamEffect(trig_e, 1.01)
    tecr_qcdMC.setParamEffect(veto_tau, 1.03)
    tecr_qcdMC.setParamEffect(qcdmu_norm, 2.0)
    tecr_qcdMC.setParamEffect(jec, 1.05)
    tecr_qcdMC.setParamEffect(id_e, 1.02)
    tecr_qcdMC.setParamEffect(reco_e, 1.02)
    btagUp=background['tecr'].integrate('process', 'QCD').integrate('systematic','btagUp')
    btagDown=background['tecr'].integrate('process', 'QCD').integrate('systematic','btagDown')
    tecr_qcdMC.setParamEffect(btag, btagUp, btagDown)

    tecr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,tecr_qcdTemplate[0].max()*2) for i,b in enumerate(tecr_qcdTemplate[0])])
    tecr_qcdBinYields = tecr_qcdBinYields * deepak15_weight['1tag']['QCD']
    tecr_qcdObservable = rl.Observable('recoil', tecr_qcdHist.axis('recoil').edges())
    tecr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, tecr_qcdObservable, tecr_qcdBinYields)

    tecr.addSample(tecr_qcd)

    ###
    # End of single electron top control region
    ###

    ###
    ###
    # Double muon control region
    ###
    ###

    ch_name = 'zmcr-'+mass+'-'+category
    zmcr = rl.Channel(ch_name)
    model.addChannel(zmcr)

    ###
    # Add data distribution to the channel
    ###

    zmcr.setObservation(template(data['zmcr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))

    zmcr_dyjetsHist = background['zmcr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    zmcr_dyjetsTemplate = template(zmcr_dyjetsHist, 'recoil')
    zmcr_dyjetsMC =  rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, zmcr_dyjetsTemplate)
    zmcr_dyjetsMC.setParamEffect(lumi, 1.027)
    zmcr_dyjetsMC.setParamEffect(trig_met, 1.01)
    zmcr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    zmcr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    zmcr_dyjetsMC.setParamEffect(jec, 1.05)
    zmcr_dyjetsMC.setParamEffect(id_mu, 1.02)
    zmcr_dyjetsMC.setParamEffect(iso_mu, 1.02)

    zmcr_dyjetsTransferFactor = zmcr_dyjetsMC.getExpectation() / sr_zjetsMC.getExpectation() * hf_fraction_weight['notag']['Z+jets']
    zmcr_dyjets = rl.TransferFactorSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, zmcr_dyjetsTransferFactor, sr_zjets)

    zmcr.addSample(zmcr_dyjets)

    ###
    # Other MC-driven processes
    ###

    zmcr_ttHist = background['zmcr'].integrate('process', 'TT').integrate('systematic','nominal')
    zmcr_ttTemplate = template(zmcr_ttHist, 'recoil')
    zmcr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, zmcr_ttTemplate)
    zmcr_ttMC.setParamEffect(lumi, 1.027)
    zmcr_ttMC.setParamEffect(trig_met, 1.01)
    zmcr_ttMC.setParamEffect(veto_tau, 1.03)
    zmcr_ttMC.setParamEffect(tt_norm, 1.4)
    zmcr_ttMC.setParamEffect(jec, 1.05)
    zmcr_ttMC.setParamEffect(id_mu, 1.02)
    zmcr_ttMC.setParamEffect(iso_mu, 1.02)

    zmcr_ttBinYields = np.array([rl.IndependentParameter(ch_name+'_tt_bin_%d' % i,b,0,zmcr_ttTemplate[0].max()*2) for i,b in enumerate(zmcr_ttTemplate[0])])
    zmcr_ttBinYields = zmcr_ttBinYields * deepak15_weight['notag']['TT']
    zmcr_ttObservable = rl.Observable('recoil', zmcr_ttHist.axis('recoil').edges())
    zmcr_tt = rl.ParametericSample(ch_name+'_tt', rl.Sample.BACKGROUND, zmcr_ttObservable, zmcr_ttBinYields)

    zmcr.addSample(zmcr_tt)

    zmcr_stHist = background['zmcr'].integrate('process', 'ST').integrate('systematic','nominal')
    zmcr_stTemplate = template(zmcr_stHist, 'recoil')
    zmcr_stMC =  rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, zmcr_stTemplate)
    zmcr_stMC.setParamEffect(lumi, 1.027)
    zmcr_stMC.setParamEffect(trig_met, 1.01)
    zmcr_stMC.setParamEffect(veto_tau, 1.03)
    zmcr_stMC.setParamEffect(st_norm, 1.2)
    zmcr_stMC.setParamEffect(jec, 1.05)
    zmcr_stMC.setParamEffect(id_mu, 1.02)
    zmcr_stMC.setParamEffect(iso_mu, 1.02)

    zmcr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,zmcr_stTemplate[0].max()*2) for i,b in enumerate(zmcr_stTemplate[0])])
    zmcr_stBinYields = zmcr_stBinYields * deepak15_weight['notag']['ST']
    zmcr_stObservable = rl.Observable('recoil', zmcr_stHist.axis('recoil').edges())
    zmcr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, zmcr_stObservable, zmcr_stBinYields)

    zmcr.addSample(zmcr_st)

    zmcr_vvHist = background['zmcr'].integrate('process', 'VV').integrate('systematic','nominal')
    zmcr_vvTemplate = template(zmcr_vvHist, 'recoil')
    zmcr_vvMC =  rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, zmcr_vvTemplate)
    zmcr_vvMC.setParamEffect(lumi, 1.027)
    zmcr_vvMC.setParamEffect(trig_met, 1.01)
    zmcr_vvMC.setParamEffect(veto_tau, 1.03)
    zmcr_vvMC.setParamEffect(vv_norm, 1.2)
    zmcr_vvMC.setParamEffect(jec, 1.05)
    zmcr_vvMC.setParamEffect(id_mu, 1.02)
    zmcr_vvMC.setParamEffect(iso_mu, 1.02)

    zmcr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,zmcr_vvTemplate[0].max()*2) for i,b in enumerate(zmcr_vvTemplate[0])])
    zmcr_vvBinYields = zmcr_vvBinYields * deepak15_weight['notag']['VV']
    zmcr_vvObservable = rl.Observable('recoil', zmcr_vvHist.axis('recoil').edges())
    zmcr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, zmcr_vvObservable, zmcr_vvBinYields)

    zmcr.addSample(zmcr_vv)

    zmcr_hbbHist = background['zmcr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    zmcr_hbbTemplate = template(zmcr_hbbHist, 'recoil')
    zmcr_hbbMC =  rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, zmcr_hbbTemplate)
    zmcr_hbbMC.setParamEffect(lumi, 1.027)
    zmcr_hbbMC.setParamEffect(trig_met, 1.01)
    zmcr_hbbMC.setParamEffect(veto_tau, 1.03)
    zmcr_hbbMC.setParamEffect(hbb_norm, 1.2)
    zmcr_hbbMC.setParamEffect(jec, 1.05)
    zmcr_hbbMC.setParamEffect(id_mu, 1.02)
    zmcr_hbbMC.setParamEffect(iso_mu, 1.02)

    zmcr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,zmcr_hbbTemplate[0].max()*2) for i,b in enumerate(zmcr_hbbTemplate[0])])
    zmcr_hbbBinYields = zmcr_hbbBinYields * deepak15_weight['notag']['Hbb']
    zmcr_hbbObservable = rl.Observable('recoil', zmcr_hbbHist.axis('recoil').edges())
    zmcr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, zmcr_hbbObservable, zmcr_hbbBinYields)

    zmcr.addSample(zmcr_hbb)

    ###
    # End of double muon control region
    ###

    ###
    ###
    # Double electron control region
    ###
    ###

    ch_name = 'zecr-'+mass+'-'+category
    zecr = rl.Channel(ch_name)
    model.addChannel(zecr)

    ###
    # Add data distribution to the channel
    ###

    zecr.setObservation(template(data['zecr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))

    zecr_dyjetsHist = background['zecr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    zecr_dyjetsTemplate = template(zecr_dyjetsHist, 'recoil')
    zecr_dyjetsMC =  rl.TemplateSample(ch_name+'_dyjetsMC', rl.Sample.BACKGROUND, zecr_dyjetsTemplate)
    zecr_dyjetsMC.setParamEffect(lumi, 1.027)
    zecr_dyjetsMC.setParamEffect(trig_e, 1.01)
    zecr_dyjetsMC.setParamEffect(veto_tau, 1.03)
    zecr_dyjetsMC.setParamEffect(dyjets_norm, 1.4)
    zecr_dyjetsMC.setParamEffect(jec, 1.05)
    zecr_dyjetsMC.setParamEffect(id_e, 1.02)
    zecr_dyjetsMC.setParamEffect(reco_e, 1.02)

    zecr_dyjetsTransferFactor = zecr_dyjetsMC.getExpectation() / sr_zjetsMC.getExpectation() * hf_fraction_weight['notag']['Z+jets']
    zecr_dyjets = rl.TransferFactorSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, zecr_dyjetsTransferFactor, sr_zjets)

    zecr.addSample(zecr_dyjets)

    ###
    # Other MC-driven processes
    ###

    zecr_ttHist = background['zecr'].integrate('process', 'TT').integrate('systematic','nominal')
    zecr_ttTemplate = template(zecr_ttHist, 'recoil')
    zecr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, zecr_ttTemplate)
    zecr_ttMC.setParamEffect(lumi, 1.027)
    zecr_ttMC.setParamEffect(trig_e, 1.01)
    zecr_ttMC.setParamEffect(veto_tau, 1.03)
    zecr_ttMC.setParamEffect(tt_norm, 1.4)
    zecr_ttMC.setParamEffect(jec, 1.05)
    zecr_ttMC.setParamEffect(id_e, 1.02)
    zecr_ttMC.setParamEffect(reco_e, 1.02)

    zecr_ttBinYields = np.array([rl.IndependentParameter(ch_name+'_tt_bin_%d' % i,b,0,zecr_ttTemplate[0].max()*2) for i,b in enumerate(zecr_ttTemplate[0])])
    zecr_ttBinYields = zecr_ttBinYields * deepak15_weight['notag']['TT']
    zecr_ttObservable = rl.Observable('recoil', zecr_ttHist.axis('recoil').edges())
    zecr_tt = rl.ParametericSample(ch_name+'_tt', rl.Sample.BACKGROUND, zecr_ttObservable, zecr_ttBinYields)

    zecr.addSample(zecr_tt)

    zecr_stHist = background['zecr'].integrate('process', 'ST').integrate('systematic','nominal')
    zecr_stTemplate = template(zecr_stHist, 'recoil')
    zecr_stMC =  rl.TemplateSample(ch_name+'_stMC', rl.Sample.BACKGROUND, zecr_stTemplate)
    zecr_stMC.setParamEffect(lumi, 1.027)
    zecr_stMC.setParamEffect(trig_e, 1.01)
    zecr_stMC.setParamEffect(veto_tau, 1.03)
    zecr_stMC.setParamEffect(st_norm, 1.2)
    zecr_stMC.setParamEffect(jec, 1.05)
    zecr_stMC.setParamEffect(id_e, 1.02)
    zecr_stMC.setParamEffect(reco_e, 1.02)

    zecr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,zecr_stTemplate[0].max()*2) for i,b in enumerate(zecr_stTemplate[0])])
    zecr_stBinYields = zecr_stBinYields * deepak15_weight['notag']['ST']
    zecr_stObservable = rl.Observable('recoil', zecr_stHist.axis('recoil').edges())
    zecr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, zecr_stObservable, zecr_stBinYields)

    zecr.addSample(zecr_st)

    zecr_vvHist = background['zecr'].integrate('process', 'VV').integrate('systematic','nominal')
    zecr_vvTemplate = template(zecr_vvHist, 'recoil')
    zecr_vvMC =  rl.TemplateSample(ch_name+'_vvMC', rl.Sample.BACKGROUND, zecr_vvTemplate)
    zecr_vvMC.setParamEffect(lumi, 1.027)
    zecr_vvMC.setParamEffect(trig_e, 1.01)
    zecr_vvMC.setParamEffect(veto_tau, 1.03)
    zecr_vvMC.setParamEffect(vv_norm, 1.2)
    zecr_vvMC.setParamEffect(jec, 1.05)
    zecr_vvMC.setParamEffect(id_e, 1.02)
    zecr_vvMC.setParamEffect(reco_e, 1.02)

    zecr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,zecr_vvTemplate[0].max()*2) for i,b in enumerate(zecr_vvTemplate[0])])
    zecr_vvBinYields = zecr_vvBinYields * deepak15_weight['notag']['VV']
    zecr_vvObservable = rl.Observable('recoil', zecr_vvHist.axis('recoil').edges())
    zecr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, zecr_vvObservable, zecr_vvBinYields)

    zecr.addSample(zecr_vv)

    zecr_hbbHist = background['zecr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    zecr_hbbTemplate = template(zecr_hbbHist, 'recoil')
    zecr_hbbMC =  rl.TemplateSample(ch_name+'_hbbMC', rl.Sample.BACKGROUND, zecr_hbbTemplate)
    zecr_hbbMC.setParamEffect(lumi, 1.027)
    zecr_hbbMC.setParamEffect(trig_e, 1.01)
    zecr_hbbMC.setParamEffect(veto_tau, 1.03)
    zecr_hbbMC.setParamEffect(hbb_norm, 1.2)
    zecr_hbbMC.setParamEffect(jec, 1.05)
    zecr_hbbMC.setParamEffect(id_e, 1.02)
    zecr_hbbMC.setParamEffect(reco_e, 1.02)

    zecr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,zecr_hbbTemplate[0].max()*2) for i,b in enumerate(zecr_hbbTemplate[0])])
    zecr_hbbBinYields = zecr_hbbBinYields * deepak15_weight['notag']['Hbb']
    zecr_hbbObservable = rl.Observable('recoil', zecr_hbbHist.axis('recoil').edges())
    zecr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, zecr_hbbObservable, zecr_hbbBinYields)

    zecr.addSample(zecr_hbb)

    ###
    # End of double electron control region
    ###

    ###
    ###
    # Single photon control region
    ###
    ###

    ch_name = 'gcr-'+mass+'-'+category
    gcr = rl.Channel(ch_name)
    model.addChannel(gcr)

    ###
    # Add data distribution to the channel
    ###

    gcr.setObservation(template(data['gcr'].integrate('process', 'SinglePhoton').integrate('systematic','nominal'), 'recoil'))

    gcr_gjetsHist = background['gcr'].integrate('process', 'G+jets').integrate('systematic','nominal')
    gcr_gjetsTemplate = template(gcr_gjetsHist, 'recoil')
    gcr_gjetsMC =  rl.TemplateSample(ch_name+'_gjetsMC', rl.Sample.BACKGROUND, gcr_gjetsTemplate)
    gcr_gjetsMC.setParamEffect(lumi, 1.027)
    gcr_gjetsMC.setParamEffect(trig_pho, 1.01)
    gcr_gjetsMC.setParamEffect(veto_tau, 1.03)
    gcr_gjetsMC.setParamEffect(gjets_norm, 1.4)
    gcr_gjetsMC.setParamEffect(jec, 1.05)
    gcr_gjetsMC.setParamEffect(id_pho, 1.02)

    gcr_gjetsTransferFactor = gcr_gjetsMC.getExpectation() / sr_zjetsMC.getExpectation() * hf_fraction_weight['notag']['G+jets']
    gcr_gjets = rl.TransferFactorSample(ch_name+'_gjets', rl.Sample.BACKGROUND, gcr_gjetsTransferFactor, sr_zjets)

    gcr.addSample(gcr_gjets)

    gcr_qcdHist = background['gcr'].integrate('process', 'QCD').integrate('systematic','nominal')
    gcr_qcdTemplate = template(gcr_qcdHist, 'recoil')
    gcr_qcdMC =  rl.TemplateSample(ch_name+'_qcdMC', rl.Sample.BACKGROUND, gcr_qcdTemplate)
    gcr_qcdMC.setParamEffect(lumi, 1.027)
    gcr_qcdMC.setParamEffect(trig_pho, 1.01)
    gcr_qcdMC.setParamEffect(veto_tau, 1.03)
    gcr_qcdMC.setParamEffect(qcdpho_norm, 2.0)
    gcr_qcdMC.setParamEffect(jec, 1.05)
    gcr_gjetsMC.setParamEffect(id_pho, 1.02)

    gcr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,gcr_qcdTemplate[0].max()*2) for i,b in enumerate(gcr_qcdTemplate[0])])
    gcr_qcdBinYields = gcr_qcdBinYields * deepak15_weight['notag']['QCD']
    gcr_qcdObservable = rl.Observable('recoil', gcr_qcdHist.axis('recoil').edges())
    gcr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, gcr_qcdObservable, gcr_qcdBinYields)

    gcr.addSample(gcr_qcd)

    with open(os.path.join(str(tmpdir), 'darkhiggsModel'+year+'.pkl'), "wb") as fout:
        pickle.dump(model, fout)

    print('Rendering')
    model.renderCombine(os.path.join(str(tmpdir), 'darkhiggsModel'+year+'/'+mass))


if __name__ == '__main__':
    if not os.path.exists('datacards'):
        os.mkdir('datacards')
    mass='mass0'
    category='monojet'
    year='2018'
    grouping=False
    darkhiggs_model('datacards',mass,category,year,grouping)
