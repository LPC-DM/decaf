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

def model(category,year,mass,grouping):
    
    def template(dictionary, process, systematic, region):
        print('Generating template for',process,'in',region)
        #print(dictionary[region].integrate('process', process).integrate('systematic',systematic).values()[()])
        output=dictionary[region].integrate('process', process).integrate('systematic',systematic).values()[()]
        binning=dictionary[region].integrate('process', process).integrate('systematic',systematic).axis('recoil').edges()
        return (output, np.arange(output.size+1), 'recoil')

    model_id='-'.join([year, category])
    if mass is not None: model_id='-'.join([year, category, mass])
    model = rl.Model('darkhiggs'+model_id)

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
        'xbb'    : rl.IndependentParameter('xbb_sf_'+year, 1., 0, 1/deepak15_pass_eff['xbb']), 
        'vqq'    : rl.IndependentParameter('vqq_sf_'+year, 1., 0, 1/deepak15_pass_eff['vqq']),
        'wcq'    : rl.IndependentParameter('wcq_sf_'+year, 1., 0, 1/deepak15_pass_eff['wcq']),
        'b'      : rl.IndependentParameter('b_sf_'+year, 1., 0, 1/deepak15_pass_eff['b']),
        'bb'     : rl.IndependentParameter('bb_sf_'+year, 1., 0, 1/deepak15_pass_eff['bb']), 
        'bc'     : rl.IndependentParameter('bc_sf_'+year, 1., 0, 1/deepak15_pass_eff['bc']),
        'c'      : rl.IndependentParameter('c_sf_'+year, 1., 0, 1/deepak15_pass_eff['c']),
        'cc'     : rl.IndependentParameter('cc_sf_'+year, 1., 0, 1/deepak15_pass_eff['cc']),
        'garbage': 1.,
        'other'  : rl.IndependentParameter('other_sf_'+year, 1., 0, 1/deepak15_pass_eff['other']),
        'tbcq'   : rl.IndependentParameter('tbcq_sf_'+year, 1., 0, 1/deepak15_pass_eff['tbcq']),
        'tbqq'   : rl.IndependentParameter('tbqq_sf_'+year, 1., 0, 1/deepak15_pass_eff['tbqq']),
        'zcc'    : rl.IndependentParameter('zcc_sf_'+year, 1., 0, 1/deepak15_pass_eff['zcc'])
    }
    
    deepak4_0tag_gentype_eff = {
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

    gentypes = {
        'Hbb': ['xbb','vqq','wcq','b','bb','bc','garbage','other','tbcq','tbqq'],
        'Z+HF': ['b','bb','c','cc','garbage','other'],
        'Z+LF': ['garbage','other'],
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
        
    with open('data/signal_fractions.json') as fin:
        signal_fractions = json.load(fin)

    with open('data/signal_fractions_no_mass.json') as fin:
        signal_fractions_no_mass = json.load(fin)

    signal_weight={}
    for process in signal_fractions.keys():
        #signal_weight[process]=np.array([])
        for gentype in signal_fractions[process].keys():
            if gentype not in gentypes[process.split('_')[0]]: continue
            print('Extracting',gentype,'fraction for',process)
            if mass is not None:
                try:
                    weight
                except:
                    weight = deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*np.array(signal_fractions[process][gentype][mass])
                else:
                    weight += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*np.array(signal_fractions[process][gentype][mass])
            else:
                try:
                    weight
                except:
                    weight = (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*np.array(signal_fractions_no_mass[process][gentype])
                else:
                    weight += (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*np.array(signal_fractions_no_mass[process][gentype])

        signal_weight[process] = weight #np.concatenate((signal_weight[process], weight))
                
    with open('data/fractions.json') as fin:
        fractions = json.load(fin)

    with open('data/fractions_no_mass.json') as fin:
        fractions_no_mass = json.load(fin)

    deepak15_weight={}
    deepak15_weight['0tag']={}
    deepak15_weight['1tag']={}
    deepak15_weight['notag']={}
    for process in ['Hbb','VV','ST','QCD','TT','Z+HF','Z+LF','W+HF','W+LF','G+HF','G+LF']:
        #deepak15_weight['0tag'][process]=np.array([])
        #deepak15_weight['1tag'][process]=np.array([])
        #deepak15_weight['notag'][process]=np.array([])
        for gentype in fractions[process].keys():
            if gentype not in gentypes[process]: continue
            print('Extracting',gentype,'fraction for',process )
            if mass is not None:
                try:
                    weight_0tag
                except:
                    weight_0tag = deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*deepak4_0tag_gentype_eff[process][gentype]*np.array(fractions[process][gentype][mass])
                else:
                    weight_0tag += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*deepak4_0tag_gentype_eff[process][gentype]*np.array(fractions[process][gentype][mass])
                try:
                    weight_1tag
                except:
                    weight_1tag = deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*(1 - deepak4_0tag_gentype_eff[process][gentype])*np.array(fractions[process][gentype][mass])
                else:
                    weight_1tag += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*(1 - deepak4_0tag_gentype_eff[process][gentype])*np.array(fractions[process][gentype][mass])
                try:
                    weight_notag
                except:
                    weight_notag = deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*np.array(fractions[process][gentype][mass])
                else:
                    weight_notag += deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype]*np.array(fractions[process][gentype][mass])
            else:
                try:
                    weight_0tag
                except:
                    weight_0tag = (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*deepak4_0tag_gentype_eff[process][gentype]*np.array(fractions_no_mass[process][gentype])
                else:
                    weight_0tag += (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*deepak4_0tag_gentype_eff[process][gentype]*np.array(fractions_no_mass[process][gentype])
                try:
                    weight_1tag
                except:
                    weight_1tag = (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*(1 - deepak4_0tag_gentype_eff[process][gentype])*np.array(fractions_no_mass[process][gentype])
                else:
                    weight_1tag += (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*(1 - deepak4_0tag_gentype_eff[process][gentype])*np.array(fractions_no_mass[process][gentype])
                try:
                    weight_notag
                except:
                    weight_notag = (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*np.array(fractions_no_mass[process][gentype])
                else:
                    weight_notag += (1 - deepak15_pass_sf[gentype]*deepak15_pass_eff[gentype])*np.array(fractions_no_mass[process][gentype])

        deepak15_weight['0tag'][process]=np.nan_to_num(weight_0tag/deepak4_0tag_process_eff[process])
        deepak15_weight['1tag'][process]=np.nan_to_num(weight_1tag/(1 - deepak4_0tag_process_eff[process]))
        deepak15_weight['notag'][process]=weight_notag

    hf_fraction_weight={}
    hf_fraction_weight['0tag']={}
    hf_fraction_weight['1tag']={}
    hf_fraction_weight['notag']={}

    hf_fraction_weight['0tag']['W+jets'] = np.nan_to_num(deepak15_weight['0tag']['W+HF']*(deepak4_0tag_process_eff['W+HF']/deepak4_0tag_process_eff['W+jets'])*whf_k*whf_fraction)
    hf_fraction_weight['0tag']['W+jets'] += np.nan_to_num(deepak15_weight['0tag']['W+LF']*(deepak4_0tag_process_eff['W+LF']/deepak4_0tag_process_eff['W+jets'])*(1 - whf_k*whf_fraction))
    hf_fraction_weight['1tag']['W+jets'] = np.nan_to_num(deepak15_weight['1tag']['W+HF']*((1-deepak4_0tag_process_eff['W+HF'])/(1-deepak4_0tag_process_eff['W+jets']))*whf_k*whf_fraction) 
    hf_fraction_weight['1tag']['W+jets'] += np.nan_to_num(deepak15_weight['1tag']['W+LF']*((1-deepak4_0tag_process_eff['W+HF'])/(1-deepak4_0tag_process_eff['W+jets']))*(1 - whf_k*whf_fraction))
    hf_fraction_weight['notag']['W+jets'] = deepak15_weight['notag']['W+HF']*whf_k*whf_fraction
    hf_fraction_weight['notag']['W+jets'] += deepak15_weight['notag']['W+LF']*(1 - whf_k*whf_fraction)

    hf_fraction_weight['0tag']['Z+jets'] = np.nan_to_num(deepak15_weight['0tag']['Z+HF']*(deepak4_0tag_process_eff['Z+HF']/deepak4_0tag_process_eff['Z+jets'])*zhf_k*zhf_fraction)
    hf_fraction_weight['0tag']['Z+jets'] += np.nan_to_num(deepak15_weight['0tag']['Z+LF']*(deepak4_0tag_process_eff['Z+LF']/deepak4_0tag_process_eff['Z+jets'])*(1 - zhf_k*zhf_fraction))
    hf_fraction_weight['1tag']['Z+jets'] = np.nan_to_num(deepak15_weight['1tag']['Z+HF']*((1-deepak4_0tag_process_eff['Z+HF'])/(1-deepak4_0tag_process_eff['Z+jets']))*zhf_k*zhf_fraction) 
    hf_fraction_weight['1tag']['Z+jets'] += np.nan_to_num(deepak15_weight['1tag']['Z+LF']*((1-deepak4_0tag_process_eff['Z+HF'])/(1-deepak4_0tag_process_eff['Z+jets']))*(1 - zhf_k*zhf_fraction))
    hf_fraction_weight['notag']['Z+jets'] = deepak15_weight['notag']['Z+HF']*zhf_k*zhf_fraction
    hf_fraction_weight['notag']['Z+jets'] += deepak15_weight['notag']['Z+LF']*(1 - zhf_k*zhf_fraction)

    hf_fraction_weight['0tag']['G+jets'] = np.nan_to_num(deepak15_weight['0tag']['G+HF']*(deepak4_0tag_process_eff['G+HF']/deepak4_0tag_process_eff['G+jets'])*ghf_k*ghf_fraction) 
    hf_fraction_weight['0tag']['G+jets'] += np.nan_to_num(deepak15_weight['0tag']['G+LF']*(deepak4_0tag_process_eff['G+LF']/deepak4_0tag_process_eff['G+jets'])*(1 - ghf_k*ghf_fraction))
    hf_fraction_weight['1tag']['G+jets'] = np.nan_to_num(deepak15_weight['1tag']['G+HF']*((1-deepak4_0tag_process_eff['G+HF'])/(1-deepak4_0tag_process_eff['G+jets']))*ghf_k*ghf_fraction) 
    hf_fraction_weight['1tag']['G+jets'] += np.nan_to_num(deepak15_weight['1tag']['G+LF']*((1-deepak4_0tag_process_eff['G+HF'])/(1-deepak4_0tag_process_eff['G+jets']))*(1 - ghf_k*ghf_fraction))
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
        'mass4': data_hists['recoil'].integrate('gentype','data').integrate('systematic','nominal').integrate('region','sr_mass4').integrate('process','MET').axis('recoil').edges(),
        'nomass': data_hists['recoil'].integrate('gentype','data').integrate('systematic','nominal').integrate('region','sr').integrate('process','MET').axis('recoil').edges()
    }    

    ###
    # Preparing histograms for fit
    ##

    data = {}
    for r in data_hists['recoil'].identifiers('region'):
        if category not in str(r): continue
        if mass is None and 'mass' in str(r): continue
        if mass is not None and mass not in str(r): continue
        m=mass
        if mass is None: m='nomass'
        print('data',r,category,m)
        data[str(r).split("_")[0]]=data_hists['recoil'].integrate('region',r).integrate('gentype').rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning[m]))

    background = {}
    for r in data.keys():
        m=mass
        if mass is None: m='nomass'
        print('bkg',r,category,m)
        background[r]=bkg_hists['recoil'].integrate('region',str(r)).integrate('gentype').rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning[m]))

    signal = {}
    for r in data.keys():
        m=mass
        if mass is None: m='nomass'
        print('sig',r,category,m)
        signal[r]=signal_hists['recoil'].integrate('region',str(r)).integrate('gentype').rebin('recoil',hist.Bin('recoil','Hadronic recoil',binning[m]))

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

    gamma_to_z_ewk = rl.NuisanceParameter('Theory_gamma_z_ewk', 'shape')

    ###
    ###
    # Signal region
    ###
    ###

    ch_name = 'sr-'+model_id
    sr = rl.Channel(ch_name)
    model.addChannel(sr)

    ###
    # Add data distribution to the channel
    ###

    #sr.setObservation(template(data['sr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))
    sr.setObservation(template(data,'MET','nominal','sr'))

    ###
    # Z(->nunu)+jets data-driven model
    ###

    #sr_zjetsHist = background['sr'].integrate('process', 'Z+jets').integrate('systematic','nominal')
    sr_zjetsTemplate = template(background,'Z+jets','nominal','sr')
    sr_zjetsMC =  rl.TemplateSample(ch_name+'_zjetsMC', rl.Sample.BACKGROUND, sr_zjetsTemplate)
    sr_zjetsMC.setParamEffect(lumi, 1.027)
    sr_zjetsMC.setParamEffect(trig_met, 1.01)
    sr_zjetsMC.setParamEffect(veto_tau, 1.03)
    sr_zjetsMC.setParamEffect(zjets_norm, 1.4)
    sr_zjetsMC.setParamEffect(jec, 1.05)
    btagUp=template(background,'Z+jets','btagUp','sr')[0]
    btagDown=template(background,'Z+jets','btagDown','sr')[0]
    sr_zjetsMC.setParamEffect(btag, btagUp, btagDown)
    sr_zjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_zjets_bin_%d' % i, b, 0, sr_zjetsTemplate[0].max()*2) for i,b in enumerate(sr_zjetsTemplate[0])]) 
    sr_zjetsBinYields = sr_zjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    sr_zjetsObservable = rl.Observable('recoil', sr_zjetsTemplate[1])
    sr_zjets = rl.ParametericSample(ch_name+'_zjets', rl.Sample.BACKGROUND, sr_zjetsObservable, sr_zjetsBinYields)

    sr.addSample(sr_zjets)

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    #sr_wjetsHist = background['sr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    sr_wjetsTemplate = template(background,'W+jets','nominal','sr')
    sr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, sr_wjetsTemplate)
    sr_wjetsMC.setParamEffect(lumi, 1.027)
    sr_wjetsMC.setParamEffect(trig_met, 1.01)
    sr_wjetsMC.setParamEffect(veto_tau, 1.03)
    sr_wjetsMC.setParamEffect(wjets_norm, 1.4)
    sr_wjetsMC.setParamEffect(jec, 1.05)
    btagUp=template(background,'W+jets','btagUp','sr')[0]
    btagDown=template(background,'W+jets','btagDown','sr')[0]
    sr_wjetsMC.setParamEffect(btag, btagUp, btagDown)
    sr_wjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_wjets_bin_%d' % i,b,0,sr_wjetsTemplate[0].max()*2) for i,b in enumerate(sr_wjetsTemplate[0])]) 
    sr_wjetsBinYields = sr_wjetsBinYields * hf_fraction_weight['0tag']['W+jets']
    sr_wjetsObservable = rl.Observable('recoil', sr_wjetsTemplate[1])
    sr_wjets = rl.ParametericSample(ch_name+'_wjets', rl.Sample.BACKGROUND, sr_wjetsObservable, sr_wjetsBinYields)

    sr.addSample(sr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    #sr_ttHist = background['sr'].integrate('process', 'TT').integrate('systematic','nominal')
    sr_ttTemplate = template(background,'TT','nominal','sr')
    sr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, sr_ttTemplate)
    sr_ttMC.setParamEffect(lumi, 1.027)
    sr_ttMC.setParamEffect(trig_met, 1.01)
    sr_ttMC.setParamEffect(veto_tau, 1.03)
    sr_ttMC.setParamEffect(tt_norm, 1.4)
    sr_ttMC.setParamEffect(jec, 1.05)
    btagUp=template(background,'TT','btagUp','sr')[0]
    btagDown=template(background,'TT','btagDown','sr')[0]
    sr_ttMC.setParamEffect(btag, btagUp, btagDown)
    sr_ttBinYields = np.array([rl.IndependentParameter(ch_name+'_tt_bin_%d' % i,b,0,sr_ttTemplate[0].max()*2) for i,b in enumerate(sr_ttTemplate[0])])
    sr_ttBinYields = sr_ttBinYields * deepak15_weight['0tag']['TT']
    sr_ttObservable = rl.Observable('recoil', sr_ttTemplate[1])
    sr_tt = rl.ParametericSample(ch_name+'_tt', rl.Sample.BACKGROUND, sr_ttObservable, sr_ttBinYields)

    sr.addSample(sr_tt)

    ###
    # Other MC-driven processes
    ###

    #sr_stHist = background['sr'].integrate('process', 'ST').integrate('systematic','nominal')
    sr_stTemplate = template(background,'ST','nominal','sr')
    sr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,sr_stTemplate[0].max()*2) for i,b in enumerate(sr_stTemplate[0])])
    sr_stBinYields = sr_stBinYields * deepak15_weight['0tag']['ST']
    sr_stObservable = rl.Observable('recoil', sr_stTemplate[1])
    sr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, sr_stObservable, sr_stBinYields)
    sr_st.setParamEffect(lumi, 1.027)
    sr_st.setParamEffect(trig_met, 1.01)
    sr_st.setParamEffect(veto_tau, 1.03)
    sr_st.setParamEffect(st_norm, 1.2)
    sr_st.setParamEffect(jec, 1.05)
    btagUp=template(background,'ST','btagUp','sr')[0]
    btagDown=template(background,'ST','btagDown','sr')[0]
    sr_st.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_st)

    #sr_dyjetsHist = background['sr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    sr_dyjetsTemplate = template(background,'DY+jets','nominal','sr')
    sr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,sr_dyjetsTemplate[0].max()*2) for i,b in enumerate(sr_dyjetsTemplate[0])])
    sr_dyjetsBinYields = sr_dyjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    sr_dyjetsObservable = rl.Observable('recoil', sr_dyjetsTemplate[1])
    sr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, sr_dyjetsObservable, sr_dyjetsBinYields)
    sr_dyjets.setParamEffect(lumi, 1.027)
    sr_dyjets.setParamEffect(trig_met, 1.01)
    sr_dyjets.setParamEffect(veto_tau, 1.03)
    sr_dyjets.setParamEffect(dyjets_norm, 1.4)
    sr_dyjets.setParamEffect(jec, 1.05)
    btagUp=template(background,'DY+jets','btagUp','sr')[0]
    btagDown=template(background,'DY+jets','btagDown','sr')[0]
    sr_dyjets.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_dyjets)

    #sr_vvHist = background['sr'].integrate('process', 'VV').integrate('systematic','nominal')
    sr_vvTemplate = template(background,'VV','nominal','sr')
    sr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,sr_vvTemplate[0].max()*2) for i,b in enumerate(sr_vvTemplate[0])])
    sr_vvBinYields = sr_vvBinYields * deepak15_weight['0tag']['VV']
    sr_vvObservable = rl.Observable('recoil', sr_vvTemplate[1])
    sr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, sr_vvObservable, sr_vvBinYields)
    sr_vv.setParamEffect(lumi, 1.027)
    sr_vv.setParamEffect(trig_met, 1.01)
    sr_vv.setParamEffect(veto_tau, 1.03)
    sr_vv.setParamEffect(vv_norm, 1.2)
    sr_vv.setParamEffect(jec, 1.05)
    btagUp=template(background,'VV','btagUp','sr')[0]
    btagDown=template(background,'VV','btagDown','sr')[0]
    sr_vv.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_vv)

    #sr_hbbHist = background['sr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    sr_hbbTemplate = template(background,'Hbb','nominal','sr')
    sr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,sr_hbbTemplate[0].max()*2) for i,b in enumerate(sr_hbbTemplate[0])])
    sr_hbbBinYields = sr_hbbBinYields * deepak15_weight['0tag']['Hbb']
    sr_hbbObservable = rl.Observable('recoil', sr_hbbTemplate[1])
    sr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, sr_hbbObservable, sr_hbbBinYields)
    sr_hbb.setParamEffect(lumi, 1.027)
    sr_hbb.setParamEffect(trig_met, 1.01)
    sr_hbb.setParamEffect(veto_tau, 1.03)
    sr_hbb.setParamEffect(hbb_norm, 1.2)
    sr_hbb.setParamEffect(jec, 1.05)
    btagUp=template(background,'Hbb','btagUp','sr')[0]
    btagDown=template(background,'Hbb','btagDown','sr')[0]
    sr_hbb.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_hbb)

    #sr_qcdHist = background['sr'].integrate('process', 'QCD').integrate('systematic','nominal')
    sr_qcdTemplate = template(background,'QCD','nominal','sr')
    sr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,sr_qcdTemplate[0].max()*2) for i,b in enumerate(sr_qcdTemplate[0])])
    sr_qcdBinYields = sr_qcdBinYields * deepak15_weight['0tag']['QCD']
    sr_qcdObservable = rl.Observable('recoil', sr_qcdTemplate[1])
    sr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, sr_qcdObservable, sr_qcdBinYields)
    sr_qcd.setParamEffect(lumi, 1.027)
    sr_qcd.setParamEffect(trig_met, 1.01)
    sr_qcd.setParamEffect(veto_tau, 1.03)
    sr_qcd.setParamEffect(qcdsig_norm, 2.0)
    sr_qcd.setParamEffect(jec, 1.05)
    btagUp=template(background,'QCD','btagUp','sr')[0]
    btagDown=template(background,'QCD','btagDown','sr')[0]
    sr_qcd.setParamEffect(btag, btagUp, btagDown)
    sr.addSample(sr_qcd)

    for s in signal['sr'].identifiers('process'):
        print(s)
        #sr_signalHist = signal['sr'].integrate('process', s).integrate('systematic','nominal')
        sr_signalTemplate = template(signal,s,'nominal','sr')
        sr_signalBinYields = np.array([rl.IndependentParameter(ch_name+'_'+str(s)+'_bin_%d' % i,b,0,sr_signalTemplate[0].max()*2) for i,b in enumerate(sr_signalTemplate[0])])
        sr_signalBinYields = sr_signalBinYields * signal_weight[str(s)]
        sr_signalObservable = rl.Observable('recoil', sr_signalTemplate[1])
        sr_signal = rl.ParametericSample(ch_name+'_'+str(s), rl.Sample.SIGNAL, sr_signalObservable, sr_signalBinYields)
        sr_signal.setParamEffect(lumi, 1.027)
        sr_signal.setParamEffect(trig_met, 1.01)
        sr_signal.setParamEffect(veto_tau, 1.03)
        sr_signal.setParamEffect(jec, 1.05)
        btagUp=template(signal, s,'btagUp','sr')[0]
        btagDown=template(signal, s,'btagDown','sr')[0]
        sr_signal.setParamEffect(btag, btagUp, btagDown)
        sr.addSample(sr_signal)


    ###
    # End of SR
    ###

    

    ###
    ###
    # Single muon W control region
    ###
    ###

    ch_name = 'wmcr-'+model_id
    wmcr = rl.Channel(ch_name)
    model.addChannel(wmcr)

    ###
    # Add data distribution to the channel
    ###

    #wmcr.setObservation(template(data['wmcr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))
    wmcr.setObservation(template(data,'MET','nominal','wmcr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    #wmcr_wjetsHist = background['wmcr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    wmcr_wjetsTemplate = template(background,'W+jets','nominal','wmcr')
    wmcr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, wmcr_wjetsTemplate)
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
    wmcr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, wmcr_wjetsTransferFactor, sr_wjets)
    wmcr.addSample(wmcr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    #wmcr_ttHist = background['wmcr'].integrate('process', 'TT').integrate('systematic','nominal')
    wmcr_ttTemplate = template(background,'TT','nominal','wmcr')
    wmcr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, wmcr_ttTemplate)
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
    wmcr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, wmcr_ttTransferFactor, sr_tt)
    wmcr.addSample(wmcr_tt)

    ###
    # Other MC-driven processes
    ###

    #wmcr_stHist = background['wmcr'].integrate('process', 'ST').integrate('systematic','nominal')
    wmcr_stTemplate = template(background,'ST','nominal','wmcr')
    wmcr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,wmcr_stTemplate[0].max()*2) for i,b in enumerate(wmcr_stTemplate[0])])
    wmcr_stBinYields = wmcr_stBinYields * deepak15_weight['0tag']['ST']
    wmcr_stObservable = rl.Observable('recoil', wmcr_stTemplate[1])
    wmcr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, wmcr_stObservable, wmcr_stBinYields)
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
    wmcr.addSample(wmcr_st)

    #wmcr_dyjetsHist = background['wmcr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    wmcr_dyjetsTemplate = template(background,'DY+jets','nominal','wmcr')
    wmcr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,wmcr_dyjetsTemplate[0].max()*2) for i,b in enumerate(wmcr_dyjetsTemplate[0])])
    wmcr_dyjetsBinYields = wmcr_dyjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    wmcr_dyjetsObservable = rl.Observable('recoil', wmcr_dyjetsTemplate[1])
    wmcr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, wmcr_dyjetsObservable, wmcr_dyjetsBinYields)
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
    wmcr.addSample(wmcr_dyjets)

    #wmcr_vvHist = background['wmcr'].integrate('process', 'VV').integrate('systematic','nominal')
    wmcr_vvTemplate = template(background,'VV','nominal','wmcr')
    wmcr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,wmcr_vvTemplate[0].max()*2) for i,b in enumerate(wmcr_vvTemplate[0])])
    wmcr_vvBinYields = wmcr_vvBinYields * deepak15_weight['0tag']['VV']
    wmcr_vvObservable = rl.Observable('recoil', wmcr_vvTemplate[1])
    wmcr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, wmcr_vvObservable, wmcr_vvBinYields)
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
    wmcr.addSample(wmcr_vv)

    #wmcr_hbbHist = background['wmcr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    wmcr_hbbTemplate = template(background,'Hbb','nominal','wmcr')
    wmcr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,wmcr_hbbTemplate[0].max()*2) for i,b in enumerate(wmcr_hbbTemplate[0])])
    wmcr_hbbBinYields = wmcr_hbbBinYields * deepak15_weight['0tag']['Hbb']
    wmcr_hbbObservable = rl.Observable('recoil', wmcr_hbbTemplate[1])
    wmcr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, wmcr_hbbObservable, wmcr_hbbBinYields)
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
    wmcr.addSample(wmcr_hbb)

    #wmcr_qcdHist = background['wmcr'].integrate('process', 'QCD').integrate('systematic','nominal')
    wmcr_qcdTemplate = template(background,'QCD','nominal','wmcr')
    wmcr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,wmcr_qcdTemplate[0].max()*2) for i,b in enumerate(wmcr_qcdTemplate[0])])
    wmcr_qcdBinYields = wmcr_qcdBinYields * deepak15_weight['0tag']['QCD']
    wmcr_qcdObservable = rl.Observable('recoil', wmcr_qcdTemplate[1])
    wmcr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, wmcr_qcdObservable, wmcr_qcdBinYields)
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
    wmcr.addSample(wmcr_qcd)

    ###
    # End of single muon W control region
    ###

    ###
    ###
    # Single muon top control region
    ###
    ###

    ch_name = 'tmcr-'+model_id
    tmcr = rl.Channel(ch_name)
    model.addChannel(tmcr)

    ###
    # Add data distribution to the channel
    ###

    #tmcr.setObservation(template(data['tmcr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))
    tmcr.setObservation(template(data,'MET','nominal','tmcr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    #tmcr_wjetsHist = background['tmcr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    tmcr_wjetsTemplate = template(background,'W+jets','nominal','tmcr')
    tmcr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, tmcr_wjetsTemplate)
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
    tmcr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, tmcr_wjetsTransferFactor, sr_wjets)
    tmcr.addSample(tmcr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    #tmcr_ttHist = background['tmcr'].integrate('process', 'TT').integrate('systematic','nominal')
    tmcr_ttTemplate = template(background,'TT','nominal','tmcr')
    tmcr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, tmcr_ttTemplate)
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
    tmcr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, tmcr_ttTransferFactor, sr_tt)
    tmcr.addSample(tmcr_tt)

    ###
    # Other MC-driven processes
    ###

    #tmcr_stHist = background['tmcr'].integrate('process', 'ST').integrate('systematic','nominal')
    tmcr_stTemplate = template(background,'ST','nominal','tmcr')
    tmcr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,tmcr_stTemplate[0].max()*2) for i,b in enumerate(tmcr_stTemplate[0])])
    tmcr_stBinYields = tmcr_stBinYields * deepak15_weight['1tag']['ST']
    tmcr_stObservable = rl.Observable('recoil', tmcr_stTemplate[1])
    tmcr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, tmcr_stObservable, tmcr_stBinYields)
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
    tmcr.addSample(tmcr_st)

    #tmcr_dyjetsHist = background['tmcr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    tmcr_dyjetsTemplate = template(background,'DY+jets','nominal','tmcr')
    tmcr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,tmcr_dyjetsTemplate[0].max()*2) for i,b in enumerate(tmcr_dyjetsTemplate[0])])
    tmcr_dyjetsBinYields = tmcr_dyjetsBinYields * hf_fraction_weight['1tag']['Z+jets']
    tmcr_dyjetsObservable = rl.Observable('recoil', tmcr_dyjetsTemplate[1])
    tmcr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, tmcr_dyjetsObservable, tmcr_dyjetsBinYields)
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
    tmcr.addSample(tmcr_dyjets)

    #tmcr_vvHist = background['tmcr'].integrate('process', 'VV').integrate('systematic','nominal')
    tmcr_vvTemplate = template(background,'VV','nominal','tmcr')
    tmcr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,tmcr_vvTemplate[0].max()*2) for i,b in enumerate(tmcr_vvTemplate[0])])
    tmcr_vvBinYields = tmcr_vvBinYields * deepak15_weight['1tag']['VV']
    tmcr_vvObservable = rl.Observable('recoil', tmcr_vvTemplate[1])
    tmcr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, tmcr_vvObservable, tmcr_vvBinYields)
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
    tmcr.addSample(tmcr_vv)

    #tmcr_hbbHist = background['tmcr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    tmcr_hbbTemplate = template(background,'Hbb','nominal','tmcr')
    tmcr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,tmcr_hbbTemplate[0].max()*2) for i,b in enumerate(tmcr_hbbTemplate[0])])
    tmcr_hbbBinYields = tmcr_hbbBinYields * deepak15_weight['1tag']['Hbb']
    tmcr_hbbObservable = rl.Observable('recoil', tmcr_hbbTemplate[1])
    tmcr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, tmcr_hbbObservable, tmcr_hbbBinYields)
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
    tmcr.addSample(tmcr_hbb)

    #tmcr_qcdHist = background['tmcr'].integrate('process', 'QCD').integrate('systematic','nominal')
    tmcr_qcdTemplate = template(background,'QCD','nominal','tmcr')
    tmcr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,tmcr_qcdTemplate[0].max()*2) for i,b in enumerate(tmcr_qcdTemplate[0])])
    tmcr_qcdBinYields = tmcr_qcdBinYields * deepak15_weight['1tag']['QCD']
    tmcr_qcdObservable = rl.Observable('recoil', tmcr_qcdTemplate[1])
    tmcr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, tmcr_qcdObservable, tmcr_qcdBinYields)
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
    tmcr.addSample(tmcr_qcd)

    ###
    # End of single muon top control region
    ###

    ###
    ###
    # Single electron W control region
    ###
    ###

    ch_name = 'wecr-'+model_id
    wecr = rl.Channel(ch_name)
    model.addChannel(wecr)

    ###
    # Add data distribution to the channel
    ###

    #wecr.setObservation(template(data['wecr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))
    wecr.setObservation(template(data,'SingleElectron','nominal','wecr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    #wecr_wjetsHist = background['wecr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    wecr_wjetsTemplate = template(background,'W+jets','nominal','wecr')
    wecr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, wecr_wjetsTemplate)
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
    wecr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, wecr_wjetsTransferFactor, sr_wjets)
    wecr.addSample(wecr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    #wecr_ttHist = background['wecr'].integrate('process', 'TT').integrate('systematic','nominal')
    wecr_ttTemplate = template(background,'TT','nominal','wecr')
    wecr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, wecr_ttTemplate)
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
    wecr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, wecr_ttTransferFactor, sr_tt)
    wecr.addSample(wecr_tt)

    ###
    # Other MC-driven processes
    ###

    #wecr_stHist = background['wecr'].integrate('process', 'ST').integrate('systematic','nominal')
    wecr_stTemplate = template(background,'ST','nominal','wecr')
    wecr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,wecr_stTemplate[0].max()*2) for i,b in enumerate(wecr_stTemplate[0])])
    wecr_stBinYields = wecr_stBinYields * deepak15_weight['0tag']['ST']
    wecr_stObservable = rl.Observable('recoil', wecr_stTemplate[1])
    wecr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, wecr_stObservable, wecr_stBinYields)
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
    wecr.addSample(wecr_st)

    #wecr_dyjetsHist = background['wecr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    wecr_dyjetsTemplate = template(background,'DY+jets','nominal','wecr')
    wecr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,wecr_dyjetsTemplate[0].max()*2) for i,b in enumerate(wecr_dyjetsTemplate[0])])
    wecr_dyjetsBinYields = wecr_dyjetsBinYields * hf_fraction_weight['0tag']['Z+jets']
    wecr_dyjetsObservable = rl.Observable('recoil',wecr_dyjetsTemplate[1])
    wecr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, wecr_dyjetsObservable, wecr_dyjetsBinYields)
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
    wecr.addSample(wecr_dyjets)

    #wecr_vvHist = background['wecr'].integrate('process', 'VV').integrate('systematic','nominal')
    wecr_vvTemplate = template(background,'VV','nominal','wecr')
    wecr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,wecr_vvTemplate[0].max()*2) for i,b in enumerate(wecr_vvTemplate[0])])
    wecr_vvBinYields = wecr_vvBinYields * deepak15_weight['0tag']['VV']
    wecr_vvObservable = rl.Observable('recoil', wecr_vvTemplate[1])
    wecr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, wecr_vvObservable, wecr_vvBinYields)
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
    wecr.addSample(wecr_vv)

    #wecr_hbbHist = background['wecr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    wecr_hbbTemplate = template(background,'Hbb','nominal','wecr')
    wecr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,wecr_hbbTemplate[0].max()*2) for i,b in enumerate(wecr_hbbTemplate[0])])
    wecr_hbbBinYields = wecr_hbbBinYields * deepak15_weight['0tag']['Hbb']
    wecr_hbbObservable = rl.Observable('recoil', wecr_hbbTemplate[1])
    wecr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, wecr_hbbObservable, wecr_hbbBinYields)
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
    wecr.addSample(wecr_hbb)

    #wecr_qcdHist = background['wecr'].integrate('process', 'QCD').integrate('systematic','nominal')
    wecr_qcdTemplate = template(background,'QCD','nominal','wecr')
    wecr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,wecr_qcdTemplate[0].max()*2) for i,b in enumerate(wecr_qcdTemplate[0])])
    wecr_qcdBinYields = wecr_qcdBinYields * deepak15_weight['0tag']['QCD']
    wecr_qcdObservable = rl.Observable('recoil', wecr_qcdTemplate[1])
    wecr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, wecr_qcdObservable, wecr_qcdBinYields)
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
    wecr.addSample(wecr_qcd)

    ###
    # End of single electron W control region
    ###

    ###
    ###
    # Single electron top control region
    ###
    ###

    ch_name = 'tecr-'+model_id
    tecr = rl.Channel(ch_name)
    model.addChannel(tecr)

    ###
    # Add data distribution to the channel
    ###

    #tecr.setObservation(template(data['tecr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))
    tecr.setObservation(template(data,'SingleElectron','nominal','tecr'))

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    #tecr_wjetsHist = background['tecr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    tecr_wjetsTemplate = template(background,'W+jets','nominal','tecr')
    tecr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, tecr_wjetsTemplate)
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
    tecr_wjets = rl.TransferFactorSample(ch_name+'_wjets', rl.Sample.BACKGROUND, tecr_wjetsTransferFactor, sr_wjets)
    tecr.addSample(tecr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    #tecr_ttHist = background['tecr'].integrate('process', 'TT').integrate('systematic','nominal')
    tecr_ttTemplate = template(background,'TT','nominal','tecr')
    tecr_ttMC =  rl.TemplateSample(ch_name+'_ttMC', rl.Sample.BACKGROUND, tecr_ttTemplate)
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
    tecr_tt = rl.TransferFactorSample(ch_name+'_tt', rl.Sample.BACKGROUND, tecr_ttTransferFactor, sr_tt)
    tecr.addSample(tecr_tt)

    ###
    # Other MC-driven processes
    ###

    #tecr_stHist = background['tecr'].integrate('process', 'ST').integrate('systematic','nominal')
    tecr_stTemplate = template(background,'ST','nominal','tecr')
    tecr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,tecr_stTemplate[0].max()*2) for i,b in enumerate(tecr_stTemplate[0])])
    tecr_stBinYields = tecr_stBinYields * deepak15_weight['1tag']['ST']
    tecr_stObservable = rl.Observable('recoil', tecr_stTemplate[1])
    tecr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, tecr_stObservable, tecr_stBinYields)
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
    tecr.addSample(tecr_st)

    #tecr_dyjetsHist = background['tecr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    tecr_dyjetsTemplate = template(background,'DY+jets','nominal','tecr')
    tecr_dyjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_dyjets_bin_%d' % i,b,0,tecr_dyjetsTemplate[0].max()*2) for i,b in enumerate(tecr_dyjetsTemplate[0])])
    tecr_dyjetsBinYields = tecr_dyjetsBinYields * hf_fraction_weight['1tag']['Z+jets']
    tecr_dyjetsObservable = rl.Observable('recoil', tecr_dyjetsTemplate[1])
    tecr_dyjets = rl.ParametericSample(ch_name+'_dyjets', rl.Sample.BACKGROUND, tecr_dyjetsObservable, tecr_dyjetsBinYields)
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
    tecr.addSample(tecr_dyjets)

    #tecr_vvHist = background['tecr'].integrate('process', 'VV').integrate('systematic','nominal')
    tecr_vvTemplate = template(background,'VV','nominal','tecr')
    tecr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,tecr_vvTemplate[0].max()*2) for i,b in enumerate(tecr_vvTemplate[0])])
    tecr_vvBinYields = tecr_vvBinYields * deepak15_weight['1tag']['VV']
    tecr_vvObservable = rl.Observable('recoil', tecr_vvTemplate[1])
    tecr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, tecr_vvObservable, tecr_vvBinYields)
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
    tecr.addSample(tecr_vv)

    #tecr_hbbHist = background['tecr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    tecr_hbbTemplate = template(background,'Hbb','nominal','tecr')
    tecr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,tecr_hbbTemplate[0].max()*2) for i,b in enumerate(tecr_hbbTemplate[0])])
    tecr_hbbBinYields = tecr_hbbBinYields * deepak15_weight['1tag']['Hbb']
    tecr_hbbObservable = rl.Observable('recoil', tecr_hbbTemplate[1])
    tecr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, tecr_hbbObservable, tecr_hbbBinYields)
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
    tecr.addSample(tecr_hbb)

    #tecr_qcdHist = background['tecr'].integrate('process', 'QCD').integrate('systematic','nominal')
    tecr_qcdTemplate = template(background,'QCD','nominal','tecr')
    tecr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,tecr_qcdTemplate[0].max()*2) for i,b in enumerate(tecr_qcdTemplate[0])])
    tecr_qcdBinYields = tecr_qcdBinYields * deepak15_weight['1tag']['QCD']
    tecr_qcdObservable = rl.Observable('recoil', tecr_qcdTemplate[1])
    tecr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, tecr_qcdObservable, tecr_qcdBinYields)
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
    tecr.addSample(tecr_qcd)

    ###
    # End of single electron top control region
    ###

    ###
    ###
    # Double muon control region
    ###
    ###

    ch_name = 'zmcr-'+model_id
    zmcr = rl.Channel(ch_name)
    model.addChannel(zmcr)

    ###
    # Add data distribution to the channel
    ###

    #zmcr.setObservation(template(data['zmcr'].integrate('process', 'MET').integrate('systematic','nominal'), 'recoil'))
    zmcr.setObservation(template(data,'MET','nominal','zmcr'))

    #zmcr_dyjetsHist = background['zmcr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    zmcr_dyjetsTemplate = template(background,'DY+jets','nominal','zmcr')
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

    #zmcr_ttHist = background['zmcr'].integrate('process', 'TT').integrate('systematic','nominal')
    zmcr_ttTemplate = template(background,'TT','nominal','zmcr')
    zmcr_ttBinYields = np.array([rl.IndependentParameter(ch_name+'_tt_bin_%d' % i,b,0,zmcr_ttTemplate[0].max()*2) for i,b in enumerate(zmcr_ttTemplate[0])])
    zmcr_ttBinYields = zmcr_ttBinYields * deepak15_weight['notag']['TT']
    zmcr_ttObservable = rl.Observable('recoil', zmcr_ttTemplate[1])
    zmcr_tt = rl.ParametericSample(ch_name+'_tt', rl.Sample.BACKGROUND, zmcr_ttObservable, zmcr_ttBinYields)
    zmcr_tt.setParamEffect(lumi, 1.027)
    zmcr_tt.setParamEffect(trig_met, 1.01)
    zmcr_tt.setParamEffect(veto_tau, 1.03)
    zmcr_tt.setParamEffect(tt_norm, 1.4)
    zmcr_tt.setParamEffect(jec, 1.05)
    zmcr_tt.setParamEffect(id_mu, 1.02)
    zmcr_tt.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_tt)

    #zmcr_stHist = background['zmcr'].integrate('process', 'ST').integrate('systematic','nominal')
    zmcr_stTemplate = template(background,'ST','nominal','zmcr')
    zmcr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,zmcr_stTemplate[0].max()*2) for i,b in enumerate(zmcr_stTemplate[0])])
    zmcr_stBinYields = zmcr_stBinYields * deepak15_weight['notag']['ST']
    zmcr_stObservable = rl.Observable('recoil', zmcr_stTemplate[1])
    zmcr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, zmcr_stObservable, zmcr_stBinYields)
    zmcr_st.setParamEffect(lumi, 1.027)
    zmcr_st.setParamEffect(trig_met, 1.01)
    zmcr_st.setParamEffect(veto_tau, 1.03)
    zmcr_st.setParamEffect(st_norm, 1.2)
    zmcr_st.setParamEffect(jec, 1.05)
    zmcr_st.setParamEffect(id_mu, 1.02)
    zmcr_st.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_st)

    #zmcr_vvHist = background['zmcr'].integrate('process', 'VV').integrate('systematic','nominal')
    zmcr_vvTemplate = template(background,'VV','nominal','zmcr')
    zmcr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,zmcr_vvTemplate[0].max()*2) for i,b in enumerate(zmcr_vvTemplate[0])])
    zmcr_vvBinYields = zmcr_vvBinYields * deepak15_weight['notag']['VV']
    zmcr_vvObservable = rl.Observable('recoil', zmcr_vvTemplate[1])
    zmcr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, zmcr_vvObservable, zmcr_vvBinYields)
    zmcr_vv.setParamEffect(lumi, 1.027)
    zmcr_vv.setParamEffect(trig_met, 1.01)
    zmcr_vv.setParamEffect(veto_tau, 1.03)
    zmcr_vv.setParamEffect(vv_norm, 1.2)
    zmcr_vv.setParamEffect(jec, 1.05)
    zmcr_vv.setParamEffect(id_mu, 1.02)
    zmcr_vv.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_vv)

    #zmcr_hbbHist = background['zmcr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    zmcr_hbbTemplate = template(background,'Hbb','nominal','zmcr')
    zmcr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,zmcr_hbbTemplate[0].max()*2) for i,b in enumerate(zmcr_hbbTemplate[0])])
    zmcr_hbbBinYields = zmcr_hbbBinYields * deepak15_weight['notag']['Hbb']
    zmcr_hbbObservable = rl.Observable('recoil', zmcr_hbbTemplate[1])
    zmcr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, zmcr_hbbObservable, zmcr_hbbBinYields)
    zmcr_hbb.setParamEffect(lumi, 1.027)
    zmcr_hbb.setParamEffect(trig_met, 1.01)
    zmcr_hbb.setParamEffect(veto_tau, 1.03)
    zmcr_hbb.setParamEffect(hbb_norm, 1.2)
    zmcr_hbb.setParamEffect(jec, 1.05)
    zmcr_hbb.setParamEffect(id_mu, 1.02)
    zmcr_hbb.setParamEffect(iso_mu, 1.02)
    zmcr.addSample(zmcr_hbb)

    ###
    # End of double muon control region
    ###

    ###
    ###
    # Double electron control region
    ###
    ###

    ch_name = 'zecr-'+model_id
    zecr = rl.Channel(ch_name)
    model.addChannel(zecr)

    ###
    # Add data distribution to the channel
    ###

    #zecr.setObservation(template(data['zecr'].integrate('process', 'SingleElectron').integrate('systematic','nominal'), 'recoil'))
    zecr.setObservation(template(data,'SingleElectron','nominal','zecr'))

    #zecr_dyjetsHist = background['zecr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    zecr_dyjetsTemplate = template(background,'DY+jets','nominal','zecr')
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

    #zecr_ttHist = background['zecr'].integrate('process', 'TT').integrate('systematic','nominal')
    zecr_ttTemplate = template(background,'TT','nominal','zecr')
    zecr_ttBinYields = np.array([rl.IndependentParameter(ch_name+'_tt_bin_%d' % i,b,0,zecr_ttTemplate[0].max()*2) for i,b in enumerate(zecr_ttTemplate[0])])
    zecr_ttBinYields = zecr_ttBinYields * deepak15_weight['notag']['TT']
    zecr_ttObservable = rl.Observable('recoil', zecr_ttTemplate[1])
    zecr_tt = rl.ParametericSample(ch_name+'_tt', rl.Sample.BACKGROUND, zecr_ttObservable, zecr_ttBinYields)
    zecr_tt.setParamEffect(lumi, 1.027)
    zecr_tt.setParamEffect(trig_e, 1.01)
    zecr_tt.setParamEffect(veto_tau, 1.03)
    zecr_tt.setParamEffect(tt_norm, 1.4)
    zecr_tt.setParamEffect(jec, 1.05)
    zecr_tt.setParamEffect(id_e, 1.02)
    zecr_tt.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_tt)

    #zecr_stHist = background['zecr'].integrate('process', 'ST').integrate('systematic','nominal')
    zecr_stTemplate = template(background,'ST','nominal','zecr')
    zecr_stBinYields = np.array([rl.IndependentParameter(ch_name+'_st_bin_%d' % i,b,0,zecr_stTemplate[0].max()*2) for i,b in enumerate(zecr_stTemplate[0])])
    zecr_stBinYields = zecr_stBinYields * deepak15_weight['notag']['ST']
    zecr_stObservable = rl.Observable('recoil', zecr_stTemplate[1])
    zecr_st = rl.ParametericSample(ch_name+'_st', rl.Sample.BACKGROUND, zecr_stObservable, zecr_stBinYields)
    zecr_st.setParamEffect(lumi, 1.027)
    zecr_st.setParamEffect(trig_e, 1.01)
    zecr_st.setParamEffect(veto_tau, 1.03)
    zecr_st.setParamEffect(st_norm, 1.2)
    zecr_st.setParamEffect(jec, 1.05)
    zecr_st.setParamEffect(id_e, 1.02)
    zecr_st.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_st)

    #zecr_vvHist = background['zecr'].integrate('process', 'VV').integrate('systematic','nominal')
    zecr_vvTemplate = template(background,'VV','nominal','zecr')
    zecr_vvBinYields = np.array([rl.IndependentParameter(ch_name+'_vv_bin_%d' % i,b,0,zecr_vvTemplate[0].max()*2) for i,b in enumerate(zecr_vvTemplate[0])])
    zecr_vvBinYields = zecr_vvBinYields * deepak15_weight['notag']['VV']
    zecr_vvObservable = rl.Observable('recoil', zecr_vvTemplate[1])
    zecr_vv = rl.ParametericSample(ch_name+'_vv', rl.Sample.BACKGROUND, zecr_vvObservable, zecr_vvBinYields)
    zecr_vv.setParamEffect(lumi, 1.027)
    zecr_vv.setParamEffect(trig_e, 1.01)
    zecr_vv.setParamEffect(veto_tau, 1.03)
    zecr_vv.setParamEffect(vv_norm, 1.2)
    zecr_vv.setParamEffect(jec, 1.05)
    zecr_vv.setParamEffect(id_e, 1.02)
    zecr_vv.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_vv)

    #zecr_hbbHist = background['zecr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    zecr_hbbTemplate = template(background,'Hbb','nominal','zecr')
    zecr_hbbBinYields = np.array([rl.IndependentParameter(ch_name+'_hbb_bin_%d' % i,b,0,zecr_hbbTemplate[0].max()*2) for i,b in enumerate(zecr_hbbTemplate[0])])
    zecr_hbbBinYields = zecr_hbbBinYields * deepak15_weight['notag']['Hbb']
    zecr_hbbObservable = rl.Observable('recoil', zecr_hbbTemplate[1])
    zecr_hbb = rl.ParametericSample(ch_name+'_hbb', rl.Sample.BACKGROUND, zecr_hbbObservable, zecr_hbbBinYields)
    zecr_hbb.setParamEffect(lumi, 1.027)
    zecr_hbb.setParamEffect(trig_e, 1.01)
    zecr_hbb.setParamEffect(veto_tau, 1.03)
    zecr_hbb.setParamEffect(hbb_norm, 1.2)
    zecr_hbb.setParamEffect(jec, 1.05)
    zecr_hbb.setParamEffect(id_e, 1.02)
    zecr_hbb.setParamEffect(reco_e, 1.02)
    zecr.addSample(zecr_hbb)

    ###
    # End of double electron control region
    ###

    ###
    ###
    # Single photon control region
    ###
    ###

    ch_name = 'gcr-'+model_id
    gcr = rl.Channel(ch_name)
    model.addChannel(gcr)

    ###
    # Add data distribution to the channel
    ###

    #gcr.setObservation(template(data['gcr'].integrate('process', 'SinglePhoton').integrate('systematic','nominal'), 'recoil'))
    gcr.setObservation(template(data,'SinglePhoton','nominal','gcr'))

    #gcr_gjetsHist = background['gcr'].integrate('process', 'G+jets').integrate('systematic','nominal')
    gcr_gjetsTemplate = template(background,'G+jets','nominal','gcr')
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

    #gcr_qcdHist = background['gcr'].integrate('process', 'QCD').integrate('systematic','nominal')
    gcr_qcdTemplate = template(background,'QCD','nominal','gcr')
    gcr_qcdBinYields = np.array([rl.IndependentParameter(ch_name+'_qcd_bin_%d' % i,b,0,gcr_qcdTemplate[0].max()*2) for i,b in enumerate(gcr_qcdTemplate[0])])
    gcr_qcdBinYields = gcr_qcdBinYields * deepak15_weight['notag']['QCD']
    gcr_qcdObservable = rl.Observable('recoil', gcr_qcdTemplate[1])
    gcr_qcd = rl.ParametericSample(ch_name+'_qcd', rl.Sample.BACKGROUND, gcr_qcdObservable, gcr_qcdBinYields)
    gcr_qcd.setParamEffect(lumi, 1.027)
    gcr_qcd.setParamEffect(trig_pho, 1.01)
    gcr_qcd.setParamEffect(veto_tau, 1.03)
    gcr_qcd.setParamEffect(qcdpho_norm, 2.0)
    gcr_qcd.setParamEffect(jec, 1.05)
    gcr_qcd.setParamEffect(id_pho, 1.02)
    gcr.addSample(gcr_qcd)

    return model

if __name__ == '__main__':
    if not os.path.exists('datacards'):
        os.mkdir('datacards')
    parser = OptionParser()
    parser.add_option('-g', '--category', help='category', dest='category', default='')
    parser.add_option('-y', '--year', help='year', dest='year', default='')
    parser.add_option('-m', '--mass', help='mass', dest='mass', default=None)
    (options, args) = parser.parse_args()
    grouping=False
    model_dict={}
    for category in ['monojet','monohs']:
        if options.category and options.category not in category: continue
        if category=='monojet':
            with open('data/darkhiggs'+options.year+'-'+category+'.model', "wb") as fout:
                pickle.dump(model(category,options.year,options.mass,grouping), fout, protocol=2)
        elif category=='monohs':
            for mass in ['mass0','mass1','mass2','mass3','mass4']:
                if options.mass and options.mass not in mass: continue
                with open('data/darkhiggs'+options.year+'-'+category+'-'+mass+'.model', "wb") as fout:
                    pickle.dump(model(category,options.year,mass,grouping), fout, protocol=2)

