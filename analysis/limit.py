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
    # Tau veto
    ###

    veto_tau = rl.NuisanceParameter('veto_tau', 'lnN')


    ###
    ###
    # Shape systematics
    ###
    ###

    ###
    # JEC/JER
    ###
    
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

    sr_zvvHist = background['sr'].integrate('process', 'Z+jets').integrate('systematic','nominal')
    sr_zvvTemplate = template(sr_zvvHist, 'recoil')
    sr_zvvMC =  rl.TemplateSample(ch_name+'_zvvMC', rl.Sample.BACKGROUND, sr_zvvTemplate)
    #sr_zvvMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_zvvHist.axis('recoil').edges(overflow='all'))-1))
    
    sr_zvvBinYields = np.array([rl.IndependentParameter(ch_name+'_zvv_bin_%d' % i, b, 0, sr_zvvTemplate[0].max()*2) for i,b in enumerate(sr_zvvTemplate[0])]) 
    sr_zvvBinYields = sr_zvvBinYields * hf_fraction_weight['0tag']['Z+jets']
    sr_zvvObservable = rl.Observable('recoil', sr_zvvHist.axis('recoil').edges())
    sr_zvv = rl.ParametericSample(ch_name+'_zvv', rl.Sample.BACKGROUND, sr_zvvObservable, sr_zvvBinYields)

    sr.addSample(sr_zvv)

    ###    
    # W(->lnu)+jets data-driven model                
    ### 

    sr_wjetsHist = background['sr'].integrate('process', 'W+jets').integrate('systematic','nominal')
    sr_wjetsTemplate = template(sr_wjetsHist, 'recoil')
    sr_wjetsMC =  rl.TemplateSample(ch_name+'_wjetsMC', rl.Sample.BACKGROUND, sr_wjetsTemplate)
    #sr_wjetsMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_wjetsHist.axis('recoil').edges(overflow='all'))-1))

    sr_wjetsBinYields = np.array([rl.IndependentParameter(ch_name+'_wjets_bin_%d' % i,b,0,sr_wjetsTemplate[0].max()*2) for i,b in enumerate(sr_wjetsTemplate[0])]) 
    sr_wjetsBinYields = sr_wjetsBinYields * hf_fraction_weight['0tag']['W+jets']
    sr_wjetsObservable = rl.Observable('recoil', sr_wjetsHist.axis('recoil').edges())
    sr_wjets = rl.ParametericSample(ch_name+'_wjets', rl.Sample.BACKGROUND, sr_wjetsObservable, sr_wjetsBinYields)
    sr.addSample(sr_wjets)

    ###    
    # top-antitop data-driven model                                                                                                                                                                  
    ### 

    sr_ttbarHist = background['sr'].integrate('process', 'TT').integrate('systematic','nominal')
    sr_ttbarTemplate = template(sr_ttbarHist, 'recoil')
    sr_ttbarMC =  rl.TemplateSample(ch_name+'_ttbarMC', rl.Sample.BACKGROUND, sr_ttbarTemplate)
    #sr_ttbarMC.setParamEffect(jec, np.random.normal(loc=1, scale=0.01, size=len(sr_ttbarHist.axis('recoil').edges(overflow='all'))-1))

    # these parameters are large, should probably log-transform them
    sr_ttbarBinYields = np.array([rl.IndependentParameter(ch_name+'_ttbar_bin_%d' % i,b,0,sr_ttbarTemplate[0].max()*2) for i,b in enumerate(sr_ttbarTemplate[0])])
    sr_ttbarBinYields = sr_ttbarBinYields * deepak15_weight['0tag']['TT']
    sr_ttbarObservable = rl.Observable('recoil', sr_ttbarHist.axis('recoil').edges())
    sr_ttbar = rl.ParametericSample(ch_name+'_ttbar', rl.Sample.BACKGROUND, sr_ttbarObservable, sr_ttbarBinYields)
    sr.addSample(sr_ttbar)

    ###
    # Other MC-driven processes
    ###

    sr_singletopHist = background['sr'].integrate('process', 'ST').integrate('systematic','nominal')
    sr_singletopTemplate = template(sr_singletopHist, 'recoil')
    sr_singletop = rl.TemplateSample(ch_name+'_singletop', rl.Sample.BACKGROUND, sr_singletopTemplate)
    sr_singletop.setParamEffect(lumi, 1.027)
    sr_singletop.setParamEffect(stop_Norm, 1.2)
    sr_singletop.setParamEffect(trig_met, 1.01)
    sr_singletop.setParamEffect(veto_tau, 1.03)
    sr.addSample(sr_singletop)

    sr_dyHist = background['sr'].integrate('process', 'DY+jets').integrate('systematic','nominal')
    sr_dyTemplate = template(sr_dyHist, 'recoil')
    sr_dy = rl.TemplateSample(ch_name+'_dy', rl.Sample.BACKGROUND, sr_dyTemplate)
    sr_dy.setParamEffect(lumi, 1.027)
    sr_dy.setParamEffect(dy_Norm, 1.2)
    sr_dy.setParamEffect(trig_met, 1.01)
    sr_dy.setParamEffect(veto_tau, 1.03)
    sr.addSample(sr_dy)

    sr_dibosonHist = background['sr'].integrate('process', 'VV').integrate('systematic','nominal')
    sr_dibosonTemplate = template(sr_dibosonHist, 'recoil')
    sr_diboson = rl.TemplateSample(ch_name+'_diboson', rl.Sample.BACKGROUND, sr_dibosonTemplate)
    sr_diboson.setParamEffect(lumi, 1.027)
    sr_diboson.setParamEffect(VV_Norm, 1.2)
    sr_diboson.setParamEffect(trig_met, 1.01)
    sr_diboson.setParamEffect(veto_tau, 1.03)
    sr.addSample(sr_diboson)

    sr_higgsHist =background['sr'].integrate('process', 'Hbb').integrate('systematic','nominal')
    sr_higgsTemplate = template(sr_higgsHist, 'recoil')
    sr_higgs = rl.TemplateSample(ch_name+'_higgs', rl.Sample.BACKGROUND, sr_higgsTemplate)
    sr_higgs.setParamEffect(lumi, 1.027)
    sr_higgs.setParamEffect(Hbb_Norm, 1.2)
    sr_higgs.setParamEffect(trig_met, 1.01)
    sr_higgs.setParamEffect(veto_tau, 1.03)
    sr.addSample(sr_higgs)

    for s in signal['sr'].identifiers('process'):
        print(s)
        sr_dmHist = signal['sr'].integrate('process', s).integrate('systematic','nominal')
        sr_dmTemplate = template(sr_dmHist, 'recoil')
        sr_dm = rl.TemplateSample(ch_name+'_'+str(s), rl.Sample.SIGNAL, sr_dmTemplate)
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
    grouping=False
    darkhiggs_model('datacards',mass,category,year,grouping)
