import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import save

def scale(folder):
    ###
    #
    # Loading MC histograms
    #
    ###

    mc_hists={}
    pd = []
    dirname = folder
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
    bkg_map["Gjets"] = ("GJets*",)
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

    return signal_hists, bkg_hists, data_hists

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    (options, args) = parser.parse_args()

    signal_hists, bkg_hists, data_hists = scale(options.folder)
    hists={}
    hists['signal']=signal_hists
    hists['bkg']=bkg_hists
    hists['data']=data_hists
    save(hists,'hists_'+options.folder+'.coffea')
