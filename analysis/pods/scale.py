import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import load, save



def scale_file(file):

    pd = []
    hists = load(file)

    for d in hists['sumw'].identifiers('dataset'):
        dataset = d.name
        print(dataset)
        if 'condor' in file:
            if dataset.split("___")[0] not in pd: pd.append(dataset.split("____")[0])
        elif 'spark' in file:
            if dataset.split("____")[0] not in pd: pd.append(dataset.split("____")[0])
    return scale(pd, hists)

def scale_directory(directory):

    hists={}
    pd = []
    for filename in os.listdir(directory):
        if '.pkl.gz' in filename:
            if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])
            with gzip.open(directory+'/'+filename) as fin:
                hin = pickle.load(fin)
                for k in hin.keys():
                    if k in hists: hists[k]+=hin[k]
                else: hists[k]=hin[k]
    
    return scale(pd, hists)

def scale(pd, hists):

    ##
    # Aggregate all the histograms that belong to a single sample
    ##

    dataset = hist.Cat("dataset", "dataset", sorting='placement')
    dataset_cats = ("dataset",)
    dataset_map = OrderedDict()
    for pdi in pd:
        #print(pdi)
        dataset_map[pdi] = (pdi+"*",)
    for key in hists.keys():
        hists[key] = hists[key].group(dataset_cats, dataset, dataset_map)

    ###
    # Rescaling MC histograms using the xsec weight
    ###

    scale={}
    for pdi in hists['sumw'].identifiers('dataset'):
        scale[pdi]=hists['sumw'].integrate('dataset', pdi).values(overflow='all')[()][1]
    for key in hists.keys():
        if key=='sumw': continue
        for pdi in hists[key].identifiers('dataset'):
            if 'MET' in pdi.name or 'SingleElectron' in pdi.name or 'SinglePhoton' in pdi.name or 'EGamma' in pdi.name: continue
            hists[key].scale({pdi:1/scale[pdi]},axis='dataset')

    ###
    # Defining 'process', to aggregate different samples into a single process
    ##

    process = hist.Cat("process", "Process", sorting='placement')
    bkg_cats = ("dataset",)
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

    signal_cats = ("dataset",)
    signal_map = OrderedDict() ### for signal samples
    signal_map["Mhs_50"] = ("*Mhs_50*",)  ## signals
    signal_map["Mhs_70"] = ("*Mhs_70*",)
    signal_map["Mhs_90"] = ("*Mhs_90*",)
    signal_map["MonoJet"] = ("MonoJet*",)  ## signals
    signal_map["MonoW"] = ("MonoW*",)    ## signals
    signal_map["MonoZ"] = ("MonoZ*",)    ## signals
    signal_hists = {}

    data_cats = ("dataset",)
    data_map = OrderedDict()
    data_map["MET"] = ("MET*", )
    data_map["SingleElectron"] = ("EGamma*", )
    data_map["SinglePhoton"] = ("EGamma*", )
    data_hists = {}

    ###
    # Storing signal and background histograms
    ###

    for key in hists.keys():
        signal_hists[key] = hists[key].group(signal_cats, process, signal_map)
        bkg_hists[key] = hists[key].group(bkg_cats, process, bkg_map)
        data_hists[key] = hists[key].group(data_cats, process, data_map)

    return signal_hists, bkg_hists, data_hists

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--file', help='file', dest='file')
    parser.add_option('-d', '--directory', help='directory', dest='directory')
    (options, args) = parser.parse_args()

    if options.directory: 
        signal_hists, bkg_hists, data_hists = scale_directory(options.directory)
        name = options.directory
    if options.file: 
        signal_hists, bkg_hists, data_hists = scale_file(options.file)
        name = options.file.split(".")[0]

    hists={}
    hists['signal']=signal_hists
    hists['bkg']=bkg_hists
    hists['data']=data_hists
    save(hists,'scaled_'+name+'.coffea')
