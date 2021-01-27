import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import load, save



def scale_file(file):

    print('Loading file:',file)    
    hists=load(file)

    pd = []
    for d in hists['sumw'].identifiers('dataset'):
        dataset = d.name
        if dataset.split("____")[0] not in pd: pd.append(dataset.split("____")[0])
    print('List of primary datasets:',pd)

    ##
    # Aggregate all the histograms that belong to a single dataset
    ##

    dataset = hist.Cat("dataset", "dataset", sorting='placement')
    dataset_cats = ("dataset",)
    dataset_map = OrderedDict()
    for pdi in pd:
        dataset_map[pdi] = (pdi+"*",)
    for key in hists.keys():
        hists[key] = hists[key].group(dataset_cats, dataset, dataset_map)
    print('Datasets aggregated')

    return scale(hists)

def scale_directory(directory):

    hists = {}
    for filename in os.listdir(directory):
        if '.merged' not in filename: continue
        print('Opening:', filename)
        hin = load(directory+'/'+filename)
        hists.update(hin)

    return scale(hists)

def scale(hists):

    ###
    # Rescaling MC histograms using the xsec weight
    ###

    scale={}
    for d in hists['sumw'].identifiers('dataset'):
        scale[d]=hists['sumw'].integrate('dataset', d).values(overflow='all')[()][1]
    print('Sumw extracted')

    for key in hists.keys():
        if key=='sumw': continue
        for d in hists[key].identifiers('dataset'):
            if 'MET' in d.name or 'SingleElectron' in d.name or 'SinglePhoton' in d.name or 'EGamma' in d.name: continue
            hists[key].scale({d:1/scale[d]},axis='dataset')
    print('Histograms scaled')


    ###
    # Defining 'process', to aggregate different samples into a single process
    ##

    process = hist.Cat("process", "Process", sorting='placement')
    cats = ("dataset",)
    sig_map = OrderedDict()
    bkg_map = OrderedDict()
    data_map = OrderedDict()
    bkg_map["Hbb"] = ("*HTo*")
    bkg_map["DY+HF"] = ("HF--DYJets*",)
    bkg_map["DY+LF"] = ("LF--DYJets*",)
    bkg_map["DY+jetsLO"] = ("lo--DYJets*",)
    bkg_map["DY+jetsNNLO"] = ("nnlo--DYJets*",)
    #bkg_map["VV"] = (["WW*","WZ*","ZZ*"],)
    bkg_map["WW"] = ("WW*", )
    bkg_map["WZ"] = ("WZ*", )
    bkg_map["ZZ"] = ("ZZ*", )
    bkg_map["ST"] = ("ST*",)
    bkg_map["TT"] = ("TT*",)
    bkg_map["W+HF"] = ("HF--WJets*",)
    bkg_map["W+LF"] = ("LF--WJets*",)
    bkg_map["W+jetsLO"] = ("lo--WJets*",)
    bkg_map["W+jetsNNLO"] = ("nnlo--WJets*",)
    bkg_map["Z+HF"] = ("HF--ZJetsToNuNu*",)
    bkg_map["Z+LF"] = ("LF--ZJetsToNuNu*",)
    bkg_map["Z+jetsLO"] = ("lo--ZJets*",)
    bkg_map["Z+jetsNNLO"] = ("nnlo--ZJets*",)
    bkg_map["G+HF"] = ("HF--GJets*",)
    bkg_map["G+LF"] = ("LF--GJets*",)
    bkg_map["G+jetsLO"] = ("lo--GJets*",)
    bkg_map["G+jetsNNLO"] = ("nnlo--GJets*",)
    bkg_map["QCD"] = ("*QCD*",)
    data_map["MET"] = ("MET*", )
    data_map["SingleElectron"] = ("SingleElectron*", )
    data_map["SinglePhoton"] = ("SinglePhoton*", )
    data_map["EGamma"] = ("EGamma*", )
    for signal in hists['sumw'].identifiers('dataset'):
        if 'mhs' not in str(signal): continue
        print(signal)
        sig_map[str(signal)] = (str(signal),)  ## signals
    print('Processes defined')
    
    ###
    # Storing signal and background histograms
    ###
    bkg_hists={}
    sig_hists={}
    data_hists={}
    for key in hists.keys():
        bkg_hists[key] = hists[key].group(cats, process, bkg_map)
        sig_hists[key] = hists[key].group(cats, process, sig_map)
        data_hists[key] = hists[key].group(cats, process, data_map)
    print('Histograms grouped')

    return bkg_hists, sig_hists, data_hists

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--file', help='file', dest='file')
    parser.add_option('-d', '--directory', help='directory', dest='directory')
    (options, args) = parser.parse_args()

    if options.directory: 
        bkg_hists, sig_hists, data_hists = scale_directory(options.directory)
        name = options.directory
    if options.file: 
        bkg_hists, sig_hists, data_hists = scale_file(options.file)
        name = options.file.split(".")[0]

    hists={
        'bkg': bkg_hists,
        'sig': sig_hists,
        'data': data_hists
    }
    save(hists,name+'.scaled')
