import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import load, save



def scale_file(file):
    return scale(load(file))

def scale_directory(directory):

    hists={}
    for filename in os.listdir(directory):
        if '.futures' not in filename: continue
        fin = directory+'/'+filename
        hin = load(fin)
        for k in hin.keys():
            if k in hists: hists[k]+=hin[k]
            else: hists[k]=hin[k]

    return scale(hists)

def scale(hists):

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
    map = OrderedDict()
    map["Hbb_merged"] = ("bb--*HToBB*")
    map["Hbb_unmerged"] = ("other--*HToBB*")
    map["DYJetsHF"] = ("HF--DYJets*",)
    map["DYJetsLF"] = ("LF--DYJets*",)
    map["VV"] = ("other--*_TuneCP5_13TeV-pythia8*",)
    map["VVbb"] = ("bb--*_TuneCP5_13TeV-pythia8*",)
    map["ST_merged"] = ("merged--ST*",)
    map["ST_unmerged"] = ("unmerged--ST*",)
    map["TTbqq"] = ("bqq--TT*",)
    map["TTqq"] = ("qq--TT*",)
    map["TTother"] = ("other--TT*",)
    map["WJetsHF"] = ("HF--WJets*",)
    map["WJetsLF"] = ("LF--WJets*",)
    map["ZJetsHF"] = ("HF--ZJetsToNuNu*",)
    map["ZJetsLF"] = ("LF--ZJetsToNuNu*",)
    map["GJetsHF"] = ("HF--GJets*",)
    map["GJetsLF"] = ("LF--GJets*",)
    map["MET"] = ("MET*", )
    map["SingleElectron"] = ("EGamma*", )
    map["SinglePhoton"] = ("EGamma*", )
    map["Monotop"] = ("Vector_Monotop*", )
    print('Processes defined')
    
    ###
    # Storing signal and background histograms
    ###

    for key in hists.keys():
        hists[key] = hists[key].group(cats, process, map)
    print('Histograms grouped')

    return hists

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--file', help='file', dest='file')
    parser.add_option('-d', '--directory', help='directory', dest='directory')
    (options, args) = parser.parse_args()

    if options.directory: 
        hists = scale_directory(options.directory)
        name = options.directory
    if options.file: 
        hists = scale_file(options.file)
        name = options.file.split(".")[0]

    save(hists,name+'.scaled')
