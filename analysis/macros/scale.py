import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import load, save

xsec = {
    ### 2018 signal, mhs = 50 GeV
    "Mz200_mhs50_Mdm100": 0.06606,
    "Mz200_mhs50_Mdm150": 0.02532,
    "Mz300_mhs50_Mdm100": 1.59,
    "Mz300_mhs50_Mdm150": 0.04397,
    "Mz500_mhs50_Mdm150": 0.9072,
    "Mz500_mhs50_Mdm250": 0.02005,
    "Mz500_mhs50_Mdm500": 0.001501,
    "Mz1000_mhs50_Mdm150": 0.5303,
    "Mz1000_mhs50_Mdm500": 0.003643,
    "Mz1000_mhs50_Mdm1000": 0.00005978,
    "Mz2000_mhs50_Mdm500": 0.02582,
    "Mz2000_mhs50_Mdm1000": 0.0002166,
    "Mz2000_mhs50_Mdm1500": 0.000002193,
    "Mz2500_mhs50_Mdm750": 0.005708,
    "Mz2500_mhs50_Mdm1250": 0.00005828,
    "Mz3000_mhs50_Mdm1000": 0.00135,
    "Mz3000_mhs50_Mdm1500": 0.00001537,

    ### 2018 signal, mhs = 70 GeV
    "Mz200_mhs70_Mdm100": 0.05611,
    "Mz200_mhs70_Mdm150": 0.02137,
    "Mz300_mhs70_Mdm100": 1.35,
    "Mz300_mhs70_Mdm150": 0.03773,
    "Mz500_mhs70_Mdm150": 0.7866,
    "Mz500_mhs70_Mdm250": 0.0176,
    "Mz500_mhs70_Mdm500": 0.001304,
    "Mz1000_mhs70_Mdm150": 0.4872,
    "Mz1000_mhs70_Mdm500": 0.003273,
    "Mz1000_mhs70_Mdm1000": 0.00005328,
    "Mz2000_mhs70_Mdm500": 0.02432,
    "Mz2000_mhs70_Mdm1000": 0.0001971,
    "Mz2000_mhs70_Mdm1500": 0.00000193,
    "Mz2500_mhs70_Mdm750": 0.005344,
    "Mz2500_mhs70_Mdm1250": 0.00005322,
    "Mz3000_mhs70_Mdm1000": 0.001265,
    "Mz3000_mhs70_Mdm1500": 0.00001412,

    ### 2018 signal, mhs = 90 GeV
    "Mz200_mhs90_Mdm100": 0.03795,
    "Mz200_mhs90_Mdm150": 0.01497,
    "Mz300_mhs90_Mdm100": 1.151,
    "Mz300_mhs90_Mdm150": 0.03218,
    "Mz500_mhs90_Mdm150": 0.6832,
    "Mz500_mhs90_Mdm250": 0.01529,
    "Mz500_mhs90_Mdm500": 0.001117,
    "Mz1000_mhs90_Mdm150": 0.4376,
    "Mz1000_mhs90_Mdm500": 0.002921,
    "Mz1000_mhs90_Mdm1000": 0.00004682,
    "Mz2000_mhs90_Mdm500": 0.02272,
    "Mz2000_mhs90_Mdm1000": 0.0001796,
    "Mz2000_mhs90_Mdm1500": 0.000001722,
    "Mz2500_mhs90_Mdm750": 0.005043,
    "Mz2500_mhs90_Mdm1250": 0.00004879,
    "Mz3000_mhs90_Mdm1000": 0.001193,
    "Mz3000_mhs90_Mdm1500": 0.00001292,

    ### 2018 signal, mhs = 110 GeV
    "Mz200_mhs110_Mdm150":   0.008799,
    "Mz300_mhs110_Mdm150":   0.02605,
    "Mz500_mhs110_Mdm150":   0.5925,
    "Mz500_mhs110_Mdm250":   0.01358,
    "Mz500_mhs110_Mdm500":   0.00104,
    "Mz1000_mhs110_Mdm150":  0.4073,
    "Mz1000_mhs110_Mdm500":  0.002749,
    "Mz1000_mhs110_Mdm1000": 0.00004688,
    "Mz2000_mhs110_Mdm500":  0.02286,
    "Mz2000_mhs110_Mdm1000": 0.0001817,
    "Mz2000_mhs110_Mdm1500": 0.000001889,
    "Mz2500_mhs110_Mdm750":  0.005186,
    "Mz2500_mhs110_Mdm1250": 0.00005072,
    "Mz3000_mhs110_Mdm1000": 0.001269,
    "Mz3000_mhs110_Mdm1500": 0.00001394,

    ### 2018 signal, mhs = 130 GeV
    "Mz200_mhs130_Mdm150":   0.003765,
    "Mz300_mhs130_Mdm150":   0.01682,
    "Mz500_mhs130_Mdm150":   0.4191,
    "Mz500_mhs130_Mdm250":   0.009572,
    "Mz500_mhs130_Mdm500":   0.0007286,
    "Mz1000_mhs130_Mdm150":  0.294,
    "Mz1000_mhs130_Mdm500":  0.001998,
    "Mz1000_mhs130_Mdm1000": 0.00003353,
    "Mz2000_mhs130_Mdm500":  0.01709,
    "Mz2000_mhs130_Mdm1000": 0.0001339,
    "Mz2000_mhs130_Mdm1500": 0.000001376,
    "Mz2500_mhs130_Mdm750":  0.003885,
    "Mz2500_mhs130_Mdm1250": 0.00003756,
    "Mz3000_mhs130_Mdm1000": 0.00095,
    "Mz3000_mhs130_Mdm1500": 0.00001024,

    ### 2018 signal, mhs = 150 GeV
    "Mz200_mhs150_Mdm150":   0.001266,
    "Mz300_mhs150_Mdm150":   0.006044,
    "Mz500_mhs150_Mdm150":   0.1564,
    "Mz500_mhs150_Mdm250":   0.003581,
    "Mz500_mhs150_Mdm500":   0.000271,
    "Mz1000_mhs150_Mdm150":  0.1131,
    "Mz1000_mhs150_Mdm500":  0.0007677,
    "Mz1000_mhs150_Mdm1000": 0.00001272,
    "Mz2000_mhs150_Mdm500":  0.006784,
    "Mz2000_mhs150_Mdm1000": 0.00005234,
    "Mz2000_mhs150_Mdm1500": 0.0000005154,
    "Mz2500_mhs150_Mdm750":  0.001556,
    "Mz2500_mhs150_Mdm1250": 0.00001474,
    "Mz3000_mhs150_Mdm1000": 0.0003786,
    "Mz3000_mhs150_Mdm1500": 0.000004007,
}

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
            if 'MET' in d.name or 'SingleElectron' in d.name or 'SinglePhoton' in d.name or 'EGamma' in d.name or 'BTagMu' in d.name: continue
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
    bkg_map['QCD-$\mu$ (l)'] = ('l--QCD*')
    bkg_map['QCD-$\mu$ (c)'] = ('c--QCD*')
    bkg_map['QCD-$\mu$ (cc)'] = ('cc--QCD*')
    bkg_map['QCD-$\mu$ (b)'] = ('b--QCD*')
    bkg_map['QCD-$\mu$ (bb)'] = ('bb--QCD*')
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
    bkg_map["QCD"] = ("QCD*",)
    data_map["MET"] = ("MET*", )
    data_map["SingleElectron"] = ("SingleElectron*", )
    data_map["SinglePhoton"] = ("SinglePhoton*", )
    data_map["EGamma"] = ("EGamma*", )
    data_map["BTagMu"] = ("BTagMu*", )
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
        data_hists[key] = hists[key].group(cats, process, data_map)
        sig_hists[key] = hists[key].group(cats, process, sig_map)
        for signal in sig_hists[key].identifiers('process'):
            print('Scaling '+str(signal)+' by xsec '+str(xsec[str(signal)]))
            sig_hists[key].scale({signal:xsec[str(signal)]},axis='process')
        
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
