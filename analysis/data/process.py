#!/usr/bin/env python
'''@package docstring
Just a giant list of processes and properties
'''

processes =    {

    #data
    'MET':('MET','Data',1),
    'EGamma':('EGamma','Data',1),
    #'MuonEG':('MuonEG','Data',1),
    #'SingleMuon':('SingleMuon','Data',1),
    #'DoubleMuon':('DoubleMuon','Data',1),
    #'DoubleEG':('DoubleEG','Data',1),
    'SingleElectron':('SingleElectron','Data',1),
    'SinglePhoton':('SinglePhoton','Data',1),
    #'JetHT':('JetHT','Data',1),
    
    # inclusive NLO V+jets 
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_nlo','MC',6025.2),
    'DYJetsToNuNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZtoNuNu_nlo','MC',11433.),
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WJets_nlo','MC',61527.),
    'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_lo_incl','MC',50400.),
    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_lo_incl_CP5','MC',50400.),

    # LO Z->nunu
    'ZJetsToNuNu_HT-100To200_13TeV-madgraph':('ZtoNuNu_ht100to200','MC',303.4),
    'ZJetsToNuNu_HT-200To400_13TeV-madgraph':('ZtoNuNu_ht200to400','MC',91.71),
    'ZJetsToNuNu_HT-400To600_13TeV-madgraph':('ZtoNuNu_ht400to600','MC',13.10),
    'ZJetsToNuNu_HT-600To800_13TeV-madgraph':('ZtoNuNu_ht600to800','MC',3.248),
    'ZJetsToNuNu_HT-800To1200_13TeV-madgraph':('ZtoNuNu_ht800to1200','MC',1.496),
    'ZJetsToNuNu_HT-1200To2500_13TeV-madgraph':('ZtoNuNu_ht1200to2500','MC',0.3425),
    'ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph':('ZtoNuNu_ht2500toinf','MC',0.005268),

    # LO Z->ll
    #2018
    #'DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht100to200_CP5','MC',160.8),
    'DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht100to200_CP5','MC',148.9),
    #'DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht200to400_CP5','MC',48.63),
    'DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht200to400_CP5','MC',40.90),
    #'DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht400to600_CP5','MC',6.975),
    'DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht400to600_CP5','MC',5.676),
    #'DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht600to800_CP5','MC',1.756),
    'DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht600to800_CP5','MC',1.359),
    #'DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht800to1200_CP5','MC',0.8094),
    'DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht800to1200_CP5','MC',0.6236),
    #'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht1200to2500_CP5','MC',0.1931),
    'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht1200to2500_CP5','MC',0.1512),
    #'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht2500toinf_CP5','MC',0.003513),
    'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8':('ZJets_ht2500toinf_CP5','MC',0.003659),
    #2017
    'DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht100to200_CP5','MC',160.8),    
    #'DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht100to200_CP5','MC',148.9),
    'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht200to400_CP5','MC',48.63),
    #'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht200to400_CP5','MC',40.90),   
    'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht400to600_CP5','MC',6.986),
    #'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht400to600_CP5','MC',5.676),
    'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht600to800_CP5','MC',1.743),
    #'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht600to800_CP5','MC',1.359),
    'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht800to1200_CP5','MC',0.8044),
    #'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht800to1200_CP5','MC',0.6236),
    'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht1200to2500_CP5','MC',0.1927),
    #'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht1200to2500_CP5','MC',0.1512),
    'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht2500toinf_CP5','MC',0.003476),
    #'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8':('ZJets_ht2500toinf_CP5','MC',0.003659),
    #2016
    'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_ht100to200','MC',148.9),
    'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_ht200to400','MC',40.90),
    'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_ht400to600','MC',5.676), 
    'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_ht600to800','MC',1.359),
    'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_ht800to1200','MC',0.6236),
    'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_ht1200to2500','MC',0.1512),
    'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_ht2500toinf','MC',0.003659),
    
    # LO W->lnu
    #2018/2017
    'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_ht100to200_CP5','MC',1392),
    'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_ht200to400_CP5','MC',410.3),
    'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_ht400to600_CP5','MC',57.85),
    'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_ht600to800_CP5','MC',12.95),
    'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_ht800to1200_CP5','MC',5.451),
    'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_ht1200to2500_CP5','MC',1.084),
    'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8':('WJets_ht2500toinf_CP5','MC',0.008061),
    #2016
    'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_ht100to200','MC',1346),
    'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_ht200to400','MC',359.1),
    'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_ht400to600','MC',48.96),
    'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_ht600to800','MC',12.05),
    'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_ht800to1200','MC',5.487),
    'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_ht1200to2500','MC',1.329),
    'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_ht2500toinf','MC',0.03215),

    # LO gamma
    #2018/2017
    'GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8':('GJets_ht100to200','MC',8604),
    'GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8':('GJets_ht200to400','MC',2194),
    'GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8':('GJets_ht400to600','MC',257.8),
    'GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8':('GJets_ht600toinf','MC',84.92),
    #2016
    #'GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('GJets_ht40to100','MC',23080.0),
    'GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('GJets_ht100to200','MC',9239),
    'GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('GJets_ht200to400','MC',2303),
    'GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('GJets_ht400to600','MC',274.5),
    'GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('GJets_ht600toinf','MC',93.53),

    # NLO W->lnu
    'WJetsToLNu_Wpt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WJets_pt50to100','MC',3258),
    'WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WJets_pt100to250','MC',677.82),
    'WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WJets_pt250to400','MC',24.083),
    'WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WJets_pt400to600','MC',3.0563),
    'WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WJets_pt600toinf','MC',0.4602),
    'W1JetsToLNu_LHEWpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8':('W1JetsToLNu_WpT50to150_CP5','MC',2496),
    'W1JetsToLNu_LHEWpT_100-150_TuneCP5_13TeV-amcnloFXFX-pythia8':('W1JetsToLNu_WpT100to150_CP5','MC',292.4),
    'W1JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8':('W1JetsToLNu_WpT150to250_CP5','MC',75.69),
    'W1JetsToLNu_LHEWpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8':('W1JetsToLNu_WpT250to400_CP5','MC',7.69),
    'W1JetsToLNu_LHEWpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8':('W1JetsToLNu_WpT400toinf_CP5','MC',0.84),
    'W2JetsToLNu_LHEWpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8':('W2JetsToLNu_WpT50to150_CP5','MC',1275),
    'W2JetsToLNu_LHEWpT_100-150_TuneCP5_13TeV-amcnloFXFX-pythia8':('W2JetsToLNu_WpT100to150_CP5','MC',303.7),
    'W2JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8':('W2JetsToLNu_WpT150to250_CP5','MC',106.4),
    'W2JetsToLNu_LHEWpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8':('W2JetsToLNu_WpT250to400_CP5','MC',10.98),
    'W2JetsToLNu_LHEWpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8':('W2JetsToLNu_WpT400toinf_CP5','MC',3.19),

    # NLO Z->ll
    'DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_pt50to100','MC',374.6800),
    'DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_pt100to250','MC',86.5200),
    'DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_pt250to400','MC',3.3247),
    'DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_pt400to650','MC',0.4491),
    'DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_pt650toinf','MC',0.0422),
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZJets_m10','MC',18610),
    'DY1JetsToLL_M-50_LHEZpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z1Jets_ZpT50to150_CP5','MC',312.4),
    'DY1JetsToLL_M-50_LHEZpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z1Jets_ZpT150to250_CP5','MC',9.462),
    'DY1JetsToLL_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z1Jets_ZpT250to400_CP5','MC',1.087),
    'DY1JetsToLL_M-50_LHEZpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z1Jets_ZpT400toinf_CP5','MC',0.1202),
    'DY2JetsToLL_M-50_LHEZpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z2Jets_ZpT50to150_CP5','MC',168.1),
    'DY2JetsToLL_M-50_LHEZpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z2Jets_ZpT150to250_CP5','MC',15.62),
    'DY2JetsToLL_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z2Jets_ZpT250to400_CP5','MC',2.718),
    'DY2JetsToLL_M-50_LHEZpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z2Jets_ZpT400toinf_CP5','MC',0.447),
    'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8':('ZJets_inclNLO_CP5','MC',6473),

    # NLO Z->nunu
    'DYJetsToNuNu_PtZ-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZtoNuNu_pt50to100','MC',593.9),
    'DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZtoNuNu_pt100to250','MC',3*54.8229),
    'DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZtoNuNu_pt250to400','MC',3*2.0705),
    'DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZtoNuNu_pt400to650','MC',3*0.2779),
    'DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZtoNuNu_pt650toinf','MC',3*0.0261),
    'Z1JetsToNuNu_M-50_LHEZpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z1JetsToNuNu_ZpT150to250_CP5','MC',18.44),
    'Z1JetsToNuNu_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z1JetsToNuNu_ZpT250to400_CP5','MC',2.06),
    'Z1JetsToNuNu_M-50_LHEZpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z1JetsToNuNu_ZpT400toinf_CP5','MC',0.2154),
    'Z2JetsToNuNu_M-50_LHEZpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z2JetsToNuNu_ZpT150to250_CP5','MC',32.93),
    'Z2JetsToNuNu_M-50_LHEZpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z2JetsToNuNu_ZpT250to400_CP5','MC',5.225),
    'Z2JetsToNuNU_M-50_LHEZpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8':('Z2JetsToNuNu_ZpT400toinf_CP5','MC',0.8596),
    

        # QCD
    #'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht100to200','MC',27990000),
    'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht200to300','MC',1735000),
    'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht300to500','MC',366800),
    'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht500to700','MC',29370),
    'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht700to1000','MC',6524),
    'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht1000to1500','MC',1064),
    'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht1500to2000','MC',121.5),
    'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('QCD_ht2000toinf','MC',25.42),

        #'QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht100to200','MC',27990000),
    'QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht200to300','MC',1735000),
    'QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht300to500','MC',366800),
    'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht500to700','MC',29370),
    'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht700to1000','MC',6524),
    'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht1000to1500','MC',1064),
    'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht1500to2000','MC',121.5),
    'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8':('QCD_ht2000toinf','MC',25.42),


        #'QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht100to200_CP5','MC',27990000),
    'QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht200to300_CP5','MC',1735000),
    'QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht300to500_CP5','MC',366800),
    'QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht500to700_CP5','MC',29370),
    'QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht700to1000_CP5','MC',6524),
    'QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht1000to1500_CP5','MC',1064),
    'QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht1500to2000_CP5','MC',121.5),
    'QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8':('QCD_ht2000toinf_CP5','MC',25.42),

        # Single tops

        'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8':('SingleTop_tTbar_lep','MC',26.22),
    'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8':('SingleTop_tT_lep','MC',44.07),
    'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1':('SingleTop_tW','MC',35.85),
    'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1':('SingleTop_tbarW','MC',35.85),
    'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1':('SingleTop_tTbar','MC',80.95),
    'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1':('SingleTop_tT','MC',136.02),
    'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1':('SingleTop_s_lep','MC',10.11),
    'ST_t-channel_5f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1':('SingleTop_t_lep','MC',216.99),
    'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1':('SingleTop_tchannel','MC',44.3),
    
        'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':('SingleTop_tbarW_CP5','MC', 35.6),
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':('SingleTop_tW_CP5','MC', 35.6),
    'ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':('SingleTop_tTbar_CP5','MC',80.95),
    'ST_t-channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV-powhegV2-madspin':('SingleTop_tTbar','MC',80.95),
    'ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':('SingleTop_tT_CP5','MC',136.02),
    'ST_t-channel_top_4f_inclusiveDecays_TuneCUETP8M2T4_13TeV-powhegV2-madspin':('SingleTop_tT','MC',136.02),

        'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8':('SingleTop_tTbar_CP518','MC',80.95),
    'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8':('SingleTop_tT_CP518','MC',136.02),
    'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8':('SingleTop_tbarW_noHad','MC',35.6*(1-0.68*0.68)),
    'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1':('SingleTop_tbarW_noHad','MC',35.6*(1-0.68*0.68)),
    'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8':('SingleTop_tW_noHad','MC',35.6*(1-0.68*0.68)),
    'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1':('SingleTop_tW_noHad','MC',35.6*(1-0.68*0.68)),


        # ttbar
    'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8':           ('TTbar_2L_Powheg','MC',87.31),
    'TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8':      ('TTbar_1L_Powheg','MC',364.36),
    'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':                    ('TTbar_MLM','MC',831.76),
    'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':             ('TTbar_2L','MC',831.76*(1-0.68)*(1-0.68)),
    'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':    ('TTbar_1LT','MC',831.76*0.68*(1-0.68)),
    'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8': ('TTbar_1LTbar','MC',831.76*0.68*(1-0.68)),
    'TTJets_SingleLeptFromT_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8': ('TTbar_FXFX_1LT','MC',831.76*0.68*(1-0.68)),
    'TTJets_SingleLeptFromTbar_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8': ('TTbar_FXFX_1LTbar','MC',831.76*0.68*(1-0.68)),
    'TT_TuneEE5C_13TeV-powheg-herwigpp':                                ('TTbar_Herwig','MC',831.76),
    'TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8':                 ('TTbar_FXFX','MC',831.76),
    'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8':                        ('TTbar_FXFX','MC',831.76),
    'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8':                           ('TTbar_Powheg','MC',831.76),
    #'TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8':                   ('TTbar_PowhegISRDown','MC',831.76),
    #'TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8':                     ('TTbar_PowhegISRUp','MC',831.76),
    #'TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8':                       ('TTbar_PowhegTuneDown','MC',831.76),
    #'TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8':                         ('TTbar_PowhegTuneUp','MC',831.76),
    'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8':('TTTo2L2Nu','MC',88.288),
    'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8':('TTZToLLNuNu','MC',0.2529),
    'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8':('TTZToQQ','MC',0.5297),
    'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8':('TTWJetsToLNu','MC',0.2043),
    'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8':('TTWJetsToQQ','MC',0.4062),

        'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8':('TTTo2L2Nu_CP5','MC',88.288),
    'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8':('TTTo2L2Nu_CP5_PSweights','MC',88.288),
    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8':('TTToSemiLeptonic_CP5','MC',1),
    'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8':('TTToSemiLeptonic_CP5_PSweights','MC',1),
    'TTToHadronic_TuneCP5_13TeV-powheg-pythia8':('TTToHadronic_CP5','MC',1),
    'TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8':('TTToHadronic_CP5_PSweights','MC',1),


        # exotic top
    'tZq_ll_4f_13TeV-amcatnlo-pythia8':('SingleTop_tZll','MC',0.0758),
    'tZq_nunu_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1':('SingleTop_tZnunu','MC',0.1379),
    'TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8':('TTbar_GJets','MC',3.786),
    'TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8':('SingleTop_tG','MC',2.967),

        # EWK V+jets
    #'EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8':('ZtoNuNu_EWK','MC',10.04),
    #'EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8':('ZJets_EWK','MC',3.99),
    #'EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8':('WJets_EWKWPlus','MC',25.81),
    #'EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8':('WJets_EWKWMinus','MC',20.35),
    #'EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8':('ZtoNuNu_EWK','MC',10.04),
    #'EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8':('ZJets_EWK','MC',3.99),
    #'EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8':('WJets_EWKWPlus','MC',25.81),
    #'EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8':('WJets_EWKWMinus','MC',20.35),

        # regular dibosons
    'WW_TuneCUETP8M1_13TeV-pythia8':('Diboson_ww','MC',118.7),
    'WZ_TuneCUETP8M1_13TeV-pythia8':('Diboson_wz','MC',47.13),
    'ZZ_TuneCUETP8M1_13TeV-pythia8':('Diboson_zz','MC',16.523),
    'WW_TuneCP5_13TeV-pythia8':('Diboson_ww_CP5','MC',118.7),
    'WZ_TuneCP5_13TeV-pythia8':('Diboson_wz_CP5','MC',47.13),
    'ZZ_TuneCP5_13TeV-pythia8':('Diboson_zz_CP5','MC',16.523),
    
        # fancy dibosons
    'WWTo2L2Nu_13TeV-powheg':('WWTo2L2Nu','MC',(118.7-3.974)*0.1086*0.1086*9), #12.178
    'WWTo4Q_13TeV-powheg':('WWTo4Q','MC',51.723),
    'WWToLNuQQ_13TeV-powheg':('WWToLNuQQ','MC',49.997),
    'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8':('WZTo1L1Nu2Q','MC',10.71),
    'WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8':('WZTo1L3Nu','MC',3.033),
    'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8':('WZTo2L2Q','MC',5.595),
    'WZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8':('WZTo2Q2Nu','MC',6.317),
    'WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('WZTo3LNu','MC',4.692),
    'ZZTo2L2Nu_13TeV_powheg_pythia8':('ZZTo2L2Nu','MC',0.5644),
    'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8':('ZZTo2L2Q','MC',3.22),
    'ZZTo4L_13TeV_powheg_pythia8':('ZZTo4L','MC',1.212),
    'ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8':('ZZTo2Q2Nu_AMC','MC',4.072),
    'ZZTo2Q2Nu_13TeV_powheg_pythia8':('ZZTo2Q2Nu_Powheg','MC',3.824),
    'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8':('ZZTo4Q','MC',6.842),
    #'GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8':('ggZZTo2e2nu','MC',0.003956),
    #'GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8':('ggZZTo2mu2tau','MC',0.003956),
    #'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8':('ggZZTo2e2mu','MC',0.0073462),
    #'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8':('ggZZTo2e2tau','MC',0.0073462),
    #'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8':('ggZZTo2mu2tau','MC',0.0073462),
    #'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8':('ggZZTo4e','MC',0.0036478),
    #'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8':('ggZZTo4mu','MC',0.0036478),
    #'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8':('ggZZTo4tau','MC',0.0036478),
    #'WGstarToLNuEE_012Jets_13TeV-madgraph':('WGstarToLNuEE_012Jets_13TeV-madgraph','MC',3.526),
    #'WGstarToLNuMuMu_012Jets_13TeV-madgraph':('WGstarToLNuMuMu_012Jets_13TeV-madgraph','MC',2.793),
    #'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','MC',489.0),
    #'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8':('ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8','MC',123.9*1.06),        

        # Triboson
    #'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8':('WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8','MC',0.2086),
    #'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8':('WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8','MC',0.16510),
    #'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8':('WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8','MC',0.05565),
    #'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8':('ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8','MC',0.01398),
    #'WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8':('WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8','MC',0.2147),

        # Higgs->bb
    'ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8':('ZvvHbb_mH125','MC',0.08912),
    'ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8':('ZllHbb_mH125','MC',0.04865),
    'ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8':('ggZvvHbb_mH125','MC',0.014366),
    'ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8':('ggZllHbb_mH125','MC',0.007842),
    'WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8':('WmLNuHbb','MC',0.100),
    'WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8':('WpLNuHbb','MC',0.159),
    'ttHTobb_M125_13TeV_powheg_pythia8':('ttHbb','MC',0.506*0.5824),
    'ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8':('ttHbb','MC',0.506*0.5824),
    'GluGluHToBB_M125_13TeV_powheg_pythia8':('ggHbb','MC',48.48*0.5824),
    'VBFHToBB_M125_13TeV_amcatnlo_pythia8':('VBFHbb','MC',3.782*0.5824),
    'VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix':('VBFHbb','MC',3.782*0.5824),
    
        # B enriched V+jets samples
    #'DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_bQuarks','MC',71.77),
    #'DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_bQuarks_pt100to200','MC',3.027),
    #'DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_bQuarks_pt200toinf','MC',0.297),
    #'WBJetsToLNu_Wpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_bQuarks_pt100to200','MC',6.004),
    #'WBJetsToLNu_Wpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_bQuarks_pt200toinf','MC',0.8524),
    #'ZBJetsToNuNu_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZtoNuNu_bQuarks_pt100to200','MC',6.03),
    #'ZBJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZtoNuNu_bQuarks_pt200toinf','MC',0.5946),
    #'DYJetsToLL_BGenFilter_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_bHadrons','MC',229.7),
    #'DYJetsToLL_BGenFilter_Zpt-100to200_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_bHadrons_pt100to200','MC',3.506),
    #'DYJetsToLL_BGenFilter_Zpt-200toInf_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZJets_bHadrons_pt200toinf','MC',0.5083),
    #'WJetsToLNu_BGenFilter_Wpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_bHadrons_pt100to200','MC',26.1),
    #'WJetsToLNu_BGenFilter_Wpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('WJets_bHadrons_pt200toinf','MC',3.545),
    #'ZJetsToNuNu_BGenFilter_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZtoNuNu_bHadrons_pt100to200','MC',2.099),
    #'ZJetsToNuNu_BGenFilter_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8':('ZtoNuNu_bHadrons_pt200toinf','MC',0.3083),
    
        # Old Signals
    #'Mjet1DM_LO_MZprime_1000_Mhs_50_Mchi_100_gSM_0p25_gDM_1p0_th_0p01_13TeV_madgraph':('Mjet1DM_LO_MZprime_1000_Mhs_50_Mchi_100_gSM_0p25_gDM_1p0_th_0p01_13TeV_madgraph','MC',1)

        # New Signals

        'MonoHs_MZprime_500_Mhs_50_Mchi_150_2018':('MonoHs_MZprime_500_Mhs_50_Mchi_150','MC',0.5302),
    'MonoHs_MZprime_500_Mhs_70_Mchi_150_2018':('MonoHs_MZprime_500_Mhs_70_Mchi_150','MC',0.4331),                                                                                         
    'MonoHs_MZprime_500_Mhs_90_Mchi_150_2018':('MonoHs_MZprime_500_Mhs_90_Mchi_150','MC',0.3656),
    'MonoJet_MZprime_500_Mchi_150_2018':('MonoJet_MZprime_500_Mchi_150','MC',2.621),
    'MonoW_MZprime_500_Mchi_150_2018':('MonoW_MZprime_500_Mchi_150','MC',0.05872),
    'MonoZ_MZprime_500_Mchi_150_2018':('MonoZ_MZprime_500_Mchi_150','MC',0.01820),
}
