###
# From https://github.com/sidnarayanan/PandaAnalysis/blob/0e411591911ae83dcbdc428496e7700445e844ae/Flat/src/CommonOps.cc#L38-L155
###

met_trigger_paths = {}
met_trigger_paths["2016"] = ["HLT_PFMET170_NoiseCleaned",
                             "HLT_PFMET170_HBHECleaned",
                             "HLT_PFMET170_JetIdCleaned",
                             "HLT_PFMET170_NotCleaned",
                             #"HLT_PFMET170_HBHE_BeamHaloCleaned",
                             #"HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight",
                             #"HLT_PFMETNoMu110_NoiseCleaned_PFMHTNoMu110_IDTight",
                             #"HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight",
                             "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
                             "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
                             "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
                             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"]
met_trigger_paths["2017"] = ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
                             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                             "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
                             "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight"]
met_trigger_paths["2018"] = ["HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60",
                             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
                             "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
                             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                             "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
                             "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight"]

singleele_trigger_paths = {}
singleele_trigger_paths["2016"] = ["HLT_Ele27_WP85_Gsf",
                                   "HLT_Ele27_WPLoose_Gsf",
                                   "HLT_Ele105_CaloIdVT_GsfTrkIdT",
                                   "HLT_Ele27_WPTight_Gsf",
                                   "HLT_Ele30_WPTight_Gsf",
                                   "HLT_Ele27_eta2p1_WPTight_Gsf",
                                   "HLT_Ele32_eta2p1_WPTight_Gsf",
                                   "HLT_Ele35_WPLoose_Gsf",
                                   "HLT_ECALHT800"]
singleele_trigger_paths["2017"] = ["HLT_Ele35_WPTight_Gsf",
                                   "HLT_Ele38_WPTight_Gsf",
                                   "HLT_Ele40_WPTight_Gsf"]
singleele_trigger_paths["2018"] = ["HLT_Ele35_WPTight_Gsf",
                                   "HLT_Ele38_WPTight_Gsf",
                                   "HLT_Ele40_WPTight_Gsf"]

singlepho_trigger_paths = {}
singlepho_trigger_paths["2016"] = ["HLT_Photon175",
                                   "HLT_Photon200",
                                   "HLT_Photon165_HE10",
                                   "HLT_Photon36_R9Id90_HE10_IsoM",
                                   "HLT_Photon50_R9Id90_HE10_IsoM",
                                   "HLT_Photon75_R9Id90_HE10_IsoM",
                                   "HLT_Photon90_R9Id90_HE10_IsoM",
                                   "HLT_Photon120_R9Id90_HE10_IsoM",
                                   "HLT_Photon165_R9Id90_HE10_IsoM",
                                   "HLT_Photon300_NoHE",
                                   "HLT_ECALHT800",
                                   "HLT_CaloJet500_NoJetID"]
singlepho_trigger_paths["2017"] = singlepho_trigger_paths["2016"]
singlepho_trigger_paths["2018"] = singlepho_trigger_paths["2016"]
