###
# From https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
###
from coffea.util import save

met_filter_flags = {}
met_filter_flags["2016"] = ["Flag_goodVertices",
                            "Flag_globalSuperTightHalo2016Filter",
                            "Flag_HBHENoiseFilter",
                            "Flag_HBHENoiseIsoFilter",
                            "Flag_EcalDeadCellTriggerPrimitiveFilter",
                            ]

met_filter_flags["2017"] = ["Flag_goodVertices",
                            "Flag_globalSuperTightHalo2016Filter",
                            "Flag_HBHENoiseFilter",
                            "Flag_HBHENoiseIsoFilter",
                            "Flag_EcalDeadCellTriggerPrimitiveFilter",
                            "Flag_BadPFMuonFilter",
                            ]

met_filter_flags["2018"] = ["Flag_goodVertices",
                            "Flag_globalSuperTightHalo2016Filter",
                            "Flag_HBHENoiseFilter",
                            "Flag_HBHENoiseIsoFilter",
                            "Flag_EcalDeadCellTriggerPrimitiveFilter",
                            "Flag_BadPFMuonFilter",
                            ]
metfilters = {}
metfilters['met_filter_flags'] = met_filter_flags
save(metfilters, 'metfilters.coffea')
