#!/usr/bin/env python

###
#
# Python script that "compiles" corrections, ids, ecc. 
# In this context compiling means generating pickle (with .coffea) extension files.
#
###

import json
import gzip
import uproot
import numexpr
import numpy as np
from coffea import hist, lookup_tools
from coffea.util import load, save
from coffea.hist import plot

from utils.triggers import met_trigger_paths, singleele_trigger_paths, singlepho_trigger_paths
triggers = {}
triggers['met_trigger_paths']       = met_trigger_paths
triggers['singleele_trigger_paths'] = singleele_trigger_paths
triggers['singlepho_trigger_paths'] = singlepho_trigger_paths
save(triggers, 'triggers.coffea')

from utils.corrections import get_ttbar_weight, get_nlo_weight, get_adhoc_weight, get_pu_weight
from utils.corrections import get_met_trig_weight, get_met_zmm_trig_weight, get_ele_trig_weight, get_pho_trig_weight
from utils.corrections import get_ecal_bad_calib
corrections = {}
corrections['get_ttbar_weight']        = get_ttbar_weight
corrections['get_nlo_weight']          = get_nlo_weight
corrections['get_adhoc_weight']        = get_adhoc_weight
corrections['get_pu_weight']           = get_pu_weight
corrections['get_met_trig_weight']     = get_met_trig_weight
corrections['get_met_zmm_trig_weight'] = get_met_zmm_trig_weight
corrections['get_ele_trig_weight']     = get_ele_trig_weight
corrections['get_pho_trig_weight']     = get_pho_trig_weight
corrections['get_ecal_bad_calib']      = get_ecal_bad_calib
save(corrections, 'corrections.coffea')

from utils.ids import isLooseElectron, isTightElectron
from utils.ids import isLooseMuon, isTightMuon
from utils.ids import isLooseTau
from utils.ids import isLoosePhoton, isTightPhoton
from utils.ids import isGoodJet, isGoodFatJet, isHEMJet
ids = {}
ids['isLooseElectron'] = isLooseElectron
ids['isTightElectron'] = isTightElectron
ids['isLooseMuon']     = isLooseMuon
ids['isTightMuon']     = isTightMuon
ids['isLooseTau']      = isLooseTau
ids['isLoosePhoton']   = isLoosePhoton
ids['isTightPhoton']   = isTightPhoton
ids['isGoodJet']       = isGoodJet
ids['isGoodFatJet']    = isGoodFatJet
ids['isHEMJet']        = isHEMJet
save(ids, 'ids.coffea')

from utils.metfilters import met_filter_flags
metfilters = {}
metfilters['met_filter_flags'] = met_filter_flags
save(metfilters, 'metfilters.coffea')

