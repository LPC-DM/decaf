#!/usr/bin/env python

#plotDir = '/cms/scratch/jongho/CMSSW_10_2_13/src/decaf/analysis/datacards/darkhiggs/pulls/plots/'
#scansDir = '/data/t3home000/snarayan/store/panda/v_8026_0_5_slim/fitting/scans/' # resonant
#scansDir = '/data/t3home000/snarayan/store/panda/v_003_0/fitting/scans2_couplings2/' # coupling scans
#scansDir = '/data/t3home000/snarayan/store/panda/v_003_0/fitting/scans2/' # mass-mass and mass-g-g and resonant
#plotDir = '/cms/scratch/matteoc/CMSSW_10_2_13/src/decaf/analysis/plots/darkhiggs2017/postfit_nounc/'
#plotDir = '/cms/scratch/matteoc/CMSSW_10_2_13/src/decaf/analysis/plots/darkhiggs2018/' 
#plotDir = '/cms/scratch/matteoc/CMSSW_10_2_13/src/decaf/analysis/plots/darkhiggs2018/pulls/'
plotDir = '/dummy/'
lumi = 59.74

from os import system
system('mkdir -p '+plotDir)
