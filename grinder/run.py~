#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import numexpr

import uproot, uproot_methods
import numpy as np
from fnal_column_analysis_tools import hist
from saiyan import Builder

with open("data/coffeabeans2016.json") as fin:
    datadef = json.load(fin)

# [pb]
dataset_xs = {k: v['xs'] for k,v in datadef.items()}
lumi = 1000.  # [1/pb]

dataset = hist.Cat("dataset", "Primary dataset")


jetpt = hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", [450, 500, 550, 600, 675, 800, 1000])
jetpt_coarse = hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", [450, 800])
jetmass = hist.Bin("AK8Puppijet0_msd", "Jet $m_{sd}$", 23, 40, 201)
jetmass_coarse = hist.Bin("AK8Puppijet0_msd", "Jet $m_{sd}$", [40, 100, 140, 200])
jetrho = hist.Bin("jetrho", r"Jet $\rho$", 13, -6, -2.1)
doubleb = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", 20, 0., 1)
doublec = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", 20, 0., 1.)
doublecvb = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", 20, 0., 1.)
doubleb_coarse = [1., 0.93, 0.92, 0.89, 0.85, 0.7]
doubleb_coarse = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", doubleb_coarse[::-1])
doublec_coarse = [0.87, 0.84, 0.83, 0.79, 0.69, 0.58]
doublec_coarse = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", doublec_coarse[::-1])
doublecvb_coarse = [0.93, 0.91, 0.86, 0.76, 0.6, 0.17, 0.12]
doublecvb_coarse = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", doublecvb_coarse[::-1])
n2ddt_coarse = hist.Bin("AK8Puppijet0_N2sdb1_ddt", "N2 DDT", [0.])


hists = {}
hists['sumw'] = hist.Hist("sumw", dataset, hist.Bin("sumw", "Weight value", [0.]))
hists['hjetpt'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
hists['hjetpt_SR'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
#hists['htagtensor'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, n2ddt_coarse, jetmass_coarse, doubleb, doublec, doublecvb, dtype='f')
hists['hsculpt'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
hists['hsculpt_SR'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')

hists['pfmet_nminus1_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("pfmet", r"PF $p_{T}^{miss}$", 40, 0, 200))
hists['opposite_ak8_n3sdb1_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0.5, 3))
hists['opposite_ak8_tau32_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 1))
hists['opposite_ak8_msd_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_msd", r"Jet $\m_{sd}$", 40, 50, 200))
hists['opposite_ak4_leadingDeepCSV_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak4_leadingDeepCSV", "Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1))

for h in hists.values(): h.clear()
nevents = defaultdict(lambda: 0.)

def clean(val, default):
    val[np.isnan(val)|(val==-999.)] = default
    return val

def processfile(dataset, file):
    # Many 'invalid value encountered in ...' due to pt and msd sometimes being zero
    # This will just fill some NaN bins in the histogram, which is fine
    tree = uproot.open(file)["Events"]
    e = Initialize({'pt':tree.array("Electron_pt"),
                    'eta':tree.array("Electron_eta"),
                    'phi':tree.array("Electron_phi"),
                    'mass':tree.array("Electron_mass"),
                    'iso':tree.array('Electron_pfRelIso03_all'),
                    'dxy':tree.array('Electron_dxy'),
                    'dz':tree.array('Electron_dz'),
                    'loose_id':tree.array('Electron_mvaSpring16GP_WP90'),
                    'tight_id':tree.array('Electron_mvaSpring16GP_WP80')})
    e['isloose'] = (e.pt>7)&(abs(e.eta)<2.4)&(abs(e.dxy)<0.05)&(abs(e.dz)<0.2)&(e.iso<0.4)&(e.loose_id)
    e['istight'] = (e.pt>30)&(abs(e.eta)<2.4)&(abs(e.dxy)<0.05)&(abs(e.dz)<0.2)&(e.tight_id)&(e.iso<0.06)
    e_loose = e[e.isloose]
    e_tight = e[e.istight]
    e_ntot = e.counts
    e_nloose = e_loose.counts
    e_ntight = e[e.istight].counts

    mu = Initialize({'pt':tree.array("Muon_pt"),
                     'eta':tree.array("Muon_eta"),
                     'phi':tree.array("Muon_phi"),
                     'mass':tree.array("Muon_mass"),
                     'iso':tree.array('Muon_pfRelIso04_all'),
                     'dxy':tree.array('Muon_dxy'),
                     'dz':tree.array('Muon_dz')})
    mu['isloose']=(mu.counts>0)&(mu.pt>5)&(abs(mu.eta)<2.4)&(abs(mu.dxy)<0.5)&(abs(mu.dz)<1.0)&(mu.iso<0.4)
    m_loose=mu[mu.isloose]
    mu_ntot = mu.counts
    mu_nloose = m_loose.counts

    tau = Initialize({'pt':tree.array('Tau_pt'),
                      'eta':tree.array('Tau_eta'),
                      'phi':tree.array('Tau_phi'),
                      'mass':tree.array('Tau_mass'),
                      'decayMode':tree.array('Tau_idDecayMode'),
                      'decayModeNew':tree.array('Tau_idDecayModeNewDMs'),
                      'id':tree.array('Tau_idMVAnew')})
    tau['isloose']=(tau.counts>0)&(tau.pt>18)&(abs(tau.eta)<2.3)&(tau.decayMode)&((tau.id&2)!=0)
    tau_loose=tau[tau.isloose]
    tau_ntot=tau.counts
    tau_nloose=tau_loose.counts

    pho = Initialize({'pt':tree.array('Photon_pt'),
                      'eta':tree.array('Photon_eta'),
                      'phi':tree.array('Photon_phi'),
                      'mass':tree.array('Photon_mass')})
    pho['isloose'] = (pho.counts>0)&(pho.pt>15)*(abs(pho.eta)<2.5)
    pho_loose=pho[pho.isloose]
    pho_ntot=pho.counts
    pho_nloose=pho_loose.counts

    fj = Initialize({'pt':tree.array('AK15Puppi_pt'),
                     'eta':tree.array('AK15Puppi_eta'),
                     'phi':tree.array('AK15Puppi_phi'),
                     'mass':tree.array('AK15Puppi_mass'),
                     'jetId':tree.array('AK15Puppi_jetId')})
    fj['isgood'] = (fj.pt > 200)&(abs(fj.eta)<2.4)&(fj.jetId > 0)
    fj['isclean'] =~fj.match(pho,1.5)&~fj.match(mu,1.5)&~fj.match(e,1.5)&fj.isgood
    fj_good=fj[fj.isgood]
    fj_clean=fj[fj.isclean]
    fj_ntot=fj.counts
    fj_ngood=fj_good.counts
    fj_nclean=fj_clean.counts

    j = Initialize({'pt':tree.array('Jet_pt'),
                    'eta':tree.array('Jet_eta'),
                    'phi':tree.array('Jet_phi'),
                    'mass':tree.array('Jet_mass'),
                    'id':tree.array('Jet_jetId')})
    j['isgood'] = (j.pt>25)&(abs(j.eta)<4.5)&((j.id&2)!=0)
    j['isclean'] = ~j.match(e,0.4)&~j.match(mu,0.4)&~j.match(pho,0.4)&j.isgood
    j['isiso'] =  ~(j.match(fj,1.5))&j.isclean
    j_good = j[j.isgood]
    j_clean = j[j.isclean]
    j_iso = j[j.isiso]
    j_ntot=j.counts
    j_ngood=j_good.counts
    j_nclean=j_clean.counts
    j_niso=j_iso.counts

    met = Initialize({'pt':tree.array("MET_pt"),
                      'eta':0,
                      'phi':tree.array("MET_phi"),
                      'mass':0})

    diele = e_loose.distincts().i0+e_loose.distincts().i1
    dimu = m_loose.distincts().i0+m_loose.distincts().i1
    uwm = met+m_loose
    uwe = met+e_loose
    uzmm = met+dimu
    uzee = met+diele
    upho = met+pho_loose

    skinny = (j_nclean>0)
    loose = (fj_nclean>0)
    zeroL = (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)
    oneM = (e_nloose==0)&(mu_nloose==1)&(tau_nloose==0)&(pho_nloose==0)
    oneE = (e_nloose==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)
    twoM = (e_nloose==0)&(mu_nloose==2)&(tau_nloose==0)&(pho_nloose==0)
    twoE = (e_nloose==2)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0)
    oneA = (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==1)
    
    sr = (j_clean[zeroL&skinny,0].pt>100)&(met[zeroL&skinny].pt>200)&(met[zeroL&skinny].delta_phi(met[zeroL&skinny].closest(j_clean[zeroL&skinny]))>0.5)

    hout = {}
    for k in hists.keys():
        h = hists[k].copy(content=False)
        if k == 'sumw':
            h.fill(dataset=dataset, sumw=genW)
        else:
            h.fill(dataset=dataset, **arrays, weight=weight)
        hout[k] = h
    return dataset, tree.numentries, hout


nworkers = 15
fileslice = slice(None)
with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = set()
    for dataset, info in datadef.items():
        futures.update(executor.submit(processfile, dataset, file) for file in info['files'][fileslice])
    try:
        total = len(futures)
        processed = 0
        while len(futures) > 0:
            finished = set(job for job in futures if job.done())
            for job in finished:
                dataset, nentries, hout = job.result()
                nevents[dataset] += nentries
                for k in hout.keys():
                    hists[k] += hout[k]
                processed += 1
                print("Processing: done with % 4d / % 4d files" % (processed, total))
            futures -= finished
        del finished
    except KeyboardInterrupt:
        print("Ok quitter")
        for job in futures: job.cancel()
    except:
        for job in futures: job.cancel()
        raise


sumw = hists.pop('sumw')
scale = {}
for ds in nevents.keys():
    ds_sumw = sumw.values(overflow='all')[(ds,)]
    print(ds, nevents[ds], ds_sumw)
    scale[ds] = lumi*dataset_xs[ds] / (ds_sumw[1]-ds_sumw[0])

for h in hists.values(): h.scale(scale, axis="dataset")

dt = time.time() - tstart
print("%.2f us*cpu/event" % (1e6*dt*nworkers/sum(nevents.values()), ))
nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in hists.values())
nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in hists.values())
print("Processed %.1fM events" % (sum(nevents.values())/1e6, ))
print("Filled %.1fM bins" % (nbins/1e6, ))
print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

# Pickle is not very fast or memory efficient, will be replaced by something better soon
with gzip.open("hists.pkl.gz", "wb") as fout:
    pickle.dump(hists, fout)

