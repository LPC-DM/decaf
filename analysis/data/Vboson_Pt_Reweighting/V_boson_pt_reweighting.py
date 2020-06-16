# imports
from __future__ import print_function
import ROOT
import sys
from array import array
from DataFormats.FWLite import Events, Handle
from math import *
from sample_info import sample_dict

# ROOT.gSystem.Load('libGenVector')

# not needed at the moment
def FindAllMothers(particle):
    mother_ids = []
    print ("particle id ", particle.pdgId())
    print ("# mothers ", particle.numberOfMothers())
    for i in range(particle.numberOfMothers()):
        print ("mother id ", particle.mother(i).pdgId())
        mother_ids.append(particle.mother(i).pdgId())
        next_mothers_ids = FindAllMothers(particle.mother(i))
        for next_mother_id in next_mothers_ids:
            mother_ids.append(next_mother_id)
    return mother_ids


# also use muR and muF variations
scales = {
    "Weight_scale_variation_muR_1p0_muF_1p0": 1,
    "Weight_scale_variation_muR_1p0_muF_2p0": 2,
    "Weight_scale_variation_muR_1p0_muF_0p5": 3,
    "Weight_scale_variation_muR_2p0_muF_1p0": 4,
    "Weight_scale_variation_muR_2p0_muF_2p0": 5,
    "Weight_scale_variation_muR_2p0_muF_0p5": 6,
    "Weight_scale_variation_muR_0p5_muF_1p0": 7,
    "Weight_scale_variation_muR_0p5_muF_2p0": 8,
    "Weight_scale_variation_muR_0p5_muF_0p5": 9,
}

# inputs
era = str(sys.argv[1])
boson = str(sys.argv[2])
postfix = str(sys.argv[3])
filenames = sys.argv[4:]

# product labels and handles
handlePruned = Handle("std::vector<reco::GenParticle>")
handlePacked = Handle("std::vector<pat::PackedGenParticle>")
eventinfo = Handle("GenEventInfoProduct")
lheinfo = Handle("LHEEventProduct")
labelPruned = "prunedGenParticles"
labelPacked = "packedGenParticles"
labelWeight = "generator"
labelLHE = "externalLHEProducer"

# binning according to https://arxiv.org/pdf/1705.04664.pdf
binning = [
    30,
    40,
    50,
    60,
    70,
    80,
    90,
    100,
    110,
    120,
    130,
    140,
    150,
    200,
    250,
    300,
    350,
    400,
    450,
    500,
    550,
    600,
    650,
    700,
    750,
    800,
    850,
    900,
    950,
    1000,
    1100,
    1200,
    1300,
    1400,
    1600,
    1800,
    2000,
    2200,
    2400,
    2600,
    2800,
    3000,
    6500,
]

# nominal histogram
v_boson_pt_hist = ROOT.TH1D(boson + "_boson_pt", boson + "_boson_pt", len(binning) - 1, array("d", binning))
v_boson_pt_hist.Sumw2()
file_ = ROOT.TFile(boson + "_boson_pt_" + era + "_" + postfix + ".root", "RECREATE")

# list of histograms for scale variations
v_boson_pt_hists = {}
for scale in scales:
    v_boson_pt_hists[scale] = ROOT.TH1D(boson + "_boson_pt_" + scale, boson + "_boson_pt_" + scale, len(binning) - 1, array("d", binning))
    v_boson_pt_hists[scale].Sumw2()

count = 0

# loop over files
for filename in filenames:
    print(filename)
    # cross section weights
    weight_xs = 1.0
    if boson == "Zvv":
        print ("looking at Zvv events")
    elif boson == "Zll":
        print ("looking at Zll events")
    elif boson == "W":
        print ("looking at W events")
    elif boson == "G":
        print ("looking at Photon events")
    else:
        print ("only W or Z boson or Photon allowed")
        exit()
    # information to calculate cross section weight
    subdict = sample_dict.get(era, None).get(boson, None)
    for key in subdict:
        if key in filename.lower():
            subsubdict = subdict.get(key, None)
            weight_xs = subsubdict.get("sigma", None) / (subsubdict.get("X", None) * subsubdict.get("N_gen", None))
    weight_xs *= 1000.0
    print ("weight_xs = ", weight_xs)

    # loop over events
    events = Events(filename)
    for event in events:
        count += 1
        if count % 10000 == 0:
            print (count)
        event.getByLabel(labelPruned, handlePruned)
        event.getByLabel(labelWeight, eventinfo)
        event.getByLabel(labelLHE, lheinfo)
        # get the products (prunedGenParticles collection, GenEventInfoProduct and LHEEventProduct)
        pruned = handlePruned.product()
        weight = eventinfo.product().weight()
        lhe_weight = lheinfo.product().originalXWGTUP()
        # list for decay products of W or Z boson
        decay_prods = []
        # list for photons either for photon+jets events or for radiated photons to add back to charged leptons
        photons = []
        # list for hadrons needed for photon isolation
        hadrons = []
        # list for isolated photons
        isolated_photons = []

        if boson == "Zvv":
            for p in pruned:
                if p.isPromptFinalState() and (abs(p.pdgId()) == 12 or abs(p.pdgId()) == 14 or abs(p.pdgId()) == 16):
                    decay_prods.append(p)
                    # print("found neutrino")
        elif boson == "Zll":
            for p in pruned:
                # need to save stable photons to calculate dressed leptons later
                if abs(p.pdgId()) == 22 and p.status() == 1 and p.statusFlags().isPrompt():
                    photons.append(p)
                    continue
                # check for prompt final state charged leptons
                if ((abs(p.pdgId()) == 11 or abs(p.pdgId()) == 13) and p.isPromptFinalState()) or (
                    abs(p.pdgId()) == 15 and p.isPromptDecayed()
                ):
                    decay_prods.append(p)
                    # print("found charged lepton")
        elif boson == "W":
            for p in pruned:
                # need to save stable photons to calculate dressed leptons later
                if abs(p.pdgId()) == 22 and p.status() == 1 and p.statusFlags().isPrompt():
                    photons.append(p)
                    continue
                # check for prompt final state charged leptons and neutrinos
                if (
                    (abs(p.pdgId()) == 11 or abs(p.pdgId()) == 12 or abs(p.pdgId()) == 13 or abs(p.pdgId()) == 14 or abs(p.pdgId()) == 16)
                    and p.isPromptFinalState()
                ) or (abs(p.pdgId()) == 15 and p.isPromptDecayed()):
                    decay_prods.append(p)
                    # print("found neutrino/charged lepton")
        elif boson == "G":
            # need packedGenParticles collection for final state hadrons
            event.getByLabel(labelPacked, handlePacked)
            packed = handlePacked.product()
            for p in pruned:
                # need to save stable photons
                if abs(p.pdgId()) == 22 and p.status() == 1 and p.statusFlags().isPrompt():
                    photons.append(p)
            for p in packed:
                # exclude final state photons and leptons from final state gen particle collection -> final state hadrons
                if p.status() == 1 and not (abs(p.pdgId()) == 22 or (abs(p.pdgId()) >= 11 and abs(p.pdgId()) <= 16)):
                    hadrons.append(p)
        else:
            print ("only W or Z boson or Photon allowed")
            exit()

        # fail-safe: check if the number of found daughters is exactly 2 as one would expect for W/Z
        if boson == "Zvv" or boson == "Zll" or boson == "W":
            # if there are more than 2 final state leptons, e.g. from additional hadron decays, take the leading two leptons (this works in 99.999% of all cases)
            if len(decay_prods) > 2:
                decay_prods.sort(key=lambda dp: dp.pt(), reverse=True)
            # if there are less than 2 final state leptons, do not consider this event (happens in O(0.01%) of all cases so basically never)
            elif len(decay_prods) < 2:
                continue
        # fail-safe: check if there is at least one final state photon for gamma+jets events
        elif boson == "G":
            if len(photons) < 1:
                continue
        else:
            print ("only W or Z boson or Photon allowed")
            exit()

        if boson == "Zvv":
            # fail-safe: check if the daughters of the Z boson are particle and anti-particle as well as same lepton flavor
            if decay_prods[0].pdgId() + decay_prods[1].pdgId() != 0:
                continue

        elif boson == "Zll":
            # fail-safe: check if the daughters of the Z boson are particle and anti-particle as well as same lepton flavor
            if decay_prods[0].pdgId() + decay_prods[1].pdgId() != 0:
                continue
            # add radiated photons back to leptons
            for decay_prod in decay_prods:
                if abs(decay_prod.pdgId()) == 11 or abs(decay_prod.pdgId()) == 13 or abs(decay_prod.pdgId()) == 15:
                    for photon in photons:
                        if sqrt(ROOT.Math.VectorUtil.DeltaR2(decay_prod.p4(), photon.p4())) < 0.1:
                            decay_prod.setP4(decay_prod.p4() + photon.p4())

        elif boson == "W":
            # fail-safe: check if the daughters of the W boson are particle and anti-particle as well as same lepton flavor adapted for W decays
            if decay_prods[0].pdgId() * decay_prods[1].pdgId() >= 0 or abs(abs(decay_prods[0].pdgId()) - abs(decay_prods[1].pdgId())) != 1:
                continue
            # add radiated photons back to lepton
            for decay_prod in decay_prods:
                if abs(decay_prod.pdgId()) == 11 or abs(decay_prod.pdgId()) == 13 or abs(decay_prod.pdgId()) == 15:
                    for photon in photons:
                        if sqrt(ROOT.Math.VectorUtil.DeltaR2(decay_prod.p4(), photon.p4())) < 0.1:
                            decay_prod.setP4(decay_prod.p4() + photon.p4())
        # photons are more complicated on theory level
        # one has to find isolated photons using the isolation prescription in https://arxiv.org/pdf/1705.04664.pdf
        elif boson == "G":
            epsilon_0_dyn = 0.1
            n_dyn = 1
            iterations = 5.0
            for photon in photons:
                isolated = True
                R_dyn = 91.1876 / (photon.pt() * sqrt(epsilon_0_dyn))
                R_0_dyn = min(1.0, R_dyn)
                for R in [R_0_dyn / iterations * i for i in range(1, int(iterations) + 1)]:
                    isolation = 0.0
                    for hadron in hadrons:
                        if sqrt(ROOT.Math.VectorUtil.DeltaR2(hadron.p4(), photon.p4())) <= R:
                            isolation += hadron.pt()
                    if isolation > (epsilon_0_dyn * photon.pt() * pow((1 - cos(R)) / (1 - cos(R_0_dyn)), n_dyn)):
                        isolated = False
                        break
                if isolated:
                    isolated_photons.append(photon)
            # if there is more than 1 isolated photon, take the leading one
            if len(isolated_photons) > 1:
                isolated_photons.sort(key=lambda dp: dp.pt(), reverse=True)
            elif len(isolated_photons) < 1:
                # print ("no isolated photon found")
                continue

        # reconstruct vector boson from the two decay products
        if boson == "Zvv" or boson == "Zll" or boson == "W":
            v_boson = decay_prods[0].p4() + decay_prods[1].p4()
        elif boson == "G":
            v_boson = isolated_photons[0].p4()
        else:
            print ("only W or Z boson or Photon allowed")
            exit()
        v_boson_pt = v_boson.pt()
        # print (v_boson_pt)
        # fill the vector boson pt
        v_boson_pt_hist.Fill(v_boson_pt, weight * weight_xs / 1000.0)
        # fill histograms for scale variations
        for scale in scales:
            # print(scale)
            # print(scales[scale])
            for i in range(lheinfo.product().weights().size()):
                # print(lheinfo.product().weights().at(i).id)
                if int(lheinfo.product().weights().at(i).id) == int(scales[scale]):
                    scale_weight = lheinfo.product().weights().at(i).wgt
                    # print(scale_weight)
                    v_boson_pt_hists[scale].Fill(v_boson_pt, weight * weight_xs / 1000.0 * scale_weight / lhe_weight)
                    break

# write all to a file
file_.WriteTObject(v_boson_pt_hist)
for scale in scales:
    file_.WriteTObject(v_boson_pt_hists[scale])
file_.Close()
print ("finished")
