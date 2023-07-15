#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
import pprint
import numpy as np
import math
import awkward
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea.arrays import Initialize
from coffea import hist, processor
from coffea.util import load, save
from optparse import OptionParser
from uproot_methods import TVector2Array, TLorentzVectorArray

class AnalysisProcessor(processor.ProcessorABC):

    lumis = { #Values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable                                                      
        '2016': 36.31,
        '2017': 41.53,
        '2018': 59.74
    }

    met_filter_flags = {
        '2016': ['goodVertices',
                 'globalSuperTightHalo2016Filter',
                 'HBHENoiseFilter',
                 'HBHENoiseIsoFilter',
                 'EcalDeadCellTriggerPrimitiveFilter',
                 'BadPFMuonFilter'
             ],
        '2017': ['goodVertices',
                 'globalSuperTightHalo2016Filter',
                 'HBHENoiseFilter',
                 'HBHENoiseIsoFilter',
                 'EcalDeadCellTriggerPrimitiveFilter',
                 'BadPFMuonFilter',
                 'ecalBadCalibFilterV2'
             ],
        '2018': ['goodVertices',
                 'globalSuperTightHalo2016Filter',
                 'HBHENoiseFilter',
                 'HBHENoiseIsoFilter',
                 'EcalDeadCellTriggerPrimitiveFilter',
                 'BadPFMuonFilter',
                 'ecalBadCalibFilterV2'
             ]
    }

    def __init__(self, year, xsec, corrections, ids, common):

        self._year = year
        self._lumi = 1000.*float(AnalysisProcessor.lumis[year])
        self._xsec = xsec

        self._gentype_map = {
            'bb':       0,
            'b':        1,
            'cc' :      2,
            'c':        3,
            'other':    4,
            'hs':       5,
        }

        self._ZHbbvsQCDwp = {
            '2016': 0.53,
            '2017': 0.61,
            '2018': 0.65
        }

        self._btagmu_triggers = {
            '2016': [
                'BTagMu_AK4Jet300_Mu5',
                'BTagMu_AK8Jet300_Mu5',
                'BTagMu_AK4DiJet170_Mu5'
                ],
            '2017': [
                'BTagMu_AK4Jet300_Mu5',
                'BTagMu_AK8Jet300_Mu5',
                'BTagMu_AK4DiJet170_Mu5'
                ],
            '2018': [
                'BTagMu_AK4Jet300_Mu5',
                'BTagMu_AK8Jet300_Mu5',
                'BTagMu_AK4DiJet170_Mu5'
                ]
        }

        self._corrections = corrections
        self._ids = ids
        self._common = common

        self._accumulator = processor.dict_accumulator({
            'sumw': hist.Hist(
                'sumw',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('sumw', 'Weight value', [0.])
            ),
            'svtemplate': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6]),
                hist.Bin('svmass','Leading Secondary Vertices (SV) mass',[-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0.,   0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9, 2.,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,  3.,   3.1,  3.2]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
            'ZHbbvsQCD': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD',15,0,1),
                hist.Bin('fjmass','AK15 Jet Mass',52,40,300),
            ),
            'fj1pt': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6]),
                hist.Bin('fj1pt','AK15 Leading SoftDrop Jet Pt',[160.0, 250.0, 280.0, 310.0, 340.0, 370.0, 400.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]),
                hist.Bin('fjmass','AK15 Jet Mass',52,40,300),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD', [0, self._ZHbbvsQCDwp[self._year], 1])
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    def process(self, events):

        dataset = events.metadata['dataset']

        isData = 'genWeight' not in events.columns
        selection = processor.PackedSelection()
        hout = self.accumulator.identity()

        ###
        #Getting ids from .coffea files
        ###

        get_msd_weight  = self._corrections['get_msd_weight']
        get_pu_weight   = self._corrections['get_pu_weight'][self._year]  
        isLooseMuon     = self._ids['isLooseMuon']
        isTightMuon     = self._ids['isTightMuon']
        isGoodFatJet    = self._ids['isGoodFatJet']
        isHEMJet        = self._ids['isHEMJet']  

        match = self._common['match']

        ###
        #Initialize physics objects
        ###

        mu = events.Muon
        leading_mu = mu[mu.pt.argmax()]

        j = events.Jet
        j['isHEM'] = isHEMJet(j.pt, j.eta, j.phi)
        j_HEM = j[j.isHEM.astype(np.bool)]
        j_nHEM = j_HEM.counts

        fj = events.AK15Puppi
        fj['sd'] = fj.subjets.sum()
        fj['isgood'] = isGoodFatJet(fj.sd.pt, fj.sd.eta, fj.jetId)
        fj['T'] = TVector2Array.from_polar(fj.pt, fj.phi)
        fj['msd_raw'] = (fj.subjets * (1 - fj.subjets.rawFactor)).sum().mass
        fj['msd_corr'] = fj.msd_raw * awkward.JaggedArray.fromoffsets(fj.array.offsets, np.maximum(1e-5, get_msd_weight(fj.sd.pt.flatten(),fj.sd.eta.flatten())))
        probQCD=fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
        probZHbb=fj.probZbb+fj.probHbb
        fj['ZHbbvsQCD'] = probZHbb/(probZHbb+probQCD)
        fj['tau21'] = fj.tau2/fj.tau1

        SV = events.SV

        ###
        # Calculating weights
        ###
        if not isData:

            gen = events.GenPart
            
            #####
            ###
            # Fat-jet dark H->bb matching at decay level
            ###
            #####
            
            ##### Get b's from dark Higgs decay
            bFromHs = gen[
                (abs(gen.pdgId) == 5) &
                gen.hasFlags(['fromHardProcess', 'isFirstCopy']) &
                (gen.distinctParent.pdgId == 54)
            ]
            
            ##### Mix fatjet subjets and gen level b's
            ##### axis=1 option to remove boundaries between fat-jets
            jetgenb = fj.subjets.flatten(axis=1).cross(bFromHs, nested=True)
            
            ##### Match subjets to b's
            mask = (bFromHs.counts>0) & ((jetgenb.i0.delta_r(jetgenb.i1) < 0.4).sum() == 1)
    
            ##### Three steps to match the jaggedness of the mask array to the fj.subjets array #####
            ##### Using the offset function to copy contents not the type of the array #####
            step1 = fj.subjets.flatten()
            step2 = awkward.JaggedArray.fromoffsets(step1.offsets, mask.content)
            step2 = step2.pad(1).fillna(0) ##### Fill None for empty arrays and convert None to False
            step3 = awkward.JaggedArray.fromoffsets(fj.subjets.offsets, step2)
    
            ##### fatjet with two subjets matched to dark Higgs b's
            fj['isHsbb'] = (step3.sum() == 2)

            
            #####
            ###
            # Fat-jet matching to one or two b's
            ###
            #####

            ##### Get gen b's
            gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])

            ##### Mix fatjet subjets and gen level b's
            ##### axis=1 option to remove boundaries between fat-jets
            jetgenb = fj.subjets.flatten(axis=1).cross(gen['isb'], nested=True)

            ##### Match subjets to b's
            mask = (gen[gen.isb].counts>0) & ((jetgenb.i0.delta_r(jetgenb.i1) < 0.4).sum() == 1)
            
            ##### Three steps to match the jaggedness of the mask array to the fj.subjets array #####
            ##### Using the offset function to copy contents not the type of the array #####
            step1 = fj.subjets.flatten()
            step2 = awkward.JaggedArray.fromoffsets(step1.offsets, mask.content)
            step2 = step2.pad(1).fillna(0) ##### Fill None for empty arrays and convert None to False
            step3 = awkward.JaggedArray.fromoffsets(fj.subjets.offsets, step2)
    
            ##### fatjet with two subjets matched to dark Higgs b's
            fj['isb']  = (step3.sum() == 1)
            fj['isbb']  = (step3.sum() == 2)


            #####
            ###
            # Fat-jet matching to one or two c's
            ###
            #####

            ##### Get gen b's
            gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])

            ##### Mix fatjet subjets and gen level c's
            ##### axis=1 option to remove boundaries between fat-jets
            jetgenc = fj.subjets.flatten(axis=1).cross(gen[gen.isc], nested=True)

            ##### Match subjets to b's
            mask = (gen[gen.isc].counts>0) & ((jetgenc.i0.delta_r(jetgenc.i1) < 0.4).sum() == 1)
            
            ##### Three steps to match the jaggedness of the mask array to the fj.subjets array #####
            ##### Using the offset function to copy contents not the type of the array #####
            step1 = fj.subjets.flatten()
            step2 = awkward.JaggedArray.fromoffsets(step1.offsets, mask.content)
            step2 = step2.pad(1).fillna(0) ##### Fill None for empty arrays and convert None to False
            step3 = awkward.JaggedArray.fromoffsets(fj.subjets.offsets, step2)
    
            ##### fatjet with two subjets matched to dark Higgs b's
            fj['isc']  = (step3.sum() == 1)
            fj['iscc']  = (step3.sum() == 2)
            

            #####
            ###
            # Calculate PU weight and systematic variations
            ###
            #####

            pu = get_pu_weight(events.Pileup.nTrueInt)

        ##### axis=1 option to remove boundaries between fat-jets #####
        ##### copy (match jaggedness and shape of array) the contents of crossed array into the fat-jet subjets #####
        ##### we're not use copy since it keeps the original array type #####
        ##### fj.subjets is a TLorentzVectorArray #####
        mu = mu[mu.isGlobal] ## Use a global muon for QCD events
        jetmu = fj.subjets.flatten(axis=1).cross(mu, nested=True)
        mask = (mu.counts>0) & ((jetmu.i0.delta_r(jetmu.i1) < 0.4) & ((jetmu.i1.pt/jetmu.i0.pt) < 0.7) & (jetmu.i1.pt > 7)).sum() == 1

        ##### Three steps to match the jaggedness of the mask array to the fj.subjets array #####
        ##### Using the offset function to copy contents not the type of the array #####
        step1 = fj.subjets.flatten()
        step2 = awkward.JaggedArray.fromoffsets(step1.offsets, mask.content)
        step2 = step2.pad(1).fillna(0) ##### Fill None for empty arrays and convert None to False
        step3 = awkward.JaggedArray.fromoffsets(fj.subjets.offsets, step2)

        ##### fatjet with two subjets matched with muons
        fj['withmu'] = step3.sum() == 2

        ###
        # Selections
        ###

        #### trigger selection ####
        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._btagmu_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('btagmu_triggers', triggers)

        #### MET filters ####
        met_filters =  np.ones(events.size, dtype=np.bool)
        if isData:
            met_filters = met_filters & events.Flag['eeBadScFilter'] #this filter is recommended for data only
        for flag in AnalysisProcessor.met_filter_flags[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters',met_filters)

        #### ak15 jet selection ####
        leading_fj = fj[fj.sd.pt.argmax()]
        leading_fj = leading_fj[leading_fj.isgood.astype(np.bool)]
        #leading_fj = leading_fj[(leading_fj.msd_corr > 40)]
        #leading_fj = leading_fj[leading_fj.withmu.astype(np.bool)]

        #### SV selection for matched with leading ak15 jet ####
        SV['ismatched'] = match(SV, leading_fj, 1.5)
        leading_SV = SV[SV.dxySig.argmax()]
        leading_SV = leading_SV[leading_SV.ismatched.astype(np.bool)]

        #fj_good = fj[fj.isgood.astype(np.bool)]
        #fj_withmu = fj_good[fj_good.withmu.astype(np.bool)]
        #fj_nwithmu = fj_withmu.counts

        noHEMj = np.ones(events.size, dtype=np.bool)
        if self._year=='2018': noHEMj = (j_nHEM==0)

        selection.add('noHEMj', noHEMj)
        selection.add('fj_pt', (leading_fj.sd.pt.max() > 250) )
        selection.add('fj_mass', (leading_fj.msd_corr.sum() > 40) ) ## optionally also <130
        selection.add('withmu', leading_fj.withmu.sum().astype(np.bool))
        #selection.add('fj_tau21', (leading_fj.tau21.sum() < 0.3) )

        isFilled = False
        if isData:
            if not isFilled:
                hout['sumw'].fill(dataset=dataset, sumw=1, weight=1)
                isFilled=True

            cut = selection.all(*selection.names)
            ##### template for bb SF #####
            hout['svtemplate'].fill(dataset=dataset,
                                    gentype=np.zeros(events.size, dtype=np.int),
                                    svmass=np.log(leading_SV.mass.sum()),
                                    ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                    weight=np.ones(events.size)*cut)
            hout['ZHbbvsQCD'].fill(dataset=dataset,
                                   gentype=np.zeros(events.size, dtype=np.int),
                                   ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                   fjmass=leading_fj.msd_corr.sum(),
                                   weight=np.ones(events.size)*cut)
            hout['fj1pt'].fill(dataset=dataset,
                                   gentype=np.zeros(events.size, dtype=np.int),
                                   fj1pt=leading_fj.sd.pt.sum(),
                                   fjmass=leading_fj.msd_corr.sum(),
                                   ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                   weight=np.ones(events.size)*cut)

        else:
            weights = processor.Weights(len(events))
            if 'L1PreFiringWeight' in events.columns: weights.add('prefiring',events.L1PreFiringWeight.Nom)
            weights.add('genw',events.genWeight)
            weights.add('pileup',pu)

            wgentype = {
                'bb' : (leading_fj.isbb).sum(),
                'b'  : ( ~leading_fj.isbb & leading_fj.isb ).sum(),
                'cc' : ( ~leading_fj.isbb & ~leading_fj.isb & leading_fj.iscc ).sum(),
                'c'  : ( ~leading_fj.isbb & ~leading_fj.isb & ~leading_fj.iscc & leading_fj.isc ).sum(),
                'other' : ( ~leading_fj.isbb & ~leading_fj.isb & ~leading_fj.iscc & ~leading_fj.isc ).sum(),
                'hs' : (leading_fj.isHsbb).sum(),
            }
            vgentype=np.zeros(events.size, dtype=np.int)
            for gentype in self._gentype_map.keys():
                vgentype += self._gentype_map[gentype]*wgentype[gentype]

            if not isFilled:
                hout['sumw'].fill(dataset=dataset, sumw=1, weight=events.genWeight.sum())
                isFilled=True

            cut = selection.all(*selection.names)
            if 'QCD' in dataset:
                ##### template for bb SF #####
                hout['svtemplate'].fill(dataset=dataset,
                                        gentype=vgentype,
                                        svmass=np.log(leading_SV.mass.sum()),
                                        ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                        weight=weights.weight()*cut)
                hout['ZHbbvsQCD'].fill(dataset=dataset,
                                       gentype=vgentype,
                                       ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                       fjmass=leading_fj.msd_corr.sum(),
                                       weight=weights.weight()*cut)
                hout['fj1pt'].fill(dataset=dataset,
                                       gentype=vgentype,
                                       fj1pt=leading_fj.sd.pt.sum(),
                                       fjmass=leading_fj.msd_corr.sum(),
                                       ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                       weight=weights.weight()*cut)
            else:
                ##### template for bb SF #####
                hout['svtemplate'].fill(dataset=dataset,
                                        gentype=vgentype,
                                        svmass=np.log(leading_SV.mass.sum()),
                                        ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                        weight=weights.weight())
                hout['ZHbbvsQCD'].fill(dataset=dataset,
                                       gentype=vgentype,
                                       ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                       fjmass=leading_fj.msd_corr.sum(),
                                       weight=weights.weight())
                hout['fj1pt'].fill(dataset=dataset,
                                       gentype=vgentype,
                                       fj1pt=leading_fj.sd.pt.sum(),
                                       fjmass=leading_fj.msd_corr.sum(),
                                       ZHbbvsQCD=leading_fj.ZHbbvsQCD.sum(),
                                       weight=weights.weight())

        return hout

    def postprocess(self, accumulator):
        scale = {}
        for d in accumulator['sumw'].identifiers('dataset'):
            print('Scaling:',d.name)
            dataset = d.name
            if '--' in dataset: dataset = dataset.split('--')[1]
            print('Cross section:',self._xsec[dataset])
            if self._xsec[dataset]!= -1: scale[d.name] = self._lumi*self._xsec[dataset]
            else: scale[d.name] = 1

        for histname, h in accumulator.items():
            if histname == 'sumw': continue
            if isinstance(h, hist.Hist):
                h.scale(scale, axis='dataset')

        return accumulator

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
    (options, args) = parser.parse_args()


    with open('metadata/'+options.metadata+'.json') as fin:
        samplefiles = json.load(fin)
        xsec = {k: v['xs'] for k,v in samplefiles.items()}

    corrections = load('data/corrections.coffea')
    ids         = load('data/ids.coffea')
    common      = load('data/common.coffea')

    processor_instance=AnalysisProcessor(year=options.year,
                                         xsec=xsec,
                                         corrections=corrections,
                                         ids=ids,
                                         common=common)

    save(processor_instance, 'data/doublebsf'+options.metadata+'.processor')
