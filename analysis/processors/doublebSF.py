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
        '2016': 35.92,
        '2017': 41.53,
        '2018': 59.74
    }


    def __init__(self, year, xsec, corrections, ids, common):

        self._year = year
        self._lumi = 1000.*float(AnalysisProcessor.lumis[year])
        self._xsec = xsec

        self._gentype_map = {
            'bb':       1,
            'bc':       2,
            'b':        3,
            'cc' :      4,
            'c':        5,
            'other':    6
        }

        self._ZHbbvsQCDwp = {
            '2016': 0.53,
            '2017': 0.61,
            '2018': 0.65
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
            'ZHbbvsQCD': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('gentype', 'Gen Type', [0, 1, 2, 3, 4, 5, 6]),
                hist.Bin('ZHbbvsQCD','ZHbbvsQCD',15,0,1)
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
        isLooseMuon     = self._ids['isLooseMuon']
        isTightMuon     = self._ids['isTightMuon']
        isGoodFatJet    = self._ids['isGoodFatJet']

        match = self._common['match']

        ###
        #Initialize physics objects
        ###

        mu = events.Muon
        mu['isloose'] = isLooseMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.looseId,self._year)
        mu['istight'] = isTightMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.tightId,self._year)
        mu['T'] = TVector2Array.from_polar(mu.pt, mu.phi)
        mu_loose=mu[mu.isloose.astype(np.bool)]
        mu_tight=mu[mu.istight.astype(np.bool)]
        mu_ntot = mu.counts
        mu_nloose = mu_loose.counts
        mu_ntight = mu_tight.counts
        leading_mu = mu[mu.pt.argmax()]
        leading_mu = leading_mu[leading_mu.istight.astype(np.bool)]

        fj = events.AK15Puppi
        fj['sd'] = fj.subjets.sum()
        #fj['isclean'] =~match(fj.sd,pho_loose,1.5)&~match(fj.sd,mu_loose,1.5)&~match(fj.sd,e_loose,1.5)
        fj['isclean'] =match(fj.sd,mu_loose,0.4)
        fj['isgood'] = isGoodFatJet(fj.sd.pt, fj.sd.eta, fj.jetId)
        fj['T'] = TVector2Array.from_polar(fj.pt, fj.phi)
        fj['msd_raw'] = (fj.subjets * (1 - fj.subjets.rawFactor)).sum().mass
        fj['msd_corr'] = fj.msd_raw * awkward.JaggedArray.fromoffsets(fj.array.offsets, np.maximum(1e-5, get_msd_weight(fj.sd.pt.flatten(),fj.sd.eta.flatten())))
        fj['rho'] = 2 * np.log(fj.msd_corr / fj.sd.pt)
        probQCD=fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
        probZHbb=fj.probZbb+fj.probHbb
        fj['ZHbbvsQCD'] = probZHbb/(probZHbb+probQCD)
        probT=fj.probTbcq+fj.probTbqq
        fj['TvsQCD'] = probT/(probT+probQCD)
        probV=fj.probWcq+fj.probWqq+fj.probZbb+fj.probZcc+fj.probZqq
        probX=probZHbb+probV
        fj['XvsQCD'] = probX/(probX+probQCD)
        fj['tau21'] = fj.tau2/fj.tau1
        fj_good = fj[fj.isgood.astype(np.bool)]
        fj_clean = fj_good[fj_good.isclean.astype(np.bool)]
        fj_ntot = fj.counts
        fj_ngood = fj_good.counts
        fj_nclean = fj_clean.counts

        ###
        #Calculating weights
        ###
        if not isData:

            gen = events.GenPart

            gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            jetgenb = fj.sd.cross(gen[gen.isb], nested=True)
            bmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < 1.5).sum()==1)&(gen[gen.isb].counts>0)
            fj['isb']  = bmatch

            bmatch = ((jetgenb.i0.delta_r(jetgenb.i1) < 1.5).sum()==2)&(gen[gen.isb].counts>0)
            fj['isbb']  = bmatch

            gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            jetgenc = fj.sd.cross(gen[gen.isc], nested=True)
            cmatch = ((jetgenc.i0.delta_r(jetgenc.i1) < 1.5).sum()==1)&(gen[gen.isc].counts>0)
            fj['isc']  = cmatch

            cmatch = ((jetgenc.i0.delta_r(jetgenc.i1) < 1.5).sum()==2)&(gen[gen.isc].counts>0)
            fj['iscc']  = cmatch



        ###
        # Selections
        ###

        leading_fj = fj[fj.sd.pt.argmax()]
        leading_fj = leading_fj[leading_fj.isgood.astype(np.bool)]
        #leading_fj = leading_fj[leading_fj.isclean.astype(np.bool)]
        selection.add('fj_pt', (leading_fj.pt.max() > 350) )
        selection.add('fj_mass', (leading_fj.msd_corr.sum() < 80) ) ## optionally also <130
        selection.add('fj_tau21', (leading_fj.tau21.sum() < 0.3) )

        selection.add('mu_pt', (leading_mu.pt.max() > 7) )
        selection.add('pt_ratio', (leading_mu.pt.max()/leading_fj.pt.max() < 0.7) )

        isFilled = False

        variables = {
            'ZHbbvsQCD':              leading_fj.ZHbbvsQCD,
        }
        print('Variables:',variables.keys())

        def fill(dataset, gentype, weight, cut):
            flat_variables = {k: v[cut].flatten() for k, v in variables.items()}
            flat_gentype = {k: (~np.isnan(v[cut])*gentype[cut]).flatten() for k, v in variables.items()}
            flat_weight = {k: (~np.isnan(v[cut])*weight[cut]).flatten() for k, v in variables.items()}

            print('variables:', flat_variables)
            print('gentype:', flat_gentype)
            print('weight:', flat_weight)
            print()
            for histname, h in hout.items():
                if not isinstance(h, hist.Hist):
                    continue
                if histname not in variables:
                    continue
                elif histname == 'sumw':
                    continue
                else:
                    flat_variable = {histname: flat_variables[histname]}
                    h.fill(dataset=dataset, gentype=flat_gentype[histname], **flat_variable, weight=flat_weight[histname])

        if isData:
            if not isFilled:
                hout['sumw'].fill(dataset=dataset, sumw=1, weight=1)
                isFilled=True
            cut = selection.all()
            #fill(dataset, np.zeros(events.size, dtype=np.int), np.ones(events.size), cut)
            fill(dataset, np.zeros(events.size, dtype=np.int), np.ones(events.size), np.ones(events.size, dtype=np.int))
        else:
            weights = processor.Weights(len(events))

            wgentype = {
                'bb' : ( leading_fj.isbb).sum(),
                #'bc' : ( ~leading_fj.isbb & (leading_fj.isb & leading_fj.isc) ).sum(),
                'b'  : ( ~leading_fj.isbb & ~(leading_fj.isb & leading_fj.isc) & leading_fj.isb ).sum(),
                'cc' : ( ~leading_fj.isbb & ~(leading_fj.isb & leading_fj.isc) & ~leading_fj.isb & leading_fj.iscc ).sum(),
                'c'  : ( ~leading_fj.isbb & ~(leading_fj.isb & leading_fj.isc) & ~leading_fj.isb & ~leading_fj.iscc & leading_fj.isc ).sum(),
                'other' : ( ~leading_fj.isbb & ~(leading_fj.isb & leading_fj.isc) & ~leading_fj.isb & ~leading_fj.iscc & ~leading_fj.isc ).sum(),
            }
            vgentype=np.zeros(events.size, dtype=np.int)
            for gentype in self._gentype_map.keys():
                vgentype += self._gentype_map[gentype]*wgentype[gentype]

            if not isFilled:
                hout['sumw'].fill(dataset=dataset, sumw=1, weight=events.genWeight.sum())
                isFilled=True

            cut = selection.all()
            #fill(dataset, vgentype, weights.weight(), cut)
            fill(dataset, vgentype, weights.weight(), np.ones(events.size, dtype=np.int))

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

    save(processor_instance, 'data/doublebSF'+options.year+'.processor')
