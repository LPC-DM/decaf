import os
import numpy as np
import json
from coffea import processor, hist, util
from coffea.util import save, load
from optparse import OptionParser
from coffea.lookup_tools.dense_lookup import dense_lookup


# VID map definition
# 
# 0  1  MinPtCut
# 2  3  PhoSCEtaMultiRangeCut
# 4  5  PhoSingleTowerHadOverEmCut
# 6  7  PhoFull5x5SigmaIEtaIEtaCut
# 8  9  PhoAnyPFIsoWithEACut
# 10 11 PhoAnyPFIsoWithEAAndQuadScalingCut
# 12 13 PhoAnyPFIsoWithEACut

def medium_id_no_sieie(photons):
    # VID bitmask encodes 7 cuts with each two bits.
    # for each cut, the more significant bit indicates
    # whether the medium WP of the cut is passed.
    #
    # To get medium or tighter, require every second bit:
    #   --> 1X-1X-1X-1X-1X-1X-1X
    #   where X = 0 or 1 (we dont care)
    #
    # To ignore the sieie requirement, also ignore bit 7
    #   --> 1X-1X-1X-XX-1X-1X-1X
    # 
    # We are going to check the logical AND of the mask
    # and the bitmap, so set all X to 0
    mask = int('10101000101010',2)
    #return (photons.vid & mask) == mask
    return (photons.vidNestedWPBitmap & mask) == mask

def medium_id_no_sieie_inv_iso(photons):
    # Same as medium_id_no_sieie, but in first step,
    # ignore isolation bits
    # --> XX-XX-XX-XX-1X-1X-1X
    mask1 = int('00000000101010',2)
    #medium_id_no_iso = (photons.vid & mask1) == mask1
    medium_id_no_iso = (photons.vidNestedWPBitmap & mask1) == mask1

    # In second step, require that at least one
    # of the most significant isolation bits
    # fails, which we achieve by using a mask
    # that would *pass* and then requiring that
    # (vid & mask) != mask
    mask2 = int('10101000000000',2)
    #inv_iso = (photons.vid & mask2) != mask2
    inv_iso = (photons.vidNestedWPBitmap & mask2) != mask2

    return medium_id_no_iso & inv_iso

def weight_shape(values, weight):
    #Broadcasts weight array to right shape for given values
    return (~np.isnan(values) * weight).flatten()

class PhotonPurity(processor.ProcessorABC):

    lumis = {
        '2016': 35.92,
        '2017': 40.66,
        '2018':'59.74'
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
                'BadPFMuonFilter'
                ],
        '2018': ['goodVertices',
                'globalSuperTightHalo2016Filter',
                'HBHENoiseFilter',
                'HBHENoiseIsoFilter',
                'EcalDeadCellTriggerPrimitiveFilter',
                'BadPFMuonFilter'
                ]
        }


    def __init__(self, year, ids, xsec, common):
        self._year = year

        self._lumi = 1000.*float(PhotonPurity.lumis[year])
        self._xsec = xsec

        self._accumulator = processor.dict_accumulator({
            'sumw': hist.Hist(
                'sumw',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('sumw', 'Weight value', [0.])
                ),
            'count': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('cat', 'Cat'),
                hist.Bin('pt', 'Photon pT', 50, 200, 1200),
                hist.Bin('sieie', 'sieie', 100, 0, 0.02)
                )
            })

        self._singlephoton_triggers = {
            '2016': [
                'Photon175',
                'Photon165_HE10'
            ],
            '2017': [
                'Photon200'
            ],
            '2018': [
                'Photon200'
            ]
        }

        self._ids = ids
        self._common = common

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):

        dataset = events.metadata['dataset']
        isData = 'genWeight' not in events.columns
        selection = processor.PackedSelection()
        hout = self.accumulator.identity()
        match = self._common['match']

        isLooseElectron = self._ids['isLooseElectron']
        isLooseMuon     = self._ids['isLooseMuon']
        isLoosePhoton   = self._ids['isLoosePhoton']
        isTightPhoton   = self._ids['isTightPhoton']
        isGoodJet       = self._ids['isGoodJet']

        #### Select loose muon and electron to select clean photon
        mu = events.Muon
        mu['isloose'] = isLooseMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.looseId,self._year)
        mu_loose=mu[mu.isloose.astype(np.bool)]

        e = events.Electron
        e['isclean'] = ~match(e,mu_loose,0.3)
        e['isloose'] = isLooseElectron(e.pt,e.eta+e.deltaEtaSC,e.dxy,e.dz,e.cutBased,self._year)
        e_clean = e[e.isclean.astype(np.bool)]
        e_loose = e_clean[e_clean.isloose.astype(np.bool)]

        #### Consider clean and tight photon for purity measurement
        pho = events.Photon
        pho['isclean']=~match(pho,mu_loose,0.5)&~match(pho,e_loose,0.5)

        _id = 'cutBasedBitmap'
        if self._year=='2016':
           _id = 'cutBased'

        def isPurityPhoton(pt, medium_id):
            mask = ~(pt == np.nan)
            if self._year == '2016':
                mask = (pt>200)&(medium_id>=2)
            else:
                mask = (pt>200)&((medium_id&2)==2)
            return mask

        pho['isloose']=isLoosePhoton(pho.pt,pho.eta,pho[_id],self._year)&(pho.electronVeto)
        pho['ispurity'] = isPurityPhoton(pho.pt, pho[_id])&(pho.isScEtaEB)&(pho.electronVeto)
        pho_clean=pho[pho.isclean.astype(np.bool)]
        pho_loose=pho_clean[pho_clean.isloose.astype(np.bool)]
        pho_purity=pho_clean[pho_clean.ispurity.astype(np.bool)]
        pho_nosieie = pho_clean[(pho_clean.pt > 200) & (pho_clean.isScEtaEB) & (pho_clean.electronVeto) & medium_id_no_sieie(pho_clean)]
        pho_nosieie_inv_iso = pho_clean[(pho_clean.pt > 200) & (pho_clean.isScEtaEB) & (pho_clean.electronVeto) & medium_id_no_sieie_inv_iso(pho_clean)]

        #### Consider AK4 jet 
        def isPurityJet(pt, eta, jet_id):
            mask = (pt > 30) & (abs(eta) < 2.4) & ((jet_id&2)==2)
            return mask

        j = events.Jet
        #30 GeV cut on jet pT, we need to check later
        #j['isgood'] = isGoodJet(j.pt, j.eta, j.jetId, j.neHEF, j.neEmEF, j.chHEF, j.chEmEF)
        j['ispurity'] = isPurityJet(j.pt, j.eta, j.jetId)
        j['isclean'] = ~match(j,e_loose,0.4)&~match(j,mu_loose,0.4)&~match(j,pho_loose,0.4)
        j_purity = j[j.ispurity.astype(np.bool)]
        j_clean = j_purity[j_purity.isclean.astype(np.bool)]
        j_nclean = j_clean.counts

        met = events.MET

        #### Genweights
        weights = processor.Weights(len(events), storeIndividual=True)

        if isData:
            weights.add('genw', np.ones(events.size))
        else:
            weights.add('genw', events.genWeight)

        #### MET filter & single photon trigger
        met_filters =  np.ones(events.size, dtype=np.bool)
        if isData: met_filters = met_filters & events.Flag['eeBadScFilter']
        for flag in PhotonPurity.met_filter_flags[self._year]:
            met_filters = met_filters & events.Flag[flag]
        #selection.add('met_filters',met_filters)

        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._singlephoton_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        #selection.add('singlephoton_triggers', triggers)

        #selection.add('jet_cut', (j_nclean>0))
        #selection.add('met60', (met.pt<60))

        event_mask = met_filters & triggers & (met.pt<60) & (j_nclean>0)

        hout['count'].fill(dataset = dataset, cat = 'medium',
                            sieie = pho_purity.sieie[event_mask].flatten(),
                            pt = pho_purity.pt[event_mask].flatten(),
                            weight=weight_shape(pho_purity.sieie[event_mask], weights.weight()[event_mask])
                            )

        hout['count'].fill(dataset = dataset, cat = 'medium_nosieie',
                            sieie = pho_nosieie.sieie[event_mask].flatten(),
                            pt = pho_nosieie.pt[event_mask].flatten(),
                            weight=weight_shape(pho_nosieie.sieie[event_mask], weights.weight()[event_mask])
                            )

        hout['count'].fill(dataset = dataset, cat = 'medium_nosieie_invertiso',
                            sieie = pho_nosieie_inv_iso.sieie[event_mask].flatten(),
                            pt = pho_nosieie_inv_iso.pt[event_mask].flatten(),
                            weight=weight_shape(pho_nosieie_inv_iso.sieie[event_mask], weights.weight()[event_mask])
                            )

        if isData:
            hout['sumw'].fill(dataset = dataset, sumw=1, weight=1)
        else:
            hout['sumw'].fill(dataset = dataset, sumw=1, weight=events.genWeight.sum())

        return hout

    def postprocess(self, accumulator):
        scale = {}
        for d in accumulator['sumw'].identifiers('dataset'):
            print('Scaling:',d.name)
            dataset = d.name
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
    (options, args) = parser.parse_args()

    with open('metadata/'+options.year+'.json') as fin:
        samplefiles = json.load(fin)
        xsec = {k: v['xs'] for k,v in samplefiles.items()}

    ids         = load('data/ids.coffea')
    common      = load('data/common.coffea')

    processor_instance=PhotonPurity(year=options.year,
                                    ids=ids,
                                    xsec=xsec,
                                    common=common)

    save(processor_instance, 'data/purity'+options.year+'.processor')
