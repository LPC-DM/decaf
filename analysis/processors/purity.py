import os
import numpy
import json
from coffea import processor, hist, util
from coffea.util import save, load
from optparse import OptionParser
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.lumi_tools import LumiMask


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
    return (photons.vid & mask) == mask

def medium_id_no_sieie_inv_iso(photons):
    # Same as medium_id_no_sieie, but in first step,
    # ignore isolation bits
    # --> XX-XX-XX-XX-1X-1X-1X
    mask1 = int('00000000101010',2)
    medium_id_no_iso = (photons.vid & mask1) == mask1

    # In second step, require that at least one
    # of the most significant isolation bits
    # fails, which we achieve by using a mask
    # that would *pass* and then requiring that
    # (vid & mask) != mask
    mask2 = int('10101000000000',2)
    inv_iso = (photons.vid & mask2) != mask2

    return medium_id_no_iso & inv_iso

class PhotonPurity(processor.ProcessorABC):

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


    def __init__(self, year, ids, common):
        self._year = year
        self._accumulator = processor.dict_accumulator({
            'sumw': hist.Hist(
                'sumw',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('sumw', 'Weight value', [0.])
            ),
            'sumw2': hist.Hist(
                'sumw2',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('sumw2', 'Weight value', [0.])
            ),
            'count':
            hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('cat', 'Cat'),
                hist.Bin('pt', 'Photon pT', 50, 200, 1200),
                hist.Bin('sieie', 'sieie', 100, 0, 0.02),
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
        match = self._common['match']

        isLooseElectron = self._ids['isLooseElectron']
        isLooseMuon     = self._ids['isLooseMuon']
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
            if year == '2016':
                mask = (pt>200)&(medium_id>=2)
            else:
                mask = (pt>200)&((medium_id&2)==2)
            return mask

        pho['ispurity'] = isPurityPhoton(pho.pt, pho[_id])&(pho.isScEtaEB)&(pho.electronVeto)
        pho_clean=pho[pho.isclean.astype(np.bool)]
        pho_purity=pho_clean[pho_clean.ispurity.astype(np.bool)]

        #### Consider AK4 jet 
        def isPurityJet(pt, eta, jet_id):
            mask = (pt > 100) & (abs(eta) < 2.4) & ((jet_id&2)==2)
            return mask

        j = events.Jet
        j['ispurity'] = isPurityJet(j.pt, j.eta, j,jetId)
        j['isclean'] = ~match(j,e_loose,0.4)&~match(j,mu_loose,0.4)&~match(j,pho_loose,0.4)
        j_purity = j[j.ispurity.astype(np.bool)]
        j_clean = j_purity[j_purity.isclean.astype(np.bool)]
        j_nclean = j_clean.counts

        met = events.MET

        '''
        #### Lumi mask
        if isData:
            if self._year == '2016':
                json = #filepath
            elif self._year == '2017':
                json = #filepath
            elif self._year == '2018':
                json = #filepath
            lumi_mask = LumiMask(json)(events['run'], events['luminosityBlock'])
        else
            lumi_mask = np.ones(events.size)==1
        '''

        #### Genweights
        weights = processor.Weights(len(events), storeIndividual=True)

        if isData:
            weights.add('genw', np.ones(events.size))
        else:
            weights.add('genw', events.genWeight)

        #### MET filter & single photon trigger
        met_filters =  np.ones(events.size, dtype=np.bool)
        if isData: met_filters = met_filters & events.Flag['eeBadScFilter']
        for flag in AnalysisProcessor.met_filter_flags[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters',met_filters)

        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._singlephoton_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('singlephoton_triggers', triggers)

        selection.add('jet_cut', (j_nclean>0))
        selection.add('met60', (met.pt<60))

        #{'met_filters', 'singlephoton_triggers', 'jet_cut', 'met60'}

    def postprocess(self, accumulator):
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
                                    common=common)

    save(processor_instance, 'data/purity'+options.year+'.processor')
