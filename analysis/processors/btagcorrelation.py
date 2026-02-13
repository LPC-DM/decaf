import os
import numpy
import json
from coffea import processor, hist, util
from coffea.util import save, load
from optparse import OptionParser
import numpy as np
from coffea.lookup_tools.dense_lookup import dense_lookup

class BTagEfficiency(processor.ProcessorABC):

    def __init__(self, year, wp, ids):
        self._year = year
        self._btagWPs = wp
        self._ids = ids
        self._common = common
        self._ZHbbvsQCDwp = {
            '2016': 0.53,
            '2017': 0.61,
            '2018': 0.65
        }
        self._accumulator = processor.dict_accumulator({
            'deepflav' :
            hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('wp', 'Working point'),
                hist.Cat('btag', 'BTag WP pass/fail'),
                hist.Cat('doublebtag', 'DBTag WP pass/fail'),
                hist.Bin('flavor', 'Jet hadronFlavour', [0, 4, 5, 6]),
                hist.Bin('pt', 'Jet pT', [20, 30, 50, 70, 100, 140, 200, 300, 600, 1000]),
                hist.Bin('abseta', 'Jet abseta', [0, 1.4, 2.0, 2.5]),
            ),
            'deepcsv' :
            hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('wp', 'Working point'),
                hist.Cat('btag', 'BTag WP pass/fail'),
                hist.Cat('doublebtag', 'DBTag WP pass/fail'),
                hist.Bin('flavor', 'Jet hadronFlavour', [0, 4, 5, 6]),
                hist.Bin('pt', 'Jet pT', [20, 30, 50, 70, 100, 140, 200, 300, 600, 1000]),
                hist.Bin('abseta', 'Jet abseta', [0, 1.4, 2.0, 2.5]),
            )
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        
        dataset = events.metadata['dataset']
        selection = processor.PackedSelection()
        isGoodJet = self._ids['isGoodJet']
        isGoodFatJet    = self._ids['isGoodFatJet']   
        match = self._common['match']

        fj = events.AK15Puppi
        fj['sd'] = fj.subjets.sum()
        fj['isgood'] = isGoodFatJet(fj.sd.pt, fj.sd.eta, fj.jetId)
        probQCD=fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
        probZHbb=fj.probZbb+fj.probHbb
        fj['ZHbbvsQCD'] = probZHbb/(probZHbb+probQCD)
        fj_good = fj[fj.isgood.astype(np.bool)]
        fj_ngood = fj_good.counts
        leading_fj = fj[fj.sd.pt.argmax()]
        leading_fj = leading_fj[leading_fj.isgood.astype(np.bool)]

        j = events.Jet
        j['isgood'] = isGoodJet(j.pt, j.eta, j.jetId, j.puId, j.neHEF, j.chHEF)
        j['isiso'] = ~match(j,leading_fj,1.5)
        j_good = j[j.isgood.astype(np.bool)]
        j_iso = j_good[j_good.isiso.astype(np.bool)]
        j_niso = j_iso.counts
        leading_j = j[j.pt.argmax()]
        leading_j = leading_j[leading_j.isgood.astype(np.bool)]
        leading_j = leading_j[leading_j.isiso.astype(np.bool)]

        name = {}
        name['deepflav']= 'btagDeepFlavB'
        name['deepcsv']= 'btagDeepB'

        out = self.accumulator.identity()

        selection.add('passdoublebtag',(leading_fj.ZHbbvsQCD.sum() > self._ZHbbvsQCDwp[self._year]))
        selection.add('faildoublebtag',~(leading_fj.ZHbbvsQCD.sum() > self._ZHbbvsQCDwp[self._year]))
        selection.add('fatjet', (fj_ngood>0))
        selection.add('isojet', (j_niso>0))

        regions = {
            'pass':['fatjet','passdoublebtag'],
            'fail':['fatjet','faildoublebtag']
        }

        for region, cuts in regions.items():
            cut = selection.all(*regions[region])
            for wp in ['loose','medium','tight']:
                for tagger in ['deepflav','deepcsv']:
                    passbtag = leading_j[name[tagger]] > self._btagWPs[tagger][self._year][wp]
                    out[tagger].fill(
                        dataset=dataset,
                        wp=wp,
                        btag='pass',
                        doublebtag=region,
                        flavor=leading_j[passbtag].hadronFlavour.sum(),
                        pt=leading_j[passbtag].pt.sum(),
                        abseta=abs(leading_j[passbtag].eta.sum()),
                        weight=np.ones(events.size)*cut
                    )
                    out[tagger].fill(
                        dataset=dataset,
                        wp=wp,
                        btag='fail',
                        doublebtag=region,
                        flavor=leading_j[~passbtag].hadronFlavour.sum(),
                        pt=leading_j[~passbtag].pt.sum(),
                        abseta=abs(leading_j[~passbtag].eta.sum()),
                        weight=np.ones(events.size)*cut
                    )
                
        return out

    def postprocess(self, a):
        return a

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
    parser.add_option('-n', '--name', help='name', dest='name')
    (options, args) = parser.parse_args()


    with open('metadata/'+options.metadata+'.json') as fin:
        samplefiles = json.load(fin)

    common = load('data/common.coffea')
    ids    = load('data/ids.coffea')
    processor_instance=BTagEfficiency(year=options.year,wp=common['btagWPs'],ids=ids)

    save(processor_instance, 'data/btagcorrelation'+options.name+'.processor')
