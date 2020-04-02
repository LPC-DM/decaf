import os
import numpy
import json
from coffea import processor, hist, util
from coffea.util import save, load
from optparse import OptionParser
from coffea.lookup_tools.dense_lookup import dense_lookup

class BTagEfficiency(processor.ProcessorABC):

    def __init__(self, year,wp):
        self._year = year
        self._btagWPs = wp
        self._accumulator = processor.dict_accumulator({
            'deepflav' :
            hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('wp', 'Working point'),
                hist.Cat('btag', 'BTag WP pass/fail'),
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

        jets = events.Jet[
            (events.Jet.pt > 30.)
            & (abs(events.Jet.eta) < 2.5)
            & (events.Jet.jetId & 2)  # tight id
        ]
        
        name = {}
        name['deepflav']= 'btagDeepFlavB'
        name['deepcsv']= 'btagDeepB'

        out = self.accumulator.identity()

        for wp in ['loose','medium','tight']:
            for tagger in ['deepflav','deepcsv']:
                passbtag = jets[name[tagger]] > self._btagWPs[tagger][self._year][wp]
                out[tagger].fill(
                    dataset=dataset,
                    wp=wp,
                    btag='pass',
                    flavor=jets[passbtag].hadronFlavour.flatten(),
                    pt=jets[passbtag].pt.flatten(),
                    abseta=abs(jets[passbtag].eta.flatten()),
                )
                out[tagger].fill(
                    dataset=dataset,
                    wp=wp,
                    btag='fail',
                    flavor=jets[~passbtag].hadronFlavour.flatten(),
                    pt=jets[~passbtag].pt.flatten(),
                    abseta=abs(jets[~passbtag].eta.flatten()),
                )
        return out

    def postprocess(self, a):
        return a

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    (options, args) = parser.parse_args()


    with open('metadata/'+options.year+'.json') as fin:
        samplefiles = json.load(fin)

    common = load('data/common.coffea')
    processor_instance=BTagEfficiency(year=options.year,wp=common['btagWPs'])

    save(processor_instance, 'data/btag'+options.year+'.processor')
