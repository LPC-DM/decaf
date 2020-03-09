import os
import numpy
from optparse import OptionParser
from coffea import processor, hist, util
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.btag_tools import BTagScaleFactor
from coffea.util import save, load

class BTagCorrector:
    deepflavWPs = {
        '2016': {
            'loose' : 0.0494,
            'medium': 0.2770,
            'tight' : 0.7264
        },
        '2017': {
            'loose' : 0.0494,
            'medium': 0.2770,
            'tight' : 0.7264
        },
        '2018': {
            'loose' : 0.0494,
            'medium': 0.2770,
            'tight' : 0.7264        
        },
    }
    deepcsvWPs = {
        '2016': {
            'loose' : 0.1241,
            'medium': 0.4184,
            'tight' : 0.7527
        },
        '2017': {
            'loose' : 0.1241,
            'medium': 0.4184,
            'tight' : 0.7527
        },
        '2018': {
            'loose' : 0.1241,
            'medium': 0.4184,
            'tight' : 0.7527
        },
    }

    btagWPs = {
        'deepflav': deepflavWPs,
        'deepcsv' : deepcsvWPs
    }

    def __init__(self, tagger, year, workingpoint):
        self._year = year
        self._wp = BTagCorrector.btagWPs[tagger][year][workingpoint]
        files = {
            'deepflav': {
                '2016': 'DeepJet_102XSF_V1.csv',
                '2017': 'DeepJet_102XSF_V1.csv',
                '2018': 'DeepJet_102XSF_V1.csv',
                },
            'deepcsv': {
                '2016': 'DeepCSV_102XSF_V1.csv',
                '2017': 'DeepCSV_102XSF_V1.csv',
                '2018': 'DeepCSV_102XSF_V1.csv',
                }
        }
        filename = 'data/'+files[tagger][year]
        print(filename)
        self.sf = BTagScaleFactor(filename, workingpoint)
        files = {
            '2016': 'btag2018.merged',
            '2017': 'btag2018.merged',
            '2018': 'btag2018.merged',
        }
        filename = 'hists/'+files[year]
        btag = load(filename)
        bpass = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag', 'pass').values()[()]
        ball = btag[tagger].integrate('dataset').integrate('wp',workingpoint).integrate('btag').values()[()]
        nom = bpass / numpy.maximum(ball, 1.)
        self.eff = dense_lookup(nom, [ax.edges() for ax in btag[tagger].axes()[3:]])

    def btag_weight(self, pt, eta, flavor, tag):
        abseta = abs(eta)
        
        #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1b_Event_reweighting_using_scale
        def zerotag(eff, sf):
            return ((1 - sf*eff) / (1 - eff)).prod()

        eff = self.eff(flavor, pt, abseta)
        sf_nom = self.sf.eval('central', flavor, abseta, pt)
        sf_up = self.sf.eval('up', flavor, abseta, pt)
        sf_down = self.sf.eval('down', flavor, abseta, pt)

        nom = zerotag(eff, sf_nom)
        up = zerotag(eff, sf_up)
        down = zerotag(eff, sf_down)
        if '-1' in tag: 
            nom = 1 - zerotag(eff, sf_nom)
            up = 1 - zerotag(eff, sf_up)
            down = 1 - zerotag(eff, sf_down)
        return nom, up, down

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-w', '--workingpoints', action='store_true', dest='workingpoints')
    (options, args) = parser.parse_args()    

    get_btag_weight = {
        'deepflav': {
            '2016': {
                'loose'  : BTagCorrector('deepflav','2016','loose').btag_weight,
                'medium' : BTagCorrector('deepflav','2016','medium').btag_weight,
                'tight'  : BTagCorrector('deepflav','2016','tight').btag_weight
            },
            '2017': {
                'loose'  : BTagCorrector('deepflav','2017','loose').btag_weight,
                'medium' : BTagCorrector('deepflav','2017','medium').btag_weight,
                'tight'  : BTagCorrector('deepflav','2017','tight').btag_weight
            },
            '2018': {
                'loose'  : BTagCorrector('deepflav','2018','loose').btag_weight,
                'medium' : BTagCorrector('deepflav','2018','medium').btag_weight,
                'tight'  : BTagCorrector('deepflav','2018','tight').btag_weight
            }
        },
        'deepcsv' : {
            '2016': {
                'loose'  : BTagCorrector('deepcsv','2016','loose').btag_weight,
                'medium' : BTagCorrector('deepcsv','2016','medium').btag_weight,
                'tight'  : BTagCorrector('deepcsv','2016','tight').btag_weight
            },
            '2017': {
                'loose'  : BTagCorrector('deepcsv','2017','loose').btag_weight,
                'medium' : BTagCorrector('deepcsv','2017','medium').btag_weight,
                'tight'  : BTagCorrector('deepcsv','2017','tight').btag_weight
            },
            '2018': {
                'loose'  : BTagCorrector('deepcsv','2018','loose').btag_weight,
                'medium' : BTagCorrector('deepcsv','2018','medium').btag_weight,
                'tight'  : BTagCorrector('deepcsv','2018','tight').btag_weight
            }
        }
    }

    if options.workingpoints: save(BTagCorrector.btagWPs, 'data/wp.coffea')
    else: save(get_btag_weight,'data/btag.coffea')
