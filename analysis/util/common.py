from coffea.util import save
import awkward
import uproot, uproot_methods
import numpy as np

def match(a, b, val):
    combinations = a.cross(b, nested=True)
    return (combinations.i0.delta_r(combinations.i1)<val).any()

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
    
common = {}
common['match'] = match
common['btagWPs'] = btagWPs
save(common, 'data/common.coffea')
