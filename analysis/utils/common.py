from coffea.util import save
import awkward
import uproot, uproot_methods
import numpy as np

def match(a, b, val):
    combinations = a.cross(b, nested=True)
    return (combinations.i0.delta_r(combinations.i1)<val).any()

def sigmoid(x,a,b,c,d):
    """
    Sigmoid function for trigger turn-on fits.
    f(x) = c + (d-c) / (1 + np.exp(-a * (x-b)))
    """
    return c + (d-c) / (1 + np.exp(-a * (x-b)))

deepflavWPs = {
    '2016': {
        'loose' : 0.0614,
        'medium': 0.3093,
        'tight' : 0.7221
    },
    '2017': {
        'loose' : 0.0521,
        'medium': 0.3033,
        'tight' : 0.7489
    },
    '2018': {
        'loose' : 0.0494,
        'medium': 0.2770,
        'tight' : 0.7264        
    },
}
deepcsvWPs = {
    '2016': {
        'loose' : 0.2217,
        'medium': 0.6321,
        'tight' : 0.8953
    },
    '2017': {
        'loose' : 0.1522,
        'medium': 0.4941,
        'tight' : 0.8001
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
common['sigmoid'] = sigmoid
common['btagWPs'] = btagWPs
save(common, 'data/common.coffea')
