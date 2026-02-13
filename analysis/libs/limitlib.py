import os
import pickle
import re

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import uproot
from scipy.interpolate import CloughTocher2DInterpolator, Rbf, interp1d
from scipy.optimize import minimize
pjoin = os.path.join


def interpolate_rbf(x,y,z,epsilon=5, function="multiquadric", smooth=3, maxval=2500, nval=100):
    ti = np.linspace(0,maxval, nval)
    xi, yi = np.meshgrid(ti, ti)

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    try:
        assert(len(x)==len(y) and len(x)==len(z))
    except AssertionError:
        print("ERROR: Inconsistent input lengths x,y,z: {0}, {1}, {2}".format(len(x),len(y),len(z)))

    # for tup in zip(x,y,z):
    #     try:
    #         assert(all([isinstance(tmp,float) for tmp in tup]))
    #     except AssertionError:
    #         print(tup)

    valid_indices = z>0
    x = x[valid_indices]
    y = y[valid_indices]
    z = z[valid_indices]


    rbf = Rbf(x, y, np.log(z), epsilon=epsilon,function=function,smooth=smooth)
    zi = np.exp(rbf(xi, yi))

    return xi, yi, zi


def find_input_files(indir):
    infiles = []
    for root, _, files in os.walk(indir):
        for file in files:
            if file.startswith("higg") and file.endswith(".root"):
                infiles.append(pjoin(root, file))
    return infiles


def load_limit(infile):
    m = re.match(
        "higgsCombinemonojet_monov_nominal_tight_(.*).Asymptotic.mH120.root",
        os.path.basename(infile)
    )
    if not m:
        return None


    f = uproot.open(infile)
    try:
        raw = f['limit']['limit'].array()
    except KeyError:
        return None
    if len(raw) != 6:
        return None

    df = pd.DataFrame()
    df['m2s']=[raw[0]]
    df['m1s']=[raw[1]]
    df['exp']=[raw[2]]
    df['p1s']=[raw[3]]
    df['p2s']=[raw[4]]
    df['obs']=[raw[5]]
    
    tag = m.groups()[0]

    df['tag']=[tag]
    if tag.startswith('80'):
        df['par1'] = [float(tag[3:7])]
        df['par2'] = [float(tag[7:112])]
    elif tag.startswith('slq'):
        m = re.match("slq_mlq(\d+)_ylq([\d|p]+)", tag)

        mlq, ylq = m.groups()
        df['par1'] = [float(mlq)]
        df['par2'] = [float(ylq.replace("p","."))]

    return df

def fill_dummy_values(df):
    mmeds = np.unique(df.par1)
    for mmed in mmeds:
        if mmed < 200:
            continue
        mdmmin = np.min(df.par2[df.par1==mmed])

        if mdmmin > mmed / 3:
            continue
        exp = np.min(df.exp[(df.par1==mmed)&(df.par2==mdmmin)])
        if np.isnan(exp):
            continue

        for m in np.linspace(1,mmed/4,10):
            df = df.append({
                "par1" : mmed,
                "par2" : m,
                "exp" : exp,
                "tag" : "dummy"
            },
        ignore_index=True
        )
    return df

def load_limits(infiles):
    return pd.concat(filter(
        lambda x: x is not None,
        map(load_limit, infiles)
        ))


def load_directory(indir):
    infiles = find_input_files(indir)

    limitfile = pjoin(indir, "limits.pkl")
    if os.path.exists(limitfile):
        with open(limitfile, "rb") as f:
            limits = pickle.loads(f.read())
    else:
        limits = load_limits(infiles)
        if len(limits):
            with open(limitfile, "wb") as f:
                pickle.dump(limits, f)
    return limits


def find_intersection(x, y, value):
    f = interp1d(x,np.log(y), fill_value="extrapolate")

    minfun = lambda tmp : (f(tmp)-np.log(value))**2

    result = minimize(minfun, x0=np.mean(x))

    return result.x


def dump_contour_to_txt(contour, outfile):
    """Write contour line points into text files

    :param contour: Contour to dump
    :type contour: Matplotlib contour
    :param outfile: Name of output file. Should contain '{}' to allow numbering
    :type outfile: str
    """
    for i, path in enumerate(contour.collections[0].get_paths()):
        with open(outfile.format(str(i)),"w") as of:
            for point in path.vertices:
                of.write("{} {}\n".format(*point))


brazilgreen = "green"
brazilyellow = "orange"
