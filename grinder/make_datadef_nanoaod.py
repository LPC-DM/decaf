#!/usr/bin/env python
import uproot
import json
import pyxrootd.client
import fnmatch
import numpy as np
import numexpr
import subprocess
import concurrent.futures
import warnings
import os
import difflib

uproot_xrootd_opts = dict(chunkbytes=30*1024, limitbytes=20*(1024**2))
fnaleos = "root://cmsxrootd.fnal.gov/"
#coffeabeans2016 = "/eos/uscms/store/group/lpccoffea/coffeabeans/nano_2016"
coffeabeans2016 = "/eos/uscms/store/group/lpcstop/noreplica/NanoTuples/v1a"
coffeabeans2017 = "/eos/uscms/store/group/lpccoffea/coffeabeans/nano_2017"

def parse_xsec(cfgfile):
    xsec_dict = {}
    with open(cfgfile) as f:
        for l in f:
            l = l.strip()
            if l.startswith('#'):
                continue
            pieces = l.split()
            samp = None
            xsec = None
            isData = False
            for s in pieces:
                if 'AOD' in s:
                    samp = s.split('/')[1]
                    if 'AODSIM' not in s:
                        isData = True
                        break
                else:
                    try:
                        xsec = float(s)
                    except ValueError:
                        try:
                            import numexpr
                            xsec = numexpr.evaluate(s).item()
                        except:
                            pass
            if samp is None:
                print('Ignore line:\n%s' % l)
            elif not isData and xsec is None:
                print('Cannot find cross section:\n%s' % l)
            else:
                xsec_dict[samp] = xsec
    return xsec_dict

xsections = parse_xsec("data/xsec.conf")
datadef = {}
for dataset in xsections.keys():
    print("looking into",coffeabeans2016+"/"+dataset)
    try:
        os.system("find "+coffeabeans2016+"/"+dataset+" -name *.root > "+dataset+".txt")
    except:
        print(dataset, "not found")
        continue
    else:
        print(dataset,"found")
        flist = open(dataset+".txt")
        urllist = [fnaleos+path.strip() for path in flist]
        xs = xsections[dataset]
        datadef[dataset] = {
            'files': urllist,
            'xs': xs,
        }
with open("data/coffeabeans2016.json", "w") as fout:
    json.dump(datadef, fout, indent=4)
os.system("rm *.txt")
