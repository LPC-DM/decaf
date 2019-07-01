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
from optparse import OptionParser
from process import *

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-p', '--pack', help='pack', dest='pack')
(options, args) = parser.parse_args()
fnaleos = "root://cmsxrootd.fnal.gov/"

beans={}
#beans['2016'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/nano_2016",
#                 "/eos/uscms/store/group/lpcstop/noreplica/NanoTuples/v1a",
#                 "/eos/uscms/store/user/hqu/NanoTuples/v1d"]
#beans['2017'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/nano_2017"]
beans['2016'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/102X/nano_2016"]
beans['2017'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/102X/nano_2017"]
beans['2018'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/102X/nano_2018"]

def split(arr, size):
     arrs = []
     while len(arr) > size:
         pice = arr[:size]
         arrs.append(pice)
         arr   = arr[size:]
     arrs.append(arr)
     return arrs

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

#xsections = parse_xsec("data/xsec.conf")
xsections={}
for k,v in processes.items():
    if v[1]=='MC':
        xsections[k] = v[2]
    else:
        xsections[k] = -1

datadef = {}
for folder in beans[options.year]:
    print("Opening",folder)
    for dataset in xsections.keys():
        if options.dataset and options.dataset not in dataset: continue
        print("Looking into",folder+"/"+dataset)
        os.system("find "+folder+"/"+dataset+" -name \'*.root\' > "+dataset+".txt")
        flist = open(dataset+".txt")
        urllist = []
        #print('file lenght:',len(flist.readlines()))
        xs = xsections[dataset]
        #sumw = 0
        for path in flist:
            s = path.strip().split('/')
            eospath = fnaleos
            for i in range (3,len(s)): eospath=eospath+'/'+s[i]
            #if xs != -1:
            #     run_tree = uproot.open(eospath)["Runs"]
            #     sumw += run_tree.array("genEventSumw")[0]
            if (not ('failed' in eospath)): urllist.append(eospath)
        print('list lenght:',len(urllist))
        urllists = split(urllist, int(options.pack))
        print(len(urllists))
        if urllist:
            for i in range(0,len(urllists)) :
                 datadef[dataset+"____"+str(i)] = {
                      'files': urllists[i],
                      'xs': xs,
                      #'sumw': sumw,
                      }
        os.system("rm "+dataset+".txt")

os.system("mkdir -p beans")
with open("beans/"+options.year+".json", "w") as fout:
    json.dump(datadef, fout, indent=4)
