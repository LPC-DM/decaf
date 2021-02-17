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
from data.process import *

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
parser.add_option('-p', '--pack', help='pack', dest='pack')
parser.add_option('-k', '--keep', action="store_true", dest="keep")
parser.add_option('-s', '--special', help='special', dest='special')
(options, args) = parser.parse_args()
fnaleos = "root://cmseos.fnal.gov/"
#fnaleos = "root://cmsxrootd.fnal.gov/"

beans={}
beans['2016'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/NanoAODv6/nano_2016","/eos/uscms/store/group/lpccoffea/coffeabeans/NanoAODv6/nano_2016/Signals/monohs"]
beans['2017'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/NanoAODv6/nano_2017","/eos/uscms/store/group/lpccoffea/coffeabeans/NanoAODv6/nano_2017/Sandeep","/eos/uscms/store/group/lpccoffea/coffeabeans/NanoAODv6/nano_2017/Signals/monohs"]
beans['2018'] = ["/eos/uscms/store/group/lpccoffea/coffeabeans/NanoAODv6/nano_2018","/eos/uscms/store/group/lpccoffea/coffeabeans/NanoAODv6/nano_2018/Signals/monohs"]
                 

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
          if not isinstance(k, str):
               print(k)
               print(options.year,k[1])
               if options.year!=str(k[1]): continue
               xsections[k[0]] = v[2]
          else: 
               xsections[k] = v[2]
     else:
          xsections[k] = -1
print(xsections)
datadef = {}
for folder in beans[options.year]:
    print("Opening",folder)
    for dataset in xsections.keys():
        if options.dataset and options.dataset not in dataset: continue
        print("Looking into",folder+"/"+dataset)
        filenames = folder+"/"+dataset+" -name \'*nano*.root\'"
        exist=False
        for filename in os.listdir('metadata'):
             if dataset+".txt" not in filename: continue
             exist=True
        if not exist:
             os.system("find "+filenames+" > metadata/"+dataset+".txt")
        with open("metadata/"+dataset+".txt") as flist:
             new_content=flist.read().replace('/eos/uscms',fnaleos)
        with open("metadata/"+dataset+".txt", 'w') as flist:
             flist.write(new_content)
        if options.keep and open("metadata/"+dataset+".txt").read(1): 
             os.system("mkdir -p metadata/"+options.year)
             os.system("cp metadata/"+dataset+".txt metadata/"+options.year)
        urllist = []
        xs = xsections[dataset]
        for path in open("metadata/"+dataset+".txt"):
            eospath = path.strip()
            if (not ('failed' in eospath)): urllist.append(eospath)
        print('list lenght:',len(urllist))
        if options.special:
             sdataset, spack = options.special.split(':')
             if sdataset in dataset:
                  urllists = split(urllist, int(spack))
             else:
                  urllists = split(urllist, int(options.pack))
        else:
             urllists = split(urllist, int(options.pack))
        print(len(urllists))
        if urllist:
            for i in range(0,len(urllists)) :
                 datadef[dataset+"____"+str(i)+"_"] = {
                      'files': urllists[i],
                      'xs': xs,
                      }
        os.system("rm metadata/"+dataset+".txt")

folder = "metadata/"+options.metadata+".json"
with open(folder, "w") as fout:
    json.dump(datadef, fout, indent=4)
