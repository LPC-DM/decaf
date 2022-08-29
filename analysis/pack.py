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

def find(_list):
    files=[]
    print('Looking into',_list)
    for path in _list:
        command='xrdfs '+fnaleos+' ls '+path
        results=os.popen(command).read()
        files.extend(results.split())
    if not any('.root' in _file for _file in files):
        files=find(files)
    return files
  
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
        xs = xsections[dataset]
        path=folder+'/'+dataset
        urllist = find([path])
        for path in urllist:
            if 'failed' in urllist: urllist.remove(path)
            if '.root' not in urllist: urllist.remove(path)
            if 'nano' not in urllist: urllist.remove(path)
            path.replace('/store/',fnaleos+'store/')
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
        
folder = "metadata/"+options.metadata+".json"
with open(folder, "w") as fout:
    json.dump(datadef, fout, indent=4)
