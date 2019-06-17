#!/usr/bin/env python
import json
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
parser.add_option('-y', '--year', help='year', dest='year')
parser.add_option('-f', '--fuse', action="store_true", dest="fuse")
(options, args) = parser.parse_args()

fnaleos = "root://cmsxrootd.fnal.gov/"
fuse = '/eos/uscms/'

with open("beans/"+options.year+".json") as fin:
    samplefiles = json.load(fin)

fileslice = slice(None)    
for dataset, info in samplefiles.items():
    if options.dataset and options.dataset not in dataset: continue
    output = dataset+".root"
    hadd = "hadd -f "+output+" "
    haddfuse = hadd
    count = 0
    for file in info['files'][fileslice]:
        count = count +1
        piece = file.strip().split('//')[2]
        fusefile = fuse+piece
        hadd=hadd+file+" "
        haddfuse=haddfuse+fusefile+" "
        if count > 1: break
    pack = hadd
    if options.fuse: pack = haddfuse
    #os.system(pack)
    print(pack)
    eospath = info['files'][fileslice][0].strip().split(dataset.strip().split('____')[0])[0]+dataset.strip().split('____')[0]
    eospath = eospath.replace('coffeabeans','packedbeans')
    fusepath = fuse+eospath.strip().split('//')[2]
    #os.system('mkdir -p '+fusepath)
    copy = 'xrdcp -f '+output+' '+eospath+'/'+output
    if options.fuse: copy= 'cp '+output+' '+fusepath+'/'+output
    #os.system(copy)
    #os.system('rm '+output)
