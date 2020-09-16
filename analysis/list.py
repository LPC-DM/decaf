#!/usr/bin/env python
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-y', '--year', help='year', dest='year')
(options, args) = parser.parse_args()

globalredirect = "root://xrootd-cms.infn.it/"

campaigns ={}
campaigns['2016'] = 'RunIISummer16MiniAODv3'
campaigns['2017'] = 'RunIIFall17MiniAODv2'
campaigns['2018'] = 'RunIIAutumn18NanoAODv6'

slist = ['DarkHiggs_MonoHs_LO_TuneCP5_13TeV-madgraph-pythia8']
datasets = []
for sample in slist:
    query="dasgoclient --query=\"dataset dataset=/"+sample+"*/"+campaigns[options.year]+"*/NANOAODSIM\""
    print(query)
    os.system(query+" >"+sample+".txt")
    datasets.append(sample+".txt")

filenames = []
for dname in datasets:                                                                                                                                                                 
    with open(dname) as infile:
        for line in infile:
            dataset = line.strip()
            query="dasgoclient --query=\"file dataset="+dataset+"\""
            print(query)
            os.system(query+" > "+dataset.split('/')[1]+"____"+dataset.split('/')[2]+".txt")
            filenames.append(dataset.split('/')[1]+"____"+dataset.split('/')[2]+".txt")
    os.system('rm '+dname)

for sample in slist:
    with open('metadata/'+sample+'.txt', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(globalredirect+line)
            os.system('rm '+fname)
