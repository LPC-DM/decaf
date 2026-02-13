#!/usr/bin/env python
import json
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-j', '--jsons', help='jsons', dest='jsons', default='')
parser.add_option('-o', '--output', help='output', dest='output', default='')
(options, args) = parser.parse_args()

dictionary={}
for j in options.jsons.split(','):
    print(j)
    with open(j) as fin:
        dictionary.update(json.load(fin))

with open(options.output, "w") as fout:
    json.dump(dictionary, fout, indent=4)
#print(dictionary)
