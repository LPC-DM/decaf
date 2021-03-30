#!/usr/bin/env python
from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option('-a', '--analysis', help='analysis', dest='analysis', default='darkhiggs')
(options, args) = parser.parse_args()

folder_list = []
for mass_point in os.listdir('datacards/'):
    if options.analysis not in mass_point: continue
    if mass_point.split('-')[0] not in folder_list: folder_list.append(mass_point.split('-')[0])

for mass_point in folder_list:
    #print('python macros/combine_cards.py -a '+mass_point)
    os.system('python macros/combine_cards.py -a '+mass_point)
