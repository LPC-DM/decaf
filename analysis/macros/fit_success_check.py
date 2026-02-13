from glob import glob
from os import path, system
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', '--fit', action='store_true', dest='fit')
parser.add_option('-a', '--limit', action='store_true', dest='limit')
parser.add_option('-d', '--directory', dest='dirt', default='datacards')
(options, args) = parser.parse_args()

names = []
names2 = []
cfit = 0

for ifile in glob(options.dirt+'/*'):
    if options.fit:
        import ROOT as rt
        try:
            print "Which mass point?", ifile.split('/')[1].split('.')[2]
            fin = rt.TFile.Open(ifile)
            success = fin.Get('fit_s')

            if success == None:
                cfit += 1
                names.append(ifile.split('/')[1].split('.')[2])
                fin.Close()
                continue
            else:
                fin.Close()
        except:
            pass

    if options.limit:
        import ROOT as rt
        try:
            print "Which mass point?", ifile.split('/')[1].split('.')[2]
            fin = rt.TFile.Open(ifile)
            success = fin.Get('limit')

            if success.GetEntries() != 6:
                cfit += 1
                names.append(ifile.split('/')[1].split('.')[2])
                fin.Close()
                continue
            else:
                names2.append(ifile.split('/')[1].split('.')[2])
                fin.Close()
        except:
            pass

if options.fit:
    print 'Failed mass points:', names, len(names)
    print 'Success mass points:', names2, len(names2)

if options.limit:
    print 'Failed mass points:', sorted(names), len(names)
    print 'Success mass points:', sorted(names2), len(names2)
