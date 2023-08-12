#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit
from optparse import OptionParser
import collections
import os

# tool to compare fitted nuisance parameters to prefit values.
#
# Also used to check for potential problems in RooFit workspaces to be used with combine
# (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsPAGPreapprovalChecks)

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
hasHelp = False
for X in ("-h", "-?", "--help"):
    if X in argv:
        hasHelp = True
        argv.remove(X)
argv.append( '-b-' )

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
argv.remove( '-b-' )
if hasHelp: argv.append("-h")

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("-p", "--poi",      dest="poi",    default="r",    type="string",  help="Name of signal strength parameter (default is 'r' as per text2workspace.py)")
parser.add_option("-o", "--outdir", dest="outdir", default='./output', type="string", help="directory to save outputs.")
parser.add_option("-w", "--writeText", dest="writeText", default='params.txt', type="string", help="save parameters in text file.")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

file = ROOT.TFile(args[0])
if file == None: raise RuntimeError, "Cannot open file %s" % args[0]
fit_s  = file.Get("fit_s")
fit_b  = file.Get("fit_b")
prefit = file.Get("nuisances_prefit")
if fit_s == None or fit_s.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the signal fit 'fit_s'"     % args[0]
if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % args[0]
if prefit == None or prefit.ClassName() != "RooArgSet":    raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'"  % args[0]

isFlagged = {}

# maps from nuisance parameter name to the row to be printed in the table
table = {}

# get the fitted parameters
fpf_b = fit_b.floatParsFinal()
fpf_s = fit_s.floatParsFinal()

pulls = []

nuis_p_i=0
data_fitb = {}
data_fits = {}
data_prefit = {}

# loop over all fitted parameters
plotsDir = options.outdir
if not os.path.exists(plotsDir):
    os.mkdir(plotsDir)

if options.writeText:
    fout = open(options.outdir+'/'+options.writeText, 'w')
    my_list = ["name", "val", "err", "flag"]
    result_string = ",".join(my_list) + "\n"
    fout.write(result_string)

for i in range(fpf_s.getSize()):

    nuis_s = fpf_s.at(i)
    name   = nuis_s.GetName();
    nuis_b = fpf_b.find(name)
    nuis_p = prefit.find(name)

    if options.poi.replace('_r','') in name:
        continue
    elif 'mcstat' in name:
        continue
    elif 'doublebtag' in name:
        continue

    mean_p, sigma_p = 0,0
    if nuis_p == None:
        # nuisance parameter NOT present in the prefit result
        pass
    else:
        # get best-fit value and uncertainty at prefit for this 
        # nuisance parameter
        mean_p, sigma_p = (nuis_p.getVal(), nuis_p.getError())

        if not sigma_p > 0: sigma_p = 1. #(nuis_p.getMax()-nuis_p.getMin())/2

    for fit_name, nuis_x in [('b', nuis_b), ('s',nuis_s)]:
        if nuis_p != None:
            if fit_name=='b':
                data_fitb[name] = {'val':nuis_x.getVal(), 'err':nuis_x.getError()}
                if options.writeText:
                    my_list = [name, str(nuis_x.getVal()), str(nuis_x.getError()), "fitb"]
                    result_string = ",".join(my_list) + "\n"
                    fout.write(result_string)
            if fit_name=='s':
                if options.writeText:
                    my_list = [name, str(nuis_x.getVal()), str(nuis_x.getError()), "fits"]
                    result_string = ",".join(my_list) + "\n"
                    fout.write(result_string)
                data_prefit[name] = {'val':mean_p, 'err':sigma_p}

            if sigma_p>0:
                # calculate the difference of the nuisance parameter
                # w.r.t to the prefit value in terms of the uncertainty
                # on the prefit value
                valShift = (nuis_x.getVal() - mean_p)/sigma_p #sqrt( (sigma_p*sigma_p) - (nuis_x.getError()*nuis_x.getError() )

                # ratio of the nuisance parameter's uncertainty
                # w.r.t the prefit uncertainty
                sigShift = nuis_x.getError()/sigma_p

            else :
                print "No definition for prefit uncertainty %s. Printing absolute shifts"%(nuis_p.GetName())
                valShift = (nuis_x.getVal() - mean_p)
                sigShift = nuis_x.getError()

if options.writeText:
    fout.close()

