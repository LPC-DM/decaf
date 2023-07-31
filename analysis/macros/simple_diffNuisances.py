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
data_prefit = {}

# loop over all fitted parameters
plotsDir = options.outdir
if not os.path.exists(plotsDir):
    os.mkdir(plotsDir)

if options.writeText:
    fout = open(options.outdir+'/'+options.writeText, 'w')
    my_list = ["name", "val", "err"]
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

        if not sigma_p > 0: sigma_p = (nuis_p.getMax()-nuis_p.getMin())/2

    for fit_name, nuis_x in [('b', nuis_b), ('s',nuis_s)]:
        if nuis_p != None:
            if fit_name=='b':
                data_fitb[name] = {'val':nuis_x.getVal(), 'err':nuis_x.getError()}
                if options.writeText:
                    my_list = [name, str(nuis_x.getVal()), str(nuis_x.getError())]
                    result_string = ",".join(my_list) + "\n"
                    fout.write(result_string)
            if fit_name=='s':
                data_prefit[name] = {'val':mean_p, 'err':sigma_p}

            if sigma_p>0:
                # calculate the difference of the nuisance parameter
                # w.r.t to the prefit value in terms of the uncertainty
                # on the prefit value
                valShift = (nuis_x.getVal() - mean_p)/sigma_p

                # ratio of the nuisance parameter's uncertainty
                # w.r.t the prefit uncertainty
                sigShift = nuis_x.getError()/sigma_p

            else :
                print "No definition for prefit uncertainty %s. Printing absolute shifts"%(nuis_p.GetName())
                valShift = (nuis_x.getVal() - mean_p)
                sigShift = nuis_x.getError()

if options.writeText:
    fout.close()

ndata = len(data_prefit.keys())
# Also make histograms for pull distributions:
hist_fit_b  = ROOT.TH1F("prefit_fit_b"   ,"B-only fit Nuisances;;#theta ",ndata,0,ndata)
hist_empty  = ROOT.TH1F("empty"   ,"empty ",ndata,0,ndata)
hist_fit_s  = ROOT.TH1F("prefit_fit_s"   ,"S+B fit Nuisances   ;;#theta ",ndata,0,ndata)
hist_prefit = ROOT.TH1F("prefit_nuisancs","Prefit Nuisances    ;;#theta ",ndata,0,ndata)

sorted_data_prefit = collections.OrderedDict(sorted(data_prefit.items()))

for i, key in enumerate(sorted_data_prefit.keys()):
    hist_empty.GetXaxis().SetBinLabel(i+1, key)
    hist_fit_b.SetBinContent(i+1, data_fitb[key]['val'])
    hist_fit_b.SetBinError(i+1, data_fitb[key]['err'])
    hist_fit_b.GetXaxis().SetBinLabel(i+1, key)
    hist_prefit.SetBinContent(i+1, data_prefit[key]['val'])
    hist_prefit.SetBinError(i+1, data_prefit[key]['err'])
    hist_prefit.GetXaxis().SetBinLabel(i+1, key)

def getGraph(hist,shift):
    gr = ROOT.TGraphErrors()
    gr.SetName(hist.GetName())
    for j in range(hist.GetNbinsX()):
        x = hist.GetBinCenter(j+1)+shift
        y = hist.GetBinContent(j+1)
        e = hist.GetBinError(j+1)
        gr.SetPoint(j,x,y)
        gr.SetPointError(j,float(abs(shift))*0.8,e)

    return gr

fname = plotsDir+'/plots.root'
fout = ROOT.TFile(fname,"RECREATE")
ROOT.gROOT.SetStyle("Plain")

canvas_nuis = ROOT.TCanvas("nuisances", "nuisances", 2500, 750)
gr_fit_b = getGraph(hist_fit_b, 0.1)
gr_fit_b.SetLineColor(ROOT.kBlue)
gr_fit_b.SetMarkerColor(ROOT.kBlue)
gr_fit_b.SetMarkerStyle(20)
gr_fit_b.SetMarkerSize(1.0)
gr_fit_b.SetLineWidth(2)

hist_empty.SetTitle("")
hist_empty.SetNdivisions(-550,"x")
hist_empty.SetStats(0)
hist_empty.GetYaxis().SetRangeUser(-3, 3)
hist_empty.Draw("histsame")

hist_prefit.SetLineWidth(2)
hist_prefit.SetTitle("")
hist_prefit.SetLineColor(ROOT.kBlack)
hist_prefit.SetFillColor(ROOT.kGray)
hist_prefit.SetMaximum(1)
hist_prefit.SetMinimum(-1)
hist_prefit.SetStats(0)
hist_prefit.SetNdivisions(-550,"x")
hist_prefit.Draw("AE2same")
hist_prefit.Draw("Ahistsame")
gr_fit_b.Draw("EPsame")

canvas_nuis.SetTopMargin(0.02)
canvas_nuis.SetBottomMargin(0.3)
canvas_nuis.SetLeftMargin(0.07)
canvas_nuis.SetRightMargin(0.03)
canvas_nuis.SetGridx()
canvas_nuis.RedrawAxis()
canvas_nuis.RedrawAxis('g')
leg=ROOT.TLegend(0.8,0.87,0.96,0.97)
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.AddEntry(hist_prefit,"Prefit","FL")
leg.AddEntry(gr_fit_b,"B-only fit","EPL")
leg.Draw()
fout.WriteTObject(canvas_nuis)

canvas_nuis.SaveAs(plotsDir+"/pulls.pdf")
canvas_nuis.SaveAs(plotsDir+"/pulls.png")
