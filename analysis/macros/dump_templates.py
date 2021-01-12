#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
import re
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

def _RooAbsCollection__iter__(self):
    it = self.iterator()
    obj = it.Next()
    while obj != None:
        yield obj
        obj = it.Next()

ROOT.RooAbsCollection.__iter__ = _RooAbsCollection__iter__

# RooAbsCollection decided to overload operator= for assignment :(
def _RooAbsCollection_assign(self, other):
    if self == other:
        return
    for item in self:
        ritem = other.find(item)
        if not ritem:
            continue
        item.setVal(ritem.getVal())
        item.setError(ritem.getError())
        item.setAsymError(ritem.getErrorLo(), ritem.getErrorHi())
        item.setAttribute("Constant", ritem.isConstant())

ROOT.RooAbsCollection.assign = _RooAbsCollection_assign


def dump_templates(args):
    fws, wsname = args.workspace.split(':')
    fin = ROOT.TFile.Open(fws)
    w = fin.Get(wsname)
    model = w.obj(args.model)
    pdf = model.GetPdf()
    observable = model.GetObservables()[args.observable]
    params = pdf.getParameters(model.GetObservables())

    if args.fit:
        ffit, fitname = args.fit.split(':')
        fin2 = ROOT.TFile.Open(ffit)
        fit = fin2.Get(fitname)
        params.assign(fit.constPars())
        params.assign(fit.floatParsFinal())

    param_override = {p[0]:float(p[1]) for p in (s.split("=") for s in args.setparams.split(",")) if len(p)==2}
    for k,v in param_override.items():
        par = params.find(k)
        if not par:
            raise ValueError("No such floating parameter '%s' in this model" % k)
        par.setVal(v)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    indexCat = pdf.indexCat()
    for i in range(indexCat.numBins(args.rangeName)):
        indexCat.setBin(i, args.rangeName)
        cat_pdf = pdf.getPdf(indexCat.getLabel())

        for shape_pdf in cat_pdf.pdfList():
            if shape_pdf.dependsOn(observable) and isinstance(shape_pdf, ROOT.RooAddPdf):
                for norm, shape in zip(shape_pdf.coefList(), shape_pdf.pdfList()):
                    if args.filter is None or args.filter.match(shape.GetName()):
                        # Who needs comments, right? Anyway...
                        # Shape here has to be some kind of PDF.
                        # Notice that createHistogram() is *NOT* a method of RooAbsPdf...
                        # What we can do instead is generatedBinned() a RooDataHist and 
                        # and createHistogram() from that instead
                        #hist = shape.createHistogram("hist_"+shape.GetName(), observable)

                        # We generate 1E6 events, set expectedData = true (i.e. the Asimov dataset)
                        # and extended = false (i.e. we don't fluctuate bin by bin).
                        # The 1E6 should be immaterial here, since we normalise the TH1 back to
                        # the normalisation variable.
                        ras = ROOT.RooArgSet(observable)
                        rdh = shape.generateBinned(ras, 1E6, True, False);
                        hist = rdh.createHistogram("hist_"+shape.GetName(), observable)
                        hist.Scale(1.0/hist.Integral()) # scale back to unit area
                        hist.Scale(norm.getVal()) # and normalise to whatever we want
                        # Great, now continue as it was before
                        
                        c = ROOT.TCanvas(shape.GetName())
                        hist.Draw("hist")
                        c.Print("%s/%s.pdf" % (args.output, shape.GetName()))
                        fout = ROOT.TFile.Open("%s/%s.root" % (args.output, shape.GetName()), "recreate")
                        hist.Write()
                        fout.Close()
            elif shape_pdf.dependsOn(observable) and isinstance(shape_pdf, ROOT.RooRealSumPdf):
                raise NotImplementedError


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Dump all templates from workspace")
    parser.add_argument("-w", "--workspace", metavar="ROOTFILE:WORKSPACE", help="Workspace to load (e.g. output of text2workspace.py)", required=True)
    parser.add_argument("-f", "--fit", metavar="ROOTFILE:FIT_NAME", help="Fit result to load")
    parser.add_argument("-m", "--model", help="Model to load (usually 'ModelConfig' for signal+bkg, 'ModelConfig_bonly' for just background)", default="ModelConfig")
    parser.add_argument("-o", "--output", help="Output directory", default="plots")
    parser.add_argument("--observable", help="Name of the observable on which the shape templates depend (default 'x' in combine)", default="x")
    parser.add_argument("--rangeName", help="Category binning range name (usually blank)", default="")
    parser.add_argument("--filter", help="Regular-expression filter on output shapes", default=None, type=re.compile)
    parser.add_argument("--setparams", help="Set comma-separated list of param=value to override any floating parameters", default="", type=str)

    args = parser.parse_args()
    ret = dump_templates(args)

