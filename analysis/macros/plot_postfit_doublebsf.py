from ROOT import *
from collections import defaultdict
from os import getenv
from array import array
from tdrStyle import *
from optparse import OptionParser
from os import system
import sys
import plotConfig
setTDRStyle()

plotextralabel = ''
PRELIM = True
new_dic = defaultdict(dict)

def getInt(h):
    nbins = h.GetNbinsX()
    total=0.
    for iB in range(1,nbins+1):
        total += h.GetBinContent(iB)*h.GetBinWidth(iB)
    return total

def plotPreFitPostFit(year, signalflag, blind=False):

    #### Define year, and category ####
    years = {
            "2016": "2016",
            "2017": "2017",
            "2018": "2018"
    }

    lumi = {
            "2016": 35.92,
            "2017": 41.53,
            "2018": 59.74
    }

    signalprocess = {
            "pass": "pass",
            "fail": "fail"
    }

    #extralabels = {
    #        "": "",
    #        "": "",
    #        "": "",
    #        "": "",
    #        "": "",
    #}

    #extralabel = extralabels[region]+' '+signalprocess[signalflag]+' '+recoilbin[recoil]

    ### Open fit result file ###
    f_mlfit = TFile(sys.argv[1], 'READ')

    ### Make printing which region, year, signal, and recoil are considered ###
    print('-----------------------------------------------------------')
    print('Considering year is:', years[year])
    print('Signal process is:', signalprocess[signalflag])
    print('Opening', "svmass"+years[year]+signalprocess[signalflag])
    print('-----------------------------------------------------------')

    ### Check directory ###
    h_dir = f_mlfit.Get("svmass"+years[year]+signalprocess[signalflag]+"_postfit")

    if h_dir == None:
        raise RuntimeError

    ### Call data ###
    h_data = f_mlfit.Get("svmass"+years[year]+signalprocess[signalflag]+"_postfit/data_obs")

    #### Call signal pre/post fit histograms ###
    #h_postfit_sig = f_mlfit.Get("shapes_fit_b/"+darkhiggs_regions['signal']+years[year]+signalprocess[signalflag]+recoilbin[recoil]+"/total_background")
    #h_prefit_sig = f_mlfit.Get("shapes_prefit/"+darkhiggs_regions['signal']+years[year]+signalprocess[signalflag]+recoilbin[recoil]+"/total_background")

    ### Define processes per region ###
    processes = [
            'genb',
            'genbb',
            'genc',
            'gencc',
            'genother',
    ]

    order = [
            'b#bar{b}',
            'b',
            'c#bar{c}',
            'c',
            'Data',
    ]

    processNames = {
            'genb': 'b',
            'genbb': 'b#bar{b}',
            'genc': 'c',
            'gencc': 'c#bar{c}',
            'genother': 'light flavor',
    }

    colors = {
        'genbb':TColor.GetColor(166, 86, 40),
        'genb':TColor.GetColor(152, 78, 163),
        'gencc':TColor.GetColor(247, 129, 191),
        'genc':TColor.GetColor(255, 127, 0),
        'genother':TColor.GetColor(77, 154, 74),
    }

    binLowE = []

    # Pre-Fit
    h_prefit = {}
    h_prefit['total'] = f_mlfit.Get("svmass"+years[year]+signalprocess[signalflag]+"_prefit/TotalProcs")
    for i in range(1,h_prefit['total'].GetNbinsX()+2):
        binLowE.append(h_prefit['total'].GetBinLowEdge(i))

    h_all_prefit = TH1F("h_all_prefit","h_all_prefit",len(binLowE)-1,array('d',binLowE))
    h_other_prefit = TH1F("h_other_prefit","h_other_prefit",len(binLowE)-1,array('d',binLowE))
    h_stack_prefit = THStack("h_stack_prefit","h_stack_prefit")

    for process in processes:
        print('which process in prefit?', process)
        h_prefit[process] = f_mlfit.Get("svmass"+years[year]+signalprocess[signalflag]+"_prefit/"+process)
        if not h_prefit[process]:
            continue
        if str(h_prefit[process].Integral()) == "nan":
            continue
        h_prefit[process].SetLineColor(kBlack)
        h_prefit[process].SetFillColor(colors[process])
        h_all_prefit.Add(h_prefit[process])
        h_stack_prefit.Add(h_prefit[process])

    # Post-Fit
    h_postfit = {}
    h_postfit['total'] = f_mlfit.Get("svmass"+years[year]+signalprocess[signalflag]+"_postfit/TotalProcs")
    h_all_postfit = TH1F("h_all_postfit","h_all_postfit",len(binLowE)-1,array('d',binLowE))
    h_other_postfit = TH1F("h_other_postfit","h_other_postfit",len(binLowE)-1,array('d',binLowE))
    h_stack_postfit = THStack("h_stack_postfit","h_stack_postfit")

    h_postfit['totalv2'] = f_mlfit.Get("svmass"+years[year]+signalprocess[signalflag]+"_postfit/TotalBkg")

    for i in range(1, h_postfit['totalv2'].GetNbinsX()+1):
        error = h_postfit['totalv2'].GetBinError(i)
        content = h_postfit['totalv2'].GetBinContent(i)

    for process in processes:
        print('which process in postfit?', process)
        h_postfit[process] = f_mlfit.Get("svmass"+years[year]+signalprocess[signalflag]+"_postfit/"+process)
        if not h_postfit[process]:
            continue
        if str(h_postfit[process].Integral()) == "nan":
            continue
        h_postfit[process].SetLineColor(kBlack)
        h_postfit[process].SetFillColor(colors[process])
        h_all_postfit.Add(h_postfit[process])
        h_stack_postfit.Add(h_postfit[process])

    gStyle.SetOptStat(0)

    c = TCanvas("c","c",600,700)
    SetOwnership(c,False)
    c.cd()
    c.SetLogy()
    c.SetBottomMargin(0.3)
    c.SetRightMargin(0.06)

    dummy = h_all_prefit.Clone("dummy")
    dummy.SetFillColor(0)
    dummy.SetLineColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.SetMarkerColor(0)
    dummy.GetYaxis().SetTitle("Events / GeV")
    dummy.GetXaxis().SetTitle("")
    dummy.GetXaxis().SetTitleSize(0)
    dummy.GetXaxis().SetLabelSize(0)
    dummy.SetMaximum(100*dummy.GetMaximum())
    dummy.SetMinimum(1.)
    dummy.GetYaxis().SetTitleOffset(1.15)
    dummy.Draw()

    h_stack_postfit.Draw("hist same")

    h_all_prefit.SetLineColor(2)
    h_all_prefit.SetLineStyle(2)
    h_all_prefit.SetLineWidth(3)
    h_all_prefit.Draw("histsame")

    h_all_postfit.SetLineColor(4)
    h_all_postfit.SetLineWidth(3)
    h_all_postfit.Draw("histsame")

    h_all_prefit.SetLineWidth(3)
    h_all_postfit.SetLineWidth(3)

    if not blind:
        h_data.SetMarkerStyle(20)
        h_data.SetMarkerSize(1.2)
        h_data.Draw("epsame")

    legend = TLegend(.72,.74,.95,.92)
    legend.SetTextSize(0.018)
    yields = {}

    if not blind:
        legend.AddEntry(h_data,"Data","elp")
        yields['Data'] = getInt(h_data)
    legend.AddEntry(h_all_prefit, "SM total (pre-fit)", "l")
    legend.AddEntry(h_all_postfit, "SM total (post-fit)", "l")

    for process in reversed(processes):
        try:
            hist = h_postfit[process]
            if not h_postfit[process]:
                continue
            if str(h_postfit[process].Integral()) == "nan":
                continue
            legend.AddEntry(hist,processNames[process],"f")
            yields[processNames[process]] = getInt(hist)
        except KeyError:
            pass
    legend.SetShadowColor(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.Draw("same")

    #l1=region+' & '
    #for o in order:
    #    if o in yields:
    #        y = yields[o]
    #        if o=='Data':
    #            l1 += ' $%i$ & '%(int(y))
    #        else:
    #            l1 += ' $%.3g$ & '%y
    #    else:
    #        l1 += ' $-$ & '
    #print(l1)

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(42)
    #latex2.SetTextSize(0.6*c.GetTopMargin())
    latex2.SetTextSize(0.45*c.GetTopMargin())
    #latex2.DrawLatex(0.16, 0.83,extralabel)
    #if 'loose' in combinecat:
    #    latex2.DrawLatex(0.16,0.78,"0.1 < BDT < 0.45")
    #elif 'tight' in combinecat:
    #    latex2.DrawLatex(0.16,0.78,"BDT > 0.45")
    latex2.SetTextAlign(31) # align right
    latex2.SetTextSize(0.5*c.GetTopMargin())
    #latex2.DrawLatex(0.94, 0.94,"%.0f fb^{-1} (13 TeV)"%(plotConfig.lumi[year]))
    latex2.DrawLatex(0.94, 0.94,"%.0f fb^{-1} (13 TeV)"%(lumi[year]))
    #latex2.DrawLatex(0.9, 0.94,"2.32 fb^{-1} (13 TeV)")
    latex2.SetTextSize(0.6*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align left
    #latex2.DrawLatex(0.16, 0.87, "CMS")
    latex2.DrawLatex(0.16, 0.94, "CMS")
    latex2.DrawLatex(0.16, 0.87, signalflag)
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)

    if PRELIM:
        #latex2.DrawLatex(0.21, 0.94, "Preliminary")
        latex2.DrawLatex(0.26, 0.94, "Preliminary")

    gPad.RedrawAxis()

    pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 0.9)
    SetOwnership(pad,False)

    pad.SetTopMargin(0.7)
    pad.SetRightMargin(0.06)
    #pad.SetLeftMargin(0.18)
    pad.SetFillColor(0)
    pad.SetGridy(0)
    pad.SetFillStyle(0)
    pad.Draw()
    pad.cd(0)

    met = []; dmet = [];
    ratio_pre = []; ratio_pre_hi = []; ratio_pre_lo = [];
    ratio_post = []; ratio_post_hi = []; ratio_post_lo = [];

    for i in range(1,h_all_prefit.GetNbinsX()+1):
        #ndata = array("d", [0.0])
        #metave = array("d",[0.0])
        #h_data.GetPoint(i-1, metave[0], ndata[0])

        #ndata = h_data.GetY()[i-1]

        if not blind:
            ndata = h_data.GetBinContent(i)
        else:
            ndata = 0
        #print ndata

        if (ndata>0.0 and not blind):
            e_data_hi = h_data.GetBinError(i)/ndata
            e_data_lo = h_data.GetBinError(i)/ndata
        else:
            e_data_hi = 0.0
            e_data_lo = 0.0

        n_all_pre = h_all_prefit.GetBinContent(i)
        n_all_post = h_all_postfit.GetBinContent(i)

        met.append(h_all_prefit.GetBinCenter(i))
        dmet.append((h_all_prefit.GetBinLowEdge(i+1)-h_all_prefit.GetBinLowEdge(i))/2)

        if (n_all_pre>0.0):
            ratio_pre.append(ndata/n_all_pre)
            ratio_pre_hi.append(ndata*e_data_hi/n_all_pre)
            ratio_pre_lo.append(ndata*e_data_lo/n_all_pre)
        else:
            ratio_pre.append(0.0)
            ratio_pre_hi.append(0.0)
            ratio_pre_lo.append(0.0)

        if (n_all_post>0.0):
            ratio_post.append(ndata/n_all_post)
            ratio_post_hi.append(ndata*e_data_hi/n_all_post)
            ratio_post_lo.append(ndata*e_data_lo/n_all_post)
        else:
            ratio_post.append(0.0)
            ratio_post_hi.append(0.0)
            ratio_post_lo.append(0.0)

    a_met = array("d", met)
    v_met = TVectorD(len(a_met),a_met)

    a_dmet = array("d", dmet)
    v_dmet = TVectorD(len(a_dmet),a_dmet)

    a_ratio_pre = array("d", ratio_pre)
    a_ratio_pre_hi = array("d", ratio_pre_hi)
    a_ratio_pre_lo = array("d", ratio_pre_lo)

    v_ratio_pre = TVectorD(len(a_ratio_pre),a_ratio_pre)
    v_ratio_pre_hi = TVectorD(len(a_ratio_pre_hi),a_ratio_pre_hi)
    v_ratio_pre_lo = TVectorD(len(a_ratio_pre_lo),a_ratio_pre_lo)

    a_ratio_post = array("d", ratio_post)
    a_ratio_post_hi = array("d", ratio_post_hi)
    a_ratio_post_lo = array("d", ratio_post_lo)

    v_ratio_post = TVectorD(len(a_ratio_post),a_ratio_post)
    v_ratio_post_hi = TVectorD(len(a_ratio_post_hi),a_ratio_post_hi)
    v_ratio_post_lo = TVectorD(len(a_ratio_post_lo),a_ratio_post_lo)

    g_ratio_pre = TGraphAsymmErrors(v_met,v_ratio_pre,v_dmet,v_dmet,v_ratio_pre_lo,v_ratio_pre_hi)
    g_ratio_pre.SetLineColor(2)
    g_ratio_pre.SetMarkerColor(2)
    g_ratio_pre.SetLineStyle(2)
    g_ratio_pre.SetMarkerStyle(24)

    g_ratio_post = TGraphAsymmErrors(v_met,v_ratio_post,v_dmet,v_dmet,v_ratio_post_lo,v_ratio_post_hi)
    g_ratio_post.SetLineColor(4)
    g_ratio_post.SetMarkerColor(4)
    g_ratio_post.SetMarkerStyle(20)

    ratiosys = h_postfit['totalv2'].Clone();
    for hbin in range(0,ratiosys.GetNbinsX()+1):
        ratiosys.SetBinContent(hbin+1,1.0)
        if (h_postfit['totalv2'].GetBinContent(hbin+1)>0):
            ratiosys.SetBinError(hbin+1,h_postfit['totalv2'].GetBinError(hbin+1)/h_postfit['totalv2'].GetBinContent(hbin+1))
        else:
            ratiosys.SetBinError(hbin+1,0)

    dummy2 = TH1F("dummy2","dummy2",len(binLowE)-1,array('d',binLowE))
    for i in range(1,dummy2.GetNbinsX()):
        dummy2.SetBinContent(i,1.0)
        dummy2.GetYaxis().SetTitle("Data / Pred.")
        #dummy2.GetXaxis().SetTitle("b-tag jet probability")
        dummy2.GetXaxis().SetTitle("log(m_{SV1, d_{xy}sig max}/GeV)")
        dummy2.GetXaxis().SetTitleOffset(1.3)
        dummy2.GetXaxis().SetTitleSize(0.04)

    dummy2.SetLineColor(0)
    dummy2.SetMarkerColor(0)
    dummy2.SetLineWidth(0)
    dummy2.SetMarkerSize(0)
    dummy2.GetYaxis().SetLabelSize(0.03)
    dummy2.GetYaxis().SetNdivisions(5);
    dummy2.GetXaxis().SetNdivisions(510)
    dummy2.GetYaxis().CenterTitle()
    dummy2.GetYaxis().SetTitleSize(0.04)
    dummy2.GetYaxis().SetTitleOffset(1.5)
    dummy2.SetMaximum(2)
    dummy2.SetMinimum(0)
    dummy2.Draw("hist")

    ratiosys.SetFillColor(kGray) #SetFillColor(ROOT.kYellow)
    ratiosys.SetLineColor(kGray) #SetLineColor(1)
    ratiosys.SetLineWidth(1)
    ratiosys.SetMarkerSize(0)
    ratiosys.Draw("e2same")

    f1 = TF1("f1","1",-5000,5000)
    f1.SetLineColor(1)
    f1.SetLineStyle(2)
    f1.SetLineWidth(2)
    f1.Draw("same")

    if not blind:
        g_ratio_pre.Draw("epsame")
        g_ratio_post.Draw("epsame")
        legend2 = TLegend(.65,.25,.8,.29)
        legend3 = TLegend(.8,.25,.95,.29)
        legend2.AddEntry(g_ratio_pre,"pre-fit","elp")
        legend3.AddEntry(g_ratio_post,"post-fit","elp")
        for l in [legend2,legend3]:
            l.SetShadowColor(0)
            l.SetFillColor(0)
            l.SetFillStyle(0)
            l.SetBorderSize(0)
            l.SetLineColor(0)
            l.Draw()


    plotextralabel_ = plotextralabel
    if PRELIM:
        plotextralabel_ += '_prelim'
    #for ext in ['pdf','png','C']:
    for ext in ['png']:
        plotDir = 'plots/doublebsf'+years[year]+'/postfit/'+signalprocess[signalflag]
        print('Directory is',plotDir)
        system('mkdir -p '+plotDir)
        c.SaveAs(plotDir+'/stack_svmass.'+ext)

    #c.SaveAs("test.pdf")
    #del c
    #del process
    #del colors
    #del h_prefit

#### Call Functions ####
if __name__ == "__main__":
    sigs = ['pass', 'fail']

    parser = OptionParser()
    parser.add_option("-y", "--year", help="year", dest="year", default="")
    (options, args) = parser.parse_args()
    year = options.year

    for isig in sigs:
        try:
            #plotPreFitPostFit('2018', isig)
            plotPreFitPostFit(year, isig)
        except:
            print("Directory do not exist in %s file! \n" % (sys.argv[1]))
            pass
