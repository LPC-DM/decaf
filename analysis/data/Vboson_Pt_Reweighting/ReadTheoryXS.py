import sys
import os
import ROOT
from array import array

filename = sys.argv[1]

binning = [
    30,
    40,
    50,
    60,
    70,
    80,
    90,
    100,
    110,
    120,
    130,
    140,
    150,
    200,
    250,
    300,
    350,
    400,
    450,
    500,
    550,
    600,
    650,
    700,
    750,
    800,
    850,
    900,
    950,
    1000,
    1100,
    1200,
    1300,
    1400,
    1600,
    1800,
    2000,
    2200,
    2400,
    2600,
    2800,
    3000,
    6500,
]
# v_boson_pt_hist = ROOT.TH1D(boson+"_boson_pt",boson+"_boson_pt",len(binning)-1,array('d',binning))
rootfile = ROOT.TFile(filename.replace(".dat", "") + ".root", "RECREATE")
histo = None
with open(filename) as f:
    for line in f:
        line_as_list = line.split(" ")
        if "BEGIN HISTO1D" in line:
            # line_as_list = line.split(' ')
            histname = line_as_list[3].replace("\n", "")
            histo = ROOT.TH1D(histname, histname, len(binning) - 1, array("d", binning))
            histo.Sumw2()
            continue
        if "END HISTO1D" in line:
            print (histo.Integral())
            rootfile.WriteTObject(histo)
            continue
        if "xlow" in line:
            continue
        if (
            "BEGIN DATA" in line
            or "Process" in line
            or "Proc" in line
            or "Sqrt" in line
            or "Unit" in line
            or "Timestamp" in line
            or "END DATA" in line
            or line == "\n"
        ):
            continue
        print ("line: ", line)
        lower_bin_edge = float(line_as_list[0])
        print ("lower bin edge: ", lower_bin_edge)
        bin_number = histo.FindBin(lower_bin_edge)
        print ("bin number: ", bin_number)
        bin_content = float(line_as_list[2])
        print ("bin content: ", bin_content)
        bin_error = float(line_as_list[3])
        print ("bin error: ", bin_error)
        bin_width = histo.GetBinWidth(bin_number)
        print ("bin width: ", bin_width)
        if "K_" in histname or "kappa_" in histname:
            histo.SetBinContent(bin_number, bin_content)
            # histo.SetBinError(bin_number,0.)
            histo.SetBinError(bin_number, bin_error)
        else:
            histo.SetBinContent(bin_number, bin_content * bin_width)
            # histo.SetBinError(bin_number,0.)
            histo.SetBinError(bin_number, bin_error * bin_width)
rootfile.Close()
