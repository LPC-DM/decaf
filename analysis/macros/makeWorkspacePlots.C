/// Utilities
std::vector<std::string> split(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

/// Colors
std::vector<std::string> majorNames = {"tt", "wjets", "zjets", "dyjets"};
std::vector<int> majorColors = {TColor::GetColor(43, 147, 34),
                                TColor::GetColor(164, 219, 120),
                                TColor::GetColor(37, 14, 210),
                                TColor::GetColor(251, 178, 93)};
std::vector<std::string> minorNames = {"qcdMC", "gjets", "hbbMC", "dyjetsMC", "vvMC", "stMC"};
std::vector<int> minorColors = {TColor::GetColor(150, 195, 220),
                                TColor::GetColor(193, 193, 193),
                                TColor::GetColor(190, 160, 204),
                                TColor::GetColor(251, 178, 93),
                                TColor::GetColor(218, 0, 24),
                                TColor::GetColor(63, 237, 30)};

/// Legend
void makeLegend() {
  TCanvas *cv = new TCanvas("cv", "cv", 600, 600);
  TH1F *histo = new TH1F("histo", "histo", 1, 0, 300);
  histo->SetLineWidth(0);
  histo->GetXaxis()->SetTitle("fjmass [GeV]");
  histo->GetYaxis()->SetTitle("Events");

  for (size_t i = 0; i != majorNames.size(); ++i) {
    histo->SetFillColor(majorColors.at(i));
    histo->SetName(majorNames.at(i).c_str());
    histo->SetTitle(majorNames.at(i).c_str());
    if (i == 0)
      histo->DrawClone();
    else
      histo->DrawClone("SAME");
  }

  for (size_t i = 0; i != minorNames.size(); ++i) {
    histo->SetFillColor(minorColors.at(i));
    histo->SetName(minorNames.at(i).c_str());
    histo->SetTitle(minorNames.at(i).c_str());
    histo->DrawClone("SAME");
  }

  TLegend *leg = cv->BuildLegend(0.2, 0.2, 0.9, 0.9);
  leg->SetBorderSize(0);
  cv->SaveAs("legend.pdf");
  cv->Draw();
}

/// Make one single plot from the workspace
void makeOnePlot(std::string wsname, std::string name, TCanvas *cv, int extra) {
  using namespace RooFit;
  std::string fileName = wsname + ".root";
  TFile *f = TFile::Open(fileName.c_str());
  RooWorkspace *w = (RooWorkspace *)f->Get(wsname.c_str());
  RooRealVar *fjmass = w->var("fjmass");
  RooAbsPdf *thePdf = nullptr;
  RooDataHist *theData = nullptr;
  thePdf = w->pdf(name.c_str());
  // theData = (RooDataHist*)w->data((name+"_observation").c_str());
  // thePdf->Print();
  // theData->Print();
  // thePdf->fitTo(*theData);

  RooPlot *xframe = fjmass->frame(Title("Model and data read from workspace"));

  // theData->plotOn(xframe);
  thePdf->plotOn(xframe, LineColor(kRed + extra));

  // Cosmetics
  xframe->Draw("SAME");
  // cv->SetLogy();
  cv->Draw();
  // TLatex lt; lt.DrawLatex(100,1E5,(name).c_str());
  // cv->SaveAs((name+std::string(".pdf")).c_str());
  // f->Close();
}

/// Stack all plots from the workspace
std::vector<TH1 *> makeComponentPlot(std::string fileName, std::string wsname, std::string catName) {
  using namespace RooFit;
  TFile *f = TFile::Open(fileName.c_str());
  RooWorkspace *w = (RooWorkspace *)f->Get(wsname.c_str());
  RooRealVar *fjmass = w->var("fjmass");
  RooAbsReal *theVar = nullptr;
  RooAbsPdf *thePdf = nullptr;
  RooDataHist *theData = nullptr;

  std::vector<TH1 *> allHistos;

  // First, get the list of PDFs
  const RooArgSet &listOfPDFs = w->allPdfs();
  // And the string that describe their contents
  const std::string listOfPDFNames = listOfPDFs.contentsString();
  // This should be easier...

  // First, get the minor MCs
  int colorIndex = 0;
  for (auto s = minorNames.begin(); s != minorNames.end(); ++s, ++colorIndex) {
    //std::string histoName = "shapeBkg_" + catName + "_" + *s + "_morph";
    std::string histoName = wsname + "_" + *s;

    //std::string normName = "n_exp_final_bin" + catName + "_proc_" + *s;
    std::string normName = histoName + "_norm";

    std::cout << "Getting histogram " << histoName << std::endl;

    thePdf = w->pdf(histoName.c_str());
    theVar = w->function(normName.c_str());
    if (thePdf == nullptr)
      continue;

    // thePdf->Print();
    TH1 *histo = thePdf->createHistogram(s->c_str(), *fjmass);
    histo->Scale(theVar->getValV());
    histo->SetLineWidth(1);
    histo->SetLineColor(kBlack);
    histo->SetFillColor(minorColors.at(colorIndex));
    allHistos.push_back(histo);
  }

  // Now the big  MCs

  colorIndex = 0;
  for (auto s = majorNames.begin(); s != majorNames.end(); ++s, ++colorIndex) {
    //std::string histoName = "shapeBkg_" + *s + "_" + catName;
    std::string histoName = wsname + "_" + *s;

    //std::string normName = "n_exp_final_bin" + catName + "_proc_" + *s;
    std::string normName = histoName + "_norm";
    std::cout << "Getting histogram " << histoName << std::endl;

    thePdf = w->pdf(histoName.c_str());
    if (thePdf == nullptr)
      continue;

    TH1 *histo = thePdf->createHistogram(s->c_str(), *fjmass);
    histo->Scale(w->function(normName.c_str())->getValV());
    histo->SetLineWidth(1);
    histo->SetLineColor(kBlack);
    histo->SetFillColor(majorColors.at(colorIndex));
    allHistos.push_back(histo);
  }

  return allHistos;
}

void makeWorkspacePlots(const std::string name, int histoTop, int latexY) {
  TCanvas *cv = new TCanvas("cv", "cv", 600, 600);
  TH1F *histo = new TH1F("histo", "histo", 1, 0, 300);
  histo->SetLineWidth(0);
  histo->GetXaxis()->SetTitle("fjmass [GeV]");
  histo->GetYaxis()->SetTitle("Events");
  histo->Draw();
  histo->GetYaxis()->SetRangeUser(0.1, histoTop);
  cv->SetLogy();

  //std::string name = "wmcr2018failrecoil3";
  std::string fileName = name + ".root";
  std::string pdfName = name + ".pdf";

  std::vector<TH1 *> allHistos = makeComponentPlot(fileName.c_str(), name.c_str(), "");
  THStack *hs = new THStack("hs", "Stacked 1D histograms");

  std::cout << "For this category  " << pdfName << " we have:";

  for (auto h = allHistos.begin(); h != allHistos.end(); ++h) {
    std::cout << (*h)->GetTitle() << " with integral " << (*h)->Integral() << "; ";
    (*h)->SetLineWidth(1);
    hs->Add(*h);
  }
  std::cout << std::endl;
  // makeOnePlot("wmcr2018failrecoil1","wmcr2018failrecoil1_wjets",cv,1);
  // makeOnePlot("wmcr2018failrecoil2","wmcr2018failrecoil2_wjets",cv,2);
  // makeOnePlot("wmcr2018failrecoil3","wmcr2018failrecoil3_wjets",cv,3);
  // makeOnePlot("wmcr2018failrecoil4","wmcr2018failrecoil4_wjets",cv,4);

  hs->Draw("HIST SAME");
  cv->RedrawAxis();
  TLatex lt;
  lt.DrawLatex(100, latexY, (name).c_str());

  cv->SaveAs(pdfName.c_str());
}
