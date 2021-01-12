void plotsDist()
{
//=========Macro generated from canvas: asdf/asdf
//=========  (Wed Jan 13 03:59:58 2021) by ROOT version 6.12/07
   TCanvas *asdf = new TCanvas("asdf", "asdf",0,0,800,800);
   gStyle->SetOptFit(1);
   asdf->SetHighLightColor(2);
   asdf->Range(-2.5,-0.6210817,2.5,5.589735);
   asdf->SetFillColor(0);
   asdf->SetBorderMode(0);
   asdf->SetBorderSize(2);
   asdf->SetFrameBorderMode(0);
   asdf->SetFrameBorderMode(0);
   
   TH1F *pulls__1 = new TH1F("pulls__1","",60,-2,2);
   pulls__1->SetBinContent(0,2);
   pulls__1->SetBinContent(8,2);
   pulls__1->SetBinContent(21,1);
   pulls__1->SetBinContent(22,1);
   pulls__1->SetBinContent(24,1);
   pulls__1->SetBinContent(25,1);
   pulls__1->SetBinContent(26,2);
   pulls__1->SetBinContent(27,1);
   pulls__1->SetBinContent(28,1);
   pulls__1->SetBinContent(29,3);
   pulls__1->SetBinContent(30,3);
   pulls__1->SetBinContent(32,2);
   pulls__1->SetBinContent(33,1);
   pulls__1->SetBinContent(34,2);
   pulls__1->SetBinContent(35,1);
   pulls__1->SetBinContent(36,2);
   pulls__1->SetBinContent(37,1);
   pulls__1->SetBinContent(38,1);
   pulls__1->SetBinContent(39,2);
   pulls__1->SetBinContent(40,1);
   pulls__1->SetBinContent(61,1);
   pulls__1->SetEntries(32);
   pulls__1->SetStats(0);
   
   TF1 *gaus1 = new TF1("gaus","gaus",-2,2, TF1::EAddToList::kNo);
   gaus1->SetFillColor(19);
   gaus1->SetFillStyle(0);
   gaus1->SetLineWidth(3);
   gaus1->SetChisquare(4.397062);
   gaus1->SetNDF(16);
   gaus1->SetParameter(0,1.328635);
   gaus1->SetParError(0,0.6281243);
   gaus1->SetParLimits(0,0,0);
   gaus1->SetParameter(1,-1.112949);
   gaus1->SetParError(1,10.91497);
   gaus1->SetParLimits(1,0,0);
   gaus1->SetParameter(2,5.15978);
   gaus1->SetParError(2,4.323747);
   gaus1->SetParLimits(2,0,5.15978);
   gaus1->SetParent(pulls__1);
   pulls__1->GetListOfFunctions()->Add(gaus1);
   pulls__1->SetMarkerStyle(20);
   pulls__1->SetMarkerSize(2);
   pulls__1->GetXaxis()->SetTitle("pull");
   pulls__1->GetYaxis()->SetTitle("Number of nuisances");
   pulls__1->Draw("pe");
   asdf->Modified();
   asdf->cd();
   asdf->SetSelected(asdf);
}
