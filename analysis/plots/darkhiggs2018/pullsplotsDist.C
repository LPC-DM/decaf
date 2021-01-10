void pullsplotsDist()
{
//=========Macro generated from canvas: asdf/asdf
//=========  (Sun Jan 10 10:03:53 2021) by ROOT version 6.12/07
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
   pulls__1->SetBinContent(9,2);
   pulls__1->SetBinContent(18,1);
   pulls__1->SetBinContent(23,1);
   pulls__1->SetBinContent(24,1);
   pulls__1->SetBinContent(25,1);
   pulls__1->SetBinContent(27,3);
   pulls__1->SetBinContent(29,2);
   pulls__1->SetBinContent(30,2);
   pulls__1->SetBinContent(31,3);
   pulls__1->SetBinContent(32,1);
   pulls__1->SetBinContent(33,1);
   pulls__1->SetBinContent(34,1);
   pulls__1->SetBinContent(35,1);
   pulls__1->SetBinContent(36,1);
   pulls__1->SetBinContent(37,1);
   pulls__1->SetBinContent(38,3);
   pulls__1->SetBinContent(39,2);
   pulls__1->SetBinContent(44,1);
   pulls__1->SetBinContent(50,1);
   pulls__1->SetBinContent(61,1);
   pulls__1->SetEntries(32);
   pulls__1->SetStats(0);
   
   TF1 *gaus1 = new TF1("gaus","gaus",-2,2, TF1::EAddToList::kNo);
   gaus1->SetFillColor(19);
   gaus1->SetFillStyle(0);
   gaus1->SetLineWidth(3);
   gaus1->SetChisquare(4.861075);
   gaus1->SetNDF(16);
   gaus1->SetParameter(0,1.474605);
   gaus1->SetParError(0,1.608511);
   gaus1->SetParLimits(0,0,0);
   gaus1->SetParameter(1,-3.105426);
   gaus1->SetParError(1,11.87369);
   gaus1->SetParLimits(1,0,0);
   gaus1->SetParameter(2,5.857344);
   gaus1->SetParError(2,3.284432);
   gaus1->SetParLimits(2,0,5.857515);
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
