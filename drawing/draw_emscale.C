#include "../headers/ana.cxx"


void draw_emscale() { 
  int ipt = 5;
  int ij = 2;
  ana anaclone;
  TFile * f1 = TFile::Open("hists_og/histsData.root");
  TFile * f2 = TFile::Open("hists/histsData.root");
  TH1D * h1 = (TH1D*)f1->Get(Form("hratio_%i_%i_0_0_0",ipt,ij));
  TH1D * h2 = (TH1D*)f2->Get(Form("hratio_%i_%i_0_0_0",ipt,ij));

  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);

  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2);
  h1->Rebin(4000);
  h1->GetYaxis()->SetRangeUser(0,h1->GetMaximum()*1.2);
  h2->SetLineColor(kRed);
  h2->SetLineWidth(2);
  h2->Rebin(4000);
  h1->GetXaxis()->SetTitle("x_{J#gamma}");
  h1->GetYaxis()->SetTitle("Counts");

  TLegend * l = new TLegend(.5,.45,.85,.55);
  l->SetLineWidth(0);
  l->AddEntry(h1,"Unscaled EMCal");
  l->AddEntry(h2,"EMCal scaled 10%");

  h1->Draw();
  h2->Draw("same");
  l->Draw("same");
  
  TF1 * func1 = new TF1("func1", "gaus",(int)(anaclone.minjete[ij]/anaclone.ptBins[ipt]/0.04 + 1)*0.04,2);
  TF1 * func2 = new TF1("func2", "gaus",(int)(anaclone.minjete[ij]/anaclone.ptBins[ipt]/0.04 + 1)*0.04,2);
  func1->SetLineColor(kBlue+1);
  func2->SetLineColor(kRed+1);
  h1->Fit(func1,"RQI");
  h2->Fit(func2,"RQI");
  
  anaclone.drawAll({
      "Data",
      "Uncalibrated Jets"
      },{
      Form("Jet R = %0.1f",anaclone.JetRs[ij]),
      Form("%0.0f GeV < cluster p_{T} < %0.0f GeV",anaclone.ptBins[ipt],anaclone.ptBins[ipt+1]),
      "Analysis Cuts",
      Form("unscaled mean = %1.3f",func1->GetParameter(1)),
      Form("scaled mean = %1.3f",func2->GetParameter(1))
      }, 0.5, 0.85, 20, 0);

  cout << (func1->GetParameter(1) - func2->GetParameter(1))/func1->GetParameter(1) << endl;


}
