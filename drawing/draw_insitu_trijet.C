#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_insitu_trijet() {
  drawer d;
  TFile * fgj = TFile::Open("/home/samson72/sphnx/gammajet/hists/insitu.root");
  TFile * ftj = TFile::Open("/home/samson72/sphnx/gammajet/hists/outputFileFinal.root");
  TH1D * h = (TH1D*)fgj->Get("hdivide1");
  //h->GetXaxis()->SetRangeUser(10.0,60.0);
  h->GetYaxis()->SetRangeUser(0.9,1.3);
  TGraphErrors * insitu_gj = new TGraphErrors(h);
  TGraphErrors * insitu_tj;
  
  ftj->GetObject("FORsam/TrijetRATIOXJ",insitu_tj);

  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetLeftMargin(.15);
  TLegend * l = new TLegend(.2,.6,.5,.85);
  l->AddEntry(insitu_gj, "#gamma-Jet #it{in situ}");
  l->AddEntry(insitu_tj, "trijet #it{in situ}");
  l->SetLineWidth(0);

  insitu_gj->SetLineColor(kBlack);
  insitu_gj->SetMarkerColor(kBlack);
  insitu_tj->SetLineColor(  kRed);
  insitu_tj->SetMarkerColor(kRed);

  insitu_tj->GetXaxis()->SetLimits(5.0, 75.0);;
  insitu_tj->GetYaxis()->SetRangeUser(0.9,1.3);
  

  int n = 16;
  double x[n];
  double y[n];
  double ex[n];
  double ey[n];
  for (int i = 0; i < insitu_gj->GetN(); i++) {
    Double_t ix, iy;
    insitu_gj->GetPoint(i, ix, iy);
    x[i]=ix;
    y[i]=iy;
    ex[i] = insitu_gj->GetErrorX(i);
    ey[i] = insitu_gj->GetErrorY(i);
  }
  for (int i = 0; i < insitu_tj->GetN(); i++) {
    Double_t ix, iy;
    insitu_tj->GetPoint(i, ix, iy);
    x[i+9]=ix;
    y[i+9]=1.0/iy;
    ex[i+9] = insitu_tj->GetErrorX(i);
    ey[i+9] = insitu_tj->GetErrorY(i);
  }
  TGraphErrors *mg = new TGraphErrors(n,x,y,ex,ey);

  mg->SetMarkerColor(kRed);
  mg->SetLineColor(kRed);
  mg->SetMarkerStyle(20);
  mg->SetMarkerSize(1);

  TF1 * func = new TF1("func", "pol2", 13.0, 70.0);
  func->SetParameter(0, 1);
  func->SetParameter(1, 0);
  func->SetParameter(2, 0);
  TF1 * func_gj = new TF1("func_gj", "pol2", 13.0, 35.0);
  func_gj->SetParameter(0, 1);
  func_gj->SetParameter(1, 0);
  func_gj->SetParameter(2, 0);
  mg->Fit(func,"0");
  TGraphErrors *band = new TGraphErrors(n);
  for (int i = 0; i < n; i++) {
    double x = 13.0 + (70.0 - 13.0) * i / (n - 1);
    band->SetPoint(i, x, func->Eval(x));
  }
  TVirtualFitter::GetFitter()->GetConfidenceIntervals(band);
  band->SetFillColor(kGreen);
  band->SetFillStyle(3001);
  
  
  int n_gj = insitu_gj->GetN();
  insitu_gj->Fit(func_gj, "0", "", 13.0, 35.0);
  TGraphErrors *band_gj = new TGraphErrors(n_gj);
  for (int i = 0; i < n_gj; i++) {
    double x = 13.0 + (35.0 - 13.0) * i / (n_gj - 1);
    band_gj->SetPoint(i, x, func_gj->Eval(x));
  }
  TVirtualFitter::GetFitter()->GetConfidenceIntervals(band_gj);
  band_gj->SetFillColorAlpha(kBlue,0.2);
  band_gj->SetFillStyle(3001);
  func_gj->SetLineColor(kBlue);

  //insitu_tj->Draw("ap");
  mg->Draw("ap");
  insitu_gj->Draw("p same");
  func->Draw("same");
  //func_gj->Draw("same");
  //band_gj->Draw("3 same");
  band->Draw("3 same");
  mg->Draw("p same");
  insitu_gj->Draw("p same");
  

  l->Draw("same");


  d.drawAll({"MC Photon"},
      {"Jet R=0.4","|vz| < 60 cm", "|#eta| < 0.7", Form("%.3f + %.3fp_{T} + %.5fp_{T}^{2}",func->GetParameter(0),func->GetParameter(1),func->GetParameter(2))},
      .52, .85, 20, c->GetWh()); 

}
