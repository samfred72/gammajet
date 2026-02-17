#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include "/home/samson72/sphnx/gammajet/src/drawer.h"

int colors[6] = {kRed, kBlue, kGreen-2, kMagenta, kOrange, kTeal};
drawer d;

void scale(TH1D * h) {
  double x_low_value = 0;
  double x_high_value = 100;

  // Get the X-axis object
  TAxis *xaxis = h->GetXaxis();

  // Find the corresponding bin numbers
  Int_t bin_low = xaxis->FindBin(x_low_value);
  Int_t bin_high = xaxis->FindBin(x_high_value);

  // Sum the bin contents within the range
  double entries_in_range = 0;
  for (Int_t bin = bin_low; bin <= bin_high; ++bin) {
    entries_in_range += h->GetBinContent(bin);
  }
  h->Scale(1.0/entries_in_range);
}

void draw_isobdt() {
  gStyle->SetOptStat(0);
  TFile * f = TFile::Open(   "/home/samson72/sphnx/gammajet/hists/histsData.root");
  TFile * f05p = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsPhoton5.root");
  TFile * f10p = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsPhoton10.root");
  TFile * f20p = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsPhoton20.root");
  TFile * f05j = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsJet5.root");
  TFile * f10j = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsJet10.root");
  TFile * f20j = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsJet20.root");
  TFile * f30j = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsJet30.root");
  TFile * f50j = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsJet50.root");
  TFile * f70j = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsJet70.root");
  
  float drawx = .45;
  float drawy = .55;
  float fontsize = 30;
  int ptbin = 5;

  const char * histname = Form("hisobdt%i",ptbin);
  TH2D * hdata = (TH2D*)f->Get(histname);
  TH2D * hmcp = d.combineMC2d(histname,1);
 
  TLine * l1 = new TLine(ana::isoBins[0], ana::bdtBins[0], ana::isoBins[0],1);
  TLine * l2 = new TLine(-1,1,ana::isoBins[0],1);
  TLine * l3 = new TLine(-1,1,-1,ana::bdtBins[0]);
  TLine * l4 = new TLine(-1,ana::bdtBins[0],ana::isoBins[0],ana::bdtBins[0]);
  
  TLine * l5 = new TLine(-1,ana::bdtCuts[0],ana::isoBins[0],ana::bdtCuts[0]);
  TLine * l6 = new TLine(ana::isoBins[0],ana::bdtCuts[0],ana::isoBins[0],ana::bdtCutsHigh[0]);
  TLine * l7 = new TLine(ana::isoBins[0],ana::bdtCutsHigh[0],-1,ana::bdtCutsHigh[0]);
  TLine * l8 = new TLine(-1,ana::bdtCutsHigh[0],-1,ana::bdtCuts[0]);
  
  TLine * l9 = new TLine(ana::isoBinsHigh[0],ana::bdtCuts[0],20,ana::bdtCuts[0]);
  TLine * l10 = new TLine(20,ana::bdtCuts[0],20,ana::bdtCutsHigh[0]);
  TLine * l11 = new TLine(20,ana::bdtCutsHigh[0],ana::isoBinsHigh[0],ana::bdtCutsHigh[0]);
  TLine * l12 = new TLine(ana::isoBinsHigh[0],ana::bdtCutsHigh[0],ana::isoBinsHigh[0],ana::bdtCuts[0]);
  
  TLine * l13 = new TLine(ana::isoBinsHigh[0],ana::bdtBins[0],20,ana::bdtBins[0]);
  TLine * l14 = new TLine(20,ana::bdtBins[0],20,1);
  TLine * l15 = new TLine(20,1,ana::isoBinsHigh[0],1);
  TLine * l16 = new TLine(ana::isoBinsHigh[0],1,ana::isoBinsHigh[0],ana::bdtBins[0]);

  const int nlines = 16;
  TLine * lines[nlines] = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16};
  for (int i = 0; i < nlines; i++) {
    lines[i]->SetLineColor(kRed);
    lines[i]->SetLineWidth(3);
  }
  TCanvas * c = new TCanvas("c","",1400,700);
  c->Divide(2,1,0,0);
  c->cd(1);
  gPad->SetLogz();
  hdata->Scale(1.0/hdata->Integral());
  hdata->GetZaxis()->SetRangeUser(1e-5,1e-2);
  hdata->Draw("col");
  for (int i = 0; i < nlines; i++) {
    lines[i]->Draw("same");
  }
  d.drawAll({"Data"},{Form("%0.0f GeV < p_{T} < %0.0f GeV",ana::ptBins[ptbin],ana::ptBins[ptbin+1]),"Paired clusters"},drawx,drawy,fontsize,700);
  d.drawText("A",.15,.90, kBlack, 25);
  d.drawText("B",.15,ana::bdtCutsHigh[0]-.05, kBlack, 25);
  d.drawText("C",.35,.90, kBlack, 25);
  d.drawText("D",.35,ana::bdtCutsHigh[0]-.05, kBlack, 25);
  
  c->cd(2);
  gPad->SetRightMargin(.15);
  gPad->SetLogz();
  hmcp->Scale(1.0/hmcp->Integral());
  hmcp->GetZaxis()->SetRangeUser(1e-5,1e-2);
  hmcp->Draw("colz");
  for (int i = 0; i < nlines; i++) {
    lines[i]->Draw("same");
  }
  //d.drawAll({"MC Photon"},{Form("%0.0f GeV < p_{T} < %0.0f GeV",ana::ptBins[ptbin],ana::ptBins[ptbin+1]),"Paired clusters"},drawx-0.2,drawy,fontsize,700);
  d.drawText("MC Photon",drawx,drawy-0.05,kBlack,(int)(fontsize*1.25));
  d.drawText("A",.05,.90, kBlack, 25);
  d.drawText("B",.05,ana::bdtCutsHigh[0]-.05, kBlack, 25);
  d.drawText("C",.25,.90, kBlack, 25);
  d.drawText("D",.25,ana::bdtCutsHigh[0]-.05, kBlack, 25);

  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/isobdt.pdf");

}
