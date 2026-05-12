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
  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/hists_Data.root");
  
  float drawx = .45;
  float drawy = .55;
  float fontsize = 30;
  int ptbin = 5;

  const char * histname = Form("hisobdt%i",ptbin);
  TH2D * hdata = (TH2D*)f->Get(histname);
  TH2D * hmcp = d.combineMC2d(histname,1);

  TBox * bdata = new TBox(-1, ana::bdtGoodLow[0], ana::isoBins[0], ana::bdtGoodHigh[0]);
  bdata->SetLineColor(kRed);
  bdata->SetFillStyle(0);
  bdata->SetLineWidth(4);


  float offset = 0.15; 
  TLine * l1 = new TLine(ana::isoBins[0], ana::bdtGoodLow[0], ana::isoBins[0], ana::bdtGoodHigh[0]);
  TLine * l2 = new TLine(-1-offset,ana::bdtGoodHigh[0],ana::isoBins[0]+offset,ana::bdtGoodHigh[0]);
  TLine * l3 = new TLine(-1,ana::bdtGoodHigh[0],-1,ana::bdtGoodLow[0]);
  TLine * l4 = new TLine(-1-offset,ana::bdtGoodLow[0],ana::isoBins[0]+offset,ana::bdtGoodLow[0]);
  
  TLine * l5 = new TLine(-1,ana::bdtBadLow[0],ana::isoBins[0],ana::bdtBadLow[0]);
  TLine * l6 = new TLine(ana::isoBins[0],ana::bdtBadLow[0],ana::isoBins[0],ana::bdtBadHigh[0]);
  TLine * l7 = new TLine(ana::isoBins[0],ana::bdtBadHigh[0],-1,ana::bdtBadHigh[0]);
  TLine * l8 = new TLine(-1,ana::bdtBadHigh[0],-1,ana::bdtBadLow[0]);
  
  TLine * l9 = new TLine(ana::isoBinsHigh[0],ana::bdtBadLow[0],20,ana::bdtBadLow[0]);
  TLine * l10 = new TLine(20,ana::bdtBadLow[0],20,ana::bdtBadHigh[0]);
  TLine * l11 = new TLine(20,ana::bdtBadHigh[0],ana::isoBinsHigh[0],ana::bdtBadHigh[0]);
  TLine * l12 = new TLine(ana::isoBinsHigh[0],ana::bdtBadHigh[0],ana::isoBinsHigh[0],ana::bdtBadLow[0]);
  
  TLine * l13 = new TLine(ana::isoBinsHigh[0],ana::bdtGoodLow[0],20,ana::bdtGoodLow[0]);
  TLine * l14 = new TLine(20,ana::bdtGoodLow[0],20,ana::bdtGoodHigh[0]);
  TLine * l15 = new TLine(20,ana::bdtGoodHigh[0],ana::isoBinsHigh[0],ana::bdtGoodHigh[0]);
  TLine * l16 = new TLine(ana::isoBinsHigh[0],ana::bdtGoodHigh[0],ana::isoBinsHigh[0],ana::bdtGoodLow[0]);

  const int nlines = 16;
  TLine * lines[nlines] = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16};
  for (int i = 0; i < nlines; i++) {
    lines[i]->SetLineColor(kRed);
    lines[i]->SetLineWidth(5);
    if (i == 3) break;
  }
  TCanvas * c = new TCanvas("c","",1400,700);
  c->Divide(2,1,0,0);
  c->cd(1);
  gPad->SetLogz();
  gPad->SetTopMargin(.05);
  gPad->SetLeftMargin(.25);
  gPad->SetFillStyle(4000);
  hdata->Scale(1.0/hdata->Integral());
  hdata->GetZaxis()->SetRangeUser(1e-5,1e-2);
  hdata->Draw("col");
  hdata->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
  hdata->GetYaxis()->SetTitle("BDT Score");
  bdata->Draw();
  //for (int i = 0; i < nlines; i++) {
  //  lines[i]->Draw("same");
  //  if (i == 3) break;
  //}
  
  TBox * box = new TBox(4, .25, 19, .63);
  TBox * box2 = new TBox(4, .25, 19, .63);
  box->SetFillColor(kWhite);
  box->SetFillStyle(1001);
  box2->SetLineWidth(2);
  box2->SetLineColor(kBlack);
  box2->SetFillStyle(0);
  box->Draw();
  box2->Draw();
  d.drawAll({"p+p #sqrt{s}=200 GeV"},{Form("%0.0f GeV < p_{T} < %0.0f GeV",ana::ptBins[ptbin],ana::ptBins[ptbin+1]),"Paired clusters"},drawx,drawy,fontsize,700);
  //d.drawText("A",.15,.90, kBlack, 25);
  //d.drawText("B",.15,ana::bdtBadHigh[0]-.05, kBlack, 25);
  //d.drawText("C",.35,.90, kBlack, 25);
  //d.drawText("D",.35,ana::bdtBadHigh[0]-.05, kBlack, 25);
  
  c->cd(2);
  gPad->SetTopMargin(.05);
  gPad->SetRightMargin(.15);
  gPad->SetLogz();
  gPad->SetFillStyle(4000);
  hmcp->Scale(1.0/hmcp->Integral());
  hmcp->GetZaxis()->SetRangeUser(1e-5,1e-2);
  hmcp->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
  hmcp->Draw("colz");
  bdata->Draw();
  //for (int i = 0; i < nlines; i++) {
  //  lines[i]->Draw("same");
  //  if (i == 3) break;
  //}
  //d.drawAll({"MC Photon"},{Form("%0.0f GeV < p_{T} < %0.0f GeV",ana::ptBins[ptbin],ana::ptBins[ptbin+1]),"Paired clusters"},drawx-0.2,drawy,fontsize,700);
  TBox * box3 = new TBox(9, .43, 19, .56);
  TBox * box4 = new TBox(9, .43, 19, .56);
  box3->SetFillColor(kWhite);
  box3->SetFillStyle(1001);
  box4->SetLineWidth(2);
  box4->SetLineColor(kBlack);
  box4->SetFillStyle(0);
  box3->Draw();
  box4->Draw();
  
  d.drawText("Pythia8 #gamma+jet",drawx,drawy-0.05,kBlack,(int)(fontsize*1.25));
  //d.drawText("A",.05,.90, kBlack, 25);
  //d.drawText("B",.05,ana::bdtBadHigh[0]-.05, kBlack, 25);
  //d.drawText("C",.25,.90, kBlack, 25);
  //d.drawText("D",.25,ana::bdtBadHigh[0]-.05, kBlack, 25);

  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/isobdt.pdf");

}
