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
  float drawy = .65;
  float fontsize = 30;
  int ptbin = 5;

  const char * histname = Form("hisobdt%i",ptbin);
  TH2D * hdata = (TH2D*)f->Get(histname);
  TH2D * hmcp = d.combineMC2d(histname,1);

  TBox * box1 = new TBox(-1, ana::bdtGoodLow[0], ana::isoBins[0], ana::bdtGoodHigh[0]);
  TBox * box2 = new TBox(-1, ana::bdtBadLow[0], ana::isoBins[0], ana::bdtBadHigh[0]);
  TBox * box3 = new TBox(ana::isoBinsHigh[0], ana::bdtGoodLow[0], 20, ana::bdtGoodHigh[0]);
  TBox * box4 = new TBox(ana::isoBinsHigh[0], ana::bdtBadLow[0], 20, ana::bdtBadHigh[0]);
  TBox * boxes[4] = {box1,box2,box3,box4};
  for (int i = 0; i < 4; i++) {
    boxes[i]->SetLineColor(kBlack);
    boxes[i]->SetFillStyle(0);
    boxes[i]->SetLineWidth(4);
  }
  
  TCanvas * c = new TCanvas("c","",1400,700);
  c->Divide(2,1,0,0);
  c->cd(1);
  gPad->SetLogz();
  gPad->SetTopMargin(.05);
  gPad->SetLeftMargin(.25);
  gPad->SetBottomMargin(0.15);
  gPad->SetFillStyle(4000);
  gPad->SetTicks(1,1);
  hdata->Scale(1.0/hdata->Integral());
  hdata->GetZaxis()->SetRangeUser(1e-5,1e-2);
  hdata->Draw("col");
  hdata->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
  hdata->GetXaxis()->ChangeLabel(-1,-1,0);
  hdata->GetYaxis()->SetTitle("BDT Score");
  hdata->GetXaxis()->SetLabelSize(.05);
  hdata->GetYaxis()->SetLabelSize(.05);
  hdata->GetXaxis()->SetTitleSize(.05);
  hdata->GetYaxis()->SetTitleSize(.05);
  for (int i = 0; i <4; i++) {
    boxes[i]->Draw();
  }
  
  TBox * lbox = new TBox(7, .34, 11, .47);
  TBox * lbox2 = new TBox(7, .34, 11, .47);
  lbox->SetFillColor(kWhite);
  lbox->SetFillStyle(1001);
  lbox2->SetLineWidth(2);
  lbox2->SetLineColor(kBlack);
  lbox2->SetFillStyle(0);
  d.drawText("A",.28,.90, kBlack, 35);
  d.drawText("B",.28,ana::bdtBadHigh[0]-.05, kBlack, 35);
  d.drawText("C",.48,.90, kBlack, 35);
  d.drawText("D",.48,ana::bdtBadHigh[0]-.05, kBlack, 35);
  lbox->Draw();
  lbox2->Draw();
  d.drawText("Data",0.55,.45,kBlack,(int)(fontsize*1.25));
  //d.drawAll({"p+p #sqrt{s}=200 GeV"},{Form("%0.0f GeV < p_{T}^{#gamma} < %0.0f GeV",ana::ptBins[ptbin],ana::ptBins[ptbin+1]),"Jet R=0.4","|#eta_{jet}| < 0.7","#Delta#phi > 7#pi/8"},drawx,drawy,fontsize,700);
  
  c->cd(2);
  gPad->SetTopMargin(.05);
  gPad->SetRightMargin(.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetLogz();
  gPad->SetFillStyle(4000);
  gPad->SetTicks(1,1);
  hmcp->Scale(1.0/hmcp->Integral());
  hmcp->GetZaxis()->SetRangeUser(1e-5,1e-2);
  hmcp->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
  hmcp->GetXaxis()->SetLabelSize(.05);
  hmcp->GetYaxis()->SetLabelSize(.05);
  hmcp->GetZaxis()->SetLabelSize(.05);
  hmcp->GetXaxis()->SetTitleSize(.05);
  hmcp->GetYaxis()->SetTitleSize(.05);
  hmcp->Draw("colz");
  for (int i = 0; i <4; i++) {
    boxes[i]->Draw();
  }
  TBox * lbox3 = new TBox(9, .34, 19, .47);
  TBox * lbox4 = new TBox(9, .34, 19, .47);
  lbox3->SetFillColor(kWhite);
  lbox3->SetFillStyle(1001);
  lbox4->SetLineWidth(2);
  lbox4->SetLineColor(kBlack);
  lbox4->SetFillStyle(0);
  
  d.drawText("A",.05,.90, kBlack, 35);
  d.drawText("B",.05,ana::bdtBadHigh[0]-.05, kBlack, 35);
  d.drawText("C",.25,.90, kBlack, 35);
  d.drawText("D",.25,ana::bdtBadHigh[0]-.05, kBlack, 35);
  
  lbox3->Draw();
  lbox4->Draw();
  d.drawText("Pythia8 #gamma+jet",drawx,.45,kBlack,(int)(fontsize*1.25));

  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/isobdt.pdf");

}
