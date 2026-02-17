#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
void draw_timing() {
  drawer d;
  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsData.root");
  TH1D * h[ana::nJetR];
  for (int ir = 0; ir < ana::nJetR; ir++) {
    h[ir] = (TH1D*)f->Get(Form("hctminusjt%i",ir));
  }
  TH1D * hc = (TH1D*)f->Get(Form("hmtminusct"));
  TH1D * hj = (TH1D*)f->Get(Form("hmtminusjt1"));
  int csize = 700;
  TCanvas * c = new TCanvas("c","",csize*2,csize);
  c->Divide(2,1);
  c->cd(1);
  gStyle->SetOptStat(0);
  gPad->SetTicks(1,1);
  int colors[ana::nJetR] = {kBlue, kBlack, kRed, kGreen};
  TLegend * l = new TLegend(.2,.6,.4,.8);
  l->SetLineWidth(0);
  for (int ir = 0; ir < ana::nJetR; ir++) {
    h[ir]->SetLineColor(colors[ir]);
    h[ir]->SetLineWidth(2);
    d.scale(h[ir]);
    h[ir]->GetYaxis()->SetRangeUser(0,h[ir]->GetMaximum()*1.2);
    h[ir]->Draw("hist same");
    l->AddEntry(h[ir],Form("R = %1.1f",ana::JetRs[ir]));
  }
  l->Draw("same");
  d.drawAll({"pp #sqrt{s}=200 GeV"},{Form("|vz| < %0.0f",ana::vzcut), "#DeltaR_{cluster,jet} > R_{jet}", "All jets and clusters"}, .55, .8, 20, csize);  

  c->cd(2);
  TPad * p = new TPad("p","",0,0,1,1);
  p->Draw();
  p->Divide(1,2);
  p->cd(1);
  hc->Draw();
  p->cd(2);
  hj->Draw();

  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/timing.pdf");
}
