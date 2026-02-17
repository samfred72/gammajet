#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
void draw_deltaphi() {
  drawer d;
  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsData.root");
  int ptbins[4] = {1,3,5,7};
  TCanvas * c = new TCanvas("c","",1000,1000);
  gStyle->SetOptStat(0);
  TPad * p = new TPad("p","",0,0,1,1);
  p->Divide(2,2,0,0);
  p->Draw();
  TLegend * l = new TLegend(.1,.65,.4,.75);
  l->SetLineWidth(0);
  TH1D * hd[4];
  TH1D * hp[4];
  int ir = 1;
  for (int i = 0; i < 4; i++) {
    p->cd(i+1);
    gPad->SetTicks(1,1);
    int ir = 1;
    hd[i] = (TH1D*)f->Get(Form("hdeltaphi%i_%i",ptbins[i],ir));
    hp[i] = d.combineMC(Form("hdeltaphi%i_%i",ptbins[i],ir),1);

    hd[i]->SetLineColor(kBlue);
    hp[i]->SetLineColor(kMagenta+1);
    hd[i]->SetLineWidth(2);
    hp[i]->SetLineWidth(2);
    d.scale(hd[i],0,M_PI);
    d.scale(hp[i],0,M_PI);
    hd[i]->GetYaxis()->SetRangeUser(0,hd[i]->GetMaximum()*1.2);
    if (i == 0 || i == 2) gPad->SetLeftMargin(0.15);
    if (i == 1 || i == 3) gPad->SetRightMargin(0.02);
    if (i == 2) {
      hd[i]->GetYaxis()->SetTitle("");
      hd[i]->GetXaxis()->SetTitle("");
    }

    hd[i]->Draw("hist same");
    hp[i]->Draw("hist same");
    d.drawText(Form("%1.1f GeV < p_{T}^{cluster} < %1.1f GeV", ana::ptBins[ptbins[i]], ana::ptBins[ptbins[i]+1]),.2,.15,1,20);
  }
  l->AddEntry(hd[0],"data");
  l->AddEntry(hp[0],"MC photon");

  c->cd();
  l->Draw("same");
  d.drawAll({"pp #sqrt{s}=200 GeV","MC Photon run28"},{Form("|vz| < %0.0f cm",ana::vzcut), "pair candidates"}, .1, .9, 20, 1000);  

  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/deltaphi.pdf");
}
