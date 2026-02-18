#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
void draw_xjbdt() {
  drawer d;
  gStyle->SetOptStat(0);
 
  int dcolors[ana::nBdtBins] = {kOrange, kRed, kBlue, kGreen};
  int pcolors[ana::nBdtBins] = {kOrange, kRed, kBlue, kGreen};

  TH1D * hd[ana::nBdtBins];
  TH1D * hp[ana::nBdtBins];
  TCanvas * c = new TCanvas("c","",700,700);
  TLegend * l = new TLegend(.55,.5,.85,.85);
  l->SetLineWidth(0);
  for (int i = 0; i < ana::nBdtBins; i++) {
    float lowcluster = 16;
    float highcluster = 25;
    float minjet = ana::jet_calib_pt_cut[1];
    float lowrange = (int)(minjet/lowcluster/0.08 + 1)*0.08;
    hd[i] = d.get(Form("hxjbdt%i_1",i),0); 
    hd[i]->SetLineColor(dcolors[i]);
    hd[i]->SetMarkerColor(dcolors[i]);
    hd[i]->SetMarkerStyle(20);
    hd[i]->Rebin(4);
    TF1 * f = d.fit(hd[i],lowrange,2,"RLQI0");
    d.scale(hd[i]);
    f->SetParameter(0,f->GetParameter(0)/hd[i]->GetEntries());
    hd[i]->GetYaxis()->SetRangeUser(0,hd[i]->GetMaximum()*1.5);
    l->AddEntry(hd[i],Form("%0.1f < BDT < %0.1f",ana::bdtBins[i],ana::bdtBins[i+1]));
    l->AddEntry((TObject*)0,Form("#mu = %2.3f #pm %2.3f",f->GetParameter(1),f->GetParError(1)), "");
    hd[i]->Draw("p same");
    f->SetLineColor(dcolors[i]);
    f->Draw("same");
  } 
  l->Draw();
  d.drawAll({"Data"},{"Analysis cuts", "R=0.4","16 GeV < p_{T}^{cluster} < 25 GeV"},.15,.85,15,700);
  
  TCanvas * c2 = new TCanvas("c2","",700,700);
  TLegend * l2 = new TLegend(.55,.5,.85,.85);
  l2->SetLineWidth(0);
  for (int i = 0; i < ana::nBdtBins; i++) {
    float lowcluster = 16;
    float highcluster = 25;
    float minjet = ana::jet_calib_pt_cut[1];
    float lowrange = (int)(minjet/lowcluster/0.08 + 1)*0.08;
    hp[i] = d.get(Form("hxjbdt%i_1",i),1);
    hp[i]->SetLineColor(pcolors[i]);
    hp[i]->SetMarkerColor(pcolors[i]);
    hp[i]->SetMarkerStyle(24);
    hp[i]->Rebin(4);
    d.scale(hp[i]);
    TF1 * f = d.fit(hp[i],lowrange,2,"RMQI0");
    hp[i]->GetYaxis()->SetRangeUser(0,hp[i]->GetMaximum()*1.5);
    l2->AddEntry(hp[i],Form("%0.1f < BDT < %0.1f",ana::bdtBins[i],ana::bdtBins[i+1]));
    l2->AddEntry((TObject*)0,Form("#mu = %2.3f #pm %2.3f",f->GetParameter(1),f->GetParError(1)), "");
    hp[i]->Draw("p same");
    f->SetLineColor(pcolors[i]);
    f->Draw("same");
  } 
  l2->Draw();
  d.drawAll({"MC Photon"},{"Analysis cuts","R=0.4","16 GeV < p_{T}^{cluster} < 25 GeV"},.15,.85,15,700);

}
