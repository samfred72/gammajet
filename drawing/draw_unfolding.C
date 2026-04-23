#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
    
void draw_unfolding() {
  drawer d(1,"pythia");
  gStyle->SetOptStat(0);

  TCanvas * c = new TCanvas("c","",700,1000);
  TPad * p1 = new TPad("p1","",0,.4,1,1);
  TPad * p2 = new TPad("p2","",0,0,1,.4);
  p1->Draw();
  p2->Draw();
  p1->cd();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);
  gPad->SetLogy();
  gPad->SetTicks(1,1);
  int itype = 1; // 1: photon sample, 2: jet sample
  int isample = -1; // -1: all, 0: photon5, 1: photon10, etc.
  const char * obj = "jet";
  const char * ir = "1";
  const char * half = "";
  TH1D * ht = d.get(Form("htruth%spt%s%s" , obj, half, ir), itype,isample);
  TH1D * hr = d.get(Form("hreco%spt%s%s"  , obj, half, ir), itype,isample);
  TH1D * hu = d.get(Form("hunfold%spt%s%s", obj, half, ir), itype,isample);

  vector<TH1D*> hists = {ht,hr,hu};
  vector<string> labels = {"Truth Jet", "Reco Jet", "Unfolded Jet"};

  int manycolors[3] = {kBlack,kBlue,kRed};
  int manystyles[3] = {20, 24, 20};
  TLegend * l = new TLegend(.55,.45,.85,.6);
  gPad->SetTicks(1,1);
  for (int i = 0; i < hists.size(); i++) {
    d.scale(hists.at(i),0,100);
    hists.at(i)->GetXaxis()->SetRangeUser(0,100);
    hists.at(i)->SetMarkerColor(manycolors[i]);
    hists.at(i)->SetMarkerStyle(manystyles[i]);
    hists.at(i)->SetMarkerSize(1);
    //hists.at(i)->GetYaxis()->SetTitleSize(0.08);
    //hists.at(i)->GetYaxis()->SetLabelSize(0.08);
    //hists.at(i)->GetYaxis()->SetTitleOffset(0.8);
    //hists.at(i)->GetXaxis()->SetTitleSize(0.08);
    //hists.at(i)->GetXaxis()->SetLabelSize(0.08);
    hists.at(i)->GetYaxis()->SetRangeUser(2e-9,hists.at(i)->GetMaximum()*10);
    hists.at(i)->GetXaxis()->SetTitle("p_{T}^{paired jet}");
    hists.at(i)->Draw("p same");
    l->AddEntry(hists.at(i),labels.at(i).c_str());
  }
  l->SetLineWidth(0);
  l->Draw();
  d.drawAll({"MC Photon combined"},{"|vz| < 60","Paired Jets R=0.4"},0.55,0.75,20,700);
  gPad->SetLeftMargin(.15);

  p2->cd();
  gPad->SetLeftMargin(.15);
  gPad->SetLogy(0);
  gPad->SetTicks(1,1);
  TH1D * hdivide = (TH1D*)ht->Clone();
  hdivide->Reset("ICES");
  hdivide->Divide(hu,ht);
  hdivide->SetLineColor(kBlack);
  hdivide->SetMarkerColor(kBlack);
  hdivide->SetMarkerSize(1);
  hdivide->SetMarkerStyle(20);
  hdivide->GetYaxis()->SetRangeUser(0.5,1.5);
  hdivide->Draw();

  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/unfolded_pt.pdf");
}
