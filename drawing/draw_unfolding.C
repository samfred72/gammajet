#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void drawData(TFile * f, string histname, string label) {
  TH1D * hist = (TH1D*)f->Get(histname.c_str());
  TLegend * l = new TLegend(.7,.45,.85,.7);
  TCanvas * c = new TCanvas(Form("c%s",histname.c_str()),"",700,700);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
  gPad->SetTicks(1,1);
  
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  hist->GetYaxis()->SetTitle("Scaled counts");
  hist->GetXaxis()->SetTitle(label.c_str());
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->Draw("p same");
}
    
void draw_unfolding() {
  drawer d(1);
  gStyle->SetOptStat(0);

  TCanvas * c = new TCanvas("c","",700,700);
  //TPad * ppt = new TPad("ppt","",0,.5,1,1);
  //TPad * ptdiv = new TPad("ptdiv","",0,0,1,.5);
  //ppt->Draw();
  //ptdiv->Draw();
  //ppt->cd();
  //ppt->SetBottomMargin(0);
  //ptdiv->SetTopMargin(0);
  gPad->SetLogy();
  TH1D * ht = d.get("htruthjetpt1",1);
  TH1D * hr = d.get("hjetpt1",1);
  TH1D * hu = d.get("hjetreco1",1);

  vector<TH1D*> hists = {ht,hr,hu};
  vector<string> labels = {"Truth Jet", "Reco Jet", "Unfolded Jet"};

  int manycolors[3] = {kBlack,kBlue,kRed};
  int manystyles[3] = {20, 24, 20};
  TLegend * l = new TLegend(.55,.45,.85,.6);
  gPad->SetTicks(1,1);
  for (int i = 0; i < hists.size(); i++) {
    d.scale(hists.at(i),0,100);
    hists.at(i)->GetXaxis()->SetRangeUser(0,60);
    hists.at(i)->SetMarkerColor(manycolors[i]);
    hists.at(i)->SetMarkerStyle(manystyles[i]);
    hists.at(i)->SetMarkerSize(1);
    //hists.at(i)->GetYaxis()->SetTitleSize(0.08);
    //hists.at(i)->GetYaxis()->SetLabelSize(0.08);
    //hists.at(i)->GetYaxis()->SetTitleOffset(0.8);
    //hists.at(i)->GetXaxis()->SetTitleSize(0.08);
    //hists.at(i)->GetXaxis()->SetLabelSize(0.08);
    hists.at(i)->GetYaxis()->SetRangeUser(2e-6,hists.at(i)->GetMaximum()*10);
    hists.at(i)->GetXaxis()->SetTitle("p_{T}^{paired jet}");
    hists.at(i)->Draw("p same");
    l->AddEntry(hists.at(i),labels.at(i).c_str());
  }
  l->SetLineWidth(0);
  l->Draw();
  d.drawAll({"MC Photon combined"},{"|vz| < 60","Paired Jets R=0.4"},0.55,0.75,20,700);
  gPad->SetLeftMargin(.15);

  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/unfolded_pt.pdf");
}
