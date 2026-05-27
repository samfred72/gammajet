#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
drawer d;

TF1 * func = new TF1("func","expo",12,35);
double piecewise_func(double *x, double *par) {
  if (x[0] >= 10 && x[0] < 11) {
    return 1.5627383;
  }
  else if (x[0] >= 11 && x[0] < 12) {
    return 1.0147688;
  }
  else return func->Eval(x[0]);
}
    
void draw_pt_formarzia() {
  drawer d("pythia");
  gStyle->SetOptStat(0);

  TCanvas * cpt = new TCanvas("cpt","",700,700);
  TPad * ppt = new TPad("ppt","",0,.5,1,1);
  TPad * ptdiv = new TPad("ptdiv","",0,0,1,.5);
  ppt->Draw();
  ptdiv->Draw();
  ppt->cd();
  ppt->SetBottomMargin(0);
  ptdiv->SetTopMargin(0);
  gPad->SetLogy();
  TH1D * hptd0 = d.get("hjetpt_formarzia2_1_0",1);
  TH1D * hptd1 = d.get("hjetpt_formarzia2_1_1",1);
  vector<TH1D*> hists = {hptd0,hptd1};
  //vector<string> labels = {"pre #Delta#phi cut", "post #Delta#phi cut"};
  vector<string> labels = {"pre BDT/Isolation cut", "post BDT/Isolation cut"};

  int manycolors[3] = {kBlack,kBlue,kRed};
  int manystyles[3] = {20, 24, 24};
  TLegend * l = new TLegend(.5,.30,.85,.50);
  gPad->SetTicks(1,1);
  for (int i = 0; i < hists.size(); i++) {
    hists.at(i)->GetXaxis()->SetRangeUser(7,40);
    hists.at(i)->SetMarkerColor(manycolors[i]);
    hists.at(i)->SetMarkerStyle(manystyles[i]);
    hists.at(i)->SetMarkerSize(1);
    hists.at(i)->GetYaxis()->SetTitleSize(0.08);
    hists.at(i)->GetYaxis()->SetLabelSize(0.08);
    hists.at(i)->GetYaxis()->SetTitleOffset(0.8);
    hists.at(i)->GetXaxis()->SetTitleSize(0.08);
    hists.at(i)->GetXaxis()->SetLabelSize(0.08);
    //hists.at(i)->GetYaxis()->SetRangeUser(0.5,hists.at(i)->GetMaximum()*10);
    hists.at(i)->GetXaxis()->SetTitle("Leading p_{T}^{jet}");
    hists.at(i)->Draw("p same");
    l->AddEntry(hists.at(i),labels.at(i).c_str());
  }
  l->SetLineWidth(0);
  l->Draw();
  d.drawAll({"MC Photon"},{"|vz| < 60","Passing all other cuts"},0.55,0.75,20,700*0.5);
  gPad->SetLeftMargin(.15);
  ptdiv->cd();
  gPad->SetLogy(0);
  gPad->SetTicks(1,1);
  gPad->SetBottomMargin(.25);
  gPad->SetLeftMargin(.15);
  TH1D * hdivide = (TH1D*)hptd0->Clone();
  hdivide->Divide(hptd0,hptd1);
  hdivide->SetMarkerColor(kBlack);
  hdivide->GetYaxis()->SetTitle("Ratio");
  hdivide->GetYaxis()->SetTitleSize(0.08);
  hdivide->GetYaxis()->SetLabelSize(0.08);
  hdivide->GetXaxis()->SetTitleSize(0.08);
  hdivide->GetXaxis()->SetLabelSize(0.08);
  hdivide->GetYaxis()->SetRangeUser(0,5);
  hdivide->Draw("p");
  d.drawLine(7,1,40,1);
  
} 
