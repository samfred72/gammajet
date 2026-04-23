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
    
void draw_kinematics() {
  drawer d("herwig");
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
  TH1D * hptd = d.get("hclusterpt",0);
  TH1D * hptp = d.get("hclusterpt",1);
  vector<TH1D*> hists = {hptd,hptp};
  vector<string> labels = {"Data", "MC Photon"};

  int manycolors[3] = {kBlack,kBlue,kRed};
  int manystyles[3] = {20, 24, 24};
  TLegend * l = new TLegend(.7,.35,.85,.55);
  gPad->SetTicks(1,1);
  for (int i = 0; i < hists.size(); i++) {
    d.scale(hists.at(i),0,100);
    hists.at(i)->GetXaxis()->SetRangeUser(7,40);
    hists.at(i)->SetMarkerColor(manycolors[i]);
    hists.at(i)->SetMarkerStyle(manystyles[i]);
    hists.at(i)->SetMarkerSize(1);
    hists.at(i)->GetYaxis()->SetTitleSize(0.08);
    hists.at(i)->GetYaxis()->SetLabelSize(0.08);
    hists.at(i)->GetYaxis()->SetTitleOffset(0.8);
    hists.at(i)->GetXaxis()->SetTitleSize(0.08);
    hists.at(i)->GetXaxis()->SetLabelSize(0.08);
    hists.at(i)->GetYaxis()->SetRangeUser(2e-6,hists.at(i)->GetMaximum()*10);
    hists.at(i)->GetXaxis()->SetTitle("Leading p_{T}^{cluster}");
    hists.at(i)->Draw("p same");
    l->AddEntry(hists.at(i),labels.at(i).c_str());
  }
  l->SetLineWidth(0);
  l->Draw();
  d.drawAll({},{"|vz| < 60","Paired Clusters"},0.55,0.75,20,700*0.5);
  gPad->SetLeftMargin(.15);
  ptdiv->cd();
  gPad->SetLogy(0);
  gPad->SetTicks(1,1);
  gPad->SetBottomMargin(.25);
  gPad->SetLeftMargin(.15);
  TH1D * hdivide = (TH1D*)hptd->Clone();
  hdivide->Divide(hptd,hptp);
  hdivide->SetMarkerColor(kBlack);
  hdivide->GetYaxis()->SetTitle("Ratio");
  hdivide->GetYaxis()->SetTitleSize(0.08);
  hdivide->GetYaxis()->SetLabelSize(0.08);
  hdivide->GetXaxis()->SetTitleSize(0.08);
  hdivide->GetXaxis()->SetLabelSize(0.08);
  hdivide->GetYaxis()->SetRangeUser(0,2);
  hdivide->Draw("p");
  d.drawLine(7,1,40,1);
  
  hdivide->Fit(func,"RQ");
  
  TF1 * pfunc = new TF1("pfunc", piecewise_func, 10, 35, 0);
  pfunc->SetNpx(10000);
  TFile * ftemp = TFile::Open("/home/samson72/sphnx/gammajet/hists/reweighting.root","RECREATE");
  hdivide->Write();
  pfunc->Write();
  TH1D * hvzd = d.get("hvz",0);
  TH1D * hvzp = d.get("hvz",1);
  hvzd->Scale(1.0/hvzd->Integral());
  hvzp->Scale(1.0/hvzp->Integral());
  TH1D * hvzr = (TH1D*)hvzd->Clone("hvz_reweight");
  hvzr->Reset("ICES");
  hvzr->Divide(hvzd,hvzp);
  hvzd->Write();
  hvzp->Write();
  hvzr->Write();
  
  ftemp->Close();

  cpt->SaveAs("/home/samson72/sphnx/gammajet/pdfs/cluster_pt.pdf");
  
  
  TCanvas * cetaphi = new TCanvas("cetaphi","",1500,700);
  TPad * pleft = new TPad("pleft","",0,0,0.4,1);
  pleft->SetRightMargin(0.01);
  pleft->Draw();
  pleft->Divide(1,2,0,0);
  pleft->cd(1);
  gPad->SetLogy(0);
  TH1D * hetad = d.get("hclustereta",0);
  TH1D * hetap = d.get("hclustereta",1);
  vector<TH1D*> histseta = {hetad,hetap};

  TLegend * leta = new TLegend(.45,.25,.75,.45);
  gPad->SetTicks(1,1);
  for (int i = 0; i < hists.size(); i++) {
    d.scale(histseta.at(i),0,100);
    histseta.at(i)->GetXaxis()->SetRangeUser(7,40);
    histseta.at(i)->SetMarkerColor(manycolors[i]);
    histseta.at(i)->SetMarkerStyle(manystyles[i]);
    histseta.at(i)->SetMarkerSize(1);
    histseta.at(i)->GetYaxis()->SetTitleSize(0.08);
    histseta.at(i)->GetYaxis()->SetLabelSize(0.08);
    histseta.at(i)->GetYaxis()->SetTitleOffset(0.95);
    histseta.at(i)->GetXaxis()->SetTitleSize(0.08);
    histseta.at(i)->GetXaxis()->SetLabelSize(0.08);
    histseta.at(i)->GetXaxis()->SetTitle("Paired cluster #eta");
    histseta.at(i)->Draw("p same");
    leta->AddEntry(histseta.at(i),labels.at(i).c_str());
  }
  leta->SetLineWidth(0);
  leta->Draw();
  d.drawAll({},{"|vz| < 60","Paired Clusters"},0.45,0.65,20,350);
  gPad->SetLeftMargin(.15);
  pleft->cd(2);
  gPad->SetLogy(0);
  gPad->SetTicks(1,1);
  gPad->SetBottomMargin(.25);
  gPad->SetLeftMargin(.15);
  hdivide = (TH1D*)hetad->Clone();
  hdivide->Divide(hetad,hetap);
  hdivide->SetMarkerColor(kBlack);
  hdivide->GetYaxis()->SetTitle("Ratio");
  hdivide->GetYaxis()->SetRangeUser(0,2);
  hdivide->GetYaxis()->SetTitleSize(0.08);
  hdivide->GetYaxis()->SetLabelSize(0.08);
  hdivide->GetXaxis()->SetTitleSize(0.08);
  hdivide->GetXaxis()->SetLabelSize(0.08);
  hdivide->Draw("p");
  d.drawLine(-1.1,1,1.1,1);
  
  cetaphi->cd();
  TPad * pright = new TPad("pright","",0.4,0,1,1);
  pright->Draw();
  pright->SetLeftMargin(0.12);
  pright->Divide(2,1,0,0);
  TH2D * hepd = d.get2d("hclusteretaphi",0);
  TH2D * hepp = d.get2d("hclusteretaphi",1);
  hepd->GetXaxis()->SetTitle("");
  hepd->GetYaxis()->SetTitle("Paired cluster #phi");
  hepp->GetXaxis()->SetTitle("Paired cluster #eta");
  pright->cd(1);
  hepd->Draw("col");
  d.drawText("Data",0.5,0.8,1,20);
  hepd->GetYaxis()->SetTitleSize(0.05);
  hepd->GetYaxis()->SetLabelSize(0.05);
  hepd->GetXaxis()->SetLabelSize(0.05);
  pright->cd(2);
  hepp->Draw("col");
  hepp->GetXaxis()->SetTitleSize(0.05);
  hepp->GetXaxis()->SetLabelSize(0.05);
  d.drawText("MC Photon",0.5,0.8,1,20);

  cetaphi->SaveAs("/home/samson72/sphnx/gammajet/pdfs/cluster_etaphi.pdf");
  
  
  TCanvas * jetaphi = new TCanvas("jetaphi","",1500,700);
  TPad * jleft = new TPad("jleft","",0,0,0.4,1);
  jleft->SetRightMargin(0.01);
  jleft->Draw();
  jleft->Divide(1,2,0,0);
  jleft->cd(1);
  gPad->SetLogy(0);
  TH1D * hjetad = d.get("hjeteta1",0);
  TH1D * hjetap = d.get("hjeteta1",1);
  vector<TH1D*> histsjeta = {hjetad,hjetap};

  TLegend * ljeta = new TLegend(.45,.20,.75,.40);
  gPad->SetTicks(1,1);
  for (int i = 0; i < hists.size(); i++) {
    d.scale(histsjeta.at(i),0,100);
    histsjeta.at(i)->GetXaxis()->SetRangeUser(7,40);
    histsjeta.at(i)->SetMarkerColor(manycolors[i]);
    histsjeta.at(i)->SetMarkerStyle(manystyles[i]);
    histsjeta.at(i)->SetMarkerSize(1);
    histsjeta.at(i)->GetXaxis()->SetTitle("Leading jet #eta");
    histsjeta.at(i)->GetYaxis()->SetTitleSize(0.08);
    histsjeta.at(i)->GetYaxis()->SetLabelSize(0.08);
    histsjeta.at(i)->GetYaxis()->SetTitleOffset(0.95);
    histsjeta.at(i)->GetXaxis()->SetTitleSize(0.08);
    histsjeta.at(i)->GetXaxis()->SetLabelSize(0.08);
    histsjeta.at(i)->Draw("p same");
    ljeta->AddEntry(histsjeta.at(i),labels.at(i).c_str());
  }
  ljeta->SetLineWidth(0);
  ljeta->Draw();
  d.drawAll({},{"|vz| < 60","Paired jets","R=0.4"},0.40,0.65,20,350);
  gPad->SetLeftMargin(.15);
  jleft->cd(2);
  gPad->SetLogy(0);
  gPad->SetTicks(1,1);
  gPad->SetBottomMargin(.25);
  gPad->SetLeftMargin(.15);
  hdivide = (TH1D*)hjetad->Clone();
  hdivide->Divide(hjetad,hjetap);
  hdivide->SetMarkerColor(kBlack);
  hdivide->GetYaxis()->SetTitle("Ratio");
  hdivide->GetYaxis()->SetRangeUser(0,2);
  hdivide->GetYaxis()->SetTitleSize(0.08);
  hdivide->GetYaxis()->SetLabelSize(0.08);
  hdivide->GetXaxis()->SetTitleSize(0.08);
  hdivide->GetXaxis()->SetLabelSize(0.08);
  hdivide->Draw("p");
  d.drawLine(-1.1,1,1.1,1);
  
  jetaphi->cd();
  TPad * jright = new TPad("jright","",0.4,0,1,1);
  jright->Draw();
  jright->SetLeftMargin(0.1);
  jright->Divide(2,1,0,0);
  TH2D * hjepd = d.get2d("hjetetaphi1",0);
  TH2D * hjepp = d.get2d("hjetetaphi1",1);
  hjepd->GetXaxis()->SetTitle("");
  hjepd->GetYaxis()->SetTitle("Paired jet #phi");
  hjepp->GetXaxis()->SetTitle("Paired jet #eta");
  jright->cd(1);
  hjepd->GetYaxis()->SetTitleSize(0.05);
  hjepd->GetYaxis()->SetLabelSize(0.05);
  hjepd->GetXaxis()->SetLabelSize(0.05);
  hjepd->Draw("col");
  d.drawText("Data",0.5,0.8,1,20);
  jright->cd(2);
  hjepp->GetXaxis()->SetTitleSize(0.05);
  hjepp->GetXaxis()->SetLabelSize(0.05);
  hjepp->Draw("col");
  d.drawText("MC Photon",0.5,0.8,1,20);

  jetaphi->SaveAs("/home/samson72/sphnx/gammajet/pdfs/jetR04_etaphi.pdf");
}
