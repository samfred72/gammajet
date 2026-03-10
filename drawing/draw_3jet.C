#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
void draw_3jet() {
  drawer d;
  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsData.root");
  TFile * fp = TFile::Open("/home/samson72/sphnx/gammajet/hists/histsPhoton20.root");
  gStyle->SetOptStat(0);
  int ir = 1;
  TH2D * hd = (TH2D*)f->Get(Form("h3jetpt%i",ir));
  TH2D * hp = d.combineMC2d(Form("h3jetpt%i",ir),1);
  float drawx = .2;
  float drawy = .85;
  float fontsize = 20; 
  
  TH1D * hfracd = new TH1D("hfracd",";Leading cluster pT;fraction of pairs with 3rd jet",ana::nPtBins,ana::ptBins);
  TH1D * hfracp = new TH1D("hfracp",";Leading cluster pT;fraction of pairs with 3rd jet",ana::nPtBins,ana::ptBins);

  for (int i = 0; i < ana::nPtBins; i++) {
    TH1D * hrd = (TH1D*)f->Get(Form("hratio_%i_%i_1_0_0_0",i,ir));
    TH1D * hrp = d.combineMC(Form("hratio_%i_%i_1_0_0_0",i,ir),1);
    //TH1D * hrp = (TH1D*)fp->Get(Form("hratio_%i_1_1_0_0_0",i));//d.combineMC(Form("hratio_%i_1_1_0_0_0",i),1);
    //cout << hp->ProjectionY(Form("pp%i",i),i+1,i+1)->Integral() << " " << hrp->Integral() << endl;
    //cout << hd->ProjectionY(Form("pp%i",i),i+1,i+1)->Integral() << " " << hrd->Integral() << endl;
    if (hrd->GetEntries() > 0) hfracd->SetBinContent(i+1, hd->ProjectionY(Form("pd%i",i),i+1,i+1)->Integral()/hrd->Integral());
    if (hrp->GetEntries() > 0) hfracp->SetBinContent(i+1, hp->ProjectionY(Form("pp%i",i),i+1,i+1)->Integral()/hrp->Integral());
  }

  TCanvas * c = new TCanvas("c","",1400,700);
  c->Divide(2,1,0,0);
  c->cd(1);
  gPad->SetLogz();
  hd->Scale(1.0/hd->Integral());
  hd->GetZaxis()->SetRangeUser(1e-5,1e-1);
  hd->Draw("col");
  d.drawAll({"Data"},{Form("R=%1.1f",ana::JetRs[ir])},drawx,drawy,fontsize,700);
  
  c->cd(2);
  gPad->SetRightMargin(.15);
  gPad->SetLogz();
  hp->Scale(1.0/hp->Integral());
  hp->GetZaxis()->SetRangeUser(1e-5,1e-1);
  hp->Draw("colz");
  d.drawText("MC Photon",drawx,drawy-0.05,kBlack,(int)(fontsize*1.25));

  TCanvas * c2 = new TCanvas("c2","",700,700);
  hfracd->GetYaxis()->SetRangeUser(0,0.8);
  hfracd->SetMarkerColor(kBlue);
  hfracp->SetMarkerColor(kMagenta+1);
  hfracd->SetLineColor(kBlue);
  hfracp->SetLineColor(kMagenta+1);
  hfracd->SetMarkerStyle(20);
  hfracp->SetMarkerStyle(20);
  hfracd->Draw("hist");
  hfracp->Draw("hist same");
  hfracd->Draw("p same");
  hfracp->Draw("p same");
  //d.drawLine(10,1,35,1);

  d.drawAll({"Data","MC Photon"},{Form("R=%1.1f",ana::JetRs[ir]),"Pairs with 3rd jet"},drawx,drawy,fontsize,700);

  TLegend * l = new TLegend(drawx,drawy-.3,drawx+.2,drawy-.2);
  l->AddEntry(hfracd,"Data");
  l->AddEntry(hfracp,"MC Photon");
  l->SetLineWidth(0);
  l->Draw();
  c2->SaveAs("/home/samson72/sphnx/gammajet/pdfs/prob_3jet.pdf");

  TCanvas * c3 = new TCanvas("c3","",700,700);
  TH1D * hdrd = d.get("h3jetdeltar5_1",0);
  TH1D * hdrp = d.get("h3jetdeltar5_1",1);
  hdrd->SetLineColor(kBlue);
  hdrp->SetLineColor(kMagenta+1);
  d.scale(hdrd,0,4);
  d.scale(hdrp,0,4);
  hdrd->Draw("hist");
  hdrp->Draw("hist same");
  TLegend * ldr = new TLegend(.50,.4,.80,.6);
  ldr->AddEntry(hdrd,"Data");
  ldr->AddEntry(hdrp,"MC Photon");
  ldr->SetLineWidth(0);
  ldr->Draw();
  d.drawAll({},{Form("%0.0f GeV < p_{T}^{cluster} < %0.0f GeV",ana::ptBins[5],ana::ptBins[6]), "Jet R=0.4"},.4,.8,20,700);
  c3->SaveAs("/home/samson72/sphnx/gammajet/pdfs/deltaR_3jet.pdf");

}
