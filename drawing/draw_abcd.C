#include "../headers/commonUtility.h"
#include "../headers/Style_jaebeom.h"
#include "../headers/ana.cxx"
ana anaclone;

void scale(TH1D * h) {
  double x_low_value = 0.5;
  double x_high_value = 2;

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

void draw_all(TCanvas * c, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratio, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_p, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_j, 
    int ia, int l, int k) {

  float drawx = 0.80;
  float drawy = 0.92;
  float fontsize = 60;
  int njrebin = 1;
  int nprebin = 1;
  int offset = 0;
  int colors[6] = {kSpring + 2, kBlue, kGreen + 3, kMagenta+1, kTeal, kSpring+2};
  float isocut = anaclone.isoBins[l];
  float bdtcut = anaclone.bdtBins[l];
  string calibstring = (k == 1 ? "JES Calibrated" : "Uncalibrated Jets");
  vector<string> t1 = {"#bf{Analysis region A}:","#bf{Analysis region B}:","#bf{Analysis region C}:","#bf{Analysis region D}:", "Combined ABCD"};
  vector<string> t2 = {Form("E_{iso} < %0.0f GeV",isocut),
                       Form("E_{iso} > %0.0f GeV",isocut),
                       Form("E_{iso} < %0.0f GeV",isocut),
                       Form("E_{iso} > %0.0f GeV",isocut),
                       ""};
  vector<string> t3 = {Form("bdt score > %0.1f",bdtcut),
                       Form("bdt score > %0.1f",bdtcut),
                       Form("bdt score < %0.1f",bdtcut), 
                       Form("bdt score < %0.1f",bdtcut),
                       ""};

  TPad * p = new TPad("p","",0,0,1,1);
  p->SetLeftMargin(0.25);
  p->SetBottomMargin(0.25);
  p->Divide(anaclone.nJetR+1,anaclone.nPtBins,0,0);
  p->Draw();

  for (int i = 0; i < anaclone.nPtBins; i++) {
    for (int j = 0; j < anaclone.nJetR; j++) {
      TH1D * h1 = hratio[i][j][k][l][ia];
      TH1D * h2 = hratiomc_p[i][j][k][l][ia];
      TH1D * h3 = hratiomc_j[i][j][k][l][ia];

      int index = i*anaclone.nJetR + j + 1 + offset;
      if (index % (anaclone.nJetR+1) == 0) {index++; offset++;}
      p->cd(index);
      if ((index + 1) % (anaclone.nJetR+1) == 0) gPad->SetRightMargin(0.01);
      gPad->SetTicks(1,1);
      if (index == (anaclone.nPtBins*(anaclone.nJetR+1) - 1)) h1->GetXaxis()->SetTitle("p_{T,max}^{Jet}/p_{T,max}^{cluster}");
      else h1->GetXaxis()->SetTitle("");
      h1->GetXaxis()->SetTitleSize(0.10);
      h1->GetXaxis()->SetLabelSize(0.08);
      if (index == (anaclone.nJetR+1)*(anaclone.nPtBins-1)+1) h1->GetXaxis()->SetLabelSize(0.07);
      h1->GetXaxis()->SetNdivisions(405,false);
      h1->GetXaxis()->ChangeLabel(-1,-1,0);
      if (index == 1) h1->GetYaxis()->SetTitle("Arbitrary Units");
      else h1->GetYaxis()->SetTitle("");
      h1->GetYaxis()->SetTitleSize(0.10);
      h1->GetYaxis()->SetLabelSize(0.08);
      if (index == (anaclone.nJetR+1)*(anaclone.nPtBins-1)+1) h1->GetYaxis()->SetLabelSize(0.07);
      h1->GetYaxis()->SetLabelOffset(0.04);
      h1->GetYaxis()->SetMaxDigits(3);
      h1->GetYaxis()->SetDecimals(2);
      scale(h1);
      h1->GetYaxis()->SetRangeUser(0,h1->GetMaximum()*1.5);
      h1->SetLineColor(colors[1]);
      h1->SetLineWidth(2);
      h1->Draw("hist same");

      scale(h2);
      h2->Scale(1.0/(float)nprebin);
      h2->SetLineColor(colors[3]);
      h2->Draw("hist same");

      scale(h3);
      h3->Scale(1.0/(float)njrebin);
      h3->SetLineColor(colors[5]);
      //h3->Draw("hist same");

      float lowcluster = anaclone.ptBins[i];
      float highcluster = anaclone.ptBins[i+1];
      float minjet = (k == 1 ? anaclone.jet_calib_pt_cut[j] : anaclone.jet_pt_cut[j]);
      float drawx = 0.15;
      if ((index - 1) % 5 == 0) drawx = 0.35;
      drawText(Form("%.0f GeV < p_{T}^{cluster} < %.0f GeV",lowcluster,highcluster),drawx,0.85,1,42);
      drawText(Form("Jet R=0.%i",j*2+2),drawx,0.75,1,42);
      TLine * line = new TLine(minjet/lowcluster,0,minjet/lowcluster,1);
      line->SetLineStyle(8);
      line->Draw();
    }
  }
  p->cd();
  anaclone.drawAll({
      "run 47289-53864",
      //"MC run28 Jet",
      "MC run28 Photon"},
      {
      Form("|vz| < %.0f cm",anaclone.vzcut),
      Form("|#eta|_{jet} < %.1f - R",anaclone.etacut),
      Form("|#eta|_{cluster} < %.1f",anaclone.etacut),
      Form("#Delta#phi > %.0f#pi/%.0f",anaclone.oppnum,anaclone.oppden),
      Form("p_{T}^{jet} > %.0f GeV",(k == 1 ? anaclone.jet_calib_pt_cut[0] : anaclone.jet_pt_cut[0])),
      Form("|#DeltaT|_{jet,cluster} < %.0f",anaclone.tcut),
      calibstring,
      t1[ia].c_str(), t2[ia].c_str(), t3[ia].c_str()},
      drawx,drawy,fontsize,c->GetWh());
  TLegend * l2 = new TLegend(drawx,drawy-.7,0.99,drawy-.55);
  l2->SetLineWidth(0);
  l2->AddEntry(hratio    [0][0][k][l][ia],("data"));
  l2->AddEntry(hratiomc_p[0][0][k][l][ia],("MC Photon reco"));
  //l2->AddEntry(hratiomc_j[0][0][l][ia],("MC Jet reco"));
  l2->Draw();
  return;
}

void draw_one(TCanvas * c, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratio, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_p, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_j, 
    int ipt, int ir, int ia, int k) {
  float drawx = 0.55;
  float drawy = 0.85;
  int i = ipt;
  int j = ir;
  string calibstring = (k == 1 ? "JES Calibrated" : "Uncalibrated Jets");
  vector<string> t1 = {"#bf{Analysis region A}","#bf{Analysis region B}","#bf{Analysis region C}","#bf{Analysis region D}", "Combined ABCD"};
  int colors[6] = {kSpring + 2, kBlue, kGreen + 3, kMagenta+1, kTeal, kSpring+2};
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.2);
  gPad->SetBottomMargin(.15);
  
  TH1D * h1 = hratio    [ipt][ir][k][0][ia];
  TH1D * h2 = hratiomc_p[ipt][ir][k][0][ia];
  TH1D * h3 = hratiomc_j[ipt][ir][k][0][ia];
  h1->GetXaxis()->SetTitle("p_{T,max}^{Jet}/p_{T,max}^{cluster}");
  h1->GetXaxis()->SetTitleSize(0.04);
  h1->GetXaxis()->SetTitleOffset(1.5);
  h1->GetXaxis()->SetLabelSize(0.04);

  h1->GetYaxis()->SetTitle("Arbitrary Units");
  h1->GetYaxis()->SetTitleSize(0.04);
  h1->GetYaxis()->SetTitleOffset(2);
  h1->GetYaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetMaxDigits(3);
  h1->GetYaxis()->SetDecimals(2);
  h1->Draw("hist e");

  h2->Draw("hist e same");

  //h3->Draw("hist e same");

  float lowcluster = anaclone.ptBins[i];
  float highcluster = anaclone.ptBins[i+1];
  float minjet = (k == 1 ? anaclone.jet_calib_pt_cut[j] : anaclone.jet_pt_cut[j]);
  TLine * line = new TLine(minjet/lowcluster,0,minjet/lowcluster,h1->GetMaximum());
  line->SetLineStyle(8);
  line->Draw();
  TLegend * l = new TLegend(.25,.75,.55,.87);
  l->SetLineWidth(0);
  l->SetFillStyle(0);
  l->AddEntry(h1,("data"));
  l->AddEntry(h2,("MC Photon reco"));
  //l->AddEntry(h3,("MC Jet reco"));
  l->Draw();
  if (ia == 0) {
    TF1 * func1 = new TF1("func1","gaus",(int)(minjet/lowcluster/0.04 + 1)*0.04,2);
    TF1 * func2 = new TF1("func2","gaus",(int)(minjet/lowcluster/0.04 + 1)*0.04,2);
    h1->Fit(func1,"RQI");
    h2->Fit(func2,"RQI");
    func1->SetLineColor(h1->GetLineColor());
    func2->SetLineColor(h2->GetLineColor());
    func1->Draw("same");
    func2->Draw("same");
  }
  
  anaclone.drawAll({
      "run 47289-53864",
      //"MC run28 Jet",
      "MC run28 Photon"},
      {Form("%0.0f GeV < p_{T}^{cluster} < %0.0f GeV",anaclone.ptBins[i],anaclone.ptBins[i+1]),
      "p_{T}^{jet} > 3 GeV", 
      Form("Jet R=%0.1f",anaclone.JetRs[j]), 
      "Analysis cuts",
      calibstring,
      t1[ia].c_str()},
    drawx,drawy,15,c->GetWh());
  return;
}


TH1D * combine_hists(TH1D * A, TH1D * B, TH1D * C, TH1D * D, int bin, string name) {
  TH1D * h = (TH1D*)A->Clone(name.c_str());
  h->Reset("ICES");
  float p = anaclone.getPurity(anaclone.ptBins[bin],anaclone.ptBins[bin+1]);
  h->Add(A,C,1/p,-(1-p)/p);
  return h;
}

void draw_abcd() {
  gStyle->SetOptStat(0);
  TFile * f = TFile::Open("hists/histsData.root");
  TFile * f05_p = TFile::Open(Form("hists/hists%s.root","Photon5"));
  TFile * f10_p = TFile::Open(Form("hists/hists%s.root","Photon10"));
  TFile * f20_p = TFile::Open(Form("hists/hists%s.root","Photon20"));
  TFile * f05_j = TFile::Open(Form("hists/hists%s.root","Jet5"));
  TFile * f10_j = TFile::Open(Form("hists/hists%s.root","Jet10"));
  TFile * f20_j = TFile::Open(Form("hists/hists%s.root","Jet20"));
  TFile * f30_j = TFile::Open(Form("hists/hists%s.root","Jet30"));
  TFile * f50_j = TFile::Open(Form("hists/hists%s.root","Jet50"));
  TFile * f70_j = TFile::Open(Form("hists/hists%s.root","Jet70"));

  const char * histname = "hratio3jet";
  cout << "collecting hists..." << endl;
  vector<vector<vector<vector<vector<TH1D*>>>>> hratio     = anaclone.collect_hists({f},{},histname,anaclone.nPtBins,anaclone.nJetR);
  vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_p = anaclone.collect_hists({f05_p,f10_p,f20_p},{5,10,20},histname,anaclone.nPtBins,anaclone.nJetR);
  vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_j = anaclone.collect_hists({f05_j,f10_j,f20_j,f30_j,f50_j,f70_j},{5,10,20,30,50,70},histname,anaclone.nPtBins,anaclone.nJetR);

  cout << "Combining ABDC..." << endl; 
  for (int i = 0; i < anaclone.nPtBins; i++) {
    for (int j = 0; j < anaclone.nJetR ; j++) {
      for (int k = 0; k < anaclone.nCalibBins; k++) {
        for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
          hratio    [i][j][k][l].push_back(combine_hists(hratio    [i][j][k][l][0],hratio    [i][j][k][l][1],hratio    [i][j][k][l][2],hratio    [i][j][k][l][3],i,Form("%s%s%i",hratio[i][j][k][l][0]->GetName(),"data",5)));
          hratiomc_p[i][j][k][l].push_back(combine_hists(hratiomc_p[i][j][k][l][0],hratiomc_p[i][j][k][l][1],hratiomc_p[i][j][k][l][2],hratiomc_p[i][j][k][l][3],i,Form("%s%s%i",hratio[i][j][k][l][0]->GetName(),"photon",5)));
          hratiomc_j[i][j][k][l].push_back(combine_hists(hratiomc_j[i][j][k][l][0],hratiomc_j[i][j][k][l][1],hratiomc_j[i][j][k][l][2],hratiomc_j[i][j][k][l][3],i,Form("%s%s%i",hratio[i][j][k][l][0]->GetName(), "jet",5)));
        }
      }
    }
  }
  

  cout << "Drawing plots..." << endl;
  int csize =  700;
  for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
    TCanvas * c = new TCanvas("c","",csize*anaclone.nJetR+200,csize*anaclone.nPtBins);
    c->SaveAs(Form("pdfs/abcd%i.pdf[",l));
    // 5 for A,B,C,D, and the combined one
    for (int ia = 0; ia < 5; ia++) {
      c->Clear();
      draw_all(c, hratio,hratiomc_p,hratiomc_j,ia,l,0);
      c->SaveAs(Form("pdfs/abcd%i.pdf",l));
    }
    c->SaveAs(Form("pdfs/abcd%i.pdf]",l));
    delete c;
  }
  int i = 5; // pt 15-17 GeV
  int j = 1; // R = 0.4
  TCanvas * cone = new TCanvas("cone","",csize,csize); 
  cone->SaveAs("pdfs/oneabcd.pdf[");
  for (int ia = 0; ia < 5; ia++) {
    cone->Clear();
    draw_one(cone, hratio,hratiomc_p,hratiomc_j,i,j,ia,0);
    cone->SaveAs("pdfs/oneabcd.pdf");
  }
  cone->SaveAs("pdfs/oneabcd.pdf]");
  
  
  
  
  cout << "Drawing calib plots..." << endl;
  for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
    TCanvas * c = new TCanvas("c","",csize*anaclone.nJetR+200,csize*anaclone.nPtBins);
    c->SaveAs(Form("pdfs/abcd%i_calib.pdf[",l));
    // 5 for A,B,C,D, and the combined one
    for (int ia = 0; ia < 5; ia++) {
      c->Clear();
      draw_all(c, hratio,hratiomc_p,hratiomc_j,ia,l,1);
      c->SaveAs(Form("pdfs/abcd%i_calib.pdf",l));
    }
    c->SaveAs(Form("pdfs/abcd%i_calib.pdf]",l));
    delete c;
  }

  TCanvas * cone_calib = new TCanvas("cone_calib","",csize,csize); 
  cone_calib->SaveAs("pdfs/oneabcd_calib.pdf[");
  for (int ia = 0; ia < 5; ia++) {
    cone_calib->Clear();
    draw_one(cone_calib, hratio,hratiomc_p,hratiomc_j,i,j,ia,1);
    cone_calib->SaveAs("pdfs/oneabcd_calib.pdf");
  }
  cone_calib->SaveAs("pdfs/oneabcd_calib.pdf]");

  /*
  for (int i = 0; i < anaclone.nPtBins; i++) {
    for (int j = 0; j < anaclone.nJetR ; j++) {
      for (int k = 0; k < anaclone.nCalibBins; k++) {
        for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
          for (int m = 0; m < 5; m++) {
            delete hratio    [i][j][k][l][m];
            delete hratiomc_p[i][j][k][l][m];
            delete hratiomc_j[i][j][k][l][m];
          }
        }
      }
    }
  }
*/

  //TCanvas * ctest = new TCanvas("ctest","",csize,csize);
  //hratio[6][1][0][0][0]->Draw();
  //hratio[6][1][0][0][4]->Draw("same");

  return;
}
