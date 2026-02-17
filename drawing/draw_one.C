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


void draw_one_plot(TCanvas * c, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratio, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_p, 
    vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_j, 
    int ipt, int ir, int ia, int k) {
  int nprebin = 1;
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
  TH1D * h2 = hratiomc_p[ipt][ir][k][0][0];
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
  scale(h1);
  h1->GetYaxis()->SetRangeUser(0,h1->GetMaximum()*1.5);
  h1->SetLineColor(colors[1]);
  h1->SetLineWidth(2);
  h1->Draw("hist e");

  scale(h2);
  h2->Scale(1.0/(float)nprebin);
  h2->SetLineColor(colors[3]);
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
  TF1 * func1 = new TF1(Form("func1%i",ia),"gaus",(int)(minjet/lowcluster/0.04 + 1)*0.04,2);
  TF1 * func2 = new TF1(Form("func2%i",ia),"gaus",(int)(minjet/lowcluster/0.04 + 1)*0.04,2);
  if (ia == 0 || ia == 4) {
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


TH1D * combine_hists(TH1D * A, TH1D * B, TH1D * C, TH1D * D, int bin) {
  TH1D * h = (TH1D*)A->Clone();
  h->Clear();
  float p = anaclone.getPurity(anaclone.ptBins[bin],anaclone.ptBins[bin+1]);
  h->Add(A,C,1/p,-(1-p)/p);
  return h;
}

void draw_one() {
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

  cout << "collecting hists..." << endl;
  vector<vector<vector<vector<vector<TH1D*>>>>> hratio     = anaclone.collect_hists({f},{},"hratio",anaclone.nPtBins,anaclone.nJetR);
  vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_p = anaclone.collect_hists({f05_p,f10_p,f20_p},{5,10,20},"hratio",anaclone.nPtBins,anaclone.nJetR);
  vector<vector<vector<vector<vector<TH1D*>>>>> hratiomc_j = anaclone.collect_hists({f05_j,f10_j,f20_j,f30_j,f50_j,f70_j},{5,10,20,30,50,70},"hratio",anaclone.nPtBins,anaclone.nJetR);

  cout << "Combining ABDC..." << endl; 
  for (int i = 0; i < anaclone.nPtBins; i++) {
    for (int j = 0; j < anaclone.nJetR ; j++) {
      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
          hratio    [i][j][k][l].push_back(combine_hists(hratio    [i][j][k][l][0],hratio    [i][j][k][l][1],hratio    [i][j][k][l][2],hratio    [i][j][k][l][3],i));
          hratiomc_p[i][j][k][l].push_back(combine_hists(hratiomc_p[i][j][k][l][0],hratiomc_p[i][j][k][l][1],hratiomc_p[i][j][k][l][2],hratiomc_p[i][j][k][l][3],i));
          hratiomc_j[i][j][k][l].push_back(combine_hists(hratiomc_j[i][j][k][l][0],hratiomc_j[i][j][k][l][1],hratiomc_j[i][j][k][l][2],hratiomc_j[i][j][k][l][3],i));
        }
      }
    }
  }
  

  int csize =  700;
  int i = 5; // pt 15-17 GeV
  int j = 1; // R = 0.4
 
  /*
  TCanvas * ctest = new TCanvas("c","",700,700);
  gPad->SetTicks();
  gStyle->SetOptStat(0);
  hratio[5][1][0][0][0]->SetLineColor(kBlue);
  hratio[5][1][0][0][0]->SetLineWidth(2);
  hratio[5][1][0][0][4]->SetLineColor(kRed);
  hratio[5][1][0][0][4]->SetLineWidth(2);
  TLegend * l = new TLegend(.25,.75,.55,.87);
  l->AddEntry(hratio[5][1][0][0][0],"Region A");
  l->AddEntry(hratio[5][1][0][0][4],"Combined ABCD");
  scale(hratio[5][1][0][0][0]);
  scale(hratio[5][1][0][0][4]);
  hratio[5][1][0][0][0]->GetYaxis()->SetRangeUser(0,hratio[5][1][0][0][0]->GetMaximum()*1.5);
  hratio[5][1][0][0][0]->Draw("hist same");
  hratio[5][1][0][0][4]->Draw("hist same");
  l->SetLineWidth(0);
  l->Draw();
  float drawx = 0.55;
  float drawy = 0.85;
  anaclone.drawAll({
      "run 47289-53864"
      },
      {Form("%0.0f GeV < p_{T}^{cluster} < %0.0f GeV",anaclone.ptBins[i],anaclone.ptBins[i+1]),
      "p_{T}^{jet} > 3 GeV", 
      Form("Jet R=%0.1f",anaclone.JetRs[j]), 
      "Analysis cuts",
      "uncalibrated jets"
      },
    drawx,drawy,15,ctest->GetWh());
  */
   TCanvas * cone = new TCanvas("cone","",csize,csize); 
  cone->SaveAs("pdfs/oneabcd.pdf[");
  for (int ia = 0; ia < 5; ia++) {
    cone->Clear();
    draw_one_plot(cone, hratio,hratiomc_p,hratiomc_j,i,j,ia,0);
    cone->SaveAs("pdfs/oneabcd.pdf");
  }
  cone->SaveAs("pdfs/oneabcd.pdf]");
  
  cout << "Drawing calib plots..." << endl;
  TCanvas * cone_calib = new TCanvas("cone_calib","",csize,csize); 
  cone_calib->SaveAs("pdfs/oneabcd_calib.pdf[");
  for (int ia = 0; ia < 5; ia++) {
    cone_calib->Clear();
    draw_one_plot(cone_calib, hratio,hratiomc_p,hratiomc_j,i,j,ia,1);
    cone_calib->SaveAs("pdfs/oneabcd_calib.pdf");
  }
  cone_calib->SaveAs("pdfs/oneabcd_calib.pdf]");


  return;
}
