#include "../headers/commonUtility.h"
#include "../headers/Style_jaebeom.h"
#include "../headers/ana.cxx"

int colors[6] = {kRed, kBlue, kGreen-2, kMagenta, kOrange, kTeal};

void scale(TH1D * h) {
  double x_low_value = 0;
  double x_high_value = 100;

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

TH1D * getcombinedhist(vector<TFile*> files, vector<int> samples, string histname, bool isphoton) {
  if (files.size() != samples.size()) return nullptr;
  ana anaclone;
  vector<TH1D*> hists;
  for (int i = 0; i < files.size(); i++) {
    hists.push_back((TH1D*)files.at(i)->Get(histname.c_str()));
  }
  TH1D * thehist = anaclone.combineMC(hists,samples,isphoton);
  return thehist;
}

TH1D * combine_hists(TH1D * A, TH1D * B, TH1D * C, TH1D * D) {
  TH1D * hnum = (TH1D*)A->Clone();
  hnum->SetName("combinednum");
  TH1D * hden = (TH1D*)A->Clone();
  hden->SetName("combinedden");
  hnum->Reset("ICES");
  hden->Reset("ICES");
  hnum->Multiply(B,C);
  hden->Multiply(A,D);
  
  TH1D * h = (TH1D*)A->Clone();
  h->SetName("combined");
  h->Reset("ICES");
  h->Divide(hnum,hden);
  for (int i = 1; i < h->GetNbinsX(); i++) {
    if (h->GetBinContent(i) > 0) h->SetBinContent(i,1-h->GetBinContent(i));
  }
  return h;
}

    
void puritymaker() {
  ana anaclone;
  gStyle->SetOptStat(0);
  TFile * f = TFile::Open("hists/histsData.root");
  TFile * f05p = TFile::Open("hists/histsPhoton5.root");
  TFile * f10p = TFile::Open("hists/histsPhoton10.root");
  TFile * f20p = TFile::Open("hists/histsPhoton20.root");
  TFile * f05j = TFile::Open("hists/histsJet5.root");
  TFile * f10j = TFile::Open("hists/histsJet10.root");
  TFile * f20j = TFile::Open("hists/histsJet20.root");
  TFile * f30j = TFile::Open("hists/histsJet30.root");
  TFile * f50j = TFile::Open("hists/histsJet50.root");
  TFile * f70j = TFile::Open("hists/histsJet70.root");
  
  const char * histname = "hclusterptabcd";
  TH1D * hd[4]; // for ABCD
  TH1D * hp[4];
  for (int i = 0; i < 4; i++) {
    hd[i] = (TH1D*)f->Get(Form("%s%i",histname,i));
    hd[i]->SetName(Form("hdata%i",i));
    hp[i] = getcombinedhist({f05p,f10p,f20p},{5,10,20},Form("%s%i",histname,i),1);
    hp[i]->SetName(Form("hphoton%i",i));
    hd[i]->Sumw2();
    hp[i]->Sumw2();
  }
  TH1D * od = combine_hists(hd[0],hd[1],hd[2],hd[3]);
  TF1 * func = new TF1("func","[0]*TMath::Erf((x - [1])/[2])",8,30);
  func->SetParameter(0,1);
  func->SetParameter(1,13);
  func->SetParameter(2,5);
  od->Fit(func,"RIMQ0");
  //od->Draw();
  //TH1D * op = combine_hists(hp[0],hp[1],hp[2],hp[3]);
  TFile * fout = TFile::Open("hists/purity.root","RECREATE");
  od->Write();
  func->Write();
  for (int i = 0; i < 4; i++) {
    hd[i]->Write();
    hp[i]->Write();
  }
  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetTicks();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  hd[0]->SetLineColor(kBlack);
  hd[1]->SetLineColor(kBlue);
  hd[2]->SetLineColor(kOrange);
  hd[3]->SetLineColor(kRed);
  TLegend * l = new TLegend(.5,.4,.8,.6);
  l->SetLineWidth(0);
  string text[4] = {"Region A","Region B","Region C","Region D"};
  for (int i = 0; i < 4; i++) {
    hd[i]->SetLineWidth(2);
    hd[i]->GetYaxis()->SetRangeUser(0.5,1e5);
    hd[i]->Draw("hist same");
    l->AddEntry(hd[i],text[i].c_str());
  }
  l->Draw();
  anaclone.drawAll({"p+p Run24"},{"Paired clusters"},.5,.8,20,700);
}
