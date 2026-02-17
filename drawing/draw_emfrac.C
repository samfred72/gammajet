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

void drawMany(vector<TH1D*> hists, vector<string> labels, string axistitle) {
  if (hists.size() != labels.size()) {
    cout << "Labels don't match up with hists!" << endl;
    return;
  }
  if (hists.size() > 6) {
    cout << "too many hists!!"; 
    return;
  }
  int manycolors[3] = {kBlack,kBlue,kRed};
  int manystyles[3] = {20, 24, 24};
  TLegend * l = new TLegend(.25,.45,.4,.7);
  TCanvas * c = new TCanvas(Form("c%s",labels.at(0).c_str()),"",700,700);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  gPad->SetTicks(1,1);
  for (int i = 0; i < hists.size(); i++) {
    scale(hists.at(i));
    hists.at(i)->SetMarkerColor(manycolors[i]);
    hists.at(i)->SetMarkerStyle(manystyles[i]);
    hists.at(i)->SetMarkerSize(1);
    hists.at(i)->GetXaxis()->SetTitle(axistitle.c_str());
    hists.at(i)->Draw("p same");
    l->AddEntry(hists.at(i),labels.at(i).c_str());
  }
  l->SetLineWidth(0);
  l->Draw();
}

void drawData(TFile * f, string histname, string label) {
  ana anaclone;
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

void drawMC(vector<TFile*> files, vector<int> samples, string histname, string label, bool isphoton) {
  if (files.size() != samples.size()) return;
  ana anaclone;
  vector<TH1D*> hists;
  for (int i = 0; i < files.size(); i++) {
    hists.push_back((TH1D*)files.at(i)->Get(histname.c_str()));
    hists.at(i)->SetName(Form("%s%i%i",histname.c_str(),i,(int)isphoton));
  }
  TH1D * thehist = anaclone.combineMC(hists,samples,isphoton);
  TLegend * l = new TLegend(.7,.45,.85,.7);
  TCanvas * c = new TCanvas(Form("c%s%i",histname.c_str(),(int)isphoton),"",700,700);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
  gPad->SetTicks(1,1);
  string trigger = (isphoton ? "Photon" : "Jet");
  
  thehist->SetMarkerColor(kBlack);
  thehist->SetMarkerStyle(20);
  thehist->SetMarkerSize(1);
  thehist->GetYaxis()->SetTitle("Scaled counts");
  thehist->GetXaxis()->SetTitle(label.c_str());
  thehist->GetXaxis()->SetTitleOffset(1.2);
  thehist->Draw("p same");
  l->AddEntry(thehist,"Sum");

  for (int i = 0; i < hists.size(); i++) {
    hists.at(i)->Scale(anaclone.scalemap[isphoton][samples.at(i)]);
    hists.at(i)->SetMarkerColor(colors[i]);
    hists.at(i)->SetMarkerStyle(24);
    hists.at(i)->SetMarkerSize(1);
    l->AddEntry(hists.at(i), Form("%s %i",trigger.c_str(),samples.at(i)));
    hists.at(i)->Draw("p same");
  }
  l->SetLineWidth(0);
  l->Draw();
}
    
void draw_emfrac() {
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
  
  float drawx = .25;
  float drawy = .85;
  float fontsize = 20;

  const char * histname1 = "hjeteta1";
  TH1D * hd1 = (TH1D*)f->Get(histname1);
  TH1D * hp1 = getcombinedhist({f05p,f10p,f20p},{5,10,20},histname1,1);
  const char * histname2 = "hjetetahighem1";
  TH1D * hd2 = (TH1D*)f->Get(histname2);
  TH1D * hp2 = getcombinedhist({f05p,f10p,f20p},{5,10,20},histname2,1);
  const char * histname3 = "hjetetalowem1";
  TH1D * hd3 = (TH1D*)f->Get(histname3);
  TH1D * hp3 = getcombinedhist({f05p,f10p,f20p},{5,10,20},histname3,1);
  
  //drawMany({hd1,hp1},{"Data","MC Photon"},"Paired Jet #eta");
  //anaclone.drawAll({},{"|vz| < 60","Jet R=0.4",Form("%0.0f GeV < cluster p_{T} < %0.0f GeV",anaclone.ptBins[4],anaclone.ptBins[5]),"Analysis cuts"},drawx,drawy,fontsize,700);
  //drawMany({hd2,hp2},{"Data","MC Photon"},"Paired Jet #eta");
  //anaclone.drawAll({},{"|vz| < 60","Jet R=0.4",Form("%0.0f GeV < cluster p_{T} < %0.0f GeV",anaclone.ptBins[4],anaclone.ptBins[5]),"Analysis cuts","Jet EMfrac > 0.8"},drawx,drawy,fontsize,700);
  drawMany({hd3,hp3},{"Data","MC Photon"},"Paired Jet #eta");
  anaclone.drawAll({},{"|vz| < 60","Jet R=0.4",Form("%0.0f GeV < cluster p_{T} < %0.0f GeV",anaclone.ptBins[4],anaclone.ptBins[5]),"Analysis cuts","Jet EMfrac < 0.5"},drawx,drawy,fontsize,700);
}
