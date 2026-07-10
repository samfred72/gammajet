#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"

TH1D * bin_to_num(TH1D * hin) {
  TH1D * hout = new TH1D(Form("%s_num",hin->GetName()),";x_{j#gamma};Normalized Counts",ana::nUnfoldXjBins,ana::unfoldXjBins);
  for (int i = 0; i < ana::nUnfoldXjBins; i++) {
    hout->SetBinContent(i+1,hin->GetBinContent(i+1)/hin->GetBinWidth(i+1));
    hout->SetBinError(i+1,hin->GetBinError(i+1)/hin->GetBinWidth(i+1));
  }
  return hout;
}

vector<TH1D*> split_histo(TH1D * hist) {
  vector<TH1D*> v;
  const char * name = hist->GetName();
  for (int i = 0; i < ana::nUnfoldPtBins; i++) {
    TH1D * h = new TH1D(Form("%s_%i",name,i),"",ana::nUnfoldXjBins,ana::unfoldXjBins);
    for (int j = 0; j < ana::nUnfoldXjBins; ++j) {
      int global = i*(ana::nUnfoldXjBins+2) + (j+2);

      h->SetBinContent(j+1,hist->GetBinContent(global));
      h->SetBinError(j+1, hist->GetBinError(global));
    }
    v.push_back(h);
  }
  return v;
}

void draw_many(TCanvas * c, 
  vector<TH1D*> hmc_reco,          
  vector<TH1D*> hmc_truth,          
  vector<TH1D*> hmc_unfold,          
  vector<TH1D*> hdata,          
  vector<TH1D*> hdata_unfold) {       
  c->cd();

  drawer d;
  float drawx = 0.77;
  float drawy = 0.85;
  float fontsize = 50;
  int offset = 0;
  int ir = 2;
  vector<string> t1 = {"#bf{Analysis region A}","#bf{Analysis region B}","#bf{Analysis region C}","#bf{Analysis region D}", "Combined ABCD"};

  TPad * p = new TPad("p","",0,0,1,1);
  p->SetLeftMargin(0.25);
  p->SetBottomMargin(0.15);
  p->Divide(2,1,0,0);
  //p->Divide(4,1,0,0);
  p->Draw();

  for (int ipt = 0; ipt < 1; ipt++) { // Specifically not using nUnfoldPtBins because I only care about the 3x3 here
  //for (int ipt = 0; ipt < ana::nUnfoldPtBins-1; ipt++) { // Specifically not using nUnfoldPtBins because I only care about the 3x3 here
    TH1D * h1 = (hmc_truth[ipt+1]);
    TH1D * h2 = (hdata[ipt+1]);
    TH1D * h3 = (hdata_unfold[ipt+1]);
    TH1D * h4 = (hmc_reco[ipt+1]);

    int index = ipt + 1 + offset;
    if (index % (3+1) == 0) {index++; offset++;}
    p->cd(index);
    if ((index + 1) % (3+1) == 0) gPad->SetRightMargin(0.01);
    gPad->SetTicks(1,1);
    if (index == (3*(3+1) - 1)) h1->GetXaxis()->SetTitle("p_{T,max}^{Jet}/p_{T,max}^{cluster}");
    else h1->GetXaxis()->SetTitle("");
    h1->GetXaxis()->SetTitleSize(0.10);
    h1->GetXaxis()->SetLabelSize(0.08);
    if (index == (3+1)*(3-1)+1) h1->GetXaxis()->SetLabelSize(0.07);
    h1->GetXaxis()->SetNdivisions(405,false);
    h1->GetXaxis()->ChangeLabel(-1,-1,0);
    if (index == 1) h1->GetYaxis()->SetTitle("Arbitrary Units");
    else h1->GetYaxis()->SetTitle("");
    h1->GetYaxis()->SetTitleSize(0.10);
    h1->GetYaxis()->SetLabelSize(0.08);
    if (index == (3+1)*(ana::nUnfoldPtBins-1)+1) h1->GetYaxis()->SetLabelSize(0.07);
    h1->GetYaxis()->SetLabelOffset(0.04);
    h1->GetYaxis()->SetMaxDigits(3);
    h1->GetYaxis()->SetDecimals(2);

    h1->Draw("p");
    h2->Draw("p same");
    h3->Draw("p same");
    //h4->Draw("p same");

    float lowcluster = ana::unfoldPtBins[ipt];
    float highcluster = ana::unfoldPtBins[ipt+1];
    float drawx = 0.15;
    if ((index - 1) % 4 == 0) drawx = 0.35;
    d.drawMany({
        Form("#bf{%.0f GeV < p_{T}^{cluster} < %.0f GeV}",lowcluster,highcluster),
        },drawx,0.88,42,c->GetWh()/3.0);
  }
  p->cd();
  d.drawAll(
      {
      "p+p #sqrt{s}=200 GeV"
      },
      {
      Form("Jet R=%1.1f",ana::JetRs[ir]),
      },
      drawx,drawy,fontsize,c->GetWh());
  TLegend * l2 = new TLegend(drawx,drawy-.6,0.99,drawy-.3);
  l2->SetLineWidth(0);
  //l2->AddEntry(hmc_reco[0], "MC reco");
  l2->AddEntry(hmc_truth[0], "MC truth");
  l2->AddEntry(hdata[0], "Data");
  l2->AddEntry(hdata_unfold[0], "Data unfolded");
  l2->Draw();
  c->SaveAs("/home/samson72/sphnx/gammajet/pdfs/unfolding/full_unfold_comp.pdf");
  return;
}

void draw_simple() {
  drawer d(1,"pythia");
  gStyle->SetOptStat(0);

  int itype = 1; // 1: photon sample, 2: jet sample
  int isample = -1; // -1: all, 0: photon5, 1: photon10, etc.
  int ir = 2;
  //const char * half = "";
  const char * half = "_half";
  TH1D * ht_full = d.get(Form("htruthxj%s%i" , half, ir), itype,isample);
  TH1D * hr_full = d.get(Form("hrecoxj%s%i"  , half, ir), itype,isample);
  TH1D * hu_full = d.get(Form("hunfoldxj%s%i", half, ir), itype,isample);
  
  vector<TH1D*> ht = split_histo(ht_full);
  vector<TH1D*> hr = split_histo(hr_full);
  vector<TH1D*> hu = split_histo(hu_full);
  vector<string> labels = {"Reco x_{J#gamma}", "Truth x_{J#gamma}", "Unfolded x_{J#gamma}"};
 
  TCanvas * c = new TCanvas("c"  ,"",700,1000);
  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/unfolding/unfolded_xj%s.pdf[",half));

  for (int i = 0; i < ana::nUnfoldPtBins; i++) {
    c->Clear();
    TPad * p1 = new TPad(Form("p1%i",i) ,"",0,.4,1,1);
    TPad * p2 = new TPad(Form("p2_%i",i),"",0,0,1,.4);
    p1->Draw();
    p2->Draw();
    p1->cd();
    p1->SetBottomMargin(0);
    p2->SetTopMargin(0);
    gPad->SetTicks(1,1);
    vector<TH1D*> hists = {bin_to_num(hr[i]),bin_to_num(ht[i]),bin_to_num(hu[i])};

    int manycolors[3] = {kBlue,kBlack,kRed};
    int manystyles[3] = {20, 20, 24};
    TLegend * l = new TLegend(.57,.71,.87,.86);
    gPad->SetTicks(1,1);
    for (int i = 0; i < hists.size(); i++) {
      hists.at(i)->SetMarkerColor(manycolors[i]);
      hists.at(i)->SetMarkerStyle(manystyles[i]);
      hists.at(i)->SetMarkerSize(1);
      hists.at(i)->GetYaxis()->SetRangeUser(0,hists.at(1)->GetMaximum()*1.5);
      hists.at(i)->GetXaxis()->SetTitle("x_{J#gamma}");
      hists.at(i)->GetYaxis()->SetTitle("Normalized Counts");
      hists.at(i)->Draw("p same");
      l->AddEntry(hists.at(i),labels.at(i).c_str());
    }
    l->SetLineWidth(0);
    l->Draw();
    d.drawAll({"MC Photon"},{"|vz| < 60",Form("R=%0.1f",ana::JetRs[ir]),Form("%.1f GeV < p_{T,#gamma} < %.1f GeV",ana::unfoldPtBins[i],ana::unfoldPtBins[i+1])},0.18,0.85,20,700);
    gPad->SetLeftMargin(.15);

    p2->cd();
    gPad->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hdivide = (TH1D*)hists.at(1)->Clone();
    hdivide->Reset("ICES");
    hdivide->Divide(hists.at(2),hists.at(1));
    hdivide->SetLineColor(kBlack);
    hdivide->SetMarkerColor(kBlack);
    hdivide->SetMarkerSize(1);
    hdivide->SetMarkerStyle(20);
    hdivide->GetYaxis()->SetRangeUser(0.5,1.5);
    hdivide->GetYaxis()->SetTitle("Ratio");
    hdivide->Draw();
    TLine * line = new TLine(0,1,27,1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/unfolding/unfolded_xj%s.pdf",half));

  }
  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/unfolding/unfolded_xj%s.pdf]",half));
}

    
void draw_xj_unfold() {
  if (!gROOT->IsBatch()) {
    //cout << "Run with -b flag or else!!" << endl;
    gROOT->SetBatch(kTRUE);
    //return;
  }
  
  draw_simple();
  //return;


  drawer du(1,"pythia"); // 1 to tell it it's the unfolding version
  drawer dn;
  gStyle->SetOptStat(0);


  int itype = 1; // 1: photon sample, 2: jet sample
  int isample = -1; // -1: all, 0: photon5, 1: photon10, etc.
  int ir = 2;
  const char * half = "";
  int nbins = ana::nUnfoldPtBins*(ana::nUnfoldXjBins + 2); // +2 for underflow and overflow bins
  // Collect flattened 2D histos
  TH1D * hmc_truth_full = du.get(Form("htruthxj%s%i" , half, ir), itype,isample);
  TH1D * hmc_reco_full = du.get(Form("hrecoxj%s%i"  , half, ir), itype,isample); 
  TH1D * hmc_unfold_full = du.get(Form("hunfoldxj%s%i", half, ir), itype,isample);
  TH1D * hdata_full = du.get(Form("hrecoxj%s%i", half, ir), 0);
  //TH1D * hdata_full = new TH1D("hdata_full","",nbins,0,nbins);
  TH2D * hresponse = du.get2d(Form("hxjresponse%s%i", half, ir), itype, isample);
  
  // Unfold
  RooUnfoldResponse * response = new RooUnfoldResponse(hmc_reco_full,hmc_truth_full,hresponse);
  RooUnfoldBayes unfold(response, hdata_full, 4, 0, 1); // four iterations
  TH1D * hdata_unfold_full = (TH1D*)unfold.Hreco();
  
  // split up flattened 2D hists
  vector<TH1D*> hdata_bin = split_histo(hdata_full);
  vector<TH1D*> hdata_unfold_bin = split_histo(hdata_unfold_full);
  vector<TH1D*> hmc_reco_bin = split_histo(hmc_reco_full);
  vector<TH1D*> hmc_truth_bin = split_histo(hmc_truth_full);
  vector<TH1D*> hmc_unfold_bin = split_histo(hmc_unfold_full);
  
  vector<TH1D*> hdata       ;
  vector<TH1D*> hdata_unfold;
  vector<TH1D*> hmc_reco    ;
  vector<TH1D*> hmc_truth   ;
  vector<TH1D*> hmc_unfold  ;
  for (int i = 0; i < ana::nUnfoldPtBins; i++) {
    hdata.push_back(bin_to_num(hdata_bin[i]));
    hdata_unfold.push_back(bin_to_num(hdata_unfold_bin[i]));
    hmc_reco.push_back(bin_to_num(hmc_reco_bin[i]));
    hmc_truth.push_back(bin_to_num(hmc_truth_bin[i]));
    hmc_unfold.push_back(bin_to_num(hmc_unfold_bin[i]));
  }
  // Format
  // The list goes: 1: MC reco, 2: MC truth, 3: MC unfold, 4: data, 5: data unfold 
  int colors[5] = {kOrange,kBlack,kRed,kBlue,kGreen};
  int styles[5] = {24, 20, 24, 21, 21};
  for (int i = 0; i < ana::nUnfoldPtBins; i++) { // Not doing the extra bin because data doesn't have it
    TH1D * hists[5] = {hmc_reco[i], hmc_truth[i], hmc_unfold[i], hdata[i], hdata_unfold[i]};
    for (int j = 0; j < 5; j++) {
      hists[j]->SetLineColor(colors[j]);
      hists[j]->SetLineWidth(1);
      hists[j]->SetMarkerColor(colors[j]);
      hists[j]->SetMarkerStyle(styles[j]);
      hists[j]->SetMarkerSize(2);
      hists[j]->Scale(1.0/hists[j]->Integral("width"));
      hists[j]->GetYaxis()->SetRangeUser(0,2.7);
    }
  }
  
  //TCanvas * cblah = new TCanvas("cblah","",1000,600);
  //hmc_truth[3]->Draw("hist");
  //gPad->SetLogy();
  //cblah->SaveAs("/home/samson72/sphnx/gammajet/pdfs/unfolding/blah.pdf");
  //return;
  

  int csize = 700;
  //TCanvas * c = new TCanvas("c","",csize*4,csize);
  TCanvas * c = new TCanvas("c","",csize*2,csize);
  draw_many(c, hmc_reco, hmc_truth, hmc_unfold, hdata, hdata_unfold);
  gApplication->Terminate();

}
