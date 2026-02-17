#ifndef DRAWER_H
#define DRAWER_H

#include "drawer.h"
using namespace std;

drawer::~drawer() {}

void drawer::drawLine(float x1, float y1, float x2, float y2) {
  TLine * l = new TLine(x1,y1,x2,y2);
  l->SetLineStyle(4);
  l->Draw();
}

void drawer::drawText(const char *text, float xp, float yp, int textColor, int textSize) {
   TLatex *tex = new TLatex(xp,yp,text);
   tex->SetTextFont(43);
   //   if(bold)tex->SetTextFont(43);
   tex->SetTextSize(textSize);
   tex->SetTextColor(textColor);
   tex->SetLineWidth(1);
   tex->SetNDC();
   tex->Draw();
}

void drawer::drawMany(vector<string> features, float drawx, float drawy, int fontsize, float csize) {
  float ydiff = fontsize * 0.002 * 700.0/csize;
  for (unsigned i = 0; i < features.size(); i++) {
    drawText(features.at(i).c_str(),drawx,drawy-ydiff*i,1,fontsize);
  }
}

void drawer::drawAll(vector<string> samples, vector<string> features, float drawx, float drawy, int fontsize, float csize) {
  float titlescale = 1.5;// * 700.0/csize;
  float subtitlescale = 1.25;// * 700.0/csize;
  float ydiff = fontsize * 0.002 * 700.0/csize;
  drawText("#bf{#it{sPHENIX}} Internal",drawx,drawy,1,(int)fontsize*titlescale);
  for (unsigned i = 0; i < samples.size(); i++) {
    drawText(samples.at(i).c_str(),drawx,drawy-ydiff*subtitlescale*(i+1),1,(int)fontsize*subtitlescale);
  }
  for (unsigned i = 0; i < features.size(); i++) {
    drawText(features.at(i).c_str(),drawx,drawy-ydiff*subtitlescale*samples.size()-ydiff*(i+1),1,fontsize);
  }
}

void drawer::format(TH1D * h, int type) {
  int colors[3] = {kBlue, kMagenta+1, kSpring+2};
  h->SetLineColor(colors[type]);
  scale(h,0,2);
  h->GetYaxis()->SetRangeUser(0,h->GetMaximum()*1.5);
  if (type == 0) h->SetLineWidth(2);
  else h->SetLineWidth(1);
}
void drawer::format(TF1 * f, int type) {
  int colors[3] = {kBlue, kMagenta+1, kSpring+2};
  f->SetLineColor(colors[type]);
  f->SetLineWidth(1);
}

void drawer::scale(TH1D * h, float low, float high) {
  // Get the X-axis object
  TAxis *xaxis = h->GetXaxis();

  // Find the corresponding bin numbers
  Int_t bin_low = xaxis->FindBin(low);
  Int_t bin_high = xaxis->FindBin(high);

  // Sum the bin contents within the range
  double entries_in_range = 0;
  for (Int_t bin = bin_low; bin <= bin_high; ++bin) {
    entries_in_range += h->GetBinContent(bin);
  }
  h->Scale(1.0/entries_in_range);
}

TF1 * drawer::fit(TH1D * h, float low, float high) {
  TF1 * func = new TF1(Form("func_%s",h->GetName()), "gaus", low, high);
  func->SetParameters(0.1,0.7,0.3);
  if (h->GetEntries() > 0) h->Fit(func,"RQIM0");
  return func;
}


TH1D * drawer::combineMC(const char * histname, bool isphoton) {
  vector<TH1D*> hists;            
  vector<TFile*> files = (isphoton ? pfiles : jfiles);
  int nfiles = (isphoton ? npfiles : njfiles);
  vector<int> samples = (isphoton ? psamples : jsamples);
  
  for (int ifile = 0; ifile < nfiles; ifile++) {
    hists.push_back((TH1D*)files[ifile]->Get(histname));
    hists[ifile]->SetName(Form("%s_%i_%i",histname,ifile,isphoton));
  }

  TH1D * thehist = (TH1D*)hists.at(0)->Clone();
  thehist->Reset("ICES");
  for (unsigned i = 0; i < hists.size(); i++) {
    hists[i]->Scale(scalemap[isphoton][samples.at(i)]);
    thehist->Add(hists[i]);
    delete hists[i];
  }
  return thehist;
}
TH2D * drawer::combineMC2d(const char * histname, bool isphoton) {
  vector<TH2D*> hists;            
  vector<TFile*> files = (isphoton ? pfiles : jfiles);
  int nfiles = (isphoton ? npfiles : njfiles);
  vector<int> samples = (isphoton ? psamples : jsamples);
  
  for (int ifile = 0; ifile < nfiles; ifile++) {
    hists.push_back((TH2D*)files[ifile]->Get(histname));
    hists[ifile]->SetName(Form("%s_%i_%i",histname,ifile,isphoton));
  }

  TH2D * thehist = (TH2D*)hists.at(0)->Clone();
  thehist->Reset("ICES");
  for (unsigned i = 0; i < hists.size(); i++) {
    hists[i]->Scale(scalemap[isphoton][samples.at(i)]);
    thehist->Add(hists[i]);
    delete hists[i];
  }
  return thehist;
}

TH1D * drawer::combine_hists(TH1D * A, TH1D * B, TH1D * C, TH1D * D, int ipt, string name) { 
  TH1D * h = (TH1D*)A->Clone(name.c_str());
  h->Reset("ICES");
  float p = ana::getPurity(ana::ptBins[ipt],ana::ptBins[ipt+1]);
  h->Add(A,C,1/p,-(1-p)/p);
  return h;
}
vector<vector<vector<vector<vector<vector<TH1D*>>>>>> drawer::collect_hists(const char * histname, int type) {
  vector<vector<vector<vector<vector<vector<TH1D*>>>>>> hists(
      ana::nPtBins, 
      vector<vector<vector<vector<vector<TH1D*>>>>>(
        ana::nJetR, 
        vector<vector<vector<vector<TH1D*>>>>(
          ana::nCalibBins, 
          vector<vector<vector<TH1D*>>>(
            ana::nIsoBdtBins,
            vector<vector<TH1D*>>(
              ana::n3jetBins,
              vector<TH1D*>(
                5)))))); // 5 for ABCD and the combined one

  bool isphoton = type == 1;
  int nrebin = 4;
  int njrebin = 1;
  int nprebin = 1;
  int mcrebin = (isphoton ? nprebin : njrebin);

  for (int i = 0; i < ana::nPtBins; i++) {
    for (int j = 0; j < ana::nJetR; j++) {
      for (int k = 0; k < ana::nCalibBins; k++) {
        for (int l = 0; l < ana::nIsoBdtBins; l++) {
          for (int m = 0; m < ana::n3jetBins; m++) {
            for (int n = 0; n < 4; n++) {
              if (type) { // It's an MC sample
                hists[i][j][k][l][m][n] = combineMC(Form("%s_%i_%i_%i_%i_%i_%i",histname,i,j,k,l,m,n),isphoton);
                hists[i][j][k][l][m][n]->Rebin(mcrebin*nrebin);
              }
              else { // it's data
                hists[i][j][k][l][m][n] = (TH1D*)dfiles[0]->Get(Form("%s_%i_%i_%i_%i_%i_%i", histname,i,j,k,l,m,n));
                hists[i][j][k][l][m][n]->Rebin(nrebin); 
              }
              hists[i][j][k][l][m][n]->SetName(Form("%s_%i_%i_%i_%i_%i_%i_%i", histname,i,j,k,l,m,n,type)); 
            }
            hists[i][j][k][l][m][4] = combine_hists(
                hists[i][j][k][l][m][0],
                hists[i][j][k][l][m][1],
                hists[i][j][k][l][m][2],
                hists[i][j][k][l][m][3],
                i,Form("%s_%i_%i_%i_%i_%i_%i_%i", histname,i,j,k,l,m,4,type));
          }
        }
      }
    }
  }
  return hists;
}

#endif // DRAWER_H
