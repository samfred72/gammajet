#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
void makeforjaebeom() {
  drawer d;
  gStyle->SetOptStat(0);
  TFile * f = TFile::Open("forjaebeom.root","RECREATE");

  TH1D * hd[ana::nPtBins];
  TH1D * hp[ana::nPtBins];
  for (int i = 0; i < ana::nPtBins; i++) {
    hd[i] = d.get(Form("hratio_%i_1_1_0_0_0",i),0);
    hd[i]->SetName(Form("xj_ptbin%i_data",i));
    hp[i] = d.get(Form("hratio_%i_1_2_0_0_0",i),1);
    hp[i]->SetName(Form("xj_ptbin%i_MCphoton",i));
  }
  for (int i = 0; i < ana::nPtBins; i++) {
    hd[i]->Write();
  }
  for (int i = 0; i < ana::nPtBins; i++) {
    hp[i]->Write();
  }
}
