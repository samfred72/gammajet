#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"
R__LOAD_LIBRARY(libgammajet.so);
    
void draw_test() {
  ana anaclone;
  drawer drawer;
  gStyle->SetOptStat(0);

  float drawx = .25;
  float drawy = .85;
  float fontsize = 20;

  const char * histname = "hratio3jet";
  vector<vector<vector<vector<vector<TH1D*>>>>> hratiod = drawer.collect_hists(histname,0); // 0 for data
  vector<vector<vector<vector<vector<TH1D*>>>>> hratiop = drawer.collect_hists(histname,1); // 1 for photon
  vector<vector<vector<vector<vector<TH1D*>>>>> hratioj = drawer.collect_hists(histname,2); // 2 for jet
                                                                                            //
                                                                                            //

  for (int i = 0; i < ana::nPtBins; i++) {
    for (int j = 0; j < ana::nJetR; j++) {
      for (int k = 0; k < ana::nCalibBins; k++) {
        for (int l = 0; l < ana::nIsoBdtBins; l++) {
          for (int m = 0; m < 5; m++) {
            cout << hratiod[i][j][k][l][m]->GetName() << endl;
            cout << hratiop[i][j][k][l][m]->GetName() << endl;
            cout << hratioj[i][j][k][l][m]->GetName() << endl;
          }
        }
      }
    }
  }

}
