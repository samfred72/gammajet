#ifndef ana_h
#define ana_h

#include <TROOT.h>
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <map>
#include <TEfficiency.h>
#include <TFile.h>
#include <TF1.h>
using namespace std;

class ana {
  public :
    ana();
    ~ana() = default;

    static Bool_t   PassEtaCut(float eta, float vz); 
    static Double_t GetShiftedEta(float _vz, float _eta);
    static float    getPurity(float low, float high);
    static Int_t findPtBin(double value);
    static Int_t findabcdBin(double iso, double bdt, int bin = 0);
    static Int_t findabcdBin(double iso, int showershape, int bin = 0);
    static Int_t findxjBin(double value);
    static Int_t findBdtBin(double value);

    static constexpr float sPHENIX_posx = 0.6;
    static constexpr float sPHENIX_posy = 0.85;
    static constexpr float posy_diff = 0.05;
    static constexpr size_t MaxClusters = 10000;
    static constexpr size_t MaxJets = 10000;

    static constexpr int nDims=7;
    static constexpr int NHIST = nDims*nDims;
    static constexpr int gridSize = 7;
    static constexpr double vzcut = 60; 
    static constexpr double oppnum = 7;
    static constexpr double oppden = 8;
    static constexpr double oppcut = oppnum*M_PI/oppden;
    static constexpr double tcut = 4;
    static constexpr double tlowcut = 0;
    static constexpr double thighcut = 4;
    static constexpr double radius = 93; 
    static constexpr double etacut = 1.1;
    static constexpr double etamin = -etacut;
    static constexpr double etamax = etacut;
    static constexpr double minclustere = 7;

    static constexpr int nCalibBins = 3;
    static constexpr int nPtBins = 9;
    static constexpr int nIsoBdtBins = 3;
    static constexpr int nBdtBins = 4;
    static constexpr int nxjBins = 3;
    static constexpr int nJetR = 4;
    static constexpr int n3jetBins = 2;
    static constexpr double ptBins[nPtBins+1] = {10,11,12,13,14,15,17,19,25,35};
    static constexpr double isoBins[nIsoBdtBins] = {2,2,2};
    static constexpr double isoBinsHigh[nIsoBdtBins] = {4,4,4};
    static constexpr double bdtGoodHigh[nIsoBdtBins] = {1.0, 1.0, 1.0};
    static constexpr double bdtGoodLow[nIsoBdtBins] = {0.8, 0.9, 0.7};
    static constexpr double bdtBadHigh[nIsoBdtBins] = {0.6, 0.6, 0.6};
    static constexpr double bdtBadLow[nIsoBdtBins] = {0.2, 0.2, 0.2};
    static constexpr double bdtBins[nBdtBins+1] = {0.4,0.7,0.8,0.9,1.0};
    static constexpr double xjBins[nxjBins+1] = {0,0.3,0.7,2.0};
    static constexpr double JetRs[nJetR] = {0.2, 0.4, 0.6, 0.8};
    static constexpr double drcut[nJetR] = {0.2, 0.4, 0.6, 0.8};
    static constexpr double jet_pt_cut[nJetR] = {3,3,3,3};
    static constexpr double jet_calib_pt_cut[nJetR] = {5,5,5,5};
    
  private:
};

#endif
