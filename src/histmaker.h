#ifndef HISTMAKER_H
#define HISTMAKER_H

#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include "/home/samson72/sphnx/gammajet/src/object.h"
#include "/home/samson72/sphnx/gammajet/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet/src/jet_object.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TSystem.h>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
using namespace std;

class histmaker {
  public:
    histmaker(string trigger) {
      this->trigger = trigger;
      f = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/gammajet%s.root",trigger.c_str()),"read");
      t = (TTree*)f->Get("towerntup");
      isMC = (trigger != "Data");
      ifstream file = ifstream("/home/samson72/sphnx/gammajet/src/MbdPmt.corr");
      string line;
      cout << "Getting MBD t0 corrections..." << endl;
      while (getline(file,line)) {
        int irunnum;
        float t0;
        istringstream iss(line);
        iss >> irunnum >> t0;
        t0map[irunnum] = t0;
      }
      t0map[0] = 0;
      treesetup();

      for (int i = 0; i < ana::nPtBins; i++) {
        for (int j = 0; j < ana::nJetR; j++) {
          for (int k = 0; k < ana::nCalibBins; k++) {
            for (int l = 0; l < ana::nIsoBdtBins; l++) {
              for (int m = 0; m < ana::n3jetBins; m++) {
                for (int n = 0; n < 4; n++) {
                  hratio[i][j][k][l][m][n] = new TH1D(Form("hratio_%i_%i_%i_%i_%i_%i",i,j,k,l,m,n),";p_{T}^{jet}/p_{T}^{#gamma};normalized counts",100,0,2);
                }
              }
            }
          }
        }
        hisobdt[i] = new TH2D(Form("hisobdt%i",i),";cluster iso;bdt score",100,-1,20,100,0,1);
      }
      for (int i = 0; i < ana::nJetR; i++) {
        hjetpt[i] =            new TH1D(Form("hjetpt%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetpt[i] =       new TH1D(Form("htruthjetpt%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetptspec[i] =   new TH1D(Form("htruthjetptspec%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetptanti[i] =   new TH1D(Form("htruthjetptanti%i",i),";jet p_{T,max};counts",100,0,100);
        hjetptprecut[i] =      new TH1D(Form("hjetptprecut%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetptprecut[i] = new TH1D(Form("htruthjetptprecut%i",i),";jet p_{T,max};counts",100,0,100);
        hmtminusjt[i] =        new TH1D(Form("hmtminusjt%i",i),";t_{mbd}-t_{jet} [ns]",100,-10,10);
        hctminusjt[i] =        new TH1D(Form("hctminusjt%i",i),";t_{cluster}-t_{jet} [ns]",100,-10,10);
        hiso[i] =              new TH1D(Form("hiso%i",i),";iso E R02;counts",100,-1,50);
        hjeteta[i] =           new TH1D(Form("hjeteta%i",i),";#eta_jet;counts",100,-1.1,1.1);
        hjetetahighem[i] =     new TH1D(Form("hjetetahighem%i",i),";#eta_jet (emfrac > 0.8);counts",100,-1.1,1.1);
        hjetetalowem[i] =      new TH1D(Form("hjetetalowem%i",i),";#eta_jet (emfrac < 0.5);counts",100,-1.1,1.1);
        hdeltar[i] =           new TH1D(Form("hdeltar%i",i),";dr [eta,phi];counts",100,0,4);

        hjetetaphi[i] =        new TH2D(Form("hjetetaphi%i",i),";#eta_{jet};#phi_{jet}",100,-1.1,1.1,100,-M_PI,M_PI);
        hiso2d[i] =            new TH2D(Form("hiso2d%i",i),";cluster E;iso E R02",100,0,50,100,0,50);
        hJES[i] =              new TH2D(Form("hJES%i",i),";Uncalibrated Jet pT;Calib scale",50,0,50,100,0,3);
        hjetsmear[i] =         new TH2D(Form("hjetsmear%i",i),";Calibrated jet pT;smeared jet pT",50,0,50,50,0,50);
        h3jetpt[i] =           new TH2D(Form("h3jetpt%i",i),";Leading cluster pT; Leading 3rd Jet pT",ana::nPtBins, ana::ptBins, 50,0,50);
        for (int j = 0; j < ana::nPtBins; j++) {
          hemfrac[j][i] = new TH1D(Form("hemfrac%i_%i",j,i),";jet emfrac;counts",100,-0.1,1.1);
          hdeltaphi[j][i] = new TH1D(Form("hdeltaphi%i_%i",j,i),";|#phi_{#gamma} - #phi_{leading jet}|;counts",100,0,M_PI);
          hdeltaphiprecut[j][i] = new TH1D(Form("hdeltaphiprecut%i_%i",j,i),";|#phi_{#gamma} - #phi_{leading jet}|;counts",100,0,M_PI);
          h3jetdeltar[j][i] = new TH1D(Form("h3jetdeltar%i_%i",j,i),";#DeltaR_{max jet, 3rd jet};counts",100,0,4);
        }
        for (int j = 0; j < ana::nxjBins; j++) {
          hjetetaxj[j][i] = new TH1D(Form("hjetetaxj%i_%i",j,i),Form(";jet eta %1.1f < xJ < %1.1f;Counts",ana::xjBins[j],ana::xjBins[j+1]),100,-1.1,1.1);
        }
        for (int j = 0; j < ana::nBdtBins; j++) {
          hxjbdt[j][i] = new TH1D(Form("hxjbdt%i_%i",j,i), ";x_{J#gamma};counts", 100, 0, 2);
        } 
      }
      for (int i = 0; i < 11; i++) {
        hbdt[i] = new TH1D(Form("hbdt%i",i),";bdt score;counts",120,-0.1,1.1);
      }
      for (int i = 0; i < 4; i++) {
        hclusterptabcd[i] = new TH1D(Form("hclusterptabcd%i",i),";p_{T}^{lead cluster};Counts",100,0,100);
      }
      
    }
    void treesetup();
    void savehists(TH1D * h[], int n);
    void savehists(TH2D * h[], int n);
    void savehists(TH1D * h[][ana::nJetR], int n, int m);
    void savehists(TH2D * h[][ana::nJetR], int n, int m);
    void end();
    float getZ(pho_object pho, vector<jet_object> jets, bool isiso);
    pho_object getmaxpho(vector<pho_object> phos);
    jet_object getmaxjet(vector<jet_object> jets, pho_object pho, int ij);
    jet_object getthirdjet(pho_object maxpho, jet_object maxjet, vector<jet_object> jets, int ir, int icalib);
    bool loop(jet_object jet, vector<jet_object> jets, int jindex, pho_object pho, int icalib, bool isthirdjet = 0);
    void make_hists();
    
    template <typename T>
      float findmaxpt(const vector<T>& objs)
      {
        float maxpt = -1;
        for (const auto& obj : objs) {
          if (obj.pt > maxpt)
            maxpt = obj.pt;
        }
        return maxpt;
      }

  private:
    TRandom3 * rand = new TRandom3();
    float mbd_t0;
    float t0corr;
    map<int,float> t0map;
    bool isMC;
    TFile * f;
    TTree * t;
    string trigger;
    int count_isc = 0;
    int count_isj[ana::nJetR] = { 0 };
    int nentries = 0;

    // -1: cluster
    //  0: jet R=0.2
    //  1: jet R=0.4
    //  2: jet R=0.6
    //  3: jet R=0.8
    map<int, map<string,int>> threshmap = {
      {-1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10",12},{"Photon20", 24}}},
      { 0,{{"Jet5", 0},{"Jet10",12},{"Jet20",20},{"Jet30",31},{"Jet50",50},{"Jet70", 70},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
      { 1,{{"Jet5", 0},{"Jet10",14},{"Jet20",22},{"Jet30",35},{"Jet50",52},{"Jet70", 71},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
      { 2,{{"Jet5", 0},{"Jet10",17},{"Jet20",35},{"Jet30",45},{"Jet50",63},{"Jet70", 79},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
      { 3,{{"Jet5", 0},{"Jet10",20},{"Jet20",40},{"Jet30",50},{"Jet50",65},{"Jet70", 80},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
    };
    map<int, map<string,int>> threshmap_high = {
      {-1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5",12},{"Photon10",24},{"Photon20",100}}},
      { 0,{{"Jet5",12},{"Jet10",20},{"Jet20",31},{"Jet30",50},{"Jet50",70},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
      { 1,{{"Jet5",14},{"Jet10",22},{"Jet20",35},{"Jet30",52},{"Jet50",71},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
      { 2,{{"Jet5",17},{"Jet10",35},{"Jet20",45},{"Jet30",63},{"Jet50",79},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
      { 3,{{"Jet5",20},{"Jet10",40},{"Jet20",50},{"Jet30",65},{"Jet50",80},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
    };
    map<int, map<string,int>> reco_threshmap_high = {
      {-1,{{"Jet5",15},{"Jet10",20},{"Jet20",30},{"Jet30",40},{"Jet50",60},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
      { 0,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
      { 1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
      { 2,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
      { 3,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}}
    };

    // Define histograms
    TH1D * hratio     [ana::nPtBins][ana::nJetR][ana::nCalibBins][ana::nIsoBdtBins][ana::n3jetBins][4]; // 4 accounts for a,b,c,d regions, 2 is for noncalib vs calib jets

    TH2D * hisobdt        [ana::nPtBins];

    TH1D * hjetetaxj      [ana::nxjBins][ana::nJetR];

    TH1D * hbdt[11]; // 11 is the number of models
    TH1D * hemfrac        [ana::nPtBins][ana::nJetR];
    TH1D * hdeltaphi      [ana::nPtBins][ana::nJetR];
    TH1D * hdeltaphiprecut[ana::nPtBins][ana::nJetR];
    TH1D * h3jetdeltar    [ana::nPtBins][ana::nJetR];

    TH1D * hxjbdt        [ana::nBdtBins][ana::nJetR]; 

    TH1D * hmtminusjt           [ana::nJetR];
    TH1D * hctminusjt           [ana::nJetR];
    TH1D * hiso                 [ana::nJetR];
    TH1D * hjetpt               [ana::nJetR];
    TH1D * htruthjetpt          [ana::nJetR];
    TH1D * htruthjetptspec      [ana::nJetR];
    TH1D * htruthjetptanti      [ana::nJetR];
    TH1D * hjetptprecut         [ana::nJetR];
    TH1D * htruthjetptprecut    [ana::nJetR];
    TH1D * hjeteta              [ana::nJetR];
    TH1D * hjetetahighem        [ana::nJetR];
    TH1D * hjetetalowem         [ana::nJetR];
    TH1D * hdeltar              [ana::nJetR];
    TH1D * hclusterptabcd       [4]; // for ABCD

    TH2D * hjetetaphi           [ana::nJetR];
    TH2D * hiso2d               [ana::nJetR];
    TH2D * hJES                 [ana::nJetR];
    TH2D * hjetsmear            [ana::nJetR];
    TH2D * h3jetpt              [ana::nJetR];


    TH1D * hfrag = new TH1D("hfrag",";cluster Z without iso;counts",100,0,1);
    TH1D * hfragiso = new TH1D("hfragiso",";cluster Z with iso;counts",100,0,1);
    TH1D * hmbdt = new TH1D("hmbdt",";time [ns]; counts",100,-10,10);
    TH1D * hclustert = new TH1D("hclustert",";time [ns]; counts",100,-10,10);
    TH1D * hjett = new TH1D("hjett",";time [ns]; counts",100,-10,10);
    TH1D * hmtminusct = new TH1D("hmtminusct",";t_{mbd}-t_{cluster} [ns]",100,-10,10);
    TH1D * hclusterptprecut = new TH1D("hclusterptprecut",";cluster p_{T,max};counts",100,0,100);
    TH1D * hclusterpt = new TH1D("hclusterpt",";p_{T}^{lead cluster};Counts",100,0,100);
    TH1D * htruthclusterptprecut = new TH1D("htruthclusterptprecut",";cluster p_{T,max};counts",100,0,100);
    TH1D * htruthclusterpt = new TH1D("htruthclusterpt",";cluster p_{T,max};counts",100,0,100);
    TH1D * hclustereta = new TH1D("hclustereta",";leading cluster #eta;counts",100,-1.2,1.2);
    TH2D * hclusteretaphi = new TH2D("hclusteretaphi",";leading cluster #eta;leading cluster #phi",100,-1.2,1.2,100,-M_PI,M_PI);
    TH2D * hmct = new TH2D("hmct",";mbd time [ns]; cluster time [ns]",100,-10,10,100,-10,10);
    TH2D * hmjt = new TH2D("hmjt",";mbd time [ns]; jet time [ns]",100,-10,10,100,-10,10);
    TH2D * hcjt = new TH2D("hcjt",";cluster time [ns]; jet time [ns]",100,-10,10,100,-10,10);

    // TTree variables
    const int MaxClusters = 10000;
    const int MaxJets = 10000;
    Int_t           RunNumber;
    Float_t         vz;
    Bool_t          ScaledTriggerBit[64];
    Bool_t          LiveTriggerBit[64];
    Int_t           Scaledowns[64];

    Int_t           mbd_nhits_south;
    Int_t           mbd_nhits_north;
    Float_t         mbd_time_south;
    Float_t         mbd_time_north;

    Short_t                  nClusters;
    vector<Float_t>         *cluster_pt;
    vector<Float_t>         *cluster_e;
    vector<Float_t>         *cluster_eta;
    vector<Float_t>         *cluster_phi;
    vector<Float_t>         *cluster_time;
    vector<vector<Float_t>> *cluster_bdt_scores;
    vector<vector<Float_t>> *cluster_showershape;

    Short_t                  nTruthClusters;
    vector<Float_t>         *truth_cluster_pt;
    vector<Float_t>         *truth_cluster_e;
    vector<Float_t>         *truth_cluster_eta;
    vector<Float_t>         *truth_cluster_phi;

    Short_t                  nJets[4];
    vector<vector<Float_t>> *jet_pt;
    vector<vector<Float_t>> *jet_pt_calib;
    vector<vector<Float_t>> *jet_pt_smear;
    vector<vector<Float_t>> *jet_e;
    vector<vector<Float_t>> *jet_eta;
    vector<vector<Float_t>> *jet_phi;
    vector<vector<Float_t>> *jet_emfrac;
    vector<vector<Float_t>> *jet_ihfrac;
    vector<vector<Float_t>> *jet_ohfrac;
    vector<vector<Float_t>> *jet_time;  

    Short_t                  nTruthJets[4];
    vector<vector<Float_t>> *truth_jet_pt;
    vector<vector<Float_t>> *truth_jet_e;
    vector<vector<Float_t>> *truth_jet_eta;
    vector<vector<Float_t>> *truth_jet_phi;

    // List of branches
    TBranch        *b_RunNumber;   //!
    TBranch        *b_vz;   //!
    TBranch        *b_ScaledTriggerBit;   //!
    TBranch        *b_LiveTriggerBit;   //!
    TBranch        *b_Scaledowns;   //!
    TBranch        *b_mbd_nhits_south;   //!
    TBranch        *b_mbd_nhits_north;   //!
    TBranch        *b_mbd_time_south;   //!
    TBranch        *b_mbd_time_north;   //!

    TBranch        *b_nClusters;   //!
    TBranch        *b_cluster_pt;
    TBranch        *b_cluster_e;
    TBranch        *b_cluster_eta;
    TBranch        *b_cluster_phi;
    TBranch        *b_cluster_showershape;   //!
    TBranch        *b_cluster_time;   //!
    TBranch        *b_cluster_bdt_scores;   //!

    TBranch        *b_nTruthClusters;   //!
    TBranch        *b_truth_cluster_pt;
    TBranch        *b_truth_cluster_e;
    TBranch        *b_truth_cluster_eta;
    TBranch        *b_truth_cluster_phi;

    TBranch        *b_nJets;   //!
    TBranch        *b_jet_pt;                           //!
    TBranch        *b_jet_pt_calib;                           //!
    TBranch        *b_jet_pt_smear;                           //!
    TBranch        *b_jet_e;                           //!
    TBranch        *b_jet_eta;                           //!
    TBranch        *b_jet_phi;                           //!
    TBranch        *b_jet_emfrac;   //!
    TBranch        *b_jet_ihfrac;   //!
    TBranch        *b_jet_ohfrac;   //!
    TBranch        *b_jet_time;   //!

    TBranch        *b_nTruthJets;   //!
    TBranch        *b_truth_jet_pt;                           //!
    TBranch        *b_truth_jet_e;                           //!
    TBranch        *b_truth_jet_eta;                           //!
    TBranch        *b_truth_jet_phi;                           //!
};
#endif // HISTMAKER_H
