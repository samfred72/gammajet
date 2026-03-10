#ifndef UNFOLDER_H
#define UNFOLDER_H

#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include "/home/samson72/sphnx/gammajet/src/object.h"
#include "/home/samson72/sphnx/gammajet/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet/src/jet_object.h"
#include <string>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
using namespace std;

class unfolder {
  public:
    unfolder(string trigger) {
      this->trigger = trigger;
      f = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/gammajet%s.root",trigger.c_str()),"read");
      t = (TTree*)f->Get("towerntup");
      treesetup();
      pho_response = RooUnfoldResponse(100, 0, 100);
      for (int i = 0; i < ana::nJetR; i++) {
        jet_response[i] = RooUnfoldResponse(100, 0, 100);
        hjetpt[i] = new TH1D(Form("hjetpt%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetpt[i] = new TH1D(Form("htruthjetpt%i",i),";jet p_{T,max};counts",100,0,100);
      }
    }
    ~unfolder(); 
    void fill_matrix();
    void treesetup();
    void end();
    pho_object getmaxpho(vector<pho_object> phos);
    jet_object getmaxjet(vector<jet_object> jets, pho_object pho, int ir);
    bool check_pair(jet_object jet, vector<jet_object> jets, int ir, pho_object pho);
    bool check_match(pho_object pr, pho_object pt, jet_object jr, jet_object jt);
    
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
    TFile * f;
    TTree * t;
    string trigger;
    int count_isc = 0;
    int count_isj[ana::nJetR] = { 0 };
    int nentries = 0;
    
    RooUnfoldResponse pho_response;
    RooUnfoldResponse jet_response[ana::nJetR];
    
    TH1D * hjetpt[ana::nJetR];
    TH1D * htruthjetpt[ana::nJetR];
    TH1D * hclusterpt = new TH1D("hclusterpt",";p_{T}^{lead cluster};Counts",100,0,100);
    TH1D * htruthclusterpt = new TH1D("htruthclusterpt",";cluster p_{T,max};counts",100,0,100);
    TH1D * hclusterreco;
    TH1D * hjetreco[ana::nJetR];
    
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

#endif // UNFOLDER_H
