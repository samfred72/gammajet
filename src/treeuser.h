#ifndef TREEUSER_H
#define TREEUSER_H

#include "/home/samson72/sphnx/gammajet/src/ana.h"
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

class treeuser {
  public:
    treeuser(string trigger, string sim = "pythia") {
      this->trigger = trigger;
      this->sim = sim;
      f = (trigger == "Data" ? 
          TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/gammajet_%s.root",trigger.c_str()),"read") :
          TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/gammajet_%s_%s.root",sim.c_str(),trigger.c_str()),"read")
          );
      t = (TTree*)f->Get("towerntup");
      isMC = (trigger != "Data");
      treesetup();
      // -1: cluster
      //  0: jet R=0.2
      //  1: jet R=0.4
      //  2: jet R=0.6
      //  3: jet R=0.8
      threshmap = (sim == "pythia" ? 
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10",12},{"Photon20", 24}}},
        { 0,{{"Jet5", 0},{"Jet10",12},{"Jet20",20},{"Jet30",31},{"Jet50",50},{"Jet70", 70},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5", 0},{"Jet10",14},{"Jet20",22},{"Jet30",35},{"Jet50",52},{"Jet70", 71},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5", 0},{"Jet10",17},{"Jet20",35},{"Jet30",45},{"Jet50",63},{"Jet70", 79},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5", 0},{"Jet10",20},{"Jet20",40},{"Jet30",50},{"Jet50",65},{"Jet70", 80},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      } : 
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 20}}},
        { 0,{{"Jet5", 0},{"Jet10",12},{"Jet20",20},{"Jet30",31},{"Jet50",50},{"Jet70", 70},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5", 0},{"Jet10",14},{"Jet20",22},{"Jet30",35},{"Jet50",52},{"Jet70", 71},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5", 0},{"Jet10",17},{"Jet20",35},{"Jet30",45},{"Jet50",63},{"Jet70", 79},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5", 0},{"Jet10",20},{"Jet20",40},{"Jet30",50},{"Jet50",65},{"Jet70", 80},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      }
      ); 
      threshmap_high = (sim == "pythia" ? 
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5",12},{"Photon10",24},{"Photon20",100}}},
        { 0,{{"Jet5",12},{"Jet10",20},{"Jet20",31},{"Jet30",50},{"Jet50",70},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5",14},{"Jet10",22},{"Jet20",35},{"Jet30",52},{"Jet50",71},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5",17},{"Jet10",35},{"Jet20",45},{"Jet30",63},{"Jet50",79},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5",20},{"Jet10",40},{"Jet20",50},{"Jet30",65},{"Jet50",80},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      } :
      std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5",15},{"Photon10",20},{"Photon20",100}}},
        { 0,{{"Jet5",12},{"Jet10",20},{"Jet20",31},{"Jet30",50},{"Jet50",70},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 1,{{"Jet5",14},{"Jet10",22},{"Jet20",35},{"Jet30",52},{"Jet50",71},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 2,{{"Jet5",17},{"Jet10",35},{"Jet20",45},{"Jet30",63},{"Jet50",79},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}},
        { 3,{{"Jet5",20},{"Jet10",40},{"Jet20",50},{"Jet30",65},{"Jet50",80},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20",  0}}}
      }
      );
      reco_threshmap_high = std::map<int, std::map<std::string, int>>{
        {-1,{{"Jet5",15},{"Jet10",20},{"Jet20",30},{"Jet30",40},{"Jet50",60},{"Jet70",200},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 0,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 1,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 2,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}},
        { 3,{{"Jet5", 0},{"Jet10", 0},{"Jet20", 0},{"Jet30", 0},{"Jet50", 0},{"Jet70",  0},{"Photon5", 0},{"Photon10", 0},{"Photon20", 0}}}
      };
    }
    void treesetup();
    vector<bool> check_keep_MC(float pho_pt, float jet_pt[], string trigger);

    bool isMC;
    TFile * f;
    TTree * t;
    string trigger;
    string sim;
    const static int nJetR = ana::nJetR;
    map<int, map<string,int>> threshmap;
    map<int, map<string,int>> threshmap_high;
    map<int, map<string,int>> reco_threshmap_high;


    // TTree variables
    Int_t           RunNumber;
    Float_t         vz;
    Bool_t          ScaledTriggerBit[64];
    Bool_t          LiveTriggerBit[64];
    Int_t           Scaledowns[64];

    Float_t         mbd_time;

    Float_t         cluster_pt;
    Float_t         cluster_e;
    Float_t         cluster_eta;
    Float_t         cluster_phi;
    Float_t         cluster_time;
    Float_t         cluster_bdt_scores[11];
    Float_t         cluster_showershape[10];

    Float_t         truth_cluster_pt;
    Float_t         truth_cluster_e;
    Float_t         truth_cluster_eta;
    Float_t         truth_cluster_phi;
    Float_t         truth_cluster_iso3;
    Float_t         truth_cluster_iso4;

    Float_t         jet_pt[nJetR];
    Float_t         jet_pt_calib[nJetR];
    Float_t         jet_pt_smear[nJetR];
    Float_t         jet_pt_smear_high[nJetR];
    Float_t         jet_pt_smear_low[nJetR];
    Float_t         jet_e[nJetR];
    Float_t         jet_eta[nJetR];
    Float_t         jet_phi[nJetR];
    Float_t         jet_emfrac[nJetR];
    Float_t         jet_time[nJetR]; 

    Bool_t          hasthirdjet[nJetR];
    Float_t         thirdjet_pt[nJetR];
    Float_t         thirdjet_dr[nJetR];

    Float_t         hadron_p[nJetR];

    Float_t         truth_jet_pt[nJetR];
    Float_t         truth_jet_e[nJetR];
    Float_t         truth_jet_eta[nJetR];
    Float_t         truth_jet_phi[nJetR];

    // List of branches
    TBranch        *b_RunNumber;   //!
    TBranch        *b_vz;   //!
    TBranch        *b_ScaledTriggerBit;   //!
    TBranch        *b_LiveTriggerBit;   //!
    TBranch        *b_Scaledowns;   //!
    TBranch        *b_mbd_time;   //!

    TBranch        *b_cluster_pt;
    TBranch        *b_cluster_e;
    TBranch        *b_cluster_eta;
    TBranch        *b_cluster_phi;
    TBranch        *b_cluster_bdt_scores;   //!
    TBranch        *b_cluster_showershape;   //!
    TBranch        *b_cluster_time;   //!

    TBranch        *b_truth_cluster_pt;
    TBranch        *b_truth_cluster_e;
    TBranch        *b_truth_cluster_eta;
    TBranch        *b_truth_cluster_phi;
    TBranch        *b_truth_cluster_iso3;
    TBranch        *b_truth_cluster_iso4;

    TBranch        *b_jet_pt;                           //!
    TBranch        *b_jet_pt_calib;                           //!
    TBranch        *b_jet_pt_smear;                           //!
    TBranch        *b_jet_pt_smear_high;                           //!
    TBranch        *b_jet_pt_smear_low;                           //!
    TBranch        *b_jet_e;                           //!
    TBranch        *b_jet_eta;                           //!
    TBranch        *b_jet_phi;                           //!
    TBranch        *b_jet_emfrac;   //!
    TBranch        *b_jet_time;   //!
    
    TBranch        *b_hasthirdjet;   //!
    TBranch        *b_thirdjet_pt;   //!
    TBranch        *b_thirdjet_dr;   //!

    TBranch        *b_hadron_p;   //!

    TBranch        *b_truth_jet_pt;                           //!
    TBranch        *b_truth_jet_e;                           //!
    TBranch        *b_truth_jet_eta;                           //!
    TBranch        *b_truth_jet_phi;                           //!
};
#endif // TREEUSER_H
