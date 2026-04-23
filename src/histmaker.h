#ifndef HISTMAKER_H
#define HISTMAKER_H

#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include "/home/samson72/sphnx/gammajet/src/object.h"
#include "/home/samson72/sphnx/gammajet/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet/src/jet_object.h"
#include "/home/samson72/sphnx/gammajet/src/treeuser.h"
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

class histmaker : public treeuser {
  public:
    histmaker(string trigger, string sim = "pythia") : treeuser(trigger, sim) {
      //ifstream file = ifstream("/home/samson72/sphnx/gammajet/src/MbdPmt.corr");
      //string line;
      //cout << "Getting MBD t0 corrections..." << endl;
      //while (getline(file,line)) {
      //  int irunnum;
      //  float t0;
      //  istringstream iss(line);
      //  iss >> irunnum >> t0;
      //  t0map[irunnum] = t0;
      //}
      //t0map[0] = 0;

      cluster_smear_file = TFile::Open("/home/samson72/sphnx/gammajet/trees/function_compare_wide.root");
      cluster_energy_smear_func = (TF1*)cluster_smear_file->Get("f_energy_quadrdiff_cE0p08");
      cluster_position_smear_func = (TF1*)cluster_smear_file->Get("f_position_quadrdiff_cE0p08");

      cout << "initializing hists..." << endl;
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
        hjetphi[i] =           new TH1D(Form("hjetphi%i",i),";#phi_jet;counts",100,-3*M_PI,3*M_PI);
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
          hdeltaphiprecut[j][i] = new TH1D(Form("hdeltaphiprecut%i_%i",j,i),";|#phi_{#gamma} - #phi_{leading jet}|;counts",100,-5*M_PI,2*M_PI);
          h3jetdeltar[j][i] = new TH1D(Form("h3jetdeltar%i_%i",j,i),";#DeltaR_{max jet, 3rd jet};counts",100,0,4);
        }
        for (int j = 0; j < ana::nxjBins; j++) {
          hjetetaxj[j][i] = new TH1D(Form("hjetetaxj%i_%i",j,i),Form(";jet eta %1.1f < xJ < %1.1f;Counts",ana::xjBins[j],ana::xjBins[j+1]),100,-1.1,1.1);
        }
        for (int j = 0; j < ana::nBdtBins; j++) {
          hxjbdt[j][i] = new TH1D(Form("hxjbdt%i_%i",j,i), ";x_{J#gamma};counts", 100, 0, 2);
        } 
        for (int j = 0; j < ana::nHadronBins; j++) {
          hhadronp[j][i] = new TH1D(Form("hhadronp%i_%i",j,i), ";max hadron p_{T};counts",100,0,100);
        }
      }
      for (int i = 0; i < 11; i++) {
        hbdt[i] = new TH1D(Form("hbdt%i",i),";bdt score;counts",120,-0.1,1.1);
      }
      for (int i = 0; i < 4; i++) {
        hclusterptabcd[i] = new TH1D(Form("hclusterptabcd%i",i),";p_{T}^{lead cluster};Counts",100,0,100);
      }
      
      treefilename = (trigger == "Data" ? Form("/home/samson72/sphnx/gammajet/hists/tree_%s.root",trigger.c_str()) : Form("/home/samson72/sphnx/gammajet/hists/tree_%s_%s.root",sim.c_str(),trigger.c_str()));
      outfile = TFile::Open(treefilename,"RECREATE");
      for (int ir = 0; ir < ana::nJetR; ir++) {
        outtree[ir] = new TTree(Form("xjtree_%i",ir),"");
        outtree[ir]->Branch("pho_pt", &outtree_pho_pt[ir]);
        outtree[ir]->Branch("jet_pt", &outtree_jet_pt[ir]);
        outtree[ir]->Branch("weight", &outtree_weight[ir]);

        outtree_3jet[ir] = new TTree(Form("xjtree_3jet_%i",ir),"");
        outtree_3jet[ir]->Branch("pho_pt", &outtree_pho_pt_3jet[ir]);
        outtree_3jet[ir]->Branch("jet_pt", &outtree_jet_pt_3jet[ir]);
        outtree_3jet[ir]->Branch("weight", &outtree_weight_3jet[ir]);

        outtree_bdt[ir] = new TTree(Form("xjtree_bdt_%i",ir),"");
        outtree_bdt[ir]->Branch("pho_pt", &outtree_pho_pt_bdt[ir]);
        outtree_bdt[ir]->Branch("jet_pt", &outtree_jet_pt_bdt[ir]);
        outtree_bdt[ir]->Branch("weight", &outtree_weight_bdt[ir]);

        outtree_iso[ir] = new TTree(Form("xjtree_iso_%i",ir),"");
        outtree_iso[ir]->Branch("pho_pt", &outtree_pho_pt_iso[ir]);
        outtree_iso[ir]->Branch("jet_pt", &outtree_jet_pt_iso[ir]);
        outtree_iso[ir]->Branch("weight", &outtree_weight_iso[ir]);

        if (isMC) {
          outtree_JERhigh[ir] = new TTree(Form("xjtree_JERhigh_%i",ir),"");
          outtree_JERhigh[ir]->Branch("pho_pt", &outtree_pho_pt_JERhigh[ir]);
          outtree_JERhigh[ir]->Branch("jet_pt", &outtree_jet_pt_JERhigh[ir]);
          outtree_JERhigh[ir]->Branch("weight", &outtree_weight_JERhigh[ir]);

          outtree_JERlow[ir] = new TTree(Form("xjtree_JERlow_%i",ir),"");
          outtree_JERlow[ir]->Branch("pho_pt", &outtree_pho_pt_JERlow[ir]);
          outtree_JERlow[ir]->Branch("jet_pt", &outtree_jet_pt_JERlow[ir]);
          outtree_JERlow[ir]->Branch("weight", &outtree_weight_JERlow[ir]);
        }
      }
    }
    void savehists(TH1D * h[], int n);
    void savehists(TH2D * h[], int n);
    void savehists(TH1D * h[][ana::nJetR], int n, int m);
    void savehists(TH2D * h[][ana::nJetR], int n, int m);
    void end();
    bool loop(jet_object jet, int jindex, pho_object pho, int icalib, float weight = 1);
    void make_hists();
    float reweight(float pt, float vz) { return reweight_func_pt->Eval(pt)*reweight_hist_vz->GetBinContent(reweight_hist_vz->FindBin(vz)); }
    
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
    
    TFile * reweight_file = TFile::Open("/home/samson72/sphnx/gammajet/hists/reweighting.root");
    TF1 * reweight_func_pt = (TF1*)reweight_file->Get("pfunc");
    TH1D * reweight_hist_vz = (TH1D*)reweight_file->Get("hvz_reweight");
    
    
    TRandom3 * rand = new TRandom3();
    float mbd_t0;
    float t0corr;
    //map<int,float> t0map;
    int nentries = 0;
    const char * treefilename;
    TFile * outfile;
    TTree * outtree[ana::nJetR];
    float outtree_pho_pt[ana::nJetR];
    float outtree_jet_pt[ana::nJetR];
    float outtree_weight[ana::nJetR];
    
    TTree * outtree_3jet[ana::nJetR];
    float outtree_pho_pt_3jet[ana::nJetR];
    float outtree_jet_pt_3jet[ana::nJetR];
    float outtree_weight_3jet[ana::nJetR];
    
    TTree * outtree_bdt[ana::nJetR];
    float outtree_pho_pt_bdt[ana::nJetR];
    float outtree_jet_pt_bdt[ana::nJetR];
    float outtree_weight_bdt[ana::nJetR];

    TTree * outtree_iso[ana::nJetR];
    float outtree_pho_pt_iso[ana::nJetR];
    float outtree_jet_pt_iso[ana::nJetR];
    float outtree_weight_iso[ana::nJetR];

    TTree * outtree_JERhigh[ana::nJetR];
    float outtree_pho_pt_JERhigh[ana::nJetR];
    float outtree_jet_pt_JERhigh[ana::nJetR];
    float outtree_weight_JERhigh[ana::nJetR];

    TTree * outtree_JERlow[ana::nJetR];
    float outtree_pho_pt_JERlow[ana::nJetR];
    float outtree_jet_pt_JERlow[ana::nJetR];
    float outtree_weight_JERlow[ana::nJetR];

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

    TH1D * hhadronp [ana::nHadronBins][ana::nJetR];

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
    TH1D * hjetphi              [ana::nJetR];
    TH1D * hjetetahighem        [ana::nJetR];
    TH1D * hjetetalowem         [ana::nJetR];
    TH1D * hdeltar              [ana::nJetR];
    TH1D * hclusterptabcd       [4]; // for ABCD

    TH2D * hjetetaphi           [ana::nJetR];
    TH2D * hiso2d               [ana::nJetR];
    TH2D * hJES                 [ana::nJetR];
    TH2D * hjetsmear            [ana::nJetR];
    TH2D * h3jetpt              [ana::nJetR];

    TH1D * hvz = new TH1D("hvz",";vz [cm];counts",400,-1000,1000);

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
    TH1D * hclusterphi = new TH1D("hclusterphi",";leading cluster #phi;counts",100,-3*M_PI,3*M_PI);
    TH2D * hclusteretaphi = new TH2D("hclusteretaphi",";leading cluster #eta;leading cluster #phi",100,-1.2,1.2,100,-M_PI,M_PI);
    TH2D * hmct = new TH2D("hmct",";mbd time [ns]; cluster time [ns]",100,-10,10,100,-10,10);
    TH2D * hmjt = new TH2D("hmjt",";mbd time [ns]; jet time [ns]",100,-10,10,100,-10,10);
    TH2D * hcjt = new TH2D("hcjt",";cluster time [ns]; jet time [ns]",100,-10,10,100,-10,10);
    TH1D * hratiosingle = new TH1D("hratiosingle",";x_{J#gamma};counts",100,0,2); 

    TFile * cluster_smear_file;
    TF1 * cluster_position_smear_func; 
    TF1 * cluster_energy_smear_func; 

};
#endif // HISTMAKER_H
