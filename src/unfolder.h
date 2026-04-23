#ifndef UNFOLDER_H
#define UNFOLDER_H

#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include "/home/samson72/sphnx/gammajet/src/object.h"
#include "/home/samson72/sphnx/gammajet/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet/src/jet_object.h"
#include "/home/samson72/sphnx/gammajet/src/treeuser.h"
#include <string>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TMarker.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
using namespace std;

class unfolder : public treeuser {
  public:
    unfolder(string trigger) : treeuser(trigger) {
      treesetup();
      for (int i = 0; i < ana::nJetR; i++) {
        hphodr[i] = new TH1D(Form("hphodr%i",i),";#DeltaR_{truth,reco};counts", 100,0,0.4);
        hjetdr[i] = new TH1D(Form("hjetdr%i",i),";#DeltaR_{truth,reco};counts",100,0,0.4);
    
        hphoeffnum[i] = new TH1D(Form("hphoeffnum%i",i),";jet p_{T};counts",100,0,100);
        hphoeffden[i] = new TH1D(Form("hphoeffden%i",i),";jet p_{T};counts",100,0,100);
        hjeteffnum[i] = new TH1D(Form("hjeteffnum%i",i),";jet p_{T};counts",100,0,100);
        hjeteffden[i] = new TH1D(Form("hjeteffden%i",i),";jet p_{T};counts",100,0,100);
        
        hphopurnum[i] = new TH1D(Form("hphopurnum%i",i),";jet p_{T};counts",100,0,100);
        hphopurden[i] = new TH1D(Form("hphopurden%i",i),";jet p_{T};counts",100,0,100);
        hjetpurnum[i] = new TH1D(Form("hjetpurnum%i",i),";jet p_{T};counts",100,0,100);
        hjetpurden[i] = new TH1D(Form("hjetpurden%i",i),";jet p_{T};counts",100,0,100);
        
        hpaireffnum[i] = new TH1D(Form("hpaireffnum%i",i),";jet p_{T};counts",100,0,100);
        hpaireffden[i] = new TH1D(Form("hpaireffden%i",i),";jet p_{T};counts",100,0,100);
        hpairpurnum[i] = new TH1D(Form("hpairpurnum%i",i),";jet p_{T};counts",100,0,100);
        hpairpurden[i] = new TH1D(Form("hpairpurden%i",i),";jet p_{T};counts",100,0,100);

        hphomissfake[i] = new TH2D(Form("hphomissfake%i",i),";reco p_{T};truth p_{T}",101,0,101,101,0,101);
        hjetmissfake[i] = new TH2D(Form("hjetmissfake%i",i),";reco p_{T};truth p_{T}",101,0,101,101,0,101);
        hpairmissfake[i] = new TH2D(Form("hpairmissfake%i",i),";reco p_{T};truth p_{T}",nbins+1,0,nbins+1,nbins+1,0,nbins+1);
        
        hrecojetpt[i] = new TH1D(Form("hrecojetpt%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetpt[i] = new TH1D(Form("htruthjetpt%i",i),";jet p_{T,max};counts",100,0,100);
        jet_response[i] = new RooUnfoldResponse(hrecojetpt[i], htruthjetpt[i]);
        
        hrecophopt[i] = new TH1D(Form("hrecophopt%i",i),";pho p_{T,max};counts",100,0,100);
        htruthphopt[i] = new TH1D(Form("htruthphopt%i",i),";pho p_{T,max};counts",100,0,100);
        pho_response[i] = new RooUnfoldResponse(hrecophopt[i], htruthphopt[i]);
        
        hrecojetpt_half[i] = new TH1D(Form("hrecojetpt_half%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetpt_half[i] = new TH1D(Form("htruthjetpt_half%i",i),";jet p_{T,max};counts",100,0,100);
        jet_response_half[i] = new RooUnfoldResponse(hrecojetpt_half[i], htruthjetpt_half[i]);
        
        hrecophopt_half[i] = new TH1D(Form("hrecophopt_half%i",i),";pho p_{T,max};counts",100,0,100);
        htruthphopt_half[i] = new TH1D(Form("htruthphopt_half%i",i),";pho p_{T,max};counts",100,0,100);
        pho_response_half[i] = new RooUnfoldResponse(hrecophopt_half[i], htruthphopt_half[i]);
       
        hrecoxj[i] = new TH1D(Form("hrecoxj%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
        htruthxj[i] = new TH1D(Form("htruthxj%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
        hunfoldxj[i] = new TH1D(Form("hunfoldxj%i",i),";unfold cluster p_{T}; x_{J#gamma}",nbins,0,nbins);
        jet_response2D[i] = new RooUnfoldResponse(hrecoxj[i], htruthxj[i],Form("response_full_jetR%d",i),Form("response_%d",i));

        hrecoxj_half[i] = new TH1D(Form("hrecoxj_half%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
        htruthxj_half[i] = new TH1D(Form("htruthxj_half%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
        hunfoldxj_half[i] = new TH1D(Form("hunfoldxj_half%i",i),";unfold cluster p_{T}; x_{J#gamma}",nbins,0,nbins);
        jet_response_half2D[i] = new RooUnfoldResponse(hrecoxj_half[i], htruthxj_half[i]);

      }
    }
    ~unfolder(); 
    void fill_matrix();
    void unfold();
    void end();
    void savehists(TH1D * h[], int n);
    void savehists(TH2D * h[], int n);
    void savehists(RooUnfoldResponse * h[], int n);
    void savehists(TEfficiency * h[], int n);
    bool check_pair(jet_object jet, int ir, pho_object pho, bool isreco);
    bool check_match(pho_object p1, pho_object p2, jet_object j1, jet_object j2);
    bool check_match(pho_object p1, pho_object p2);
    bool check_match(jet_object j1, jet_object j2);
    void set_dodraw(bool draw) {dodraw = draw; }
    
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
    float dodraw = false;
    int count_isc = 0;
    int count_isj[ana::nJetR] = { 0 };
    int nentries = 0;
    int niterate = 1;
    TRandom rand;
    
    RooUnfoldResponse * pho_response[ana::nJetR];
    RooUnfoldResponse * jet_response[ana::nJetR];
    RooUnfoldResponse * pho_response_half[ana::nJetR];
    RooUnfoldResponse * jet_response_half[ana::nJetR];

    TH1D * hphodr[ana::nJetR];
    TH1D * hjetdr[ana::nJetR];

    TH1D * hphoeffnum[ana::nJetR];
    TH1D * hphoeffden[ana::nJetR];
    TEfficiency * hphoeff[ana::nJetR];
    TH1D * hjeteffnum[ana::nJetR];
    TH1D * hjeteffden[ana::nJetR];
    TEfficiency * hjeteff[ana::nJetR];
    
    TH1D * hphopurnum[ana::nJetR];
    TH1D * hphopurden[ana::nJetR];
    TEfficiency * hphopur[ana::nJetR];
    TH1D * hjetpurnum[ana::nJetR];
    TH1D * hjetpurden[ana::nJetR];
    TEfficiency * hjetpur[ana::nJetR];
  
    TH1D * hpaireffnum[ana::nJetR];
    TH1D * hpaireffden[ana::nJetR];
    TEfficiency * hpaireff[ana::nJetR];
    TH1D * hpairpurnum[ana::nJetR];
    TH1D * hpairpurden[ana::nJetR];
    TEfficiency * hpairpur[ana::nJetR];

    TH2D * hphomissfake[ana::nJetR];
    TH2D * hjetmissfake[ana::nJetR];
    TH2D * hpairmissfake[ana::nJetR];

    TH1D * hrecojetpt[ana::nJetR];
    TH1D * htruthjetpt[ana::nJetR];
    TH1D * hunfoldjetpt[ana::nJetR];
    TH2D * hjetresponse[ana::nJetR];
    TH1D * hrecophopt[ana::nJetR];
    TH1D * htruthphopt[ana::nJetR];
    TH1D * hunfoldphopt[ana::nJetR];
    TH2D * hphoresponse[ana::nJetR];
    
    TH1D * hrecojetpt_half[ana::nJetR];
    TH1D * htruthjetpt_half[ana::nJetR];
    TH1D * hunfoldjetpt_half[ana::nJetR];
    TH2D * hjetresponse_half[ana::nJetR];
    TH1D * hrecophopt_half[ana::nJetR];
    TH1D * htruthphopt_half[ana::nJetR];
    TH1D * hunfoldphopt_half[ana::nJetR];
    TH2D * hphoresponse_half[ana::nJetR];
    
   
    int nbins = ana::nPtBins * ana::nUnfoldBins; 
    RooUnfoldResponse * jet_response2D[ana::nJetR];
    RooUnfoldResponse * jet_response_half2D[ana::nJetR];
    
    TH1D * hrecoxj         [ana::nJetR];
    TH1D * htruthxj        [ana::nJetR];
    TH1D * hunfoldxj       [ana::nJetR];
    TH2D * hxjresponse     [ana::nJetR];
    
    TH1D * hrecoxj_half    [ana::nJetR];
    TH1D * htruthxj_half   [ana::nJetR];
    TH1D * hunfoldxj_half  [ana::nJetR];
    TH2D * hxjresponse_half[ana::nJetR];

};

#endif // UNFOLDER_H
