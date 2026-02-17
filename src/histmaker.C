#include <iostream>
#include <vector>
#include <TSystem.h>
#include "/home/samson72/sphnx/headers/ana.h"
#include "/home/samson72/sphnx/headers/TreeSetting.h"
#include "/home/samson72/sphnx/headers/drawer.h"
#include "/home/samson72/sphnx/headers/object.h"

// global variables
ana anaclone;
float mbd_t0;
bool isMC;

// Define histograms
TH1D * hratio     [anaclone.nPtBins][anaclone.nJetR][anaclone.nCalibBins][anaclone.nIsoBdtBins][4]; // 4 accounts for a,b,c,d regions, 2 is for noncalib vs calib jets
TH1D * hratio3jet [anaclone.nPtBins][anaclone.nJetR][anaclone.nCalibBins][anaclone.nIsoBdtBins][4]; // 4 accounts for a,b,c,d regions, 2 is for noncalib vs calib jets

TH2D * hisobdt        [anaclone.nPtBins];

TH1D * hjetetaxj      [anaclone.nxjBins][anaclone.nJetR];

TH1D * hbdt[11]; // 11 is the number of models
TH1D * hemfrac        [anaclone.nPtBins][anaclone.nJetR];
TH1D * hdeltaphi      [anaclone.nPtBins][anaclone.nJetR];
TH1D * hdeltaphiprecut[anaclone.nPtBins][anaclone.nJetR];

TH1D * hmtminusjt           [anaclone.nJetR];
TH1D * hctminusjt           [anaclone.nJetR];
TH1D * hiso                 [anaclone.nJetR];
TH1D * hjetpt               [anaclone.nJetR];
TH1D * htruthjetpt          [anaclone.nJetR];
TH1D * htruthjetptspec      [anaclone.nJetR];
TH1D * htruthjetptanti      [anaclone.nJetR];
TH1D * hjetptprecut         [anaclone.nJetR];
TH1D * htruthjetptprecut    [anaclone.nJetR];
TH1D * htruthjetptprecutspec[anaclone.nJetR];
TH1D * htruthjetptprecutanti[anaclone.nJetR];
TH1D * hjeteta              [anaclone.nJetR];
TH1D * hjetetahighem        [anaclone.nJetR];
TH1D * hjetetalowem         [anaclone.nJetR];
TH1D * hdeltar              [anaclone.nJetR];
TH1D * hclusterptabcd       [4]; // for ABCD

TH2D * hjetetaphi           [anaclone.nJetR];
TH2D * hiso2d               [anaclone.nJetR];
TH2D * hJES                 [anaclone.nJetR];

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

void savehists(TH1D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void savehists(TH2D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void savehists(TH1D * h[][anaclone.nJetR], int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      h[i][j]->Write();
    }
  }
}

// Initialize histograms as needed
void inith() {
  for (int i = 0; i < anaclone.nPtBins; i++) {
    for (int j = 0; j < anaclone.nJetR; j++) {
      for (int k = 0; k < anaclone.nCalibBins; k++) {
        for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
          for (int m = 0; m < 4; m++) {
            hratio[i][j][k][l][m] = new TH1D(Form("hratio_%i_%i_%i_%i_%i",i,j,k,l,m),";p_{T}^{jet}/p_{T}^{#gamma};normalized counts",100,0,2);
            hratio3jet[i][j][k][l][m] = new TH1D(Form("hratio3jet_%i_%i_%i_%i_%i",i,j,k,l,m),";p_{T}^{jet}/p_{T}^{#gamma};normalized counts",100,0,2);
          }
        }
      }
    }
    hisobdt[i] = new TH2D(Form("hisobdt%i",i),";cluster iso;bdt score",100,-1,20,100,0,1);
  }
  for (int i = 0; i < anaclone.nJetR; i++) {
    hjetpt[i] =                new TH1D(Form("hjetpt%i",i),";jet p_{T,max};counts",100,0,100);
    htruthjetpt[i] =           new TH1D(Form("htruthjetpt%i",i),";jet p_{T,max};counts",100,0,100);
    htruthjetptspec[i] =       new TH1D(Form("htruthjetptspec%i",i),";jet p_{T,max};counts",100,0,100);
    htruthjetptanti[i] =       new TH1D(Form("htruthjetptanti%i",i),";jet p_{T,max};counts",100,0,100);
    hjetptprecut[i] =          new TH1D(Form("hjetptprecut%i",i),";jet p_{T,max};counts",100,0,100);
    htruthjetptprecut[i] =     new TH1D(Form("htruthjetptprecut%i",i),";jet p_{T,max};counts",100,0,100);
    htruthjetptprecutspec[i] = new TH1D(Form("htruthjetptprecutspec%i",i),";jet p_{T,max};counts",100,0,100);
    htruthjetptprecutanti[i] = new TH1D(Form("htruthjetptprecutanti%i",i),";jet p_{T,max};counts",100,0,100);
    hmtminusjt[i] =            new TH1D(Form("hmtminusjt%i",i),";t_{mbd}-t_{jet} [ns]",100,-10,10);
    hctminusjt[i] =            new TH1D(Form("hctminusjt%i",i),";t_{cluster}-t_{jet} [ns]",100,-10,10);
    hiso[i] =                  new TH1D(Form("hiso%i",i),";iso E R02;counts",100,-1,50);
    hjeteta[i] =               new TH1D(Form("hjeteta%i",i),";#eta_jet;counts",100,-1.1,1.1);
    hjetetahighem[i] =         new TH1D(Form("hjetetahighem%i",i),";#eta_jet (emfrac > 0.8);counts",100,-1.1,1.1);
    hjetetalowem[i] =          new TH1D(Form("hjetetalowem%i",i),";#eta_jet (emfrac < 0.5);counts",100,-1.1,1.1);
    hdeltar[i] =               new TH1D(Form("hdeltar%i",i),";dr [eta,phi];counts",100,0,4);

    hjetetaphi[i] =            new TH2D(Form("hjetetaphi%i",i),";#eta_{jet};#phi_{jet}",100,-1.1,1.1,100,-M_PI,M_PI);
    hiso2d[i] =                new TH2D(Form("hiso2d%i",i),";cluster E;iso E R02",100,0,50,100,0,50);
    hJES[i] =                  new TH2D(Form("hJES%i",i),";Uncalibrated Jet pT;Calib scale",50,0,50,100,0,3);
    for (int j = 0; j < anaclone.nPtBins; j++) {
      hemfrac[j][i] = new TH1D(Form("hemfrac%i_%i",j,i),";jet emfrac;counts",100,-0.1,1.1);
      hdeltaphi[j][i] = new TH1D(Form("hdeltaphi%i_%i",j,i),";|#phi_{#gamma} - #phi_{leading jet}|;counts",100,0,M_PI);
      hdeltaphiprecut[j][i] = new TH1D(Form("hdeltaphiprecut%i_%i",j,i),";|#phi_{#gamma} - #phi_{leading jet}|;counts",100,0,M_PI);
    }
    for (int j = 0; j < anaclone.nxjBins; j++) {
      hjetetaxj[j][i] = new TH1D(Form("hjetetaxj%i_%i",j,i),Form(";jet eta %1.1f < xJ < %1.1f;Counts",anaclone.xjBins[j],anaclone.xjBins[j+1]),100,-1.1,1.1);
    }
  }
  for (int i = 0; i < 11; i++) {
    hbdt[i] = new TH1D(Form("hbdt%i",i),";bdt score;counts",120,-0.1,1.1);
  }
  for (int i = 0; i < 4; i++) {
    hclusterptabcd[i] = new TH1D(Form("hclusterptabcd%i",i),";p_{T}^{lead cluster};Counts",100,0,100);
  }
}

// Gets the fragmentation function of the cluster
// Only considers jets within R=0.4 of the cluster. Has option for isolation energy or not
float getZ(pho_object pho, vector<jet_object> jets, bool isiso) {
  jet_object closest;
  float closestdr = 100;
  float Z = 0;
  for (int i = 0; i < jets.size(); i++) {
    jet_object jet = jets.at(i);
    float dr = pho.deltaR(jet);
    if (dr < 0.4 && (!isiso && dr < closestdr) || (isiso && dr < closestdr && pho.iso4 < 2)) {
      closestdr = dr;
      Z = pho.pt/jet.pt;
    }
  }
  return Z;
}


// Doesn't actually loop. Just checks that the found photon and found jet are a correct match for each other.
// This is where the xJ plot is filled
bool loop(jet_object jet,int jindex, pho_object pho, int icalib, bool isthirdjet = 0) {
  bool useshowershape = 0;
  float dphi = jet.deltaPhi(pho);
  //float jetval = jet.pt * (1 + 0.1 * jet.emfrac);
  //float phoval = pho.pt * 1.1;
  float jetval = jet.pt;
  float phoval = pho.pt;
  float val = jetval/phoval;
  int ipt = anaclone.findPtBin(pho.pt);
  
  if (abs(pho.eta) > anaclone.etacut) return false;
  if (abs(jet.eta) > anaclone.etacut - anaclone.JetRs[jindex]) return false;
  if (!isMC && abs(pho.t - jet.t) > anaclone.tcut) return false;
  if (dphi < anaclone.oppcut) return false;
  if (ipt < 0) return false;
  if (useshowershape && pho.showershape == 0) return false;
  bool istight = pho.showershape == 2;

  // Get the ABCD info
  for (int iib = 0; iib < anaclone.nIsoBdtBins; iib++) {
    if (!useshowershape && (pho.iso4 > anaclone.isoBins[iib] && pho.iso4 < anaclone.isoBinsHigh[iib])) continue;
    if (pho.bdt < anaclone.bdtCuts[iib]) continue;
    
    bool isiso = pho.iso4 < anaclone.isoBins[iib];
    bool isbdt = pho.bdt > anaclone.bdtBins[iib];
    int iabcd;
    if (useshowershape) {
      iabcd = (((isiso << 0b1) | istight) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    }
    else {
      iabcd = (((isiso << 0b1) | isbdt) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    }
    hratio[ipt][jindex][icalib][iib][iabcd]->Fill(val); 
    if (isthirdjet) hratio3jet[ipt][jindex][icalib][iib][iabcd]->Fill(val);
  }
 
  int ixj = anaclone.findxjBin(val);
  if (icalib == 0 && isthirdjet == 0) { 
    hdeltaphi[ipt][jindex]->Fill(dphi);
    hjetpt[jindex]->Fill(jet.pt); 
    hjeteta[jindex]->Fill(jet.eta); 
    if (jet.emfrac > 0.8) hjetetahighem[jindex]->Fill(jet.eta); 
    if (jet.emfrac < 0.5) hjetetalowem[jindex]->Fill(jet.eta); 
    hjetetaphi[jindex]->Fill(jet.eta,jet.phi); 
    hemfrac[ipt][jindex]->Fill(jet.emfrac);
    if (ixj >= 0) hjetetaxj[ixj][jindex]->Fill(jet.eta);
  }
  return true;
}

// Finds the truth pT for clusters and jets. 
// Returns {max, max (unmatched to a photon), max (matched to a photon)}. The last two are for the jet samples
vector<float> findmaxpt(vector<pho_object> phos) {
  float maxpt = 0;
  for (int i = 0; i < phos.size(); i++) {
    pho_object pho = phos.at(i);
    
    if (pho.pt > maxpt) {
      maxpt = pho.pt;
    }
  }
  // spec is away from a photon, anti is matched to photon
  return {maxpt,0,0};
}
vector<float> findmaxpt(vector<jet_object> jets, vector<pho_object> phos) {
  float maxpt = 0;
  float maxptspec = 0;
  float maxptanti = 0;
  for (int i = 0; i < jets.size(); i++) {
    jet_object jet = jets.at(i);
    
    bool isnotphoton = true;
    for (int j = 0; j < phos.size(); j++) {
      pho_object pho = phos.at(j);
      float dr = jet.deltaR(pho);
      if (dr < 0.2) {
        isnotphoton = false;
        break;
      }
    }
    
    if (jet.pt > maxpt) {
      maxpt = jet.pt;
    }
    if (jet.pt > maxptspec && isnotphoton) {
      maxptspec = jet.pt;
    }
    if (jet.pt > maxptanti && !isnotphoton) {
      maxptanti = jet.pt;
    }

  }
  // spec is away from a photon, anti is matched to photon
  return {maxpt,maxptspec,maxptanti};
}

// Finds the max jet on the opposite side of the detector given the max cluster 
jet_object getmaxjet(vector<jet_object> jets, pho_object pho,int ij) {
  jet_object max;
  for (int i = 0; i < jets.size(); i++) {
    jet_object jet = jets.at(i);
    float dr = pho.deltaR(jet);
    hdeltar[ij]->Fill(dr);
    float dphi = jet.deltaPhi(pho);
    int ipt = anaclone.findPtBin(pho.pt);
    if (ipt >= 0) hdeltaphiprecut[ipt][ij]->Fill(dphi);
    if (!isMC && (mbd_t0 - jet.t > anaclone.thighcut || mbd_t0 - jet.t < anaclone.tlowcut)) continue;
    if (dr < anaclone.drcut) continue;
    if (jet.pt > max.pt) {
      max = jet;
    }
  }
  return max;
}
// Finds the max cluster
pho_object getmaxpho(vector<pho_object> phos) {
  pho_object max;
  for (int i = 0; i < phos.size(); i++) {
    pho_object pho = phos.at(i);
    //if (obj.iso4 > 2 || obj.bdt < 0.8) continue;
    if (!isMC && (mbd_t0 - pho.t > anaclone.thighcut || mbd_t0 - pho.t < anaclone.tlowcut)) continue;
    if (pho.pt > max.pt) {
      max = pho;
    }
  }
  return max;
}


// Checks if there's a third jet in the event
bool getthirdjet(pho_object maxpho, jet_object maxjet, vector<jet_object> jets) {
  for (int i = 0; i < jets.size(); i++) {
    jet_object jet = jets.at(i);
    if (maxpho.deltaR(jet) > anaclone.drcut && maxjet.deltaR(jet) > anaclone.drcut) return true;
  }
  return false;
}


// The called function
void histmaker(string trigger = "Data")
{
  gInterpreter->GenerateDictionary("vector<vector<float> >");
  gSystem->Load("/home/samson72/sphnx/headers/ana.o");
  inith();
  isMC = (trigger != "Data");
  map<int,float> t0map;
  int runnum = 0;
  ifstream file = ifstream("../headers/MbdPmt.corr");
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
  float t0corr = t0map[runnum];
   

  cout << "Getting files..." << endl;
  // Filename management
  TFile *f = new TFile(Form("trees/gammajet%s.root",trigger.c_str()),"read");
  TTree *t = (TTree*) f->Get("towerntup");
  treesetup(t,isMC); // Here all the branches of the ttree are set

  const char * wfilename = Form("hists/hists%s.root",trigger.c_str());
  TFile *wf = new TFile(wfilename,"recreate");
  
  // Min and max thresholds for MC samples
  string strig = trigger;
  string trig = strig.substr(0,strig.find("_"));
  if (isMC) cout << trig << endl;
  
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

  int drawcount = 0;
  int count_isc = 0;
  vector<int> count_isj(anaclone.nJetR);
  Long64_t nentries = t->GetEntriesFast();
  cout << "running..." << endl;
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (RunNumber != runnum) {
      runnum = RunNumber;
      t0corr = t0map[runnum];
    }
    if(e % 1000==0) std::cout << "entry " << e << "/" << nentries << " (" << (float)e/nentries*100. << "%)" << "\t\r" << std::flush;
    if (fabs(vz) > anaclone.vzcut) continue;
    if (!isMC && !ScaledTriggerBit[27]) continue;
    mbd_t0 = (mbd_time_south + mbd_time_north)/2.0 - t0corr;
    
    if (!isMC) {
      // Fill time histos before any cuts
      hmbdt->Fill(mbd_t0);
      for (int i = 0; i < nJets[1]; i++) {
        hjett->Fill(jet_time->at(1).at(i)*17.6);
        hmjt->Fill(mbd_t0,jet_time->at(1).at(i)*17.6);
      }
      for (int i = 0; i < anaclone.nJetR; i++) {
        for (int j = 0; j < jet_time->at(i).size(); j++ ){
          hmtminusjt[i]->Fill(mbd_t0-jet_time->at(i).at(j)*17.6);
        }
      }

      for (int i = 0; i < nClusters; i++) {
        hclustert->Fill(cluster_time->at(i)*17.6);
        hmct->Fill(mbd_t0,cluster_time->at(i)*17.6);
        hmtminusct->Fill(mbd_t0-cluster_time->at(i)*17.6);
        for (int j = 0; j < anaclone.nJetR; j++) {
          for (int k = 0; k < jet_time->at(j).size(); k++ ){
            hctminusjt[j]->Fill(cluster_time->at(i)*17.6-jet_time->at(j).at(k)*17.6);
          }
        }
      }
    }

    vector<pho_object> truth_clusters;
    vector<pho_object> clusters;
    vector<vector<jet_object>> truth_jets;
    vector<vector<jet_object>> jets;
    vector<vector<jet_object>> jets_calib;
    //vector<vector<jet_object>> jets_smear;
    
    // Check if the MC event should be kept
    bool isc = 1;
    vector<bool> isj = {1,1,1,1};
    bool isphoton = (trig == "Photon5" || trig == "Photon10" || trig == "Photon20");

    if (isMC) {
      // Check the clusters
      truth_clusters = make_clusters(*truth_cluster_pt, *truth_cluster_e, *truth_cluster_eta, *truth_cluster_phi);

      vector<float> truthcpt = findmaxpt(truth_clusters);
      isc = (isphoton ? (truthcpt[0] > threshmap[-1][trig] && truthcpt[0] < threshmap_high[-1][trig]) : 1);
      count_isc += isc;
      if (truthcpt[0] > 0) htruthclusterptprecut->Fill(truthcpt[0]);
      if (isc) htruthclusterpt->Fill(truthcpt[0]);
      
      // Check the jets
      for (int i = 0; i < anaclone.nJetR; i++) {
        truth_jets.push_back(make_jets(truth_jet_pt->at(i), truth_jet_e->at(i), truth_jet_eta->at(i), truth_jet_phi->at(i)));
        vector<float> truthjpt = findmaxpt(truth_jets[i], truth_clusters);
        isj[i] = (isphoton ? 1 : (truthjpt[0] > threshmap[i][trig] && truthjpt[0] < threshmap_high[i][trig]));
        if (truthjpt[0] > 0) htruthjetptprecut[i]->Fill(truthjpt[0]);
        if (truthjpt[1] > 0) htruthjetptprecutspec[i]->Fill(truthjpt[1]);
        if (truthjpt[2] > 0) htruthjetptprecutanti[i]->Fill(truthjpt[2]);
        if (isj[i] && isc) {
          htruthjetpt[i]->Fill(truthjpt[0]);
          htruthjetptspec[i]->Fill(truthjpt[1]);
          htruthjetptanti[i]->Fill(truthjpt[2]);
        }
        count_isj[i] += isj[i];
      }
    }

    // Store all the jets into nice arrays
    clusters = make_clusters(*cluster_pt, *cluster_e, *cluster_eta, *cluster_phi, *cluster_showershape, *cluster_time, *cluster_bdt_scores);
    for (int i = 0; i < nClusters; i++) {
      for (int j = 0; j < 11; j++) {
        hbdt[j]->Fill(cluster_bdt_scores->at(j).at(i));
      }
    }
    for (int i = 0; i < anaclone.nJetR; i++) {
      jets.push_back(      make_jets(jet_pt->at(i),       jet_e->at(i), jet_eta->at(i), jet_phi->at(i), jet_emfrac->at(i), jet_ihfrac->at(i), jet_ohfrac->at(i), jet_time->at(i)));
      jets_calib.push_back(make_jets(jet_pt_calib->at(i), jet_e->at(i), jet_eta->at(i), jet_phi->at(i), jet_emfrac->at(i), jet_ihfrac->at(i), jet_ohfrac->at(i), jet_time->at(i)));
      //if (isMC) jets_smear.push_back(make_jets(jet_pt_smear->at(i), jet_e->at(i), jet_eta->at(i), jet_phi->at(i), jet_emfrac->at(i), jet_ihfrac->at(i), jet_ohfrac->at(i), jet_time->at(i)));
    }
    for (int ir = 0; ir < anaclone.nJetR; ir++) {
      for (int ij = 0; ij < jets.at(ir).size(); ij++) {
        hJES[ir]->Fill(jets.at(ir).at(ij).pt,jets.at(ir).at(ij).pt/jets_calib.at(ir).at(ij).pt);
      }
    }
    
    // Now actually find the jet and photon objects
    pho_object maxpho = getmaxpho(clusters); 

    vector<jet_object> maxjet(anaclone.nJetR);
    vector<jet_object> maxjet_calib(anaclone.nJetR);
    //vector<jet_object> maxjet_smear(anaclone.nJetR);
    int ipt = anaclone.findPtBin(maxpho.pt);
    if (ipt < 0) continue; 

    bool ispaired[anaclone.nJetR] = { 0 };
    bool anypaired = false;
    for (int ir = 0; ir < anaclone.nJetR; ir++) {
      // one for uncalib, calibrated, and JER smeared`
      maxjet[ir] = getmaxjet(jets[ir], maxpho,ir);
      maxjet_calib[ir] = getmaxjet(jets_calib[ir], maxpho,ir);
      //if (isMC) maxjet_smear[ir] = getmaxjet(jets_smear[ir], maxpho,ir);
      bool hasthirdjet = getthirdjet(maxpho, maxjet[ir], jets[ir]);
      // Fill the xJ histograms
      if (maxjet[ir].pt > anaclone.jet_pt_cut[ir] && isc && isj[ir]) {
        ispaired[ir] =    loop(maxjet[ir], ir, maxpho, 0);
        if (!hasthirdjet) loop(maxjet[ir], ir, maxpho, 0, 1);
      }
      if (maxjet_calib[ir].pt > anaclone.jet_calib_pt_cut[ir] && isc && isj[ir]) {
                          loop(maxjet_calib[ir], ir, maxpho, 1);
        if (!hasthirdjet) loop(maxjet_calib[ir], ir, maxpho, 1, 1);
      }
      //if (isMC && maxjet_smear[ir].pt > anaclone.jet_calib_pt_cut[ir] && isc && isj[ir]) {
      //                    loop(maxjet_smear[ir], ir, maxpho, 2);
      //  if (!hasthirdjet) loop(maxjet_smear[ir], ir, maxpho, 2, 1);
      //}
      
      anypaired |= ispaired[ir];
    }
   
    if (!anypaired) continue;

    // max cluster histos
    hisobdt[ipt]->Fill(maxpho.iso4,maxpho.bdt); 
    hiso[0]->Fill(maxpho.iso3);
    hiso[1]->Fill(maxpho.iso4);
    hiso2d[0]->Fill(maxpho.pt,maxpho.iso3);
    hiso2d[1]->Fill(maxpho.pt,maxpho.iso4);
    
    hclusterpt->Fill(maxpho.pt);
    hclustereta->Fill(maxpho.eta); 
    hclusteretaphi->Fill(maxpho.eta,maxpho.phi); 
    
    bool isiso = maxpho.iso4 < anaclone.isoBins[0];
    bool isbdt = maxpho.bdt > anaclone.bdtBins[0];
    int iabcd = (((isiso << 0b1) | isbdt) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    hclusterptabcd[iabcd]->Fill(maxpho.pt); 
    
    // cluster and jet kinematic histos
    float frag =    getZ(maxpho, jets[1], 0);
    float fragiso = getZ(maxpho, jets[1], 1);
    if (frag > 0) hfrag->Fill(frag);
    if (fragiso > 0) hfragiso->Fill(fragiso);
    
    if (maxpho.pt > 0) {
      hclusterptprecut->Fill(maxpho.pt);
    }
    for (int i = 0; i < anaclone.nJetR; i++) {
      if (maxjet[i].pt > 0) {
        hjetptprecut[i]->Fill(maxjet[i].pt);
      }
    }
    
  }
  // the end
  for (int i = 0; i < anaclone.nPtBins; i++) {
    for (int j = 0; j < anaclone.nJetR; j++) {
      for (int k = 0; k < anaclone.nCalibBins; k++) {
        for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
          for (int m = 0; m < 4; m++) {
            hratio[i][j][k][l][m]->Write();
          }
        }
      }
    }
  }
  for (int i = 0; i < anaclone.nPtBins; i++) {
    for (int j = 0; j < anaclone.nJetR; j++) {
      for (int k = 0; k < anaclone.nCalibBins; k++) {
        for (int l = 0; l < anaclone.nIsoBdtBins; l++) {
          for (int m = 0; m < 4; m++) {
            hratio3jet[i][j][k][l][m]->Write();
          }
        }
      }
    }
  }
  savehists(hisobdt,anaclone.nPtBins);
  savehists(hclusterptabcd,4); // for ABCD
  
  savehists(hemfrac,anaclone.nPtBins,anaclone.nJetR);
  savehists(hdeltaphi,anaclone.nPtBins,anaclone.nJetR);
  savehists(hdeltaphiprecut,anaclone.nPtBins,anaclone.nJetR);
  savehists(hjetetaxj,anaclone.nxjBins,anaclone.nJetR);
  
  savehists(hjetpt,anaclone.nJetR);
  savehists(hjetptprecut,anaclone.nJetR);
  savehists(htruthjetpt,anaclone.nJetR);
  savehists(htruthjetptspec,anaclone.nJetR);
  savehists(htruthjetptanti,anaclone.nJetR);
  savehists(htruthjetptprecut,anaclone.nJetR);
  savehists(htruthjetptprecutspec,anaclone.nJetR);
  savehists(htruthjetptprecutanti,anaclone.nJetR);
  savehists(hiso,anaclone.nJetR);
  savehists(hiso2d,anaclone.nJetR);
  savehists(hmtminusjt,anaclone.nJetR);
  savehists(hctminusjt,anaclone.nJetR);
  savehists(hjeteta,anaclone.nJetR);
  savehists(hjetetahighem,anaclone.nJetR);
  savehists(hjetetalowem,anaclone.nJetR);
  savehists(hjetetaphi,anaclone.nJetR);
  savehists(hdeltar,anaclone.nJetR);
  savehists(hJES,anaclone.nJetR);
  
  savehists(hbdt,11);
  
  hclusterpt->Write();
  htruthclusterpt->Write();
  hclusterptprecut->Write();
  htruthclusterptprecut->Write();
  hclustereta->Write();
  hclusteretaphi->Write();
  hfrag->Write();
  hfragiso->Write();
  hmbdt->Write();
  hclustert->Write();
  hjett->Write();
  hmtminusct->Write();
  hmct->Write();
  hmjt->Write();
  hcjt->Write();

  std::cout << std::endl << "All done with " << trigger << "!" << std::endl;
  if (isMC) {
    cout << "Events with cluster: " << count_isc  << "/" << nentries << ": " << (int)((float)count_isc /(float)nentries*100) << "%" << endl;
    for (int i = 0; i < anaclone.nJetR; i++) {
      cout << Form("Events with jet R=0.%i:   ",2*(i+1)) << count_isj[i] << "/" << nentries << ": " << (int)((float)count_isj[i]/(float)nentries*100) << "%" << endl;
    }
  }
  return;
}
