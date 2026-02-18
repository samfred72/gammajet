#include "/home/samson72/sphnx/gammajet/src/histmaker.h"
using namespace std;

// Doesn't actually loop. Just checks that the found photon and found jet are a correct match for each other.
// This is where the xJ plot is filled
bool histmaker::loop(jet_object jet, vector<jet_object> jets, int ir, pho_object pho, int icalib, bool isthirdjet) {
  bool useshowershape = 0;
  float dphi = jet.deltaPhi(pho);
  float jetval = jet.pt;
  float phoval = pho.pt;
  float val = jetval/phoval;
  int ipt = ana::findPtBin(pho.pt);
  if (icalib == 0) {
    hdeltaphi[ipt][ir]->Fill(dphi);
  }
  
  if (abs(pho.eta) > ana::etacut) return false;
  if (abs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (!isMC && abs(pho.t - jet.t) > ana::tcut) return false;
  if (dphi < ana::oppcut) return false;
  if (ipt < 0) return false;
  if (useshowershape && pho.showershape == 0) return false;
  
  jet_object thirdjet = getthirdjet(pho, jet, jets, ir, icalib);
  bool hasthirdjet = thirdjet.pt > 0;

  // Get the ABCD info and fill
  int iabcd[ana::nIsoBdtBins];
  for (int iib = 0; iib < ana::nIsoBdtBins; iib++) {
    if (useshowershape) {
      iabcd[iib] = ana::findabcdBin(pho.iso4, pho.showershape, iib);
    }
    else {
      iabcd[iib] = ana::findabcdBin(pho.iso4, pho.bdt, iib);
    }
    if (iabcd[iib] == -1) continue;

    hratio[ipt][ir][icalib][iib][0][iabcd[iib]]->Fill(val); 
    if (!hasthirdjet) hratio[ipt][ir][icalib][iib][1][iabcd[iib]]->Fill(val);
  }
 
  if (icalib == 1) { 
    hjetpt[ir]->Fill(jet.pt); 
    hjeteta[ir]->Fill(jet.eta); 
    if (jet.emfrac > 0.8) hjetetahighem[ir]->Fill(jet.eta); 
    if (jet.emfrac < 0.5) hjetetalowem[ir]->Fill(jet.eta); 
    hjetetaphi[ir]->Fill(jet.eta,jet.phi); 
    hemfrac[ipt][ir]->Fill(jet.emfrac);
    int ixj = ana::findxjBin(val);
    if (ixj >= 0) hjetetaxj[ixj][ir]->Fill(jet.eta);
    if (iabcd[0] == 0 && hasthirdjet) {
      h3jetpt[ir]->Fill(pho.pt,thirdjet.pt);
      h3jetdeltar[ipt][ir]->Fill(thirdjet.deltaR(jet));
    }
    if (pho.iso4 < ana::isoBins[0] && pho.pt > 16 && pho.pt < 25) {
      int ibdt = ana::findBdtBin(pho.bdt);
      hxjbdt[ibdt][ir]->Fill(val);
    }
  }
  return true;
}

// The called function
void histmaker::make_hists()
{
  nentries = t->GetEntriesFast();
  cout << "running..." << endl;
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    t0corr = t0map[RunNumber];
    if(e % 1000==0) std::cout << "entry " << e << "/" << nentries << " (" << (float)e/nentries*100. << "%)" << "\t\r" << std::flush;
    if (fabs(vz) > ana::vzcut) continue;
    if (!isMC && !ScaledTriggerBit[27] && !ScaledTriggerBit[38]) continue;
    mbd_t0 = (mbd_time_south + mbd_time_north)/2.0 - t0corr;
    
    vector<pho_object> truth_clusters;
    vector<pho_object> clusters;
    vector<vector<jet_object>> truth_jets;
    vector<vector<jet_object>> jets;
    vector<vector<jet_object>> jets_calib;
    vector<vector<jet_object>> jets_smear;
    
    clusters = make_clusters(*cluster_pt, *cluster_e, *cluster_eta, *cluster_phi, *cluster_showershape, *cluster_time, *cluster_bdt_scores);
    for (int i = 0; i < ana::nJetR; i++) {
      jets.push_back(make_jets(jet_pt->at(i), jet_e->at(i), jet_eta->at(i), jet_phi->at(i), jet_emfrac->at(i), jet_ihfrac->at(i), jet_ohfrac->at(i), jet_time->at(i)));
      jets_calib.push_back(make_jets(jet_pt_calib->at(i), jet_e->at(i), jet_eta->at(i), jet_phi->at(i), jet_emfrac->at(i), jet_ihfrac->at(i), jet_ohfrac->at(i), jet_time->at(i)));
      if (isMC) jets_smear.push_back(make_jets(jet_pt_smear->at(i), jet_e->at(i), jet_eta->at(i), jet_phi->at(i), jet_emfrac->at(i), jet_ihfrac->at(i), jet_ohfrac->at(i), jet_time->at(i)));
    }
    
    // Fill time histos before any cuts
    if (!isMC) {
      hmbdt->Fill(mbd_t0);
      for (int i = 0; i < jets.at(1).size(); i++) {
        hjett->Fill(jets.at(1).at(i).t);
        hmjt->Fill(mbd_t0,jets.at(1).at(i).t);
      }
      for (int i = 0; i < ana::nJetR; i++) {
        for (int j = 0; j < jet_time->at(i).size(); j++ ){
          hmtminusjt[i]->Fill(mbd_t0-jets.at(i).at(j).t);
        }
      }

      for (int i = 0; i < clusters.size(); i++) {
        pho_object cluster = clusters.at(i);
        hclustert->Fill(cluster.t);
        hmct->Fill(mbd_t0,cluster.t);
        hmtminusct->Fill(mbd_t0-cluster.t);
        for (int j = 0; j < ana::nJetR; j++) {
          for (int k = 0; k < jets.at(j).size(); k++ ){
            jet_object jet = jets.at(j).at(k);
            if (jet.deltaR(cluster) > ana::JetRs[j]) hctminusjt[j]->Fill(cluster.t-jet.t);
          }
        }
      }
    }

    
    // Check if the MC event should be kept
    bool isc = 1;
    vector<bool> isj = {1,1,1,1};
    bool isphoton = (trigger == "Photon5" || trigger == "Photon10" || trigger == "Photon20");

    if (isMC) {
      // Check the clusters
      truth_clusters = make_clusters(*truth_cluster_pt, *truth_cluster_e, *truth_cluster_eta, *truth_cluster_phi);

      float truthcpt = findmaxpt(truth_clusters);
      isc = (isphoton ? (truthcpt > threshmap[-1][trigger] && truthcpt < threshmap_high[-1][trigger]) : 1);
      count_isc += isc;
      if (truthcpt > 0) htruthclusterptprecut->Fill(truthcpt);
      if (isc) htruthclusterpt->Fill(truthcpt);
      
      // Check the jets
      for (int i = 0; i < ana::nJetR; i++) {
        truth_jets.push_back(make_jets(truth_jet_pt->at(i), truth_jet_e->at(i), truth_jet_eta->at(i), truth_jet_phi->at(i)));
        float truthjpt = findmaxpt(truth_jets[i]);
        isj[i] = (isphoton ? 1 : (truthjpt > threshmap[i][trigger] && truthjpt < threshmap_high[i][trigger]));
        if (truthjpt > 0) htruthjetptprecut[i]->Fill(truthjpt);
        if (isj[i] && isc) {
          htruthjetpt[i]->Fill(truthjpt);
        }
        count_isj[i] += isj[i];
      }
    }

    // Store all the jets into nice arrays
    for (int i = 0; i < nClusters; i++) {
      for (int j = 0; j < 11; j++) {
        hbdt[j]->Fill(cluster_bdt_scores->at(j).at(i));
      }
    }
    for (int ir = 0; ir < ana::nJetR; ir++) {
      for (int ij = 0; ij < jets.at(ir).size(); ij++) {
        hJES[ir]->Fill(jets.at(ir).at(ij).pt,jets.at(ir).at(ij).pt/jets_calib.at(ir).at(ij).pt);
        if (isMC) hjetsmear[ir]->Fill(jets.at(ir).at(ij).pt,jets_smear.at(ir).at(ij).pt);
      }
    }
    
    // Now actually find the jet and photon objects
    pho_object maxpho = getmaxpho(clusters); 

    vector<jet_object> maxjet(ana::nJetR);
    vector<jet_object> maxjet_calib(ana::nJetR);
    vector<jet_object> maxjet_smear(ana::nJetR);
    int ipt = ana::findPtBin(maxpho.pt);
    if (ipt < 0) continue; 

    bool ispaired[ana::nJetR] = { 0 };
    bool anypaired = false;
    for (int ir = 0; ir < ana::nJetR; ir++) {
      // one for uncalib, calibrated, and JER smeared`
      maxjet[ir] = getmaxjet(jets[ir], maxpho,ir);
      maxjet_calib[ir] = getmaxjet(jets_calib[ir], maxpho,ir);
      if (isMC) maxjet_smear[ir] = getmaxjet(jets_smear[ir], maxpho,ir);
      // Fill the xJ histograms
      if (maxjet[ir].pt > ana::jet_pt_cut[ir] && isc && isj[ir]) {
        ispaired[ir] =    loop(maxjet[ir], jets[ir], ir, maxpho, 0);
      }
      if (maxjet_calib[ir].pt > ana::jet_calib_pt_cut[ir] && isc && isj[ir]) {
        loop(maxjet_calib[ir], jets_calib[ir], ir, maxpho, 1);
      }
      if (isMC && maxjet_smear[ir].pt > ana::jet_calib_pt_cut[ir] && isc && isj[ir]) {
        loop(maxjet_smear[ir], jets_smear[ir], ir, maxpho, 2);
      }
      
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
    
    bool isiso = maxpho.iso4 < ana::isoBins[0];
    bool isbdt = maxpho.bdt > ana::bdtBins[0];
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
    for (int i = 0; i < ana::nJetR; i++) {
      if (maxjet[i].pt > 0) {
        hjetptprecut[i]->Fill(maxjet[i].pt);
      }
    }
    
  }
  // the end
  return;
}


// Finds the max jet on the opposite side of the detector given the max cluster 
jet_object histmaker::getmaxjet(vector<jet_object> jets, pho_object pho,int ij) {
  jet_object max;
  for (int i = 0; i < jets.size(); i++) {
    jet_object jet = jets.at(i);
    float dr = pho.deltaR(jet);
    hdeltar[ij]->Fill(dr);
    float dphi = jet.deltaPhi(pho);
    int ipt = ana::findPtBin(pho.pt);
    if (ipt >= 0) hdeltaphiprecut[ipt][ij]->Fill(dphi);
    if (!isMC && (mbd_t0 - jet.t > ana::thighcut || mbd_t0 - jet.t < ana::tlowcut)) continue;
    if (dr < ana::drcut[ij]) continue;
    if (jet.pt > max.pt) {
      max = jet;
    }
  }
  return max;
}
// Finds the max cluster
pho_object histmaker::getmaxpho(vector<pho_object> phos) {
  pho_object max;
  for (int i = 0; i < phos.size(); i++) {
    pho_object pho = phos.at(i);
    if (!isMC && (mbd_t0 - pho.t > ana::thighcut || mbd_t0 - pho.t < ana::tlowcut)) continue;
    if (pho.pt > max.pt) {
      max = pho;
    }
  }
  return max;
}
// Checks if there's a third jet in the event
jet_object histmaker::getthirdjet(pho_object maxpho, jet_object maxjet, vector<jet_object> jets, int ir, int icalib) {
  jet_object max3jet;
  //if (jets.size() > 2) {
  //cout << "Jets of radius " << ana::JetRs[ir] << ": " << jets.size() << endl;
  //cout << " pho eta: " << maxpho.eta << "  pho phi: " << maxpho.phi << endl;
  //cout << " jet eta: " << maxjet.eta << "  jet phi: " << maxjet.phi << endl;
  //}
  float ptcut = (icalib ? ana::jet_calib_pt_cut[ir] : ana::jet_pt_cut[ir]);
  for (int i = 0; i < jets.size(); i++) {
    jet_object jet = jets.at(i);
    if (jet.pt < ptcut) continue;
    //if (jets.size() > 2) cout << "3jet eta: " << jet.eta << " 3jet phi: " << jet.phi << " 3jet dRp: " << maxpho.deltaR(jet) << " 3jet dRj: " << maxjet.deltaR(jet) << " 3jet pt: " << jet.pt << endl;
    if (maxpho.deltaR(jet) > ana::drcut[ir] && maxjet.deltaR(jet) > ana::drcut[ir]) {
      //cout << "3jet eta: " << jet.eta << " 3jet phi: " << jet.phi << " 3jet dRp: " << maxpho.deltaR(jet) << " 3jet dRj: " << maxjet.deltaR(jet) << " 3jet pt: " << jet.pt << endl;
      if (jet.pt > max3jet.pt) max3jet = jet;
    }
  }
  return max3jet;
}

// Gets the fragmentation function of the cluster
// Only considers jets within R=0.4 of the cluster. Has option for isolation energy or not
float histmaker::getZ(pho_object pho, vector<jet_object> jets, bool isiso) {
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

void histmaker::savehists(TH1D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void histmaker::savehists(TH2D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void histmaker::savehists(TH1D * h[][ana::nJetR], int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      h[i][j]->Write();
    }
  }
}
void histmaker::savehists(TH2D * h[][ana::nJetR], int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      h[i][j]->Write();
    }
  }
}
void histmaker::end() {
  const char * wfilename = Form("/home/samson72/sphnx/gammajet/hists/hists%s.root",trigger.c_str());
  TFile *wf = TFile::Open(wfilename,"recreate");
  
  for (int i = 0; i < ana::nPtBins; i++) {
    for (int j = 0; j < ana::nJetR; j++) {
      for (int k = 0; k < ana::nCalibBins; k++) {
        for (int l = 0; l < ana::nIsoBdtBins; l++) {
          for (int m = 0; m < ana::n3jetBins; m++) {
            for (int n = 0; n < 4; n++) {
              hratio[i][j][k][l][m][n]->Write();
            }
          }
        }
      }
    }
  }
  savehists(hisobdt,ana::nPtBins);
  savehists(hclusterptabcd,4); // for ABCD
  
  savehists(hemfrac,ana::nPtBins,ana::nJetR);
  savehists(hdeltaphi,ana::nPtBins,ana::nJetR);
  savehists(hdeltaphiprecut,ana::nPtBins,ana::nJetR);
  savehists(h3jetdeltar,ana::nPtBins, ana::nJetR);
  savehists(hjetetaxj,ana::nxjBins,ana::nJetR);
  savehists(hxjbdt,ana::nBdtBins,ana::nJetR);
  
  savehists(hjetpt,ana::nJetR);
  savehists(hjetptprecut,ana::nJetR);
  savehists(htruthjetpt,ana::nJetR);
  savehists(htruthjetptspec,ana::nJetR);
  savehists(htruthjetptanti,ana::nJetR);
  savehists(htruthjetptprecut,ana::nJetR);
  savehists(hiso,ana::nJetR);
  savehists(hiso2d,ana::nJetR);
  savehists(hmtminusjt,ana::nJetR);
  savehists(hctminusjt,ana::nJetR);
  savehists(hjeteta,ana::nJetR);
  savehists(hjetetahighem,ana::nJetR);
  savehists(hjetetalowem,ana::nJetR);
  savehists(hjetetaphi,ana::nJetR);
  savehists(hdeltar,ana::nJetR);
  savehists(hJES,ana::nJetR);
  savehists(hjetsmear,ana::nJetR);
  savehists(h3jetpt,ana::nJetR);
  
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
    for (int i = 0; i < ana::nJetR; i++) {
      cout << Form("Events with jet R=0.%i:   ",2*(i+1)) << count_isj[i] << "/" << nentries << ": " << (int)((float)count_isj[i]/(float)nentries*100) << "%" << endl;
    }
  }
}

void histmaker::treesetup() {
  cluster_pt = 0;
  cluster_e = 0;
  cluster_eta = 0;
  cluster_phi = 0;
  cluster_time = 0;
  cluster_bdt_scores = 0;
  cluster_showershape = 0;

  truth_cluster_pt = 0;
  truth_cluster_e = 0;
  truth_cluster_eta = 0;
  truth_cluster_phi = 0;

  jet_pt = 0;
  jet_pt_calib = 0;
  jet_pt_smear = 0;
  jet_e = 0;
  jet_eta = 0;
  jet_phi = 0;
  jet_emfrac = 0;
  jet_ihfrac = 0;
  jet_ohfrac = 0;
  jet_time = 0;  

  truth_jet_pt = 0;
  truth_jet_e = 0;
  truth_jet_eta = 0;
  truth_jet_phi = 0;

  // Set branch addresses and branch pointers
  if (!t) return;

  t->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  t->SetBranchAddress("vz", &vz, &b_vz);
  t->SetBranchAddress("ScaledTriggerBit", ScaledTriggerBit, &b_ScaledTriggerBit);
  t->SetBranchAddress("LiveTriggerBit", LiveTriggerBit, &b_LiveTriggerBit);
  t->SetBranchAddress("Scaledowns", Scaledowns, &b_Scaledowns);
  t->SetBranchAddress("mbd_nhits_south", &mbd_nhits_south, &b_mbd_nhits_south);
  t->SetBranchAddress("mbd_nhits_north", &mbd_nhits_north, &b_mbd_nhits_north);
  t->SetBranchAddress("mbd_time_south", &mbd_time_south, &b_mbd_time_south);
  t->SetBranchAddress("mbd_time_north", &mbd_time_north, &b_mbd_time_north);

  t->SetBranchAddress("nClusters", &nClusters, &b_nClusters);
  t->SetBranchAddress("cluster_pt" , &cluster_pt , &b_cluster_pt );
  t->SetBranchAddress("cluster_e"  , &cluster_e  , &b_cluster_e  );
  t->SetBranchAddress("cluster_eta", &cluster_eta, &b_cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi, &b_cluster_phi);
  t->SetBranchAddress("cluster_showershape", &cluster_showershape, &b_cluster_showershape);
  t->SetBranchAddress("cluster_time", &cluster_time, &b_cluster_time);
  t->SetBranchAddress("cluster_bdt_scores", &cluster_bdt_scores, &b_cluster_bdt_scores);

  t->SetBranchAddress("nJets", nJets, &b_nJets);
  t->SetBranchAddress("jet_pt" , &jet_pt , &b_jet_pt );
  t->SetBranchAddress("jet_pt_calib" , &jet_pt_calib , &b_jet_pt_calib );
  t->SetBranchAddress("jet_e"  , &jet_e  , &b_jet_e  );
  t->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
  t->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
  t->SetBranchAddress("jet_emfrac", &jet_emfrac, &b_jet_emfrac);
  t->SetBranchAddress("jet_ihfrac", &jet_ihfrac, &b_jet_ihfrac);
  t->SetBranchAddress("jet_ohfrac", &jet_ohfrac, &b_jet_ohfrac);
  t->SetBranchAddress("jet_time", &jet_time, &b_jet_time);

  if (isMC) {   
    t->SetBranchAddress("jet_pt_smear" , &jet_pt_smear , &b_jet_pt_smear );
    t->SetBranchAddress("nTruthClusters", &nTruthClusters, &b_nTruthClusters);
    t->SetBranchAddress("truth_cluster_pt" , &truth_cluster_pt , &b_truth_cluster_pt );
    t->SetBranchAddress("truth_cluster_e"  , &truth_cluster_e  , &b_truth_cluster_e  );
    t->SetBranchAddress("truth_cluster_eta", &truth_cluster_eta, &b_truth_cluster_eta);
    t->SetBranchAddress("truth_cluster_phi", &truth_cluster_phi, &b_truth_cluster_phi);
    
    t->SetBranchAddress("nTruthJets", nTruthJets, &b_nTruthJets);
    t->SetBranchAddress("truth_jet_pt" , &truth_jet_pt , &b_truth_jet_pt );
    t->SetBranchAddress("truth_jet_e"  , &truth_jet_e  , &b_truth_jet_e  );
    t->SetBranchAddress("truth_jet_eta", &truth_jet_eta, &b_truth_jet_eta);
    t->SetBranchAddress("truth_jet_phi", &truth_jet_phi, &b_truth_jet_phi);
  }
}
