#include "unfolder.h"
using namespace std;

unfolder::~unfolder() {}
pho_object unfolder::getmaxpho(vector<pho_object> phos) {
  pho_object max;
  for (int i = 0; i < phos.size(); i++) {
    pho_object pho = phos.at(i);
    if (pho.pt > max.pt) {
      max = pho;
    }
  }
  return max;
}
jet_object unfolder::getmaxjet(vector<jet_object> jets, pho_object pho, int ir) {
  jet_object max;
  for (int i = 0; i < jets.size(); i++) {
    jet_object jet = jets.at(i);
    float dr = pho.deltaR(jet);
    float dphi = jet.deltaPhi(pho);
    if (dr < ana::drcut[ir]) continue;
    if (jet.pt > max.pt) {
      max = jet;
    }
  }
  return max;
}

bool unfolder::check_pair(jet_object jet, vector<jet_object> jets, int ir, pho_object pho) {
  float dphi = jet.deltaPhi(pho);
  int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
  
  if (iabcd != 0) return false;
  if (abs(pho.eta) > ana::etacut) return false;
  if (abs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;
  if (pho.pt < 0) return false;
  
  return true;
}

bool unfolder::check_match(pho_object pr, pho_object pt, jet_object jr, jet_object jt) {
  if (pr.deltaR(pt) < 0.1 && jr.deltaR(jt) < 0.1) return true;
  else return false;
}
void unfolder::fill_matrix() {

  nentries = t->GetEntriesFast();
  cout << "running..." << endl;

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (e % 1000 == 0)
      std::cout << "entry " << e << "/" << nentries
        << " (" << (float)e/nentries*100. << "%)\t\r" << std::flush;

    // -----------------------
    // Build clusters & jets
    // -----------------------
    auto clusters = make_clusters(*cluster_pt, *cluster_e, *cluster_eta, *cluster_phi,
        *cluster_showershape, *cluster_time, *cluster_bdt_scores);
    auto truth_clusters = make_clusters(*truth_cluster_pt, *truth_cluster_e,
        *truth_cluster_eta, *truth_cluster_phi);

    vector<vector<jet_object>> jets(ana::nJetR);
    vector<vector<jet_object>> truth_jets(ana::nJetR);

    for (int i = 0; i < ana::nJetR; i++) {
      jets[i] = make_jets(jet_pt_smear->at(i), jet_e->at(i), jet_eta->at(i), jet_phi->at(i),
          jet_emfrac->at(i), jet_ihfrac->at(i), jet_ohfrac->at(i), jet_time->at(i));

      truth_jets[i] = make_jets(truth_jet_pt->at(i), truth_jet_e->at(i),
          truth_jet_eta->at(i), truth_jet_phi->at(i));
    }

    // -----------------------
    // Event selection
    // -----------------------
    bool isphoton = (trigger == "Photon5" || trigger == "Photon10" || trigger == "Photon20");

    float truth_max_pho_pt = findmaxpt(truth_clusters);
    bool isc = (isphoton ? (truth_max_pho_pt > threshmap[-1][trigger] &&
          truth_max_pho_pt < threshmap_high[-1][trigger])
        : true);
    count_isc += isc;

    vector<bool> isj(ana::nJetR, true);
    for (int i = 0; i < ana::nJetR; i++) {
      float truth_max_jet_pt = findmaxpt(truth_jets[i]);
      isj[i] = (isphoton ? true : (truth_max_jet_pt > threshmap[i][trigger] &&
            truth_max_jet_pt < threshmap_high[i][trigger]));
      count_isj[i] += isj[i];
    }

    // -----------------------
    // Leading photon & isolation
    // -----------------------
    pho_object maxpho = getmaxpho(clusters);
    pho_object maxpho_truth = getmaxpho(truth_clusters);

    float maxpho_iso4 = 0;
    // add cluster contribution
    for (const auto& c : truth_clusters) {
      float dr = maxpho_truth.deltaR(c);
      if (dr > 0.01 && dr < 0.4) maxpho_iso4 += c.pt;
    }
    // add jet contribution (only for R=1 jets here)
    for (const auto& j : truth_jets[1]) {
      float dr = maxpho_truth.deltaR(j);
      if (dr > 0.01 && dr < 0.4) maxpho_iso4 += j.pt;
    }
    maxpho_truth.iso4 = maxpho_iso4;
    maxpho_truth.bdt = 0.99;

    // -----------------------
    // Leading jets
    // -----------------------
    vector<jet_object> maxjet(ana::nJetR);
    vector<jet_object> maxjet_truth(ana::nJetR);
    vector<bool> ispaired(ana::nJetR, false);
    vector<bool> ispaired_truth(ana::nJetR, false);
    bool anypaired = false, anypaired_truth = false, isanymatch = false;

    for (int ir = 0; ir < ana::nJetR; ir++) {
      maxjet[ir] = getmaxjet(jets[ir], maxpho, ir);
      maxjet_truth[ir] = getmaxjet(truth_jets[ir], maxpho_truth, ir);

      // -----------------------
      // Pairing
      // -----------------------
      if (maxpho.pt > ana::cluster_pt_cut && maxjet[ir].pt > ana::jet_pt_cut[ir] && isc && isj[ir]) {
        ispaired[ir] = check_pair(maxjet[ir], jets[ir], ir, maxpho);
        anypaired |= ispaired[ir];
      }
      if (maxpho_truth.pt > ana::cluster_pt_cut && maxjet_truth[ir].pt > ana::jet_pt_cut[ir] && isc && isj[ir]) {
        ispaired_truth[ir] = check_pair(maxjet_truth[ir], truth_jets[ir], ir, maxpho_truth);
        anypaired_truth |= ispaired_truth[ir];
      }

      // -----------------------
      // Fill response matrix
      // -----------------------
      if (ispaired[ir] && !ispaired_truth[ir])
        jet_response[ir].Fake(maxjet[ir].pt);
      else if (!ispaired[ir] && ispaired_truth[ir])
        jet_response[ir].Miss(maxjet_truth[ir].pt);
      else if (ispaired[ir] && ispaired_truth[ir]) {
        bool ismatch = check_match(maxpho, maxpho_truth, maxjet[ir], maxjet_truth[ir]);
        isanymatch |= ismatch;
        if (ismatch)
          jet_response[ir].Fill(maxjet[ir].pt, maxjet_truth[ir].pt);
        else {
          jet_response[ir].Miss(maxjet_truth[ir].pt);
          jet_response[ir].Fake(maxjet[ir].pt);
        }
      }

      if (ispaired[ir]) hjetpt[ir]->Fill(maxjet[ir].pt);
      if (ispaired_truth[ir]) htruthjetpt[ir]->Fill(maxjet_truth[ir].pt);
    }

    // -----------------------
    // Photon histograms & response
    // -----------------------
    if (anypaired) hclusterpt->Fill(maxpho.pt);
    if (anypaired_truth) htruthclusterpt->Fill(maxpho_truth.pt);

    if (anypaired && !isanymatch)
      pho_response.Fake(maxpho.pt);
    else if (!anypaired && anypaired_truth)
      pho_response.Miss(maxpho_truth.pt);
    else if (anypaired && anypaired_truth) {
      if (isanymatch)
        pho_response.Fill(maxpho.pt, maxpho_truth.pt);
      else {
        pho_response.Miss(maxpho_truth.pt);
        pho_response.Fake(maxpho.pt);
      }
    }
  }

  // -----------------------
  // Unfold
  // -----------------------
  RooUnfoldBayes pho_unfold(&pho_response, hclusterpt, 4);
  hclusterreco = (TH1D*)pho_unfold.Hreco();
  hclusterreco->SetName("hclusterreco");

  for (int ir = 0; ir < ana::nJetR; ir++) {
    RooUnfoldBayes jet_unfold(&jet_response[ir], hjetpt[ir], 4);
    hjetreco[ir] = (TH1D*)jet_unfold.Hreco();
    hjetreco[ir]->SetName(Form("hjetreco%i", ir));
  }
}

void unfolder::end() {
  const char * wfilename = Form("/home/samson72/sphnx/gammajet/hists/%s_unfolding.root",trigger.c_str());
  cout << "Writing files to " << wfilename << endl;
  TFile::Open(wfilename, "RECREATE");
  
  hclusterpt->Write();
  htruthclusterpt->Write();
  hclusterreco->Write();
  for (int ir = 0; ir < ana::nJetR; ir++) {
    hjetpt[ir]->Write();
    htruthjetpt[ir]->Write();
    hjetreco[ir]->Write();
  }
}

void unfolder::treesetup() {
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

