#include "unfolder.h"
using namespace std;

unfolder::~unfolder() {}

bool unfolder::check_pair(jet_object jet, int ir, pho_object pho, bool isreco) {
  float dphi = jet.deltaPhi(pho);
  int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
  
  if (iabcd != 0) return false;
  if (abs(pho.eta) > ana::etacut) return false;
  if (abs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;
  
  return true;
}

bool unfolder::check_match(pho_object p1, pho_object p2, jet_object j1, jet_object j2) {
  if (p1.deltaR(p2) < 0.1 && j1.deltaR(j2) < 0.3) return true;
  else return false;
}
bool unfolder::check_match(pho_object p1, pho_object p2) {
  if (p1.deltaR(p2) < 0.1) return true;
  else return false;
}
bool unfolder::check_match(jet_object j1, jet_object j2) {
  if (j1.deltaR(j2) < 0.3) return true;
  else return false;
}
void unfolder::fill_matrix() {

  nentries = t->GetEntriesFast();
  cout << "running..." << endl;

  TCanvas * c = new TCanvas("c","",500,1000);
  gStyle->SetOptStat(0);
  if (dodraw) c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/event_displays_%s.pdf[",trigger.c_str()));
  int ndraw = 0;

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    bool use_half = rand.Integer(2) % 2;
    if (e % 1000 == 0)
      std::cout << "entry " << e << "/" << nentries
        << " (" << (float)e/nentries*100. << "%)\t\r" << std::flush;
    
    if (fabs(vz) > ana::vzcut) continue;

    // -----------------------
    // Event selection
    // -----------------------

    vector<bool> keepMC = check_keep_MC(truth_cluster_pt, truth_jet_pt, trigger);
    bool keep = 0;
    for (int i = 0; i < ana::nJetR + 1; i++) {
      keep |= keepMC.at(i);
    }
    if (!keep) continue;

    // -----------------------
    // Leading photon & isolation
    // -----------------------
    pho_object maxpho = pho_object(
        cluster_pt, 
        cluster_e, 
        cluster_eta, 
        cluster_phi, 
        cluster_showershape[8],
        cluster_showershape[9],
        cluster_time, 
        cluster_bdt_scores[9], 
        pho_object::get_showershape(cluster_showershape, cluster_pt)
    ); 
    pho_object maxpho_truth = pho_object(
        truth_cluster_pt, 
        truth_cluster_e, 
        truth_cluster_eta, 
        truth_cluster_phi, 
        //cluster_showershape[8], // TEMPORARY!!!!!
        //cluster_showershape[9],
        truth_cluster_iso3,
        truth_cluster_iso4,
        0, // no time object for truth 
        0.99, // truth photon is guaranteed a photon
        2 // truth photon is guaranteed a photon
    ); 


    vector<jet_object> maxjet(ana::nJetR);
    vector<jet_object> maxjet_truth(ana::nJetR);
    vector<bool> ispaired(ana::nJetR, false);
    vector<bool> ispaired_truth(ana::nJetR, false);

    for (int ir = 0; ir < ana::nJetR; ir++) {
    
    
      maxjet[ir] = jet_object(
          jet_pt_smear[ir], 
          jet_e[ir], 
          jet_eta[ir], 
          jet_phi[ir], 
          jet_emfrac[ir], 
          0, 
          0, 
          jet_time[ir]
      );
      maxjet_truth[ir] = jet_object(
          truth_jet_pt[ir], 
          truth_jet_e[ir], 
          truth_jet_eta[ir], 
          truth_jet_phi[ir], 
          0, 0, 0, 0);
      
      hphodr[ir]->Fill(maxpho.deltaR(maxpho_truth));
      hjetdr[ir]->Fill(maxjet[ir].deltaR(maxjet_truth[ir]));
      
      hphopurden[ir]->Fill(maxpho.pt);
      hphoeffden[ir]->Fill(maxpho_truth.pt);
      if (check_match(maxpho, maxpho_truth)) {
        hphopurnum[ir]->Fill(maxpho.pt);
        hphoeffnum[ir]->Fill(maxpho_truth.pt);
      }
      hjetpurden[ir]->Fill(maxjet[ir].pt);
      hjeteffden[ir]->Fill(maxjet_truth[ir].pt);
      if (check_match(maxjet[ir], maxjet_truth[ir])) {
        hjetpurnum[ir]->Fill(maxjet[ir].pt);
        hjeteffnum[ir]->Fill(maxjet_truth[ir].pt);
      }
      
      float xj = maxjet[ir].pt/maxpho.pt;
      float xj_truth = maxjet_truth[ir].pt/maxpho_truth.pt;
      int bin = ana::findUnfoldBin(xj,maxpho.pt);
      int bin_truth = ana::findUnfoldBin(xj_truth,maxpho_truth.pt);


      // -----------------------
      // Pairing
      // -----------------------
      if (maxpho.pt > ana::cluster_pt_cut && maxjet[ir].pt > ana::jet_calib_pt_cut[ir]) {
        ispaired[ir] = check_pair(maxjet[ir], ir, maxpho,1);
      }
      if (maxpho_truth.pt > ana::cluster_pt_cut && maxjet_truth[ir].pt > ana::jet_calib_pt_cut[ir]) {
        ispaired_truth[ir] = check_pair(maxjet_truth[ir], ir, maxpho_truth,1);
      }

      // -----------------------
      // Fill response matrix
      // -----------------------
      //if (bin < 0 && ispaired[ir]) cout << "ISSUE RECO!!! xj: " << xj << " pt: " << maxpho.pt << endl;
      //if (bin_truth < 0 && ispaired_truth[ir]) cout << "ISSUE TRUTH!!! xj: " << xj_truth << " pt: " << maxpho_truth.pt << endl;
      bool ismatch = false;
      if (ispaired[ir] && !ispaired_truth[ir] && bin >= 0) {
        jet_response[ir]->Fake(maxjet[ir].pt);
        pho_response[ir]->Fake(maxpho.pt);
        jet_response2D[ir]->Fake(bin);
        hphomissfake[ir]->Fill(maxpho.pt,100);
        hjetmissfake[ir]->Fill(maxjet[ir].pt,100);
        hpairmissfake[ir]->Fill(bin,ana::nPtBins*ana::nUnfoldBins);
        if (use_half) {
          jet_response_half[ir]->Fake(maxjet[ir].pt);
          pho_response_half[ir]->Fake(maxpho.pt);
          jet_response_half2D[ir]->Fake(bin);
        }
      }
      else if (!ispaired[ir] && ispaired_truth[ir] && bin_truth >= 0) {
        jet_response[ir]->Miss(maxjet_truth[ir].pt);
        pho_response[ir]->Miss(maxpho_truth.pt);
        jet_response2D[ir]->Miss(bin_truth);
        hphomissfake[ir]->Fill(100,maxpho_truth.pt);
        hjetmissfake[ir]->Fill(100,maxjet_truth[ir].pt);
        hpairmissfake[ir]->Fill(ana::nPtBins*ana::nUnfoldBins,bin_truth);
        if (use_half) {
          jet_response_half[ir]->Miss(maxjet_truth[ir].pt);
          pho_response_half[ir]->Miss(maxpho_truth.pt);
          jet_response_half2D[ir]->Miss(bin_truth);
        }
      }
      else if (ispaired[ir] && ispaired_truth[ir] && bin >= 0 && bin_truth >= 0) {
        ismatch = check_match(maxpho, maxpho_truth, maxjet[ir], maxjet_truth[ir]);
        if (ismatch) {
          jet_response[ir]->Fill(maxjet[ir].pt, maxjet_truth[ir].pt);
          pho_response[ir]->Fill(maxpho.pt, maxpho_truth.pt);
          jet_response2D[ir]->Fill(bin, bin_truth);
        
          hphomissfake[ir]->Fill(maxpho.pt,maxpho_truth.pt);
          hjetmissfake[ir]->Fill(maxjet[ir].pt,maxjet_truth[ir].pt);
          hpairmissfake[ir]->Fill(bin,bin_truth);
          
          if (use_half) {
            jet_response_half[ir]->Fill(maxjet[ir].pt, maxjet_truth[ir].pt);
            pho_response_half[ir]->Fill(maxpho.pt, maxpho_truth.pt);
            jet_response_half2D[ir]->Fill(bin, bin_truth); 
          }
        }
        else {
          jet_response[ir]->Miss(maxjet_truth[ir].pt);
          jet_response[ir]->Fake(maxjet[ir].pt);
          pho_response[ir]->Miss(maxpho_truth.pt);
          pho_response[ir]->Fake(maxpho.pt);
          jet_response2D[ir]->Miss(bin_truth);
          jet_response2D[ir]->Fake(bin);
        
          hphomissfake[ir]->Fill(maxpho.pt,100);
          hjetmissfake[ir]->Fill(maxjet[ir].pt,100);
          hpairmissfake[ir]->Fill(bin,ana::nPtBins*ana::nUnfoldBins);
          hphomissfake[ir]->Fill(100,maxpho_truth.pt);
          hjetmissfake[ir]->Fill(100,maxjet_truth[ir].pt);
          hpairmissfake[ir]->Fill(ana::nPtBins*ana::nUnfoldBins,bin_truth);
          
          
          if (use_half) {
            jet_response_half[ir]->Miss(maxjet_truth[ir].pt);
            jet_response_half[ir]->Fake(maxjet[ir].pt);
            pho_response_half[ir]->Miss(maxpho_truth.pt);
            pho_response_half[ir]->Fake(maxpho.pt);
            jet_response_half2D[ir]->Miss(bin_truth);
            jet_response_half2D[ir]->Fake(bin);
          }
        }
      }

      // ------------------
      // Fill Histograms
      // ------------------
     
      int ipt = ana::findPtBin(maxpho.pt); 
      if (ispaired[ir]) {
        hpairpurden[ir]->Fill(bin);
        
        hrecojetpt[ir]->Fill(maxjet[ir].pt);
        hrecophopt[ir]->Fill(maxpho.pt);
        hrecoxj[ir]->Fill(bin);
        if (!use_half) {
          hrecojetpt_half[ir]->Fill(maxjet[ir].pt);
          hrecophopt_half[ir]->Fill(maxpho.pt);
          hrecoxj_half[ir]->Fill(bin);
        }
      }

      if (ispaired_truth[ir]) {
        hpaireffden[ir]->Fill(bin_truth);

        htruthjetpt[ir]->Fill(maxjet_truth[ir].pt);
        htruthphopt[ir]->Fill(maxpho_truth.pt);
        htruthxj[ir]->Fill(bin_truth);
        if (!use_half) {
          htruthjetpt_half[ir]->Fill(maxjet_truth[ir].pt);
          htruthphopt_half[ir]->Fill(maxpho_truth.pt);
          htruthxj_half[ir]->Fill(bin_truth);
        }
      }

      if (ismatch) {
        hpairpurnum[ir]->Fill(bin);
        hpaireffnum[ir]->Fill(bin_truth);
      }


      // Drawing event displays
      if (dodraw && ndraw < 100 && ir == 1 && (!ispaired_truth[ir] && ispaired[ir])) {
        if (!ispaired[ir]) {
          float dphi = maxjet[ir].deltaPhi(maxpho);
          int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
          cout << "Event " << ndraw + 1 << " reco failed because: ";
          if (iabcd != 0) cout << endl << "abcd cut: BDT: " << maxpho.bdt << " ISO: " << maxpho.iso4;
          if (abs(maxpho.eta) >= ana::etacut) cout << endl << "photon eta: " << abs(maxpho.eta);
          if (abs(maxjet[ir].eta) >= ana::etacut - ana::JetRs[ir]) cout << endl << "jet eta: " << abs(maxjet[ir].eta);
          if (dphi <= ana::oppcut) cout << endl << "dphi: " << dphi;
          if (maxpho.pt <= ana::cluster_pt_cut) cout << endl << "cluster pt: " << maxpho.pt;
          if (maxjet[ir].pt <= ana::jet_calib_pt_cut[ir]) cout << endl << "jet pt: " << maxjet[ir].pt;
          cout << endl;
        }

        if (!ispaired_truth[ir]) {
          float dphi = maxjet_truth[ir].deltaPhi(maxpho_truth);
          int iabcd = ana::findabcdBin(maxpho_truth.iso4, maxpho_truth.bdt, 0);
          cout << "Event " << ndraw + 1 << " truth failed because: ";
          if (iabcd != 0) cout << endl << "abcd cut: BDT: " << maxpho_truth.bdt << " ISO: " << maxpho_truth.iso4;
          if (abs(maxpho_truth.eta) >= ana::etacut) cout << endl << "photon eta: " << abs(maxpho_truth.eta);
          if (abs(maxjet_truth[ir].eta) >= ana::etacut - ana::JetRs[ir]) cout << endl << "jet eta: " << abs(maxjet_truth[ir].eta);
          if (dphi <= ana::oppcut) cout << endl << "dphi: " << dphi;
          if (maxpho_truth.pt <= ana::cluster_pt_cut) cout << endl << "cluster pt: " << maxpho_truth.pt;
          if (maxjet_truth[ir].pt <= ana::jet_calib_pt_cut[ir]) cout << endl << "jet pt: " << maxjet_truth[ir].pt;
          cout << endl;
        }
        

        TH2D * h = new TH2D("heventdisplay",";eta;phi",100,-1.5,1.5,100,-M_PI,M_PI);
        h->Draw();

        TMarker * star = new TMarker(maxpho.eta, maxpho.phi, 29);
        star->SetMarkerColor(kRed);
        star->SetMarkerSize(2);
        TMarker * star_truth = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 29);
        star_truth->SetMarkerColor(kBlack);
        star_truth->SetMarkerSize(3);
        
        TMarker * circle = new TMarker(maxjet[ir].eta, maxjet[ir].phi, 20);
        circle->SetMarkerColorAlpha(kRed,0.5);
        circle->SetMarkerSize(15);
        TMarker * circle_truth = new TMarker(maxjet_truth[ir].eta, maxjet_truth[ir].phi, 20);
        circle_truth->SetMarkerColorAlpha(kBlack, 0.5);
        circle_truth->SetMarkerSize(15);
        
        // For the legend
        TMarker * dummy_truth = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 29);
        dummy_truth->SetMarkerColor(kBlack);
        dummy_truth->SetMarkerSize(2);
        TMarker * dummy_circle = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 20);
        dummy_circle->SetMarkerColorAlpha(kRed,0.5);
        dummy_circle->SetMarkerSize(2);
        TMarker * dummy_circle_truth = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 20);
        dummy_circle_truth->SetMarkerColorAlpha(kBlack,0.5);
        dummy_circle_truth->SetMarkerSize(2);

        if (maxpho_truth.pt > 0) star_truth->Draw();
        if (maxpho.pt > 0) star->Draw();
        if (maxjet[ir].pt > 0) circle->Draw();
        if (maxjet_truth[ir].pt > 0) circle_truth->Draw();

        TLegend * l = new TLegend(0,0.8,.3,1);
        l->AddEntry(star, "Reco Cluster");
        l->AddEntry(dummy_truth, "Truth Cluster");
        l->AddEntry(dummy_circle, "Reco Jet");
        l->AddEntry(dummy_circle_truth, "Truth Jet");
        l->Draw();

        TLatex latex;
        latex.SetNDC();
       
        if (!ispaired[ir]) latex.DrawLatex(0.3,0.95,"Reco Failed");
        if (!ispaired_truth[ir]) latex.DrawLatex(0.3,0.9,"Truth Failed");

        ndraw++;
        c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/event_displays_%s.pdf",trigger.c_str()));
        delete h;
        delete l;
        c->Clear();

        if (ndraw == 100) c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/event_displays_%s.pdf]",trigger.c_str()));
      }
    }
  }
}

void unfolder::unfold() {
  // -----------------------
  // Unfold
  // -----------------------

  for (int ir = 0; ir < ana::nJetR; ir++) {
    RooUnfoldBayes jet_unfold(jet_response[ir], hrecojetpt[ir], niterate, 0, 1);
    hunfoldjetpt[ir] = (TH1D*)jet_unfold.Hreco();
    hunfoldjetpt[ir]->SetName(Form("hunfoldjetpt%i", ir));
    hjetresponse[ir] = (TH2D*)jet_response[ir]->Hresponse();
    hjetresponse[ir]->SetName(Form("hjetresponse%i", ir));
  
    RooUnfoldBayes pho_unfold(pho_response[ir], hrecophopt[ir], niterate, 0, 1);
    hunfoldphopt[ir] = (TH1D*)pho_unfold.Hreco();
    hunfoldphopt[ir]->SetName(Form("hunfoldphopt%i",ir));
    hphoresponse[ir] = (TH2D*)pho_response[ir]->Hresponse();
    hphoresponse[ir]->SetName(Form("hphoresponse%i",ir));
    
    RooUnfoldBayes jet_unfold_half(jet_response_half[ir], hrecojetpt_half[ir], niterate, 0, 1);
    hunfoldjetpt_half[ir] = (TH1D*)jet_unfold_half.Hreco();
    hunfoldjetpt_half[ir]->SetName(Form("hunfoldjetpt_half%i", ir));
    hjetresponse_half[ir] = (TH2D*)jet_response_half[ir]->Hresponse();
    hjetresponse_half[ir]->SetName(Form("hjetresponse_half%i", ir));
  
    RooUnfoldBayes pho_unfold_half(pho_response_half[ir], hrecophopt_half[ir], niterate, 0, 1);
    hunfoldphopt_half[ir] = (TH1D*)pho_unfold_half.Hreco();
    hunfoldphopt_half[ir]->SetName(Form("hunfoldphopt_half%i",ir));
    hphoresponse_half[ir] = (TH2D*)pho_response_half[ir]->Hresponse();
    hphoresponse_half[ir]->SetName(Form("hphoresponse_half%i",ir));
    
    RooUnfoldBayes jet_unfold2D(jet_response2D[ir], hrecoxj[ir], niterate, 0, 1);
    hunfoldxj[ir] = (TH1D*)jet_unfold2D.Hreco();
    hunfoldxj[ir]->SetName(Form("hunfoldxj%i",  ir));
    hxjresponse[ir] = (TH2D*)jet_response2D[ir]->Hresponse();
    hxjresponse[ir]->SetName(Form("hxjresponse%i", ir));

    RooUnfoldBayes jet_unfold_half2D(jet_response_half2D[ir], hrecoxj_half[ir], niterate, 0, 1);
    hunfoldxj_half[ir] = (TH1D*)jet_unfold_half2D.Hreco();
    hunfoldxj_half[ir]->SetName(Form("hunfoldxj_half%i", ir));
    hxjresponse_half[ir] = (TH2D*)jet_response_half2D[ir]->Hresponse();
    hxjresponse_half[ir]->SetName(Form("hxjresponse_half%i", ir));
  }
}

void unfolder::savehists(TH1D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void unfolder::savehists(TH2D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void unfolder::savehists(RooUnfoldResponse * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void unfolder::savehists(TEfficiency * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}

void unfolder::end() {
  const char * wfilename = Form("/home/samson72/sphnx/gammajet/hists/%s_%s_unfolding.root",trigger.c_str(),sim.c_str());
  cout << "Writing files to " << wfilename << endl;
  TFile::Open(wfilename, "RECREATE");

  savehists(hphodr,ana::nJetR);
  savehists(hjetdr,ana::nJetR);

  savehists(hrecojetpt,ana::nJetR);
  savehists(htruthjetpt,ana::nJetR);
  savehists(hunfoldjetpt,ana::nJetR);
  savehists(hjetresponse,ana::nJetR);
  savehists(hrecophopt,ana::nJetR);
  savehists(htruthphopt,ana::nJetR);
  savehists(hunfoldphopt,ana::nJetR);
  savehists(hphoresponse,ana::nJetR);

  savehists(hrecojetpt_half,ana::nJetR);
  savehists(htruthjetpt_half,ana::nJetR);
  savehists(hunfoldjetpt_half,ana::nJetR);
  savehists(hjetresponse_half,ana::nJetR);
  savehists(hrecophopt_half,ana::nJetR);
  savehists(htruthphopt_half,ana::nJetR);
  savehists(hunfoldphopt_half,ana::nJetR);
  savehists(hphoresponse_half,ana::nJetR);

  savehists(hrecoxj,ana::nJetR);
  savehists(htruthxj,ana::nJetR);
  savehists(hunfoldxj,ana::nJetR);
  savehists(hxjresponse,ana::nJetR);

  savehists(hrecoxj_half,ana::nJetR);
  savehists(htruthxj_half,ana::nJetR);
  savehists(hunfoldxj_half,ana::nJetR);
  savehists(hxjresponse_half,ana::nJetR);

  savehists(jet_response2D,ana::nJetR);

  savehists(hphopurden,ana::nJetR);
  savehists(hphopurnum,ana::nJetR);
  savehists(hphoeffden,ana::nJetR);
  savehists(hphoeffnum,ana::nJetR);
  savehists(hjetpurden,ana::nJetR);
  savehists(hjetpurnum,ana::nJetR);
  savehists(hjeteffden,ana::nJetR);
  savehists(hjeteffnum,ana::nJetR);
  savehists(hpairpurden,ana::nJetR);
  savehists(hpairpurnum,ana::nJetR);
  savehists(hpaireffden,ana::nJetR);
  savehists(hpaireffnum,ana::nJetR);

  for (int ir = 0; ir < ana::nJetR; ir++) {
    hphoeff[ir] = new TEfficiency(*hphoeffnum[ir], *hphoeffden[ir]);
    hphopur[ir] = new TEfficiency(*hphopurnum[ir], *hphopurden[ir]);
    hjeteff[ir] = new TEfficiency(*hjeteffnum[ir], *hjeteffden[ir]);
    hjetpur[ir] = new TEfficiency(*hjetpurnum[ir], *hjetpurden[ir]);
    hpaireff[ir] = new TEfficiency(*hpaireffnum[ir], *hpaireffden[ir]);
    hpairpur[ir] = new TEfficiency(*hpairpurnum[ir], *hpairpurden[ir]);

    hphoeff [ir]->SetName(Form("hphoeff%i",ir));
    hphopur [ir]->SetName(Form("hphopur%i",ir));
    hjeteff [ir]->SetName(Form("hjeteff%i",ir));
    hjetpur [ir]->SetName(Form("hjetpur%i",ir));
    hpaireff[ir]->SetName(Form("hpaireff%i",ir));
    hpairpur[ir]->SetName(Form("hpairpur%i",ir));
  }
  savehists(hphoeff,ana::nJetR);
  savehists(hphopur,ana::nJetR);
  savehists(hjeteff,ana::nJetR);
  savehists(hjetpur,ana::nJetR);
  savehists(hpaireff,ana::nJetR);
  savehists(hpairpur,ana::nJetR);

  savehists(hphomissfake,ana::nJetR);
  savehists(hjetmissfake,ana::nJetR);
  savehists(hpairmissfake,ana::nJetR);
}
