#include "/home/samson72/sphnx/gammajet/src/histmaker.h"
using namespace std;

// Doesn't actually loop. Just checks that the found photon and found jet are a correct match for each other.
// This is where the xJ plot is filled
bool histmaker::loop(jet_object jet, int ir, pho_object pho, int icalib, float weight) {
  bool useshowershape = 0;
  float dphi = jet.deltaPhi(pho);
  float jetval = jet.pt;
  float phoval = pho.pt;
  float val = jetval/phoval;
  int ipt = ana::findPtBin(pho.pt);
  bool issingle = pho.pt > ana::singleptlow && pho.pt < ana::singlepthigh;
  if (icalib == 0) {
    hdeltaphi[ipt][ir]->Fill(dphi);
  }
 
  //cout << pho.eta << " " << jet.eta << " " << dphi << " " << ipt << endl; 
  if (abs(pho.eta) > ana::etacut) return false;
  if (abs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (!isMC && abs(pho.t - jet.t) > ana::tcut) return false;
  if (dphi < ana::oppcut) return false;
  if (ipt < 0) return false;
  if (useshowershape && pho.showershape == 0) return false;
  
  // Get the ABCD info and fill
  int iabcd[ana::nIsoBdtBins];
  bool b_hasthirdjet = hasthirdjet[ir]; 
  for (int iib = 0; iib < ana::nIsoBdtBins; iib++) {
    if (useshowershape) {
      iabcd[iib] = ana::findabcdBin(pho.iso4, pho.showershape, iib);
    }
    else {
      iabcd[iib] = ana::findabcdBin(pho.iso4, pho.bdt, iib);
    }
    if (iabcd[iib] == -1) continue;
    
    hratio[ipt][ir][icalib][iib][0][iabcd[iib]]->Fill(val, weight); 
    
    if (!b_hasthirdjet) hratio[ipt][ir][icalib][iib][1][iabcd[iib]]->Fill(val);
    if (issingle && ir == 1 && icalib == 2 && iib == 0 && iabcd[iib] == 0) hratiosingle->Fill(val);
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
    if (iabcd[0] == 0 && hasthirdjet[ir]) {
      h3jetpt[ir]->Fill(pho.pt,thirdjet_pt[ir]);
      h3jetdeltar[ipt][ir]->Fill(thirdjet_dr[ir]);
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
  cout << "isMC " << isMC << endl;
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    //t0corr = t0map[RunNumber];
    if(e % 1000==0) std::cout << "entry " << e << "/" << nentries << " (" << (float)e/nentries*100. << "%)" << "\t\r" << std::flush;
    //if (RunNumber != 51154) continue;
    if (fabs(vz) > ana::vzcut) continue;
    if (!isMC && !ScaledTriggerBit[27] && !ScaledTriggerBit[38]) continue;
    mbd_t0 = mbd_time;// - t0corr;
    
    // Fill time histos before any cuts
    if (!isMC) {
      hmbdt->Fill(mbd_t0);
      if (jet_pt[1] != 0) hjett->Fill(jet_time[1]);
      if (jet_pt[1] != 0) hmjt->Fill(mbd_t0,jet_time[1]);

      for (int i = 0; i < ana::nJetR; i++) {
        if (jet_pt[i] != 0) hmtminusjt[i]->Fill(mbd_t0-jet_time[i]);
      }

      hclustert->Fill(cluster_time);
      hmct->Fill(mbd_t0,cluster_time);
      hmtminusct->Fill(mbd_t0-cluster_time);
      if (jet_pt[1] != 0) hcjt->Fill(cluster_time,jet_time[1]);
      for (int j = 0; j < ana::nJetR; j++) {
        if (jet_pt[j] != 0) hctminusjt[j]->Fill(cluster_time-jet_time[j]);
      }
    }
    
    // Check if the MC event should be kept
    vector<bool> keepMC(ana::nJetR + 1);
    for (int i = 0; i < ana::nJetR + 1; i++) {
      keepMC.at(i) = 1;
    }
    if (isMC) {
      if (truth_cluster_pt > 0) htruthclusterptprecut->Fill(truth_cluster_pt);
      for (int i = 0; i < ana::nJetR; i++) {  
        if (truth_jet_pt[i] > 0) htruthjetptprecut[i]->Fill(truth_jet_pt[i]);
      }
      
      vector<bool> keepMC = check_keep_MC(truth_cluster_pt, truth_jet_pt, trigger);
      //for( int i = 0; i < keepMC.size(); i++) {
      //  cout << keepMC[i] << " ";
      //}
      //cout << endl;
      if (!keepMC.at(keepMC.size()-1)) continue;

      htruthclusterpt->Fill(truth_cluster_pt);
      for (int i = 0; i < ana::nJetR; i++) {  
        if (keepMC.at(i)) htruthjetpt[i]->Fill(truth_jet_pt[i]);
      }
    }

    for (int i = 0; i < 11; i++) {
      hbdt[i]->Fill(cluster_bdt_scores[i]);
    }
    
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
    float newE = rand->Gaus(cluster_e, cluster_energy_smear_func->Eval(cluster_e));
    newE *= 1.007;
    float newPhi = rand->Gaus(cluster_phi, cluster_position_smear_func->Eval(cluster_e));
    float newEta = rand->Gaus(cluster_eta, cluster_position_smear_func->Eval(cluster_e));

    pho_object maxpho_smear = pho_object(
        newE/TMath::CosH(newEta),
        newE,
        newEta,
        newPhi,
        cluster_showershape[8],
        cluster_showershape[9],
        cluster_time, 
        cluster_bdt_scores[9], 
        pho_object::get_showershape(cluster_showershape, cluster_pt)
    ); 
    if (maxpho.pt > ana::cluster_pt_cut) {
      hclusterphi->Fill(maxpho.phi);
    }

    vector<jet_object> maxjet(ana::nJetR);
    vector<jet_object> maxjet_calib(ana::nJetR);
    vector<jet_object> maxjet_smear(ana::nJetR);
    vector<jet_object> maxjet_smear_high(ana::nJetR);
    vector<jet_object> maxjet_smear_low(ana::nJetR);
    int ipt = ana::findPtBin(maxpho.pt);
    if (ipt < 0) continue; 

    bool ispaired[ana::nJetR] = { 0 };
    bool anypaired = false;
    for (int ir = 0; ir < ana::nJetR; ir++) {
      if (!keepMC.at(ir)) continue;
      // one for uncalib, calibrated, and JER smeared`
      maxjet[ir] = jet_object(
          jet_pt[ir], 
          jet_e[ir],
          jet_eta[ir], 
          jet_phi[ir], 
          jet_emfrac[ir], 
          0, 
          0, 
          jet_time[ir]
      );
      maxjet_calib[ir] = jet_object(
          jet_pt_calib[ir], 
          jet_e[ir],
          jet_eta[ir], 
          jet_phi[ir], 
          jet_emfrac[ir], 
          0, 
          0, 
          jet_time[ir]
      );
      if (isMC) {
        maxjet_smear[ir] = jet_object(
          jet_pt_smear[ir], 
          jet_e[ir],
          jet_eta[ir], 
          jet_phi[ir], 
          jet_emfrac[ir], 
          0, 
          0, 
          jet_time[ir]
        );
        maxjet_smear_high[ir] = jet_object(
          jet_pt_smear_high[ir], 
          jet_e[ir],
          jet_eta[ir], 
          jet_phi[ir], 
          jet_emfrac[ir], 
          0, 
          0, 
          jet_time[ir]
        );
        maxjet_smear_low[ir] = jet_object(
          jet_pt_smear_low[ir], 
          jet_e[ir],
          jet_eta[ir], 
          jet_phi[ir], 
          jet_emfrac[ir], 
          0, 
          0, 
          jet_time[ir]
        );
      }
      else {
        maxjet_smear[ir] = maxjet_calib[ir];
        maxjet_smear_high[ir] = maxjet_calib[ir];
        maxjet_smear_low[ir] = maxjet_calib[ir];
      }
      
      if (maxjet[ir].pt > ana::jet_pt_cut[ir]) {
        hdeltar[ir]->Fill(maxjet[ir].deltaR(maxpho));
        hjetphi[ir]->Fill(maxjet[ir].phi);
        hdeltaphiprecut[ipt][ir]->Fill(maxjet[ir].deltaPhi(maxpho));
      }
      // Fill the xJ histograms
      if (maxjet[ir].pt > ana::jet_pt_cut[ir]) {
        loop(maxjet[ir],  ir, maxpho, 0);
      }
      if (maxjet_calib[ir].pt > ana::jet_calib_pt_cut[ir]) {
        loop(maxjet_calib[ir], ir, maxpho, 1);
      }
      if (maxjet_smear[ir].pt > ana::jet_calib_pt_cut[ir]) {
        //if (ir == 1) maxpho.print();
        ispaired[ir] = loop(maxjet_smear[ir], ir, maxpho, 2);
        loop(maxjet_smear[ir], ir, maxpho, 3, (isMC ? reweight(maxpho.pt,vz) : 1.0));
        loop(maxjet_smear[ir], ir, (isMC ? maxpho_smear : maxpho), 4);
      }
      if (maxjet_smear_high[ir].pt > ana::jet_calib_pt_cut[ir]) {
        loop(maxjet_smear_high[ir], ir, maxpho, 5);
      }
      if (maxjet_smear_low[ir].pt > ana::jet_calib_pt_cut[ir]) {
        loop(maxjet_smear_low[ir], ir, maxpho, 6);
      }
      
      anypaired |= ispaired[ir];


      if (isMC && ispaired[ir]) {
        int ihadron = ana::findHadronBin(hadron_p[ir]);
        hhadronp[ihadron][ir];
      }

      // Fill skimmed ttrees
      if (ispaired[ir]) {
        int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
        int iabcd_bdt = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 1);
        int iabcd_iso = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 2);
        float weight = (isMC ? reweight(maxpho.pt, vz) : 1.0);
        if (iabcd == 0) { // Nominal
          outtree_pho_pt[ir] = maxpho.pt;
          outtree_jet_pt[ir] = maxjet_smear[ir].pt;
          outtree_weight[ir] = weight;
          outtree[ir]->Fill();
          if (isMC) {
            if (jet_pt_smear_high[ir] > ana::jet_calib_pt_cut[ir]) { // JER + uncertainty
              outtree_pho_pt_JERhigh[ir] = maxpho.pt;
              outtree_jet_pt_JERhigh[ir] = maxjet_smear_high[ir].pt;
              outtree_weight_JERhigh[ir] = weight;
              outtree_JERhigh[ir]->Fill();
            }
            if (jet_pt_smear_low[ir] > ana::jet_calib_pt_cut[ir]) { // JER - uncertainty
              outtree_pho_pt_JERlow[ir] = maxpho.pt;
              outtree_jet_pt_JERlow[ir] = maxjet_smear_low[ir].pt;
              outtree_weight_JERlow[ir] = weight;
              outtree_JERlow[ir]->Fill();
            }
          }
        }
        if (iabcd == 0 && !hasthirdjet[ir]) { // 3jet
          outtree_pho_pt_3jet[ir] = maxpho.pt;
          outtree_jet_pt_3jet[ir] = maxjet_smear[ir].pt;
          outtree_weight_3jet[ir] = weight;
          outtree_3jet[ir]->Fill();
        }
        if (iabcd_bdt == 0) { // narrow BDT
          outtree_pho_pt_bdt[ir] = maxpho.pt;
          outtree_jet_pt_bdt[ir] = maxjet_smear[ir].pt;
          outtree_weight_bdt[ir] = weight;
          outtree_bdt[ir]->Fill();
        }
        if (iabcd_iso == 0) { // narrow isolation energy
          outtree_pho_pt_iso[ir] = maxpho.pt;
          outtree_jet_pt_iso[ir] = maxjet_smear[ir].pt;
          outtree_weight_iso[ir] = weight;
          outtree_iso[ir]->Fill();
        }
      }
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
    
    int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
    hclusterptabcd[iabcd]->Fill(maxpho.pt); 
    
    // cluster and jet kinematic histos
    //float frag =    getZ(maxpho, jets[1], 0);
    //float fragiso = getZ(maxpho, jets[1], 1);
    //if (frag > 0) hfrag->Fill(frag);
    //if (fragiso > 0) hfragiso->Fill(fragiso);
    
    if (maxpho.pt > 0) {
      hclusterptprecut->Fill(maxpho.pt);
    }
    for (int i = 0; i < ana::nJetR; i++) {
      if (maxjet[i].pt > 0) {
        hjetptprecut[i]->Fill(maxjet[i].pt);
      }
    } 
    hvz->Fill(vz);
  }
  // the end
  return;
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
  
  std::cout << "Writing tree to " << treefilename << std::endl;
  outfile->cd();
  for (int ir = 0; ir < ana::nJetR; ir++) {
    outtree[ir]->Write();
    outtree_3jet[ir]->Write();
    outtree_bdt[ir]->Write();
    outtree_iso[ir]->Write(); 
    if (isMC) {
      outtree_JERhigh[ir]->Write();
      outtree_JERlow[ir]->Write();
    }
  }


  const char * wfilename = (trigger == "Data" ? Form("/home/samson72/sphnx/gammajet/hists/hists_%s.root",trigger.c_str()) : Form("/home/samson72/sphnx/gammajet/hists/hists_%s_%s.root",sim.c_str(),trigger.c_str()));
  std::cout << "Writing files to " << wfilename << std::endl;
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

  hvz->Write();

  savehists(hisobdt,ana::nPtBins);
  savehists(hclusterptabcd,4); // for ABCD
  
  savehists(hemfrac,ana::nPtBins,ana::nJetR);
  savehists(hdeltaphi,ana::nPtBins,ana::nJetR);
  savehists(hdeltaphiprecut,ana::nPtBins,ana::nJetR);
  savehists(h3jetdeltar,ana::nPtBins, ana::nJetR);
  savehists(hjetetaxj,ana::nxjBins,ana::nJetR);
  savehists(hxjbdt,ana::nBdtBins,ana::nJetR);
  savehists(hhadronp,ana::nHadronBins,ana::nJetR);
  
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
  savehists(hjetphi,ana::nJetR);
  savehists(hjetetahighem,ana::nJetR);
  savehists(hjetetalowem,ana::nJetR);
  savehists(hjetetaphi,ana::nJetR);
  savehists(hdeltar,ana::nJetR);
  savehists(hJES,ana::nJetR);
  savehists(hjetsmear,ana::nJetR);
  savehists(h3jetpt,ana::nJetR);
  
  savehists(hbdt,11);
  
  hclusterpt->Write();
  hclusterphi->Write();
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
  hratiosingle->Write();

  std::cout << std::endl << "All done with " << trigger << "!" << std::endl;
}

