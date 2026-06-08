#include "/home/samson72/sphnx/gammajet/src/histmaker.h"
using namespace std;

bool histmaker::check_pair(jet_object jet, int ir, pho_object pho) {
  float dphi = jet.deltaPhi(pho);
  int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
  
  //if (iabcd != 0) return false;
  if (fabs(pho.eta) > ana::etacut) return false;
  if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;
  if (!isMC && fabs(pho.t - jet.t) > ana::tcut) return false;
  
  return true;
}

// Doesn't actually loop. Just checks that the found photon and found jet are a correct match for each other.
// This is where the xJ plot is filled
bool histmaker::loop(jet_object jet, int ir, pho_object pho, int icalib, float weight, pho_object truthpho) {
  bool useshowershape = 0;
  float dphi = jet.deltaPhi(pho);
  float jetval = jet.pt;
  float phoval = pho.pt;
  float val = jetval/phoval;
  int ipt = ana::findPtBin(pho.pt);
  float lowval = ana::jet_calib_pt_cut[ir]/ana::ptBins[ipt];
  float lowbin = ((int)(lowval/0.08 + 1))*0.08;
  if (val < lowbin) return false;
  bool issingle = pho.pt > ana::singleptlow && pho.pt < ana::singlepthigh;
  if (icalib == 0) {
    hdeltaphi[ipt][ir]->Fill(dphi);
  }

  bool passabcd = ana::findabcdBin(pho.iso4, pho.showershape, 0) == 0;
 
  //cout << pho.eta << " " << jet.eta << " " << dphi << " " << ipt << endl; 
  if (fabs(pho.eta) > ana::etacut) return false;
  if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (!isMC && fabs(pho.t - jet.t) > ana::tcut) return false;
  if (ipt < 0) return false;
  if (useshowershape && pho.showershape == 0) return false;
  if (passabcd && icalib == 3) hjetpt_formarzia[ir][0][0]->Fill(jet.pt);
  if (dphi < ana::oppcut) return false;
  if (passabcd && icalib == 3) hjetpt_formarzia[ir][0][1]->Fill(jet.pt);
  
  if (icalib == 3) hjetpt_formarzia[ir][1][0]->Fill(jet.pt);
  if (passabcd && icalib == 3) hjetpt_formarzia[ir][1][1]->Fill(jet.pt);
  // Get the ABCD info and fill
  int iabcd[ana::nIsoBdtBins];
  bool b_hasthirdjet = hasthirdjet[ir]; 
  for (int iib = 0; iib < ana::nIsoBdtBins; iib++) {
    if (useshowershape) {
      iabcd[iib] = ana::findabcdBin(pho.iso4, pho.showershape, iib);
    }
    else {
      //iabcd[iib] = ana::findabcdBin(pho.iso4, pho.bdt, pho.pt, iib);
      iabcd[iib] = ana::findabcdBin(pho.iso4, pho.bdt, iib);
    }
    if (iabcd[iib] == -1) continue;// || !passabcd) continue;
    
    hratio[ipt][ir][icalib][iib][0][iabcd[iib]]->Fill(val, weight); 
    
    if (!b_hasthirdjet) hratio[ipt][ir][icalib][iib][1][iabcd[iib]]->Fill(val);
    if (issingle && ir == 1 && icalib == 2 && iib == 0 && iabcd[iib] == 0) hratiosingle->Fill(val);


    if (iib == 0 && iabcd[iib] == 0 && icalib == 3) {
      int ief = ana::findEmfracBin(jet.emfrac);
      hratio_emfrac[ipt][ir][ief]->Fill(val,weight);
      //cout << pho.deltaR(truthpho) << endl;
      if (pho.deltaR(truthpho) < 0.05) hratio_direct[ipt][ir][0]->Fill(val,weight);
      else hratio_direct[ipt][ir][1]->Fill(val,weight);
    }
  }
 
  if (icalib == 3) { 
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
    float treepho;
    float treejet;
    
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
      
      keepMC = check_keep_MC(truth_cluster_pt, cluster_pt, truth_jet_pt, jet_pt_smear, trigger);
      //cout << truth_cluster_pt << " " << cluster_pt << " " <<  truth_jet_pt[2] << " " << jet_pt_smear[2] << " " << trigger << endl; 
      //cout << keepMC[2]  << endl;
      //for( int i = 0; i < keepMC.size(); i++) {
      //  cout << keepMC[i] << " ";
      //}
      //cout << endl;
      if (!keepMC.at(keepMC.size()-1)) continue;
      //if (keepMC[2] && cluster_pt > 14) cout << "AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" << truth_cluster_pt << " " << cluster_pt << " " <<  truth_jet_pt[2] << " " << jet_pt_smear[2] << " " << trigger << endl;

      htruthclusterpt->Fill(truth_cluster_pt);
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
    pho_object truthpho = pho_object(
        truth_cluster_pt, 
        truth_cluster_e,
        truth_cluster_eta, 
        truth_cluster_phi, 
        0,
        00,
        0, 
        .99, 
        2
    ); 
    if (!isMC) truthpho = maxpho;
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
    pho_object maxpho_scalehigh = pho_object(
        cluster_pt + (isMC ? (cluster_pt*0.011) : 0),
        cluster_e,
        cluster_eta,
        cluster_phi,
        cluster_showershape[8],
        cluster_showershape[9],
        cluster_time, 
        cluster_bdt_scores[9], 
        pho_object::get_showershape(cluster_showershape, cluster_pt)
    ); 
    pho_object maxpho_scalelow = pho_object(
        cluster_pt - (isMC ? (cluster_pt*0.011) : 0),
        cluster_e,
        cluster_eta,
        cluster_phi,
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
      //cout << "ir: " << ir << " keep: " << keepMC.at(ir) << endl;
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
      bool userecalib = rand->Integer(2);
      if (isMC) {
        maxjet_smear[ir] = jet_object(
          //(userecalib ? jet_pt_smear[ir] / 0.975 : jet_pt_smear[ir]), 
          jet_pt_smear[ir], 
          //jet_pt_recalib[ir], 
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
        //loop(maxjet[ir],  ir, maxpho, 0, 1, truthpho);
      }
      if (maxjet_calib[ir].pt > ana::jet_calib_pt_cut[ir]) {
        //loop(maxjet_calib[ir], ir, maxpho, 1, 1, truthpho);
      }
      if (maxjet_smear[ir].pt > ana::jet_calib_pt_cut[ir]) {
        loop(maxjet_smear[ir], ir, maxpho, 2,  1, truthpho); // Normal smearing
        ispaired[ir] = loop(maxjet_smear[ir], ir, maxpho, 3, (isMC ? reweight(maxpho.pt,vz,maxjet[ir].emfrac) : 1.0), truthpho); // smearing with weights
        loop(maxjet_smear[ir], ir, (isMC ? maxpho_smear : maxpho), 4,  1, truthpho); // smearing the photon, too
        loop(maxjet_smear[ir], ir, maxpho_scalelow, 0, (isMC ? reweight(maxpho.pt,vz,maxjet[ir].emfrac) : 1.0), truthpho); // photon scale
        loop(maxjet_smear[ir], ir, maxpho_scalehigh , 1, (isMC ? reweight(maxpho.pt,vz,maxjet[ir].emfrac) : 1.0), truthpho); // photon scale
      }
      if (maxjet_smear_high[ir].pt > ana::jet_calib_pt_cut[ir]) {
        loop(maxjet_smear_high[ir], ir, maxpho, 5,  (isMC ? reweight(maxpho.pt,vz,maxjet[ir].emfrac) : 1.0), truthpho);
      }
      if (maxjet_smear_low[ir].pt > ana::jet_calib_pt_cut[ir]) {
        loop(maxjet_smear_low[ir], ir, maxpho, 6,  (isMC ? reweight(maxpho.pt,vz,maxjet[ir].emfrac) : 1.0), truthpho);
      }
    
      int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
      if (iabcd == 0 && ispaired[ir] && ir == 2) {
        
        int ptbin = ana::findPtBin(maxpho.pt);
        float val = maxjet_smear[ir].pt/maxpho.pt;
        float lowval = ana::jet_calib_pt_cut[ir]/ana::ptBins[ipt];
        float lowbin = ((int)(lowval/0.08 + 1))*0.08;
        if (val > lowbin) {
          if (userecalib) {
            treepho = maxpho.pt;
            treejet = maxjet_smear[ir].pt;
            h1_forjin[ptbin]->Fill(val);
          }
          else {
            h2_forjin[ptbin]->Fill(val);
          }
        }
      }

      anypaired |= ispaired[ir];


      if (isMC && ispaired[ir]) {
        int ihadron = ana::findHadronBin(hadron_p[ir]);
        hhadronp[ihadron][ir];
      }

      // Fill skimmed ttrees
      if (check_pair(maxjet[ir], ir, maxpho)) {
        //int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, maxpho.pt, 0);
        //int iabcd_bdt = ana::findabcdBin(maxpho.iso4, maxpho.bdt, maxpho.pt, 1);
        //int iabcd_iso = ana::findabcdBin(maxpho.iso4, maxpho.bdt, maxpho.pt, 2);
        int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
        int iabcd_bdt = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 1);
        int iabcd_iso = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 2);
        float weight = (isMC ? reweight(maxpho.pt, vz, maxjet[ir].emfrac) : 1.0);
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
            if (jet_pt_smear[ir] > ana::jet_calib_pt_cut[ir]) { // photon + scale
              outtree_pho_pt_scalehigh[ir] = maxpho_scalehigh.pt;
              outtree_jet_pt_scalehigh[ir] = maxjet_smear[ir].pt;
              outtree_weight_scalehigh[ir] = weight;
              outtree_scalehigh[ir]->Fill();
            }
            if (jet_pt_smear[ir] > ana::jet_calib_pt_cut[ir]) { // photon - scale
              outtree_pho_pt_scalelow[ir] = maxpho_scalelow.pt;
              outtree_jet_pt_scalelow[ir] = maxjet_smear[ir].pt;
              outtree_weight_scalelow[ir] = weight;
              outtree_scalelow[ir]->Fill();
            }
            if (ir == 2 && userecalib && jet_pt_smear[ir] > ana::jet_calib_pt_cut[ir]) { // Jin test
              float val =  maxjet_smear[ir].pt / maxpho.pt;
              float lowval = ana::jet_calib_pt_cut[ir]/ana::ptBins[ipt];
              float lowbin = ((int)(lowval/0.08 + 1))*0.08;
              if (val > lowbin && val < 2) { 
                if (maxpho.pt != treepho || maxjet_smear[ir].pt != treejet) {
                  //cout << "BAD NUMBERS!!" << endl;
                  //cout << maxpho.pt << " " << treepho << " " << maxjet_smear[ir].pt << " " << treejet << endl;
                }
                outtree_pho_pt_forjin[ir] = maxpho.pt;
                outtree_jet_pt_forjin[ir] = maxjet_smear[ir].pt;
                outtree_weight_forjin[ir] = scalemap[trigger];
                //cout << setprecision(15) << scalemap[trigger] << endl;
                outtree_forjin[ir]->Fill();
              }
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
    for (int i = 0; i < ana::nJetR; i++) {  
      if (keepMC.at(i) && ispaired[i]) htruthjetpt[i]->Fill(truth_jet_pt[i]);
    }

    // max cluster histos
    hisobdt[ipt]->Fill(maxpho.iso4,maxpho.bdt); 
    hiso[0]->Fill(maxpho.iso3);
    hiso[1]->Fill(maxpho.iso4);
    hiso2d[0]->Fill(maxpho.pt,maxpho.iso3);
    hiso2d[1]->Fill(maxpho.pt,maxpho.iso4);
    
    hclusterpt->Fill(maxpho.pt);
    hclustereta->Fill(maxpho.eta); 
    hclusteretaphi->Fill(maxpho.eta,maxpho.phi); 
    
    //int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, maxpho.pt, 0);
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
      outtree_scalehigh[ir]->Write();
      outtree_scalelow[ir]->Write();
      outtree_forjin[ir]->Write();
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
  for (int i = 0; i < ana::nPtBins; i++) {
    for (int j = 0; j < ana::nJetR; j++) {
      for (int k = 0; k < ana::nEmfracBins; k++) {
        hratio_emfrac[i][j][k]->Write();
      }
    }
  }
  for (int i = 0; i < ana::nPtBins; i++) {
    for (int j = 0; j < ana::nJetR; j++) {
      for (int k = 0; k < 2; k++) {
        hratio_direct[i][j][k]->Write();
      }
    }
  }

  hvz->Write();

  savehists(hisobdt,ana::nPtBins);
  savehists(h1_forjin,ana::nPtBins);
  savehists(h2_forjin,ana::nPtBins);
  savehists(hclusterptabcd,4); // for ABCD
  
  savehists(hemfrac,ana::nPtBins,ana::nJetR);
  savehists(hdeltaphi,ana::nPtBins,ana::nJetR);
  savehists(hdeltaphiprecut,ana::nPtBins,ana::nJetR);
  savehists(h3jetdeltar,ana::nPtBins, ana::nJetR);
  savehists(hjetetaxj,ana::nxjBins,ana::nJetR);
  savehists(hxjbdt,ana::nBdtBins,ana::nJetR);
  savehists(hhadronp,ana::nHadronBins,ana::nJetR);
  

  for (int i = 0; i < ana::nJetR; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        hjetpt_formarzia[i][j][k]->Write();
      }
    }
  }
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

