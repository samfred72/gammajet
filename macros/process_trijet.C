#include "/home/samson72/sphnx/gammajet/src/ana.h"
R__LOAD_LIBRARY(libgammajet.so);

void process_trijet(const char * JER = "", int ir = 1) {
  const char * rname = ana::rnames[ir];
  
  TFile * f;
  if (strcmp(JER, "JERHigh") == 0 )     f = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%spythiaJERHigh.root",rname),"READ");
  else if (strcmp(JER, "JERLow") == 0 ) f = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%spythiaJERLow.root",rname),"READ");
  else if (strcmp(JER, "JERReco") == 0 ) f = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%spythiaJERReco.root",rname),"READ");
  else if (strcmp(JER, "HERWIG") == 0 ) f = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sherwig.root",rname),"READ");
  else f = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%spythia.root",rname),"READ");
  TTree * t0 = (TTree*)f->Get("ttree");
  TTree * t1 = (TTree*)f->Get("ttree8");
  TTree * t2 = (TTree*)f->Get("ttree12");
  TTree * t3 = (TTree*)f->Get("ttree20");
  TTree * t4 = (TTree*)f->Get("ttree30");
  const int ntrees = 4;
  TTree * t[ntrees] = {t1,t2,t3,t4};
  TH1D * h[ntrees][ana::nTrijetPtBins];
  TH1D * h0[ana::nTrijetPtBins];
  for (int i = 0; i < ntrees; i++) {
    for (int j = 0; j < ana::nTrijetPtBins; j++) {
      h[i][j] = new TH1D(Form("h%i_%i",i,j),"",45,0.4,2.65);
      if (i == 0) h0[j] = new TH1D(Form("hdata%i",j),"",45,0.4,2.65);
    }
  }

  for (int it = 0; it < ntrees; it++) {
    cout << "Processing tree: " << it << "..." << endl;
    //t[it]->Show(0);
    float leading_pt;
    float subleading_pt;
    float subleading_eta;
    float subleading_phi;
    float subsubleading_pt;
    float subsubleading_eta;
    float subsubleading_phi;
    float weight;
    t[it]->SetBranchAddress("leadingPT",&leading_pt);
    t[it]->SetBranchAddress("SLPT",     &subleading_pt);
    t[it]->SetBranchAddress("SLeta",        &subleading_eta);
    t[it]->SetBranchAddress("SLphi",        &subleading_phi);
    t[it]->SetBranchAddress("SSLPT",    &subsubleading_pt);
    t[it]->SetBranchAddress("SSLeta",       &subsubleading_eta);
    t[it]->SetBranchAddress("SSLphi",       &subsubleading_phi);
    t[it]->SetBranchAddress("weight", &weight);
    
    int nentries = t[it]->GetEntries();
    for (int e = 0; e < nentries; e++) {
      if(e % 1000==0) std::cout << "entry " << e << "/" << nentries << " (" << (float)e/nentries*100. << "%)" << "\t\r" << std::flush;
      t[it]->GetEntry(e);
      int ipt = ana::findTrijetPtBin(leading_pt);
      //if (it == 1) cout << leading_pt << " "  << ipt << endl;
      if (ipt < 0) continue;
      TLorentzVector v1 = TLorentzVector();
      TLorentzVector v2 = TLorentzVector();
      v1.SetPtEtaPhiE(subleading_pt, subleading_eta, subleading_phi, subleading_pt);
      v2.SetPtEtaPhiE(subsubleading_pt, subsubleading_eta, subsubleading_phi, subsubleading_pt);
      TLorentzVector vsum = v1+v2;
      float multi_pt = vsum.Pt();
      h[it][ipt]->Fill(leading_pt/multi_pt, weight);
    }
  }
    
  cout << "processing data..." << endl;
  float leading_pt;
  float subleading_pt;
  float subleading_eta;
  float subleading_phi;
  float subsubleading_pt;
  float subsubleading_eta;
  float subsubleading_phi;
  float weight;
  t0->SetBranchAddress("leadingPT",&leading_pt);
  t0->SetBranchAddress("SLPT",     &subleading_pt);
  t0->SetBranchAddress("SLeta",        &subleading_eta);
  t0->SetBranchAddress("SLphi",        &subleading_phi);
  t0->SetBranchAddress("SSLPT",    &subsubleading_pt);
  t0->SetBranchAddress("SSLeta",       &subsubleading_eta);
  t0->SetBranchAddress("SSLphi",       &subsubleading_phi);

  int nentries = t0->GetEntries();
  for (int e = 0; e < nentries; e++) {
    if(e % 1000==0) std::cout << "entry " << e << "/" << nentries << " (" << (float)e/nentries*100. << "%)" << "\t\r" << std::flush;
    t0->GetEntry(e);
    int ipt = ana::findTrijetPtBin(leading_pt);
    //if (it == 1) cout << leading_pt << " "  << ipt << endl;
    if (ipt < 0) continue;
    TLorentzVector v1 = TLorentzVector();
    TLorentzVector v2 = TLorentzVector();
    v1.SetPtEtaPhiE(subleading_pt, subleading_eta, subleading_phi, subleading_pt);
    v2.SetPtEtaPhiE(subsubleading_pt, subsubleading_eta, subsubleading_phi, subsubleading_pt);
    TLorentzVector vsum = v1+v2;
    float multi_pt = vsum.Pt();
    h0[ipt]->Fill(leading_pt/multi_pt);
  }

  TFile * wf;
  if (strcmp(JER, "JERHigh") == 0)     wf = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_JERhigh_%s.root",rname),"RECREATE");
  else if (strcmp(JER, "JERLow") == 0) wf = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_JERlow_%s.root",rname),"RECREATE");
  else if (strcmp(JER, "JERReco") == 0) wf = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_JERreco_%s.root",rname),"RECREATE");
  else if (strcmp(JER, "HERWIG") == 0) wf = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_HERWIG_%s.root",rname),"RECREATE");
  else wf = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_%s.root",rname),"RECREATE");
  cout << "writing file to: " << wf->GetName() << endl;
  for (int j = 0; j < ana::nTrijetPtBins; j++) {
    TH1D * hsum = (TH1D*)h[0][j]->Clone(Form("htrijet%i",j));
    hsum->Reset("ICES");
    hsum->Add( h[0][j],  h[1][j]);
    hsum->Add( hsum,     h[2][j]);
    hsum->Add( hsum,     h[3][j]);
    hsum->Write();
  }
  for (int j = 0; j < ana::nTrijetPtBins; j++) {
    h0[j]->Write();
  }
  cout << "All done!" << endl;
}
