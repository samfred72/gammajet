
#include "/home/samson72/sphnx/gammajet/src/ana.h"
#include "/home/samson72/sphnx/gammajet/src/drawer.h"

void make_table() {
  drawer d("pythia");
  drawer dh("herwig");
  //TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/hists

  cout << std::left
    << std::setw(10) << "radius" << " | "
    << std::setw(28) << "gammajet Herwig/Pythia in situ" << " | "
    << std::setw(28) << "trijet Herwig/Pythia in situ"
    << endl;

  // optional separator line
  cout << std::string(10, '-') << "-+-"
    << std::string(28, '-') << "-+-"
    << std::string(28, '-')
    << endl;
  for (int i = 0; i < ana::nJetR; i++) {
    // gammajet pythia
    TH1D * hnumP = new TH1D(Form("hnumP%i",i),"",ana::nPtBins, ana::ptBins);
    TH1D * hdenP = new TH1D(Form("hdenP%i",i),"",ana::nPtBins, ana::ptBins);
    for (int j = 0; j < ana::nPtBins; j++) {
      TH1D *hd = d.get(Form("hratio_%i_%i_2_0_0_0",j,i),0);
      TH1D *hp = d.get(Form("hratio_%i_%i_2_0_0_0",j,i),1);
      hd->Rebin(4);
      hp->Rebin(4);
      int lowbin = (int)(ana::jet_calib_pt_cut[i]/ana::ptBins[j]/0.08 + 1) + 1;
      hd->GetXaxis()->SetRange(lowbin,25);
      hp->GetXaxis()->SetRange(lowbin,25);
      hnumP->SetBinContent(j+1, hd->GetMean());
      hnumP->SetBinError(j+1, hd->GetMeanError());
      hdenP->SetBinContent(j+1, hp->GetMean());
      hdenP->SetBinError(j+1, hp->GetMeanError());
      delete hd;
      delete hp;
    }
    
    // gammajet herwig
    TH1D * hnumH = new TH1D(Form("hnumH%i",i),"",ana::nPtBins, ana::ptBins);
    TH1D * hdenH = new TH1D(Form("hdenH%i",i),"",ana::nPtBins, ana::ptBins);
    for (int j = 0; j < ana::nPtBins; j++) {
      TH1D *hd = dh.get(Form("hratio_%i_%i_2_0_0_0",j,i),0);
      TH1D *hp = dh.get(Form("hratio_%i_%i_2_0_0_0",j,i),1);
      hd->Rebin(4);
      hp->Rebin(4);
      int lowbin = (int)(ana::jet_calib_pt_cut[i]/ana::ptBins[j]/0.08 + 1) + 1;
      hd->GetXaxis()->SetRange(lowbin,25);
      hp->GetXaxis()->SetRange(lowbin,25);
      hnumH->SetBinContent(j+1, hd->GetMean());
      hnumH->SetBinError(j+1, hd->GetMeanError());
      hdenH->SetBinContent(j+1, hp->GetMean());
      hdenH->SetBinError(j+1, hp->GetMeanError());
    }
    
    // trijet pythia
    const char * rname = ana::rnames[i];
    TFile * fP = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_%s.root",rname),"READ");
    TFile * fH = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_HERWIG_%s.root",rname),"READ");
    TH1D * hnum3 = new TH1D(Form("hnum3%i",i),"",ana::nTrijetPtBins, ana::trijetPtBins);
    TH1D * hden3 = new TH1D(Form("hden3%i",i),"",ana::nTrijetPtBins, ana::trijetPtBins);
    for (int j = 0; j < ana::nTrijetPtBins; j++) {
      TH1D *hn = (TH1D*)fP->Get(Form("htrijet%i",j));
      TH1D *hd = (TH1D*)fH->Get(Form("htrijet%i",j));
      hnum3->SetBinContent(j+1, hn->GetMean());
      hnum3->SetBinError(j+1,   hn->GetMeanError());
      hden3->SetBinContent(j+1, hd->GetMean());
      hden3->SetBinError(j+1,   hd->GetMeanError());
      delete hn;
      delete hd;
    }
    
    TH1D * hrat = (TH1D*)hdenP->Clone();
    hrat->Divide(hdenP,hdenH);

    TH1D * hrat3 = (TH1D*)hnum3->Clone();
    hrat3->Divide(hnum3,hden3);
    
    TF1 * f = new TF1(Form("f%i",i),"pol0",10,35);
    hrat->Fit(f,"RQIM0");
    TF1 * f3 = new TF1(Form("f3%i",i),"pol0",30,70);
    hrat3->Fit(f3,"RQIM0");
    cout << std::left
     << std::setw(10) << rname << " | "
     << f->GetParameter(0) << " $\\pm$ " << f->GetParError(0) << " | "
     << f3->GetParameter(0) << " $\\pm$ " << f3->GetParError(0) 
     << endl;
    



  }
}
