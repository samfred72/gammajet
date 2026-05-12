
#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"


// ============================================================
//  MINUIT FIT + COVARIANCE + UNCERTAINTY BAND
// ============================================================

#include <vector>
#include <cmath>
#include <iostream>

#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>

using namespace std;

// ------------------------------------------------------------
// You already have these structs in your code
// ------------------------------------------------------------
struct GammaEvent { float pt, base; int bin; };
struct TrijetEvent { float lead, sl, slphi, ssl, sslphi; int bin; };

// ------------------------------------------------------------
// Helper (same as yours)
// ------------------------------------------------------------
inline float compute_pt_sum(float pt1, float phi1,
                            float pt2, float phi2)
{
  float px = pt1*cos(phi1) + pt2*cos(phi2);
  float py = pt1*sin(phi1) + pt2*sin(phi2);
  return sqrt(px*px + py*py);
}

// ------------------------------------------------------------
// Data container
// ------------------------------------------------------------
struct FitData {
  vector<GammaEvent>* gamma;
  vector<TrijetEvent>* trijet;

  float* meanp;
  float* merrp;
  float* mean3;
  float* merr3;

  int nTot;
};

// ------------------------------------------------------------
// CHI2 FUNCTION
// ------------------------------------------------------------
double compute_chi2(double pa, double pb, FitData& data) {

  //cout << "computing chi2" << endl;
  vector<float> sum(data.nTot, 0.f), sum2(data.nTot, 0.f);
  vector<int> count(data.nTot, 0);

  // gammajet
  for (auto &ev : *data.gamma) {
    float f = pa + pb*ev.pt;
    float x = ev.base / f;

    float lowval = (int)(ana::jet_calib_pt_cut[1]/ana::ptBins[ev.bin]/0.08 + 1) * 0.08;
    if (x > lowval && x <= 2) {
      sum[ev.bin] += x;
      sum2[ev.bin] += x*x;
      count[ev.bin]++;
    }
  }

  // trijet
  for (auto &ev : *data.trijet) {
    float f1 = pa + pb*ev.sl;
    float f2 = pa + pb*ev.ssl;

    float pt1 = ev.sl / f1;
    float pt2 = ev.ssl / f2;

    float multi = compute_pt_sum(pt1, ev.slphi, pt2, ev.sslphi);

    float fL = pa + pb*ev.lead;
    float x = ev.lead / multi / fL;

    int b = ana::nPtBins + ev.bin;
    if (x > 0.4 && x < 2.65) {
      sum[b] += x;
      sum2[b] += x*x;
      count[b]++;
    }
  }

  double chisq = 0;

  for (int i = 0; i < data.nTot; i++) {
    if (count[i] == 0) continue;

    float mean = sum[i]/count[i];
    float var  = sum2[i]/count[i] - mean*mean;
    float err  = sqrt(var/count[i]);

    float refm = (i < ana::nPtBins) ? data.meanp[i] : data.mean3[i-ana::nPtBins];
    float refe = (i < ana::nPtBins) ? data.merrp[i] : data.merr3[i-ana::nPtBins];

    float diff = 1 - mean/refm;
    //float errt = sqrt(err*err + refe*refe);
    float errt = sqrt((err*err) / (refm*refm) + (mean*mean) * (refe*refe) / (refm*refm*refm*refm));

    //cout << setprecision(2) << diff << " " << errt << ",   "; 

    chisq += (diff*diff)/(errt*errt);
  }
  //cout << endl;
  //std::cout << std::defaultfloat;

  //cout << pa << " " << pb << " " << chisq << endl;
  return chisq;
}

//// ------------------------------------------------------------
//// FUNCTOR FOR MINUIT
//// ------------------------------------------------------------
//struct Chi2Functor {
//  FitData& data;
//  Chi2Functor(FitData& d) : data(d) {}
//
//  double operator()(const double *par) {
//    return compute_chi2(par[0], par[1], data);
//  }
//};

// ============================================================
// MAIN ENTRY (call this after filling gamma/trijet etc.)
// ============================================================
void grid_insitu_minuit(const char * form = "nominal", int ir = 2) {

  drawer d;
  if (strcmp(form,"HERWIG") == 0) {
    d = drawer("herwig");
  }

  // -----------------------------
  // Reference (MC) values
  // -----------------------------
  float meanp[ana::nPtBins];
  float merrp[ana::nPtBins];
  float mean3[ana::nTrijetPtBins];
  float merr3[ana::nTrijetPtBins];

  float minjet = ana::jet_calib_pt_cut[ir];
  const char * rname = ana::rnames[ir];

  // Gammajet reference
  for (int i = 0; i < ana::nPtBins; i++ ) {
    float lowcluster = ana::ptBins[i];
    int lowbin = (int)(minjet/lowcluster/0.08 + 1);

    const char * histname;
    if (strcmp(form, "nominal") == 0) {
      histname = Form("hratio_%i_%i_2_0_0_0",i,ir);
    }
    else if (strcmp(form, "3jet") == 0) {
      histname = Form("hratio_%i_%i_2_0_1_0",i,ir);
    }
    else if (strcmp(form, "bdt") == 0) {
      histname = Form("hratio_%i_%i_2_1_0_0",i,ir);
    }
    else if (strcmp(form, "iso") == 0) {
      histname = Form("hratio_%i_%i_2_2_0_0",i,ir);
    }
    else if (strcmp(form, "JERhigh") == 0) {
      histname = Form("hratio_%i_%i_5_0_0_0",i,ir);
    }
    else if (strcmp(form, "JERlow") == 0) {
      histname = Form("hratio_%i_%i_6_0_0_0",i,ir);
    }
    else if (strcmp(form, "HERWIG") == 0) {
      histname = Form("hratio_%i_%i_2_0_0_0",i,ir);
    }
    else {
      cout << "Bad variation!" << endl;
      return;
    }

    TH1D* h = d.get(histname,1);
    h->Rebin(4);
    h->GetXaxis()->SetRange(lowbin+1,25);

    meanp[i] = h->GetMean();
    merrp[i] = h->GetMeanError();
  }

  // Trijet reference
  TFile * ftrijet;
  if (strcmp(form,"JERhigh") == 0) {
    ftrijet = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_JERhigh_%s.root",rname),"READ");
  }
  else if (strcmp(form,"JERlow") == 0) {
    ftrijet = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_JERlow_%s.root",rname),"READ");
  }
  else if (strcmp(form,"HERWIG") == 0) {
    ftrijet = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_HERWIG_%s.root",rname),"READ");
  }
  else {
    ftrijet = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/hists_trijet_%s.root",rname),"READ");
  }
  for (int i = 0; i < ana::nTrijetPtBins; i++ ) {
    TH1D* h = (TH1D*)ftrijet->Get(Form("htrijet%i",i));
    mean3[i] = h->GetMean();
    merr3[i] = h->GetMeanError();
  }

  // -----------------------------
  // Load trees ONCE
  // -----------------------------
  
  //gammajet data
  TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/tree_Data.root"));

  TTree * t;
  if (strcmp(form, "nominal") == 0 || strcmp(form, "JERhigh") == 0 || strcmp(form, "JERlow") == 0 || strcmp(form, "HERWIG") == 0) {
    t = (TTree*)f->Get(Form("xjtree_%i",ir));
  }
  else {
    t = (TTree*)f->Get(Form("xjtree_%s_%i",form,ir));
  }

  float pho_pt, jet_pt;
  t->SetBranchAddress("pho_pt",&pho_pt);
  t->SetBranchAddress("jet_pt",&jet_pt);
  t->SetBranchStatus("weight",0);

  // multijet data
  TFile * f3;
  if (strcmp(form,"JERhigh") == 0) {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sPYTHIAJERHigh.root",rname),"READ");
  }
  else if (strcmp(form,"JERlow") == 0) {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sPYTHIAJERLow.root",rname),"READ");
  }
  else if (strcmp(form,"HERWIG") == 0) {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sHERWIG.root",rname),"READ");
  }
  else {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sPYTHIA.root",rname),"READ");
  }
  

  TTree * t3 = (TTree*)f3->Get("ttree");

  float leading_pt, subleading_pt, subleading_phi;
  float subsubleading_pt, subsubleading_phi;

  t3->SetBranchAddress("leadingPT",&leading_pt);
  t3->SetBranchAddress("SLPT",&subleading_pt);
  t3->SetBranchAddress("SLphi",&subleading_phi);
  t3->SetBranchAddress("SSLPT",&subsubleading_pt);
  t3->SetBranchAddress("SSLphi",&subsubleading_phi);

  // -----------------------------
  // Cache events
  // -----------------------------
  struct Result { float pa, pb, chisq; };

  vector<GammaEvent> gamma;
  vector<TrijetEvent> trijet;
  vector<Result> result;

  int nentries = t->GetEntries();
  for (int e = 0; e < nentries; e++) {
    t->GetEntry(e);
    int ipt = ana::findPtBin(pho_pt);
    if (ipt < 0) continue;
    gamma.push_back({jet_pt, jet_pt/pho_pt, ipt});
  }

  int nentries3 = t3->GetEntries();
  for (int e = 0; e < nentries3; e++) {
    t3->GetEntry(e);
    int ipt = ana::findTrijetPtBin(leading_pt);
    if (ipt < 0) continue;

    trijet.push_back({leading_pt, subleading_pt, subleading_phi,
        subsubleading_pt, subsubleading_phi, ipt});
  }

  cout << "Cached events: gamma=" << gamma.size()
    << " trijet=" << trijet.size() << endl;

  // -----------------------------
  // Setup data
  // -----------------------------
  cout << mean3[0] << endl;
  FitData data;
  data.gamma  = &gamma;
  data.trijet = &trijet;
  data.meanp  = meanp;
  data.merrp  = merrp;
  data.mean3  = mean3;
  data.merr3  = merr3;
  data.nTot   = ana::nPtBins + ana::nTrijetPtBins;

  // -----------------------------
  // Minuit setup
  // -----------------------------
  cout << "Running minuit minimzer..." << endl;
  //Chi2Functor fcn(data);

  ROOT::Minuit2::Minuit2Minimizer minimizer;

  ROOT::Math::Functor fcn(
      [&](const double *p) {
      return compute_chi2(p[0], p[1], data);
      }, 2);

  minimizer.SetFunction(fcn);

  minimizer.SetLimitedVariable(0, "pa", 1.01, 0.0001, 0.8, 1.2);
  minimizer.SetLimitedVariable(1, "pb", 0, 1e-6, -0.02, 0.02);

  minimizer.Minimize();
  minimizer.Hesse();

  // -----------------------------
  // Results
  // -----------------------------
  const double *xs   = minimizer.X();
  const double *errs = minimizer.Errors();

  double pa = xs[0];
  double pb = xs[1];

  cout << "\n===== FIT RESULT =====\n";
  cout << "pa = " << pa << " ± " << errs[0] << endl;
  cout << "pb = " << pb << " ± " << errs[1] << endl;

  // -----------------------------
  // Covariance
  // -----------------------------
  double cov[4];
  minimizer.GetCovMatrix(cov);

  double var_a  = cov[0];
  double cov_ab = cov[1];
  double cov_ba = cov[2];
  double var_b  = cov[3];

  cout << "\nCovariance matrix:\n";
  cout << var_a  << "  " << cov_ab << endl;
  cout << cov_ba << "  " << var_b  << endl;

  double corr = cov_ab / sqrt(var_a * var_b);
  cout << "Correlation = " << corr << endl;
  
  // used for both bands
  const int nPoints = 200;

  // ----------------------------
  // Build chi2 band
  // ----------------------------
  vector<double> pa_scan, pb_scan, chi_scan;

  double pa0 = pa;
  double pb0 = pb;

  // scan range (tune if needed)
  double dpa = 5 * errs[0];
  double dpb = 5 * errs[1];

  int Nscan = 100;

  for (int ia = 0; ia < Nscan; ia++) {
    double pa_i = pa0 - dpa + 2*dpa * ia / (Nscan-1);

    for (int ib = 0; ib < Nscan; ib++) {
      double pb_i = pb0 - dpb + 2*dpb * ib / (Nscan-1);

      double chi = compute_chi2(pa_i, pb_i, data);

      pa_scan.push_back(pa_i);
      pb_scan.push_back(pb_i);
      chi_scan.push_back(chi);
    }
  }


  double delta = 2.30;  // 1σ joint
  double chi_min = compute_chi2(pa0, pb0, data);

  cout << "minchisq: " << chi_min << endl;

  double x[nPoints];
  double y_low_env[nPoints];
  double y_high_env[nPoints];

  for (int i = 0; i < nPoints; i++) {
    x[i] = 100.0 * i / (nPoints - 1);

    y_low_env[i]  =  1e9;
    y_high_env[i] = -1e9;
  }

  // loop over scan points
  for (size_t k = 0; k < pa_scan.size(); k++) {

    if (chi_scan[k] > chi_min + delta) continue;

    double pa_i = pa_scan[k];
    double pb_i = pb_scan[k];

    for (int i = 0; i < nPoints; i++) {
      double y = pa_i + pb_i * x[i];

      if (y < y_low_env[i])  y_low_env[i]  = y;
      if (y > y_high_env[i]) y_high_env[i] = y;
    }
  }

  TGraph *g_env = new TGraph(2*nPoints);

  for (int i = 0; i < nPoints; i++) {
    g_env->SetPoint(i, x[i], y_high_env[i]);
    g_env->SetPoint(2*nPoints-1-i, x[i], y_low_env[i]);
  }

  // -----------------------------
  // Build uncertainty band
  // -----------------------------
  //double x[nPoints], y[nPoints];
  double y[nPoints];
  double yup[nPoints], ydn[nPoints];

  for (int i = 0; i < nPoints; i++) {
    double pt = 100.0 * i / (nPoints - 1);

    double f = pa + pb * pt;

    double sigma2 =
        var_a
      + pt*pt * var_b
      + 2 * pt * cov_ab;

    double sigma = sqrt(sigma2);

    x[i]   = pt;
    y[i]   = f;
    yup[i] = f + sigma;
    ydn[i] = f - sigma;
  }
  TGraph *g_band = new TGraph(2*nPoints);
  for (int i = 0; i < nPoints; i++) {
    g_band->SetPoint(i, x[i], yup[i]);
    g_band->SetPoint(2*nPoints-1-i, x[i], ydn[i]);
  }

  // -----------------------------
  // ROOT objects
  // -----------------------------
  TGraph *g_central = new TGraph(nPoints, x, y);
  TGraph *g_up      = new TGraph(nPoints, x, yup);
  TGraph *g_dn      = new TGraph(nPoints, x, ydn);


  TF1 *fbest = new TF1("fbest","[0]+[1]*x",0,100);
  fbest->SetParameters(pa,pb);
  cout << fbest->GetParameter(0) << endl;


  TCanvas *c_compare = new TCanvas("c_compare","bands",800,600);

  // central
  fbest->GetYaxis()->SetRangeUser(0.9,1.15);
  fbest->SetLineColor(kBlack);
  fbest->Draw();

  // covariance band (your current)
  g_band->SetFillColorAlpha(kBlue, 0.3);
  g_band->Draw("F same");

  // envelope band
  g_env->SetFillColorAlpha(kRed, 0.2);
  g_env->Draw("F same");
  
  
  // -----------------------------
  // chisq grid
  // -----------------------------
  
  //int na = 100; float lowa =  0.95, higha = 1.05;
  //int nb = 100; float lowb = -0.002, highb = 0.002;
  //TH2D * hchisq = new TH2D("hchisq",";pa;pb",na,0,na,nb,0,nb);
  //for (int ia = 0; ia < na; ia++) {
  //  float pa = lowa + ia*(higha-lowa)/na;
  //  for (int ib = 0; ib < nb; ib++) {
  //    float pb = lowb + ib*(highb-lowb)/nb;

  //    float chi = compute_chi2(pa, pb, data);
  //    hchisq->Fill(ia+1,ib+1,chi);
  //  }
  //}

  // -----------------------------
  // Save
  // -----------------------------
  TFile *wf = TFile::Open("fit_with_cov.root","RECREATE");

  fbest->Write();
  g_central->Write("central");
  g_up->Write("up");
  g_dn->Write("down");
  g_band->Write("band");
  g_env->Write("envelope");
  //hchisq->Write();

  wf->Close();

  cout << "\nSaved to fit_with_cov.root\n";




}


