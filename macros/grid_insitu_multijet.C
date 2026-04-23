#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

#include <vector>
#include <cmath>
#include <cfloat>
#include <iostream>

using namespace std;

// -----------------------------
// Helper: fast 2-vector pT sum
// -----------------------------
inline float compute_pt_sum(
  float pt1, float phi1,
  float pt2, float phi2)
{
  float px = pt1 * cos(phi1) + pt2 * cos(phi2);
  float py = pt1 * sin(phi1) + pt2 * sin(phi2);
  return sqrt(px*px + py*py);
}

// -----------------------------
// Main optimized function
// -----------------------------
void grid_insitu_multijet(const char * form = "nominal", int ir = 1) {


  drawer d;
  if (strcmp(form,"HERWIG") == 0) {
    d = drawer("herwig");
  }

  // -----------------------------
  // Reference (MC) values
  // -----------------------------
  float mean3[ana::nTrijetPtBins];
  float merr3[ana::nTrijetPtBins];

  const char * rname = ana::rnames[ir];

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
  TFile * f3;
  if (strcmp(form,"JERhigh") == 0) {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sJERHigh.root",rname),"READ");
  }
  else if (strcmp(form,"JERlow") == 0) {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sJERLow.root",rname),"READ");
  }
  else if (strcmp(form,"HERWIG") == 0) {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%sHERWIG.root",rname),"READ");
  }
  else {
    f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/trees/SAMfile_%s.root",rname),"READ");
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
  struct TrijetEvent { float lead, sl, slphi, ssl, sslphi; int bin; };
  struct Result { float pa, pb, chisq; };

  vector<TrijetEvent> trijet;
  vector<Result> result;

  int nentries3 = t3->GetEntries();
  for (int e = 0; e < nentries3; e++) {
    t3->GetEntry(e);
    int ipt = ana::findTrijetPtBin(leading_pt);
    if (ipt < 0) continue;

    trijet.push_back({leading_pt, subleading_pt, subleading_phi,
        subsubleading_pt, subsubleading_phi, ipt});
  }

  cout << " trijet=" << trijet.size() << endl;

  // -----------------------------
  // Grid setup
  // -----------------------------
  int na = 100; float lowa =  0.95, higha = 1.05;
  int nb = 100; float lowb = -0.002, highb = 0.002;
  //int na = 100; float lowa =  0.5, higha = 1.5;
  //int nb = 100; float lowb = -0.02, highb = 0.02;

  int nTot = ana::nTrijetPtBins;

  vector<float> sum(nTot), sum2(nTot);
  vector<int> count(nTot);

  float minchisq = FLT_MAX;
  float minpa=0, minpb=0;

  TH1D * hchisq = new TH1D("hchisq",";#Chi^{2};counts",1000,0,10000);
  float bestmeans    [nTot];
  float bestmerrs    [nTot];
  float standardmeans[nTot];
  float standardmerrs[nTot];

  // -----------------------------
  // GRID LOOP
  // -----------------------------
  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;
    if (ia % 10 == 0) cout << ia << "/" << na << endl;

    for (int ib = 0; ib < nb; ib++) {
      float pb = lowb + ib*(highb-lowb)/nb;

      fill(sum.begin(), sum.end(), 0.f);
      fill(sum2.begin(), sum2.end(), 0.f);
      fill(count.begin(), count.end(), 0);

      // trijet
      for (auto &ev : trijet) {
        float f1 = pa + pb*ev.sl;
        float f2 = pa + pb*ev.ssl;

        float pt1 = ev.sl / f1;
        float pt2 = ev.ssl / f2;

        float multi = compute_pt_sum(pt1, ev.slphi, pt2, ev.sslphi);
        //float multi = compute_pt_sum(ev.sl, ev.slphi, ev.ssl, ev.sslphi);

        float fL = pa + pb*ev.lead;
        float x = ev.lead / multi / fL;

        int b = ev.bin;
        if ( x > 0.4 && x < 2.65) {
          sum[b] += x;
          sum2[b] += x*x;
          count[b]++;
        }
      }

      // chisq
      float chisq = 0;
      float mean;
      float err;
      for (int i = 0; i < nTot; i++) {
        if (count[i] == 0) continue;

        float mean = sum[i]/count[i];
        float var  = sum2[i]/count[i] - mean*mean;
        float err  = sqrt(var/count[i]);

        float refm = mean3[i];
        float refe = merr3[i];

        float diff = 1 - mean/refm;
        float errt = sqrt(err*err + refe*refe);

        chisq += (diff*diff)/(errt*errt);
      }

      if (chisq < minchisq) {
        minchisq = chisq;
        minpa = pa;
        minpb = pb;
        for (int i = 0; i < nTot; i++) {
          float mean = sum[i]/count[i];
          float var  = sum2[i]/count[i] - mean*mean;
          float err  = sqrt(var/count[i]);
          bestmeans[i] = mean;
          bestmerrs[i] = err;
        }

      }

      hchisq->Fill(chisq);
      result.push_back({pa,pb,chisq});
    }
  }

  cout << "\nFINAL RESULT\n";
  cout << "chi2 = " << minchisq << endl;
  cout << "f(pT) = " << minpa << " + "
       << minpb << "*pT\n";
  
  //cout << "Final pass on standard func" << endl;
  fill(sum.begin(), sum.end(), 0.f);
  fill(sum2.begin(), sum2.end(), 0.f);
  fill(count.begin(), count.end(), 0);

  // trijet
  for (auto &ev : trijet) {
    float pt1 = ev.sl;
    float pt2 = ev.ssl;

    float multi = compute_pt_sum(pt1, ev.slphi, pt2, ev.sslphi);

    float x = ev.lead / multi;

    int b = ev.bin;
    if (x > 0.4 && x < 2.65) {
      sum[b] += x;
      sum2[b] += x*x;
      count[b]++;
    }
  }
  for (int i = 0; i < nTot; i++) {
    if (count[i] == 0) continue;

    float mean = sum[i]/count[i];
    float var  = sum2[i]/count[i] - mean*mean;
    float err  = sqrt(var/count[i]);

    cout << "Standard mean for bin " << i << " is " << mean << " / " << (mean3[i]) << endl;
    standardmeans[i] = mean;
    standardmerrs[i] = err;
  }
  
  // Collect trijet points
  TH1D * hinsitu3 = new TH1D("hinsitu3", ";p_{T} [GeV];corrected <x_{j,data}>/<x_{j,MC}>", ana::nTrijetPtBins, ana::trijetPtBins);
  TH1D * hstandardinsitu3 = new TH1D("hstandardinsitu3", ";p_{T} [GeV];corrected <x_{j,data}>/<x_{j,MC}>", ana::nTrijetPtBins, ana::trijetPtBins);
  for (int i = 0; i < ana::nTrijetPtBins; i++) {
    cout << bestmeans[i]/mean3[i] << " " << bestmerrs[i]*bestmerrs[i] + merr3[i]*merr3[i] << " " << standardmeans[i]/mean3[i] << " " << standardmerrs[i]*standardmerrs[i] + merr3[i]*merr3[i] << endl;
    hinsitu3->SetBinContent(i+1, bestmeans[i]/mean3[i]);
    hinsitu3->SetBinError(i+1,  TMath::Sqrt(bestmerrs[i]*bestmerrs[i] + merr3[i]*merr3[i] ));
    hstandardinsitu3->SetBinContent(i+1, standardmeans[i]/mean3[i]);
    hstandardinsitu3->SetBinError(i+1,  TMath::Sqrt(standardmerrs[i]*standardmerrs[i] + merr3[i]*merr3[i] ));
  }

  std::sort(result.begin(), result.end(), [](const Result& a, const Result& b) {
      return a.chisq > b.chisq; // keep your ordering
      });

  float minChi = result.back().chisq;   // smallest
  int nKeep; 
  for (int i = 0; i < result.size(); i++) {
    if (result.at(i).chisq < minChi + 3.53) {
      nKeep = result.size() - i;
      break;
    }
  }
  float maxChi = result.at(result.size() - nKeep).chisq;  // largest

  const int nPoints = 1000;
  float yLow[nPoints];
  std::fill(std::begin(yLow), std::end(yLow), FLT_MAX);
  float yHigh[nPoints] = { 0 };
  float xfunc[nPoints] = { 0 };
  
  TLegend * leg2 = new TLegend(.15,.65,.45,.85);
  for (int i = result.size() - nKeep; i < result.size(); i++) {

    float red = (maxChi - result.at(i).chisq) / (maxChi - minChi);
    int color = TColor::GetColor(red, 0.0, 0.0);

    TF1 *func = new TF1(Form("tmp%i",i), "pol1", 0, 100);
    func->SetParameters(result.at(i).pa, result.at(i).pb);
    
    for (int j = 0; j < nPoints; ++j) {
      double xi = 100.0/nPoints * j;

      float y = func->Eval(xi);
      float ymin = std::min(y, yLow[j]);
      float ymax = std::max(y, yHigh[j]);

      yLow[j] = ymin;
      yHigh[j] = ymax;
      xfunc[j] = xi;
    }
    delete func;
  }
  TF1 * func = new TF1("fbest","pol1",0,100);
  func->SetParameters(minpa,minpb);
  
  //cout << "Creating min and max functions" << endl;
  TGraph *gHigh = new TGraph(nPoints, xfunc, yHigh);
  TSpline3 *sHigh = new TSpline3("sHigh", gHigh);
  TF1 * fHigh = new TF1("fHigh",
      [sHigh](double *xx, double *) {
      return sHigh->Eval(xx[0]);
      }, 0, 100, 0);

  
  TGraph *gLow = new TGraph(nPoints, xfunc, yLow);
  TSpline3 *sLow = new TSpline3("sLow", gLow);
  TF1 *fLow = new TF1("fLow",
      [sLow](double *xx, double *) {
      return sLow->Eval(xx[0]);
      }, 0, 100, 0);

  const char * wfilename = Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_multijet_%s_%s.root",form,rname);
  cout << "Writing output to " << wfilename << endl;
  TFile * wf = TFile::Open(wfilename,"RECREATE");
  func->Write();
  fHigh->Write();
  fLow->Write();
  hinsitu3->Write();
  hstandardinsitu3->Write();

  TTree * wt = new TTree("T","tree");
  wt->Branch("chisq", &minchisq, "chisq/F");
  wt->Fill();
  wt->Write();  
  wf->Close();

}
  
  
