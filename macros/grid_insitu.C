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
void grid_insitu(const char * form = "nominal", int ir = 2, int version = 0) {
  bool dogamma = (version == 0 || version == 1);
  bool domulti = (version == 0 || version == 2);
  bool use_quad = false;

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
  const char * histname;
  for (int i = 0; i < ana::nPtBins; i++ ) {
    float lowcluster = ana::ptBins[i];
    int lowbin = (int)(minjet/lowcluster/0.08 + 1);

    if (strcmp(form, "nominal") == 0) {
      histname = Form("hratio_%i_%i_3_0_0_0",i,ir);
    }
    else if (strcmp(form, "3jet") == 0) {
      histname = Form("hratio_%i_%i_3_0_1_0",i,ir);
    }
    else if (strcmp(form, "bdt") == 0) {
      histname = Form("hratio_%i_%i_3_1_0_0",i,ir);
    }
    else if (strcmp(form, "iso") == 0) {
      histname = Form("hratio_%i_%i_3_2_0_0",i,ir);
    }
    else if (strcmp(form, "JERhigh") == 0) {
      histname = Form("hratio_%i_%i_5_0_0_0",i,ir);
    }
    else if (strcmp(form, "JERlow") == 0) {
      histname = Form("hratio_%i_%i_6_0_0_0",i,ir);
    }
    else if (strcmp(form, "scalehigh") == 0) {
      histname = Form("hratio_%i_%i_0_0_0_0",i,ir);
    }
    else if (strcmp(form, "scalelow") == 0) {
      histname = Form("hratio_%i_%i_1_0_0_0",i,ir);
    }
    else if (strcmp(form, "HERWIG") == 0) {
      histname = Form("hratio_%i_%i_3_0_0_0",i,ir);
    }
    else {
      cout << "Bad variation!" << endl;
      return;
    }

    TH1D* h = d.get(histname,2);

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
  TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/tree_Data.root"));

  TTree * t;
  if (strcmp(form, "nominal") == 0 || strcmp(form, "JERhigh") == 0 || strcmp(form, "JERlow") == 0 || strcmp(form, "scalehigh") == 0 || strcmp(form, "scalelow") == 0 || strcmp(form, "HERWIG") == 0) {
    t = (TTree*)f->Get(Form("xjtree_%i",ir));
  }
  else {
    t = (TTree*)f->Get(Form("xjtree_%s_%i",form,ir));
  }

  float pho_pt, jet_pt;
  t->SetBranchAddress("pho_pt",&pho_pt);
  t->SetBranchAddress("jet_pt",&jet_pt);
  t->SetBranchStatus("weight",0);

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
  struct GammaEvent { float pt, base; int bin; };
  struct TrijetEvent { float lead, sl, slphi, ssl, sslphi; int bin; };
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
  // Grid setup
  // -----------------------------
  int na = (!domulti ? 1000 : 200); float lowa =  0.80, higha = 1.00;
  //int na = (!domulti ? 1000 : 200); float lowa =  0.95, higha = 1.05;
  int nb = (!domulti ? 1 : 100); float lowb = (!domulti? 0 : -0.002), highb = (!domulti? 0 : 0.002);
  int nc = (use_quad ? 100 : 1); float lowc = (use_quad ? -0.002 : 0), highc = (use_quad ? 0.002 : 0);

  int nTot = ana::nPtBins + ana::nTrijetPtBins;

  vector<float> sum(nTot), sum2(nTot);
  vector<int> count(nTot);

  float minchisq = FLT_MAX;
  float minpa=0, minpb=0, minpc=0;

  TH1D * hchisq = new TH1D("hchisq",";#Chi^{2};counts",1000,0,10000);
  float bestmeans    [nTot];
  float bestmerrs    [nTot];
  float standardmeans[nTot];
  float standardmerrs[nTot];

  TH2D *hchisq2d = new TH2D("hchisq2d",";pa;pb",na,lowa-0.00001,higha-0.00001,nb,lowb-0.00001,highb-0.00001);

  TH1D *hxj_MC[ana::nPtBins];
  TH1D *hxj_data[ana::nPtBins];
  TH1D *hxj_corrected[ana::nPtBins];
  for (int i = 0; i < ana::nPtBins; i++) {
    if (strcmp(form, "nominal") == 0) {
      histname = Form("hratio_%i_%i_3_0_0_0",i,ir);
    }
    else if (strcmp(form, "3jet") == 0) {
      histname = Form("hratio_%i_%i_3_0_1_0",i,ir);
    }
    else if (strcmp(form, "bdt") == 0) {
      histname = Form("hratio_%i_%i_3_1_0_0",i,ir);
    }
    else if (strcmp(form, "iso") == 0) {
      histname = Form("hratio_%i_%i_3_2_0_0",i,ir);
    }
    else if (strcmp(form, "JERhigh") == 0) {
      histname = Form("hratio_%i_%i_5_0_0_0",i,ir);
    }
    else if (strcmp(form, "JERlow") == 0) {
      histname = Form("hratio_%i_%i_6_0_0_0",i,ir);
    }
    else if (strcmp(form, "scalehigh") == 0) {
      histname = Form("hratio_%i_%i_0_0_0_0",i,ir);
    }
    else if (strcmp(form, "scalelow") == 0) {
      histname = Form("hratio_%i_%i_1_0_0_0",i,ir);
    }
    else if (strcmp(form, "HERWIG") == 0) {
      histname = Form("hratio_%i_%i_3_0_0_0",i,ir);
    }
    hxj_MC[i] = (TH1D*)d.get(histname,2,-1)->Clone(Form("hxj_MC_%i",i));
    hxj_data[i] = (TH1D*)d.get(histname,0,-1)->Clone(Form("hxj_data_%i",i));
    hxj_corrected[i] = new TH1D(Form("hxj_corrected_%i",i),";x_J;Normalized Counts",25,0,2);
  }
  // -----------------------------
  // GRID LOOP
  // -----------------------------
  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;
    if (ia % 10 == 0) cout << ia << "/" << na << endl;

    for (int ib = 0; ib < nb; ib++) {
      float pb = lowb + ib*(highb-lowb)/nb;
      for (int ic = 0; ic < nc; ic++) {
        float pc = lowc + ic*(highc-lowc)/nc;

        fill(sum.begin(), sum.end(), 0.f);
        fill(sum2.begin(), sum2.end(), 0.f);
        fill(count.begin(), count.end(), 0);

        // gammajet
        for (auto &ev : gamma) {
          float f = pa + pb*ev.pt + pc*ev.pt*ev.pt;
          //float f = pa*(1 + pb*ev.pt);
          float x = ev.base / f;

          float lowval = (int)(ana::jet_calib_pt_cut[1]/ana::ptBins[ev.bin]/0.08 + 1) * 0.08;
          if (x > lowval  && x <= 2) {
            sum[ev.bin] += x;
            sum2[ev.bin] += x*x;
            count[ev.bin]++;
          }
        }

        // trijet
        for (auto &ev : trijet) {
          float f1 = pa + pb*ev.sl  + pc*ev.sl*ev.sl;
          float f2 = pa + pb*ev.ssl + pc*ev.ssl*ev.ssl;

          float pt1 = ev.sl / f1;
          float pt2 = ev.ssl / f2;

          float multi = compute_pt_sum(pt1, ev.slphi, pt2, ev.sslphi);

          float fL = pa + pb*ev.lead + pc*ev.lead*ev.lead;
          float x = ev.lead / multi / fL;

          int b = ana::nPtBins + ev.bin;
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

          float refm = (i < ana::nPtBins) ? meanp[i] : mean3[i-ana::nPtBins];
          float refe = (i < ana::nPtBins) ? merrp[i] : merr3[i-ana::nPtBins];

          float diff = 1 - mean/refm;
          float errt = sqrt((err*err) / (refm*refm) + (mean*mean) * (refe*refe) / (refm*refm*refm*refm));

          if (dogamma && i < ana::nPtBins)  chisq += (diff*diff)/(errt*errt); // For checking gammajet
          if (domulti && i >= ana::nPtBins) chisq += (diff*diff)/(errt*errt); // For checking multijet

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
        hchisq2d->Fill(pa,pb,chisq);
        result.push_back({pa,pb,chisq});
      }
    }
  }

  cout << "\nFINAL RESULT\n";
  cout << "chi2 = " << minchisq << endl;
  cout << "f(pT) = " << minpa;
  if (domulti) cout << " + " << minpb << "*pT";
  if (use_quad) cout << " + " << minpc << "*pT^2";
  cout << endl;




  //cout << "Final pass on standard func" << endl;
  fill(sum.begin(), sum.end(), 0.f);
  fill(sum2.begin(), sum2.end(), 0.f);
  fill(count.begin(), count.end(), 0);

  // gammajet
  for (auto &ev : gamma) {
    float x = ev.base;
    double f = minpa + minpb*ev.pt + minpc*ev.pt*ev.pt;
    double x_corrected = ev.base / f;

    float lowval = (int)(ana::jet_calib_pt_cut[2]/ana::ptBins[ev.bin]/0.08 + 1) * 0.08;
    if (x > lowval  && x <= 2) {
      sum[ev.bin] += x;
      sum2[ev.bin] += x*x;
      count[ev.bin]++;
    }
    // Specifically used here for the closure test
    if (x_corrected > lowval && x < 2) hxj_corrected[ev.bin]->Fill(x_corrected);
  }

  // trijet
  for (auto &ev : trijet) {
    float pt1 = ev.sl;
    float pt2 = ev.ssl;
    float multi = compute_pt_sum(pt1, ev.slphi, pt2, ev.sslphi);

    float x = ev.lead / multi;

    int b = ana::nPtBins + ev.bin;
    if (x > 0.4 && x < 2.65) {
      sum[b] += x;
      sum2[b] += x*x;
      count[b]++;
    }
  }
  for (int i = 0; i < nTot; i++) {
    float mean = sum[i]/count[i];
    float var  = sum2[i]/count[i] - mean*mean;
    float err  = sqrt(var/count[i]);

    if (i == 0) cout << "Gammajet bins:" << endl;
    if (i == ana::nPtBins) cout << "Multijet bins:" << endl;
    cout << "Default means for bin " << i << " is " << mean << " / " << ((i < ana::nPtBins) ? meanp[i] : mean3[i-ana::nPtBins]) << endl;
    standardmeans[i] = mean;
    standardmerrs[i] = err;
  }

  // Collect gammajet points

  TH1D * hinsitu = new TH1D("hinsitu", ";p_{T} bin;corrected <x_{j,data}>/<x_{j,MC}>", ana::nPtBins, ana::ptBins);
  TH1D * hstandardinsitu = new TH1D("hstandardinsitu", ";p_{T} bin;corrected <x_{j,data}>/<x_{j,MC}>", ana::nPtBins, ana::ptBins);
  for (int i = 0; i < ana::nPtBins; i++) {
    hinsitu->SetBinContent(i+1, bestmeans[i]/meanp[i]);
    hinsitu->SetBinError(i+1,  TMath::Sqrt(bestmerrs[i]*bestmerrs[i] + merrp[i]*merrp[i] ));
    hstandardinsitu->SetBinContent(i+1, standardmeans[i]/meanp[i]);
    float err = TMath::Sqrt(standardmerrs[i]*standardmerrs[i] + merrp[i]*merrp[i] );
    //if (err < 0.01) err = 0.01;
    hstandardinsitu->SetBinError(i+1, err );
  }
  // Collect trijet points
  TH1D * hinsitu3 = new TH1D("hinsitu3", ";p_{T} [GeV];corrected <x_{j,data}>/<x_{j,MC}>", ana::nTrijetPtBins, ana::trijetPtBins);
  TH1D * hstandardinsitu3 = new TH1D("hstandardinsitu3", ";p_{T} [GeV];corrected <x_{j,data}>/<x_{j,MC}>", ana::nTrijetPtBins, ana::trijetPtBins);
  for (int i = 0; i < ana::nTrijetPtBins; i++) {
    hinsitu3->SetBinContent(i+1, bestmeans[ana::nPtBins+i]/mean3[i]);
    hinsitu3->SetBinError(i+1,  TMath::Sqrt(bestmerrs[ana::nPtBins+i]*bestmerrs[ana::nPtBins+i] + merr3[i]*merr3[i] ));
    hstandardinsitu3->SetBinContent(i+1, standardmeans[ana::nPtBins+i]/mean3[i]);
    float err = TMath::Sqrt(standardmerrs[ana::nPtBins+i]*standardmerrs[ana::nPtBins+i] + merr3[i]*merr3[i] );
    //if (err < 0.01) err = 0.01;
    hstandardinsitu3->SetBinError(i+1, err);
  }

  std::sort(result.begin(), result.end(), [](const Result& a, const Result& b) {
      return a.chisq > b.chisq; // sort highest to lowest 
      });

  float minChi = result.back().chisq;   // smallest
  int nKeep; 
  float chisq_thresh = 2.30;
  if (use_quad) chisq_thresh = 3.52;
  if (!domulti) chisq_thresh = 1.0;
  for (int i = 0; i < result.size(); i++) {
    if (result.at(i).chisq < minChi + chisq_thresh) {
      nKeep = result.size() - i;
      cout << "Found nkeep:" << nKeep << endl;
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

  if (dogamma && !domulti) cout << minpa << "+" << yHigh[0] - minpa << "/-" << minpa - yLow[0] << endl;

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

  const char * gammatext = "";
  if (dogamma && !domulti) gammatext = "_gammajet";
  const char * wfilename = Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_%s_jet_%s%s.root",form,rname,gammatext);
  cout << "Writing output to " << wfilename << endl;
  TFile * wf = TFile::Open(wfilename,"RECREATE");
  func->Write();
  fHigh->Write();
  fLow->Write();
  hinsitu->Write();
  hstandardinsitu->Write();
  hinsitu3->Write();
  hstandardinsitu3->Write();
  hchisq->Write();
  hchisq2d->Write();

  for(int i = 0; i < ana::nPtBins; i++) {
    hxj_MC[i]->Scale(1.0/hxj_MC[i]->Integral());
    hxj_data[i]->Scale(1.0/hxj_data[i]->Integral());
    hxj_corrected[i]->Scale(1.0/hxj_corrected[i]->Integral());
    hxj_MC[i]->Write(0);
    hxj_data[i]->Write(0);
    hxj_corrected[i]->Write(0);
  }

  TTree * wt = new TTree("T","tree");
  wt->Branch("chisq", &minchisq, "chisq/F");
  wt->Fill();
  wt->Write();  
  wf->Close();

}

  
