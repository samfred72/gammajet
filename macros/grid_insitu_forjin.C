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
inline double compute_pt_sum(
  double pt1, double phi1,
  double pt2, double phi2)
{
  double px = pt1 * cos(phi1) + pt2 * cos(phi2);
  double py = pt1 * sin(phi1) + pt2 * sin(phi2);
  return sqrt(px*px + py*py);
}

// -----------------------------
// Main optimized function
// -----------------------------
void grid_insitu_forjin(const char * form = "nominal", int ir = 2) {


  drawer d;
  if (strcmp(form,"HERWIG") == 0) {
    d = drawer("herwig");
  }

  // -----------------------------
  // Reference (MC) values
  // -----------------------------
  double meanp[ana::nPtBins];
  double merrp[ana::nPtBins];
  double mean3[ana::nTrijetPtBins];
  double merr3[ana::nTrijetPtBins];

  double minjet = ana::jet_calib_pt_cut[ir];
  const char * rname = ana::rnames[ir];

  // Gammajet reference
  for (int i = 0; i < ana::nPtBins; i++ ) {
    double lowcluster = ana::ptBins[i];
    int lowbin = (int)(minjet/lowcluster/0.08 + 1);

    const char * histname;
    histname = Form("h2_%i",i);

    TH1D* h = d.get(histname,1,-1);
    //h->Rebin(4);
    //h->GetXaxis()->SetRange(lowbin+1,25);

    meanp[i] = h->GetMean();
    merrp[i] = h->GetMeanError();
  }

  // -----------------------------
  // Load trees ONCE
  // -----------------------------
  TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/tree_forjin.root"));

  TTree * t;
  t = (TTree*)f->Get(Form("xjtree_forjin_%i",ir));

  double pho_pt, jet_pt;
  double weight;
  t->SetBranchAddress("pho_pt",&pho_pt);
  t->SetBranchAddress("jet_pt",&jet_pt);
  t->SetBranchAddress("weight",&weight);

  // -----------------------------
  // Cache events
  // -----------------------------
  struct GammaEvent { double pt, base; double weight; int bin; };
  struct Result { double pa, pb, chisq; };

  vector<GammaEvent> gamma;
  vector<Result> result;

  int nentries = t->GetEntries();
  for (int e = 0; e < nentries; e++) {
    t->GetEntry(e);
    int ipt = ana::findPtBin(pho_pt);
    if (ipt < 0) continue;
    //gamma.push_back({jet_pt, jet_pt/pho_pt, 1, ipt});
    gamma.push_back({jet_pt, jet_pt/pho_pt, weight, ipt});
    //cout << setprecision(15) << gamma.at(gamma.size()-1).weight << endl;
  }

  cout << "Cached events: gamma=" << gamma.size() << endl;

  // -----------------------------
  // Grid setup
  // -----------------------------
  int na = 100; double lowa =  1.00, higha = 1.05;
  int nb = 100; double lowb = -0.002, highb = 0.002;

  int nTot = ana::nPtBins;

  vector<double> sum(nTot), sum2(nTot);
  vector<double> sumW(nTot), sumW2(nTot);

  double minchisq = FLT_MAX;
  double minpa=0, minpb=0;

  TH1D * hchisq = new TH1D("hchisq",";#Chi^{2};counts",1000,0,10000);
  double bestmeans    [nTot];
  double bestmerrs    [nTot];
  double standardmeans[nTot];
  double standardmerrs[nTot];

  TH2D *hchisq2d = new TH2D("hchisq2d",";pa;pb",na,lowa-0.00001,higha-0.00001,nb,lowb-0.00001,highb-0.00001);
  TH1D *hxj_default[ana::nPtBins];
  TH1D *hxj_scaled[ana::nPtBins];
  TH1D *hxj_corrected[ana::nPtBins];
  for (int i = 0; i < ana::nPtBins; i++) {
    const char * histname = Form("h2_%i",i);
    hxj_default[i] = (TH1D*)d.get(histname,1,-1)->Clone(Form("hxj_default_%i",i));
    hxj_scaled[i] = new TH1D(Form("hxj_scaled_%i",i),";x_J;Normalized Counts",25,0,2);
    hxj_corrected[i] = new TH1D(Form("hxj_corrected_%i",i),";x_J;Normalized Counts",25,0,2);
  }
  //TH2D *hchisq2d = new TH2D("hchisq2d",";pa;pb",na,lowa-0.00001,higha-0.00001,nb,lowb-0.00001,highb-0.00001);
  

  // -----------------------------
  // GRID LOOP
  // -----------------------------
  for (int ia = 0; ia < na; ia++) {
    double pa = lowa + ia*(higha-lowa)/na;
    if (ia % 10 == 0) cout << ia << "/" << na << endl;

    for (int ib = 0; ib < nb; ib++) {
      double pb = lowb + ib*(highb-lowb)/nb;

      fill(sum.begin(), sum.end(), 0.f);
      fill(sum2.begin(), sum2.end(), 0.f);
      fill(sumW.begin(), sumW.end(), 0);
      fill(sumW2.begin(), sumW2.end(), 0);

      // gammajet
      for (auto &ev : gamma) {
        double f = pa + pb*ev.pt;
        //double f = pa*(1 + pb*ev.pt);
        double x = ev.base / f;

        double lowval = ((int)(ana::jet_calib_pt_cut[ir]/ana::ptBins[ev.bin]/0.08 + 1)) * 0.08;
        if (x > lowval  && x <= 2) {
          double w = ev.weight;
          //cout << setprecision(15) << w << endl;

          sum[ev.bin]  += w*x;
          sum2[ev.bin] += w*x*x;
          sumW[ev.bin] += w;
          sumW2[ev.bin] += w*w;
        }
      }


      // chisq
      double chisq = 0;
      double mean;
      double err;
      for (int i = 0; i < nTot; i++) {
        if (sumW[i] == 0) continue;

        double mean = sum[i]/sumW[i];
        double var  = sum2[i]/sumW[i] - mean*mean;
        double Neff = sumW[i]*sumW[i]/sumW2[i];
        double err  = sqrt(var/Neff);

        double refm = (i < ana::nPtBins) ? meanp[i] : mean3[i-ana::nPtBins];
        double refe = (i < ana::nPtBins) ? merrp[i] : merr3[i-ana::nPtBins];

        double diff = 1 - mean/refm;
        //double errt = sqrt(err*err + refe*refe);
        //if (errt < 0.01) errt = 0.01;
        double errt = sqrt((err*err) / (refm*refm) + (mean*mean) * (refe*refe) / (refm*refm*refm*refm));

        chisq += (diff*diff)/(errt*errt);
        //if (i < ana::nPtBins) chisq += (diff*diff)/(errt*errt); // For checking just gammajet
        //if (i >= ana::nPtBins) chisq += (diff*diff)/(errt*errt); // For checking just multijet
      }

      if (chisq < minchisq) {
        minchisq = chisq;
        minpa = pa;
        minpb = pb;
        for (int i = 0; i < nTot; i++) {
          double mean = sum[i]/sumW[i];
          double var  = sum2[i]/sumW[i] - mean*mean;
          double Neff = sumW[i]*sumW[i]/sumW2[i];
          double err = sqrt(var/Neff);
          bestmeans[i] = mean;
          bestmerrs[i] = err;
        }

        //cout << "Better chisq: " << chisq
        //  << "  f=" << pa << " + " << pb << "*pT\n";
      }

      hchisq->Fill(chisq);
      hchisq2d->Fill(pa,pb,chisq);
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
  fill(sumW.begin(), sumW.end(), 0);
  fill(sumW2.begin(), sumW2.end(), 0);

  // gammajet
  int count = 0;
  for (auto &ev : gamma) {
    double x = ev.base;

    double f = minpa + minpb*ev.pt;
    double x_corrected = ev.base / f;
    //cout << ev.pt << " " << x << endl;

    double lowval = ((int)(ana::jet_calib_pt_cut[ir]/ana::ptBins[ev.bin]/0.08 + 1)) * 0.08;
    if (ev.bin == 0 && (x < lowval || x > 2)) cout << "WHAT" << endl;
    if (ev.bin == 0) count++;  
    double w = ev.weight;

      sum[ev.bin] += w*x;
      sum2[ev.bin] += w*x*x;
      sumW[ev.bin] += w;
      //if (ev.bin == 0) cout << setprecision(15) << w << " " << sumW[0] << " " << count*1 << endl;
      sumW2[ev.bin] += w*w;
      //if (ev.bin == 0) cout << sum[0] << " " << sum2[0] << " " << sumW[0] << " " << sumW2[0] << endl;
    //}
    hxj_scaled[ev.bin]->Fill(x,w);
    if (x_corrected > lowval && x < 2) hxj_corrected[ev.bin]->Fill(x_corrected,w);
  }
  for (int i = 0; i < nTot; i++) {
    //if (count[i] == 0) continue;

    double mean = sum[i]/sumW[i];
    double var  = sum2[i]/sumW[i] - mean*mean;
    double Neff = sumW[i]*sumW[i]/sumW2[i];
    count += Neff;
    //cout << setprecision(15) << sumW[i] << " " << sumW2[i] << endl;
    //cout << "Neff for bin " << i << ": " << Neff << endl;
    double err = sqrt(var/Neff);

    cout << "Standard mean for bin " << i << " is " << mean << " / " << ((i < ana::nPtBins) ? meanp[i] : mean3[i-ana::nPtBins]) << endl;
    standardmeans[i] = mean;
    standardmerrs[i] = err;
  }
  
  // Collect gammajet points

  TH1D * hinsitu = new TH1D("hinsitu", ";p_{T} bin;corrected <x_{j,data}>/<x_{j,MC}>", ana::nPtBins, ana::ptBins);
  TH1D * hstandardinsitu = new TH1D("hstandardinsitu", ";p_{T} bin;corrected <x_{j,data}>/<x_{j,MC}>", ana::nPtBins, ana::ptBins);
  for (int i = 0; i < ana::nPtBins; i++) {
    hinsitu->SetBinContent(i+1, 1.0/(bestmeans[i]/meanp[i]));
    hinsitu->SetBinError(i+1,  TMath::Sqrt(bestmerrs[i]*bestmerrs[i] + merrp[i]*merrp[i] ));
    hstandardinsitu->SetBinContent(i+1, 1.0/(standardmeans[i]/meanp[i]));
    double err = TMath::Sqrt(standardmerrs[i]*standardmerrs[i] + merrp[i]*merrp[i] );
    //if (err < 0.01) err = 0.01;
    hstandardinsitu->SetBinError(i+1, err );
  }

  std::sort(result.begin(), result.end(), [](const Result& a, const Result& b) {
      return a.chisq > b.chisq; // keep your ordering
      });

  double minChi = result.back().chisq;   // smallest
  int nKeep; 
  for (int i = 0; i < result.size(); i++) {
    //if (result.at(i).chisq < minChi + 3.53) {
    if (result.at(i).chisq < minChi + 2.30) { // THIS IS THE RIGHT NUMBER FOR 2D
      nKeep = result.size() - i;
      break;
    }
  }
  double maxChi = result.at(result.size() - nKeep).chisq;  // largest

  const int nPoints = 1000;
  double yLow[nPoints];
  std::fill(std::begin(yLow), std::end(yLow), FLT_MAX);
  double yHigh[nPoints] = { 0 };
  double xfunc[nPoints] = { 0 };
  
  TLegend * leg2 = new TLegend(.15,.65,.45,.85);
  for (int i = result.size() - nKeep; i < result.size(); i++) {

    float red = (maxChi - result.at(i).chisq) / (maxChi - minChi);
    int color = TColor::GetColor(red, 0.0, 0.0);

    TF1 *func = new TF1(Form("tmp%i",i), "pol1", 0, 100);
    func->SetParameters(result.at(i).pa, result.at(i).pb);
    
    for (int j = 0; j < nPoints; ++j) {
      double xi = 100.0/nPoints * j;

      double y = func->Eval(xi);
      double ymin = std::min(y, yLow[j]);
      double ymax = std::max(y, yHigh[j]);

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

  const char * wfilename = Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_forjin.root");
  cout << "Writing output to " << wfilename << endl;
  TFile * wf = TFile::Open(wfilename,"RECREATE");
  func->Write();
  fHigh->Write();
  fLow->Write();
  hinsitu->Write();
  hstandardinsitu->Write();
  hchisq->Write();
  hchisq2d->Write();

  for(int i = 0; i < ana::nPtBins; i++) {
    hxj_default[i]->Write(0);
    hxj_scaled[i]->Write(0);
    hxj_corrected[i]->Write(0);
  }

  TTree * wt = new TTree("T","tree");
  wt->Branch("chisq", &minchisq, "chisq/F");
  wt->Fill();
  wt->Write();  
  wf->Close();

}
  
  
