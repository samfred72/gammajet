
#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_global_systematics(int ir=2, int type=3) {
  drawer d;
  
  bool pythia_default = 0;
  bool sym = 0;
  const char * filetag;

  if (type == 0) {
    pythia_default = 1;
    sym = 0;
    filetag = "_pythiadefault_asym";
  }
  if (type == 1) {
    pythia_default = 1;
    sym = 1;
    filetag = "_pythiadefault_sym_JERasym";
  }
  if (type == 2) {
    pythia_default = 0;
    sym = 0;
    filetag = "_average_asym";
  }
  if (type == 3) {
    pythia_default = 0;
    sym = 1;
    filetag = "_average_sym_JERasym";
  }


  const char * rname = ana::rnames[ir];
  rname = "R04";
  float model_vals[ana::nJetR] = {0.9767, 0.9746, 0.9990, 1.0139, 1.0327, 1.0420, 1.0513};
  float model_sys[ana::nJetR];
  for (int i = 0; i < ana::nJetR; i++) {
    float val0 = fabs(1-model_vals[2])/2.0;
    float val = fabs(1-model_vals[i])/2.0;
    model_sys[i] = TMath::Sqrt(fabs(val0*val0 - val*val));
    cout << model_sys[i] << endl;
  }
  TFile * fn  = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_nominal_%s.root"  ,rname));
  TFile * fb  = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_bdt_%s.root"      ,rname));
  TFile * fi  = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_iso_%s.root"      ,rname));
  TFile * f3  = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_3jet_%s.root"     ,rname));
  TFile * fh  = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_JERhigh_%s.root"  ,rname));
  TFile * fl  = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_JERlow_%s.root"   ,rname));
  TFile * fsh = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_scalehigh_%s.root",rname));
  TFile * fsl = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_scalelow_%s.root" ,rname));
  TFile * fH  = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_HERWIG_%s.root"   ,rname));
  rname = ana::rnames[ir];
  const int nfiles = 9;
  TFile * f[nfiles] = {fn,fb,fi,f3,fh,fl,fsh,fsl,fH};
  TF1 * func[nfiles];
  int colors[nfiles] = {kBlack, kBlue, kRed, kGreen, kOrange, kMagenta, kMagenta+2, kOrange+2, kCyan};
  string info[nfiles] = {"Nominal","Narrow BDT", "Narrow Isolation", "Third Jet Cut", "High JER Smearing", "Low JER Smearing", "High EM Scale", "Low EM Scale", "Model Dependence"};
    

  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  for (int i = 0; i < nfiles; i++) {
    func[i] = (TF1*)f[i]->Get("fbest");
    func[i]->SetName(Form("f%i",i));
    func[i]->SetLineColor(colors[i]);
    func[i]->SetLineWidth(2);
    if (i == 0) {
      func[i]->SetTitle("");
      func[i]->GetYaxis()->SetRangeUser(-0.1,0.2);
      func[i]->GetYaxis()->SetTitle("Data-to-MC JES correction");
      func[i]->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    }
  }
 
  // This is the average of the pythia and herwig results. Used as a potential actual final result
  TF1 * faverage = new TF1(Form("faverage"), 
      [=](double *x, double *) {
      double v1 = func[0]->Eval(x[0]);
      double v2 = func[nfiles-1]->Eval(x[0]); 
      return (v1+v2)/2.0;
      },0,100,0);
  // This is their difference divided by 2 (used for the Model systematic)
  TF1 * diff = new TF1(Form("diff"), 
      [=](double *x, double *) {
      double v1 = func[0]->Eval(x[0]);
      double v2 = func[nfiles-1]->Eval(x[0]);
      return (v2-v1)/2.0;
      },0,100,0);
  

  TLegend * lbottom = new TLegend(0.5754,0.575,0.8954,0.865);
  lbottom->SetFillStyle(0);
  TF1 * ratios[nfiles-1];
  TF1 * ratios_sym[nfiles-1];
  for (int i = 1; i < nfiles; i++) {
    TF1 * ratio = new TF1(Form("r%i",i), // The ratio is relative to the pythia result, because we assume the systematics to be uncorrelated to the model
        [=](double *x, double *) {
        double denom = func[0]->Eval(x[0]);
        //if (denom == 0) return 0.0;
        return func[i]->Eval(x[0]) / denom - 1.0;
      },0,100,0);
    ratios[i-1] = (i != nfiles - 1 || pythia_default ? ratio : diff); // For the model "ratio", we don't take a ratio, we take the diff function from before
    ratios[i-1]->SetLineColor(colors[i]);
    ratios[i-1]->SetLineWidth(2);
    if (i == 1) {
      ratios[i-1]->SetTitle("");
      ratios[i-1]->GetYaxis()->SetRangeUser(-0.1,0.2);
      ratios[i-1]->GetXaxis()->SetRangeUser(10,100);
      //ratios[i-1]->GetXaxis()->SetRangeUser(10,60);
      ratios[i-1]->GetYaxis()->SetTitle("Relative Uncertainty");
      ratios[i-1]->GetXaxis()->SetTitle("jet p_{T} [GeV]");
      ratios[i-1]->GetXaxis()->SetTitleOffset(1.3);
      ratios[i-1]->Draw();
    }
    else ratios[i-1]->Draw("same");

    ratios_sym[i-1] =new TF1(Form("r_sym%i",i),
        [=](double*x, double *) {
        return -ratios[i-1]->Eval(x[0]);
        }, 0, 100, 0);
    ratios_sym[i-1]->SetLineColor(colors[i]);
    ratios_sym[i-1]->SetLineWidth(2);
    ratios_sym[i-1]->SetLineStyle(7);
    
    bool isJER = (i-1 == 3 || i-1 == 4 || i-1 == 5 || i-1 == 6);
    bool isHERWIG = (i-1 == 7);
    bool forceSym = (type == 2 || type == 3) && isHERWIG;
    if ((sym && !isJER) || forceSym) ratios_sym[i-1]->Draw("same");

    lbottom->AddEntry(ratios[i-1],info[i].c_str(),"L");
  }

  TLine * l = new TLine(0,0,1,1);
  l->SetLineStyle(7);
  l->SetLineColor(kGray+1);
  l->SetLineWidth(2);
  lbottom->AddEntry(l, "Symmeterized Systematic","L");

  // The systematic uncertanty with nothing symmeterized except the model
  TF1 * systematic_up = new TF1("fsystematic_up",
      [=](double *x, double *) {
      double val = 0;
      for (int i = 0; i < nfiles - 1; i++ ) { 
        double r = ratios[i]->Eval(x[0]);
        bool isJER = (i == 3 || i == 4 || i == 5 || i == 6);
        bool isHERWIG = (i == 7);
        bool forceSym = (type == 2 || type == 3) && isHERWIG;

        if ((sym && !isJER) || forceSym ) {
          // symmetrize everything except JER
          val += r*r;
        }
        else {
          // asymmetric treatment (original behavior)
          if (r > 0) val += r*r;
        }
      }
      val += model_sys[ir]*model_sys[ir];///4.0; // From the photonjet
      return TMath::Sqrt(val);
      },0,100,0);
  // The systematic uncertanty with nothing symmeterized except the model
  TF1 * systematic_down = new TF1("fsystematic_down",
      [=](double *x, double *) {
      double val = 0;
      for (int i = 0; i < nfiles - 1; i++ ) { 
        double r = ratios[i]->Eval(x[0]);
        bool isJER = (i == 3 || i == 4 || i == 5 || i == 6);
        bool isHERWIG = (i == 7);
        bool forceSym = (type == 2 || type == 3) && isHERWIG;

        if ((sym && !isJER) || forceSym ) {
        // symmetrize everything except JER
          val += r*r;
        }
        else {
          // asymmetric treatment (original behavior)
          if (r < 0) val += r*r;
        }
      }
      
      val += model_sys[ir]*model_sys[ir];///4.0; // From the photonjet
      return -TMath::Sqrt(val);
      },0,100,0);

  systematic_up->SetLineColor(kBlack);
  systematic_up->SetLineWidth(3);
  systematic_up->SetLineStyle(9);
  systematic_down->SetLineColor(kBlack);
  systematic_down->SetLineWidth(3);
  systematic_down->SetLineStyle(9);
  
  systematic_up->Draw("same");
  systematic_down->Draw("same");
  
  d.drawAll({"p + p #sqrt{s} = 200 GeV"},{Form("|#eta_{jet}| < %1.1f", 1.1 - ana::JetRs[ir]),Form("Jet R=%1.1f",ana::JetRs[ir])},.17,.83,25, 700);
  
  lbottom->AddEntry(systematic_up,"Total Systematic","L");
  lbottom->SetLineWidth(0);
  lbottom->Draw();
  
  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/note/global_systematic_%s%s.pdf",rname,filetag));
  
  TCanvas * c2 = new TCanvas("c2","",700,700);
  // Drawing fits
  TF1 * flb = (TF1*)f[0]->Get("fbest");
  TF1 * flh = (TF1*)f[0]->Get("fHigh");
  TF1 * fll = (TF1*)f[0]->Get("fLow");
  TTree * tl = (TTree*)f[0]->Get("T");
  float chisql;
  tl->SetBranchAddress("chisq",&chisql);
  tl->GetEntry(0);

  float xmin = 10;
  float xmax = 100;
  //float xmax = 60;
  float ymin = 0.85;
  float ymax = 1.10;
  TH1F* frame2 = c2->DrawFrame(xmin, ymin, xmax, ymax);
  frame2->GetXaxis()->SetTitle("jet p_{T} [GeV]");
  frame2->GetYaxis()->SetTitle("Data-to-MC JES correction");
  frame2->GetXaxis()->SetTitleOffset(1.3);
  c2->GetFrame()->SetLineWidth(0);
  c2->Modified();
  c2->Update();
  gStyle->SetOptStat(0);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  
  
  int nPoints = 1000;
  double xMin = flb->GetXmin();
  double xMax = flb->GetXmax();

  // The statistical error band
  TGraphAsymmErrors* bandl = new TGraphAsymmErrors(nPoints);

  for (int i = 0; i < nPoints; ++i) {
    double x = xMin + (xMax - xMin) * i / (nPoints - 1);

    double y_best = flb->Eval(x) + (!pythia_default ? diff->Eval(x) : 0);
    double y_low  = fll->Eval(x) + (!pythia_default ? diff->Eval(x) : 0);
    double y_high = flh->Eval(x) + (!pythia_default ? diff->Eval(x) : 0);

    double err_low  = y_best - y_low;
    double err_high = y_high - y_best;

    bandl->SetPoint(i, x, y_best);
    bandl->SetPointError(i, 0, 0, err_low, err_high);
  }
  
  // The systematic error band
  TGraphAsymmErrors* bands = new TGraphAsymmErrors(nPoints);

  for (int i = 0; i < nPoints; ++i) {
    double x = xMin + (xMax - xMin) * i / (nPoints - 1);

    double y_best = flb->Eval(x) + (!pythia_default ? ratios[nfiles-2]->Eval(x) : 0);
    double y_low  = y_best + systematic_down->Eval(x);
    double y_high = y_best + systematic_up->Eval(x);

    double err_low  = y_best - y_low;
    double err_high = y_high - y_best;

    bands->SetPoint(i, x, y_best);
    bands->SetPointError(i, 0, 0, err_low, err_high);
  }
  
  bandl->SetFillColorAlpha(kGreen,0.6);
  //bandl->SetFillStyle(3001);
  bandl->SetLineColor(kGreen);
  bandl->SetLineWidth(1);
  bandl->Draw("L3 same");
  
  bands->SetFillColorAlpha(kBlue,0.2);
  //bands->SetFillStyle(3001);
  bands->SetLineColor(kBlue);
  bands->SetLineWidth(1);
  bands->Draw("L3 same");
 
  if (pythia_default) {
    flb->SetLineColor(kBlack);
    flb->Draw("same");
  }
  else {
    faverage->SetLineColor(kBlack);
    faverage->Draw("same");
  }

  TF1 * errhigh = new TF1(Form("errhigh"), 
      [=](double *x, double *) {
      double e1 = systematic_up->Eval(x[0]);
      double e2 = flh->Eval(x[0])-flb->Eval(x[0]);
      return TMath::Sqrt(e1*e1+e2*e2);
      },0,100,0);
  
  TF1 * errlow = new TF1(Form("errlow"), 
      [=](double *x, double *) {
      double e1 = systematic_down->Eval(x[0]);
      double e2 = fll->Eval(x[0])-flb->Eval(x[0]);
      return TMath::Sqrt(e1*e1+e2*e2);
      },0,100,0);
  
  TF1 * linehigh = new TF1(Form("fhigh"), 
      [=](double *x, double *) {
      if (pythia_default) {
        return flb->Eval(x[0])+errhigh->Eval(x[0]);
      }
      else{
        return faverage->Eval(x[0])+errhigh->Eval(x[0]);
      }
      },0,100,0);
  
  TF1 * linelow = new TF1(Form("flow"), 
      [=](double *x, double *) {
      return faverage->Eval(x[0])-errlow->Eval(x[0]);
      },0,100,0);
  
  linehigh->SetLineWidth(2);
  linehigh->SetLineStyle(9);
  linehigh->SetLineColor(kRed);
  linelow->SetLineWidth(2);
  linelow->SetLineStyle(9);
  linelow->SetLineColor(kRed);
  linehigh->Draw("same");
  linelow->Draw("same");

  TLegend * leg2 = new TLegend(0.2,0.15,0.6,0.3);
  leg2->AddEntry(bandl, "statistical uncertainty","F");
  //leg2->AddEntry(bandl, Form("#splitline{%.3f + %.1e*pT}{#chi^{2} = %.2f}", flb->GetParameter(0), flb->GetParameter(1), chisql));
  leg2->AddEntry(bands, Form("systematic uncertainty"),"F");
  leg2->AddEntry(linehigh, Form("total uncertainty"),"L");
  
  leg2->SetLineWidth(0);
  leg2->SetFillStyle(0);
  leg2->Draw();
  
  d.drawAll({"p + p #sqrt{s} = 200 GeV"},{Form("|#eta_{jet}| < %1.1f",1.1 - ana::JetRs[ir]), Form("Jet R=%1.1f",ana::JetRs[ir])},.2,.83,25, 700);
  //d.drawAll({"p + p #sqrt{s} = 200 GeV"},{"|vz| < 60 cm", Form("|#eta| < %1.1f",1.1 - ana::JetRs[ir]), Form("#it{in situ} correction: %.3f + %.1e*pT", (flb->GetParameter(0)+func[nfiles-1]->GetParameter(0))/2.0, (flb->GetParameter(1)+func[nfiles-1]->GetParameter(1))/2.0),Form("Jet R=%1.1f",ana::JetRs[ir])},.2,.83,20, 700);

  c2->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/note/global_insitu_%s%s.pdf",rname,filetag));

  TFile * ff = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/funcs_insitu_%s.root",rname),"RECREATE");
  faverage->SetName("finsitu");
  faverage->Write();
  linehigh->Write();
  linelow->Write();

}
