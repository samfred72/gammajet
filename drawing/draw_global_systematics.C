
#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_global_systematics(int ir=1) {
  drawer d;

  const char * rname = ana::rnames[ir];
  float model_sys[ana::nJetR] = {0.0079, 0, 0.0206, 0.0337};
  TFile * fn = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_nominal_%s.root","R04"));//rname));
  TFile * fb = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_bdt_%s.root"    ,"R04"));//rname));
  TFile * fi = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_iso_%s.root"    ,"R04"));//rname));
  TFile * f3 = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_3jet_%s.root"   ,"R04"));//rname));
  TFile * fh = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_JERhigh_%s.root","R04"));//rname));
  TFile * fl = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_JERlow_%s.root" ,"R04"));//rname));
  TFile * fH = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_HERWIG_%s.root" ,"R04"));//rname));
  const int nfiles = 7;
  TFile * f[nfiles] = {fn,fb,fi,f3,fh,fl,fH};
  TF1 * func[nfiles];
  int colors[nfiles] = {kBlack, kBlue, kRed, kGreen, kOrange, kMagenta, kCyan};
  string info[nfiles] = {"Nominal","Narrow BDT", "Narrow Isolation", "Third Jet Cut", "High JER smearing", "Low JER smearing", "HERWIG"};
    

  TCanvas * c = new TCanvas("c","",700,1000);
  TPad * p1 = new TPad("p1","",0,0.5,1,1);
  TPad * p2 = new TPad("p2","",0,0,1,0.5);
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);
  p1->Draw();
  p2->Draw();
  p1->SetTicks(1,1);
  p2->SetTicks(1,1);
  p1->cd();
  TLegend * ltop = new TLegend(0.15,0.50,0.5,0.86);
  for (int i = 0; i < nfiles; i++) {
    func[i] = (TF1*)f[i]->Get("fbest");
    func[i]->SetName(Form("f%i",i));
    func[i]->SetLineColor(colors[i]);
    func[i]->SetLineWidth(2);
    //func[i]->GetYaxis()->SetRangeUser(0.8,1.2);
    if (i == 0) {
      func[i]->SetTitle("");
      func[i]->GetYaxis()->SetRangeUser(0.9,1.2);
      func[i]->GetYaxis()->SetTitle("#it{in situ} correction");
      func[i]->GetXaxis()->SetTitle("p_{T}^{jet}");
      func[i]->Draw();
    }
    else func[i]->Draw("same");
    ltop->AddEntry(func[i], info[i].c_str());
  }
  func[0]->Draw("same");
  ltop->SetLineWidth(0);
  ltop->Draw();

  d.drawAll({"p+p #sqrt{s}=200 GeV"},{"Variant global #it{in situ} corrections",Form("Jet R=%1.1f",ana::JetRs[ir])},.5, .80, 20, 500);
  
  
  TF1 * faverage = new TF1(Form("faverage"), 
      [=](double *x, double *) {
      double v1 = func[0]->Eval(x[0]);
      double v2 = func[nfiles-1]->Eval(x[0]); 
      return (v1+v2)/2.0;
      },0,100,0);
  TF1 * diff = new TF1(Form("diff"), 
      [=](double *x, double *) {
      double v1 = func[0]->Eval(x[0]);
      double v2 = func[nfiles-1]->Eval(x[0]);
      return (v2-v1)/2.0;
      },0,100,0);
  
  
  p2->cd();
  TLegend * lbottom = new TLegend(0.15,0.6,0.5,0.89);
  TF1 * ratios[nfiles-1];
  for (int i = 1; i < nfiles; i++) {
    TF1 * ratio = new TF1(Form("r%i",i), 
        [=](double *x, double *) {
        double denom = func[0]->Eval(x[0]);
        //if (denom == 0) return 0.0;
        return func[i]->Eval(x[0]) / denom;
      },0,100,0);
    ratios[i-1] = (i != nfiles -1 ? ratio : diff);
    ratio->SetLineColor(colors[i]);
    ratio->SetLineWidth(2);
    //ratio->GetYaxis()->SetRangeUser(0.8,1.2);
    if (i == 1) {
      ratio->SetTitle("");
      ratio->GetYaxis()->SetRangeUser(0.9,1.2);
      ratio->GetYaxis()->SetTitle("ratio");
      ratio->GetXaxis()->SetTitle("p_{T}^{jet}");
      ratio->Draw();
    }
    else ratio->Draw("same");
    lbottom->AddEntry(ratio,info[i].c_str());
  }
  lbottom->SetLineWidth(0);
  lbottom->Draw();

  TF1 * systematic = new TF1("fsystematic",
      [=](double *x, double *) {
      double val = 0;
      for (int i = 0; i < nfiles - 1; i++ ) {
        if (i != nfiles - 2) val += (1-ratios[i]->Eval(x[0]))*(1-ratios[i]->Eval(x[0]));
        else val += ratios[i]->Eval(x[0])*ratios[i]->Eval(x[0]);
      }
      val += model_sys[ir]*model_sys[ir];
      return TMath::Sqrt(val);
      },0,100,0);
  d.drawMany({"Ratio:","Nominal to Variant"},.55, .85, 20, 500);
  
  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/global_systematic_%s.pdf",rname));
  
  
  TCanvas * c2 = new TCanvas("c2","",700,700);
  // Drawing fits
  TF1 * flb = (TF1*)f[0]->Get("fbest");
  TF1 * flh = (TF1*)f[0]->Get("fHigh");
  TF1 * fll = (TF1*)f[0]->Get("fLow");
  TTree * tl = (TTree*)f[0]->Get("T");
  float chisql;
  tl->SetBranchAddress("chisq",&chisql);
  tl->GetEntry(0);

  float xmin = 5;
  float xmax = 100;
  float ymin = 0.9;
  float ymax = 1.15;
  TH1F* frame2 = c2->DrawFrame(xmin, ymin, xmax, ymax);
  frame2->GetXaxis()->SetTitle("jet p_{T} [GeV]");
  frame2->GetYaxis()->SetTitle("#it{in situ} value");
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

    double y_best = flb->Eval(x) + diff->Eval(x);
    double y_low  = fll->Eval(x) + diff->Eval(x);
    double y_high = flh->Eval(x) + diff->Eval(x);

    double err_low  = y_best - y_low;
    double err_high = y_high - y_best;

    bandl->SetPoint(i, x, y_best);
    bandl->SetPointError(i, 0, 0, err_low, err_high);
  }
  
  bandl->SetFillColor(kGreen);
  bandl->SetFillStyle(3001);
  bandl->SetLineColor(kGreen);
  bandl->SetLineWidth(2);
  bandl->Draw("3 same");
  
  // The systematic error band
  TGraphAsymmErrors* bands = new TGraphAsymmErrors(nPoints);

  for (int i = 0; i < nPoints; ++i) {
    double x = xMin + (xMax - xMin) * i / (nPoints - 1);

    double y_best = flb->Eval(x) + diff->Eval(x);
    double y_low  = y_best + systematic->Eval(x);
    double y_high = y_best - systematic->Eval(x);

    double err_low  = y_best - y_low;
    double err_high = y_high - y_best;

    bands->SetPoint(i, x, y_best);
    bands->SetPointError(i, 0, 0, err_low, err_high);
  }
  
  bands->SetFillColorAlpha(kBlue,0.5);
  bands->SetFillStyle(3001);
  bands->SetLineWidth(0);
  bands->SetMarkerSize(0);
  bands->Draw("3 same");
  
  bandl->SetFillColorAlpha(kGreen,0.8);
  bandl->SetFillStyle(3001);
  bandl->SetLineColor(kGreen);
  bandl->SetLineWidth(2);
  //bandl->Draw("3 same");
  
  
  faverage->SetLineColor(kBlack);
  faverage->Draw("same");

  TF1 * errhigh = new TF1(Form("errhigh"), 
      [=](double *x, double *) {
      double e1 = systematic->Eval(x[0]);
      double e2 = flh->Eval(x[0])-flb->Eval(x[0]);
      return TMath::Sqrt(e1*e1+e2*e2);
      },0,100,0);
  
  TF1 * errlow = new TF1(Form("errlow"), 
      [=](double *x, double *) {
      double e1 = systematic->Eval(x[0]);
      double e2 = fll->Eval(x[0])-flb->Eval(x[0]);
      return TMath::Sqrt(e1*e1+e2*e2);
      },0,100,0);
  
  TF1 * linehigh = new TF1(Form("fhigh"), 
      [=](double *x, double *) {
      return faverage->Eval(x[0])+errhigh->Eval(x[0]);
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

  TLegend * leg2 = new TLegend(0.2,0.15,0.8,0.4);
  leg2->AddEntry(bandl, "statistical uncertainty");
  //leg2->AddEntry(bandl, Form("#splitline{%.3f + %.1e*pT}{#chi^{2} = %.2f}", flb->GetParameter(0), flb->GetParameter(1), chisql));
  leg2->AddEntry(bands, Form("systematic uncertainty"));
  leg2->AddEntry(linehigh, Form("total uncertainty"));
  
  leg2->SetLineWidth(0);
  leg2->SetFillStyle(0);
  leg2->Draw();
  
  d.drawAll({"p + p #sqrt{s} = 200 GeV"},{"|vz| < 60 cm", "|#eta| < 1.1",Form("#it{in situ} correction: %.3f + %.1e*pT", (flb->GetParameter(0)+func[nfiles-1]->GetParameter(0))/2.0, (flb->GetParameter(1)+func[nfiles-1]->GetParameter(1))/2.0),Form("Jet R=%1.1f",ana::JetRs[ir])},.2,.83,20, 700);

  c2->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/global_insitu_%s.pdf",rname));

  TFile * ff = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/funcs_insitu_%s.root",rname),"RECREATE");
  faverage->SetName("finsitu");
  faverage->Write();
  linehigh->Write();
  linelow->Write();

}
