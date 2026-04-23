
#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_insitu_fit(const char * form = "nominal", int ir = 1) {
  const char * rname = ana::rnames[ir];
  TFile * fl = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_%s_%s.root",form,rname));
  TFile * fq = TFile::Open("/home/samson72/sphnx/gammajet/hists/insitu_fit_quad.root");
  TH1D * hs = (TH1D*)fl->Get("hstandardinsitu"); //s for standard, l for linear, q for quadratic
  TH1D * hs3 = (TH1D*)fl->Get("hstandardinsitu3");
  TH1D * hl = (TH1D*)fl->Get("hinsitu");
  TH1D * hl3 = (TH1D*)fl->Get("hinsitu3");
  TH1D * hq = (TH1D*)fq->Get("hinsitu");
  TH1D * hq3 = (TH1D*)fq->Get("hinsitu3");
   
  double xmin = 10;
  double xmax = 60;
  double ymin = 0.9;
  double ymax = 1.15;
  cout << "Drawing insitu" << endl;
  TCanvas * c = new TCanvas("c","",700,700);
  TH1F* frame = c->DrawFrame(xmin, ymin, xmax, ymax);
  frame->GetXaxis()->SetTitle("cluster p_{T} [GeV]");
  frame->GetYaxis()->SetTitle("ratio data <x_{J}>/MC <x_{J}>");
  frame->GetXaxis()->SetTitleColor(kBlue);
  frame->GetXaxis()->SetLabelColor(kBlue);


  TGaxis *topAxis = new TGaxis(xmin, ymax, xmax, ymax, xmin, xmax, 510, "-");
  topAxis->SetTitle("leading jet p_{T} [GeV]");
  topAxis->SetTitleColor(kOrange+1);
  topAxis->SetLabelColor(kOrange+1);
  topAxis->SetTitleFont(42);
  topAxis->SetLabelFont(42);
  topAxis->Draw();

  c->GetFrame()->SetLineWidth(0);
  c->Modified();
  c->Update();
  gStyle->SetOptStat(0);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  
  hs->SetLineColor(kBlue);
  hs->SetMarkerColor(kBlue);
  hs->SetMarkerStyle(24);
  hs->SetMarkerSize(1);
  hs->Draw("same");
  

  hl->SetLineColor(kBlue);
  hl->SetMarkerColor(kBlue);
  hl->SetMarkerStyle(21);
  hl->SetMarkerSize(1);
  hl->Draw("same");
  
  //hq->SetLineColor(kBlue);
  //hq->SetMarkerColor(kBlue);
  //hq->SetMarkerStyle(21);
  //hq->SetMarkerSize(1);
  //hq->Draw("same");

  hs3->SetLineColor(kOrange+1);
  hs3->SetMarkerColor(kOrange+1);
  hs3->SetMarkerStyle(24);
  hs3->SetMarkerSize(1);
  hs3->Draw("same");

  hl3->SetLineColor(kOrange+1);
  hl3->SetMarkerColor(kOrange+1);
  hl3->SetMarkerStyle(21);
  hl3->SetMarkerSize(1);
  hl3->Draw("same");
  
  //hq3->SetLineColor(kOrange+10);
  //hq3->SetMarkerColor(kOrange+10);
  //hq3->SetMarkerStyle(21);
  //hq3->SetMarkerSize(1);
  //hq3->Draw("same");

  TLegend * leg = new TLegend(0.17,0.6,0.6,0.89);
  leg->AddEntry(hs, "Uncorrected #gamma-Jet ratio");
  leg->AddEntry(hl, "Linearly corrected #gamma-Jet ratio");
  //leg->AddEntry(hq, "Quadratically corrected #gamma-Jet ratio");
  leg->AddEntry(hs3, "Uncorrected multijet ratio");
  leg->AddEntry(hl3, "Linearly corrected multijet ratio");
  //leg->AddEntry(hq3, "Quadratically corrected multijet ratio");
  leg->SetLineWidth(0);
  leg->Draw("same");

  // Line for 1
  TLine * l = new TLine(10,1,60, 1);
  l->SetLineStyle(9);
  l->Draw("same");

  drawer d;
  d.drawAll({"p + p #sqrt{s} = 200 GeV"},{"|vz| < 60 cm", "|#eta| < 1.1",Form("Jet R=%1.1f",ana::JetRs[ir])},.6,.83,15, 700);

  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/corrected_insitu_%s_%s.pdf", form, rname));



  // Drawing fits
  TF1 * flb = (TF1*)fl->Get("fbest");
  TF1 * flh = (TF1*)fl->Get("fHigh");
  TF1 * fll = (TF1*)fl->Get("fLow");
  TF1 * fqb = (TF1*)fq->Get("fbest");
  TF1 * fqh = (TF1*)fq->Get("fHigh");
  TF1 * fql = (TF1*)fq->Get("fLow");
  TTree * tl = (TTree*)fl->Get("T");
  TTree * tq = (TTree*)fq->Get("T");
  float chisql;
  float chisqq;
  tl->SetBranchAddress("chisq",&chisql);
  tq->SetBranchAddress("chisq",&chisqq);
  tl->GetEntry(0);
  tq->GetEntry(0);

  TCanvas * c2 = new TCanvas("c2","",700,700);
  xmin = 5;
  xmax = 80;
  ymin = 0.9;
  ymax = 1.15;
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

  TGraphAsymmErrors* bandl = new TGraphAsymmErrors(nPoints);

  for (int i = 0; i < nPoints; ++i) {
    double x = xMin + (xMax - xMin) * i / (nPoints - 1);

    double y_best = flb->Eval(x);
    double y_low  = fll->Eval(x);
    double y_high = flh->Eval(x);

    double err_low  = y_best - y_low;
    double err_high = y_high - y_best;

    bandl->SetPoint(i, x, y_best);
    bandl->SetPointError(i, 0, 0, err_low, err_high);
  }
  
  TGraphAsymmErrors* bandq = new TGraphAsymmErrors(nPoints);

  for (int i = 0; i < nPoints; ++i) {
    double x = xMin + (xMax - xMin) * i / (nPoints - 1);

    double y_best = fqb->Eval(x);
    double y_low  = fql->Eval(x);
    double y_high = fqh->Eval(x);

    double err_low  = y_best - y_low;
    double err_high = y_high - y_best;

    bandq->SetPoint(i, x, y_best);
    bandq->SetPointError(i, 0, 0, err_low, err_high);
  }
  
  bandq->SetFillColor(kViolet);
  bandq->SetFillStyle(3001);
  bandq->SetLineColor(kBlue);
  bandq->SetLineWidth(2);
  //bandq->Draw("3 same");
  
  bandl->SetFillColor(kSpring-2);
  bandl->SetFillStyle(3001);
  bandl->SetLineColor(kGreen);
  bandl->SetLineWidth(2);
  bandl->Draw("3 same");
  
  
  
  flb->SetLineColor(kGreen);
  flb->Draw("same");
  fqb->SetLineColor(kBlue);
  //fqb->Draw("same");

  TLegend * leg2 = new TLegend(0.2,0.15,0.8,0.4);
  leg2->AddEntry(bandl, Form("#splitline{%.3f + %.1e*pT}{#chi^{2} = %.2f}", flb->GetParameter(0), flb->GetParameter(1), chisql));
  //leg2->AddEntry(bandq, Form("#splitline{%.3f + %.2e*pT + %.2e*pT^{2}}{#chi^{2} = %.2f}", fqb->GetParameter(0), fqb->GetParameter(1), fqb->GetParameter(2), chisqq));
  
  leg2->SetLineWidth(0);
  leg2->SetFillStyle(0);
  leg2->Draw();
  
  d.drawAll({"p + p #sqrt{s} = 200 GeV"},{"|vz| < 60 cm", "|#eta| < 1.1"},.2,.83,20, 700);

  //d.drawText("Best func: ", .15,.25,1,20);
  //d.drawText(Form("%.5f + %.5f*pT + %.5f*pT^{2}",minpa,minpb,minpc), .15,.20,1,20);
}
