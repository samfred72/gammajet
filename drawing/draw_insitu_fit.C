
#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_insitu_fit(const char * form = "nominal", int ir = 2) {
  const char * rname = ana::rnames[ir];
  TFile * fl = TFile::Open(Form("/home/samson72/sphnx/gammajet/hists/insitu_fit_%s_%s.root",form,rname));
  TFile * fq = TFile::Open("/home/samson72/sphnx/gammajet/hists/insitu_fit_quad.root");
  TH1D * hs = (TH1D*)fl->Get("hstandardinsitu"); //s for standard, l for linear, q for quadratic
  TH1D * hs3 = (TH1D*)fl->Get("hstandardinsitu3");
  TH1D * hl = (TH1D*)fl->Get("hinsitu");
  TH1D * hl3 = (TH1D*)fl->Get("hinsitu3");
  TH1D * hq = (TH1D*)fq->Get("hinsitu");
  TH1D * hq3 = (TH1D*)fq->Get("hinsitu3");
   
  cout << "Drawing insitu" << endl;
  TCanvas * c = new TCanvas("c","",700,700);

  gStyle->SetOptStat(0);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.15);
  
  hs->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  hs->GetYaxis()->SetTitle("<x_{J#gamma}>_{Data} / <x_{J#gamma}>_{Sim}");
  hs->GetXaxis()->SetTitleOffset(1.3);
  hs->GetYaxis()->SetRangeUser(0.9,1.15);
  hs->SetLineColor(kBlue);
  hs->SetMarkerColor(kBlue);
  hs->SetMarkerStyle(24);
  hs->SetMarkerSize(1);
  hs->GetXaxis()->SetLabelSize(hs->GetXaxis()->GetLabelSize()*1.15);
  hs->GetYaxis()->SetLabelSize(hs->GetYaxis()->GetLabelSize()*1.15);
  hs->GetXaxis()->SetTitleSize(hs->GetXaxis()->GetTitleSize()*1.15);
  hs->GetYaxis()->SetTitleSize(hs->GetYaxis()->GetTitleSize()*1.15);
  hs->Draw("same");
  
  hl->SetLineColor(kBlue);
  hl->SetMarkerColor(kBlue);
  hl->SetMarkerStyle(21);
  hl->SetMarkerSize(1);
  //hl->Draw("same");
  
  const int n = hl->GetNbinsX();

  std::vector<double> x(n), y(n), ex(n), ey(n);

  for (int i = 0; i < n; ++i) {
    x[i]  = hl->GetBinCenter(i+1) + 0.05; // small offset
    y[i]  = hl->GetBinContent(i+1);
    ex[i] = hl->GetBinWidth(i+1)*0.5;
    ey[i] = hl->GetBinError(i+1);
  }

  gStyle->SetEndErrorSize(0);
  TGraphErrors *gdivide_data =
    new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());

  gdivide_data->SetMarkerColor(kBlue);
  gdivide_data->SetLineColor(kBlue);
  gdivide_data->SetMarkerStyle(21);
  gdivide_data->SetMarkerSize(1);
  gdivide_data->Draw("P SAME"); 


  
  TLegend * leg = new TLegend(0.17,0.53,0.6,0.68);
  leg->AddEntry(hs, "Uncorrected #gamma-Jet ratio");
  leg->AddEntry(hl, "Corrected #gamma-Jet ratio");
  leg->SetLineWidth(0);
  leg->Draw();
  
  // Line for 1
  TLine * l = new TLine(10,1,35, 1);
  l->SetLineStyle(9);
  l->Draw("same");
  
  drawer d;
  d.drawAll({"p + p #sqrt{s} = 200 GeV"},{
      //"|vz| < 60 cm", 
      Form("|#eta_{jet}| < %1.1f",1.1-ana::JetRs[ir]),
      Form("Jet R=%1.1f", ana::JetRs[ir]),
      },.20,.83,21, 700);

  c->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/note/corrected_insitu_gammajet_%s_%s.pdf", form, rname));
  
  TCanvas * c3 = new TCanvas("c3","",700,700);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.15);
  
  hs3->GetXaxis()->SetTitle("p_{T,1}^{jet} [GeV]");
  hs3->GetYaxis()->SetTitle("<x_{J}^{-1}>_{Data} / <x_{J}^{-1}>_{Sim}");
  hs3->GetXaxis()->SetTitleOffset(1.3);
  hs3->GetYaxis()->SetRangeUser(0.9,1.15);
  hs3->SetLineColor(kRed+2);
  hs3->SetMarkerColor(kRed+2);
  hs3->SetMarkerStyle(24);
  hs3->SetMarkerSize(1);
  hs3->GetXaxis()->SetLabelSize(hs3->GetXaxis()->GetLabelSize()*1.15);
  hs3->GetYaxis()->SetLabelSize(hs3->GetYaxis()->GetLabelSize()*1.15);
  hs3->GetXaxis()->SetTitleSize(hs3->GetXaxis()->GetTitleSize()*1.15);
  hs3->GetYaxis()->SetTitleSize(hs3->GetYaxis()->GetTitleSize()*1.15);
  hs3->Draw("same");

  hl3->SetLineColor(kRed+2);
  hl3->SetMarkerColor(kRed+2);
  hl3->SetMarkerStyle(21);
  hl3->SetMarkerSize(1);
  //hl3->Draw("same");
  
  const int n3 = hl3->GetNbinsX();

  std::vector<double> x3(n), y3(n), ex3(n), ey3(n);

  for (int i = 0; i < n3; ++i) {
    x3[i]  = hl3->GetBinCenter(i+1) + 0.08; // small offset
    y3[i]  = hl3->GetBinContent(i+1);
    ex3[i] = hl3->GetBinWidth(i+1)*0.5;
    ey3[i] = hl3->GetBinError(i+1);
  }

  gStyle->SetEndErrorSize(0);
  TGraphErrors *gdivide_data3 =
    new TGraphErrors(n3, x3.data(), y3.data(), ex3.data(), ey3.data());

  gdivide_data3->SetMarkerColor(kRed+2);
  gdivide_data3->SetLineColor(kRed+2);
  gdivide_data3->SetMarkerStyle(21);
  gdivide_data3->SetMarkerSize(1);
  gdivide_data3->Draw("P SAME"); 
  
  TLegend * leg3 = new TLegend(0.17,0.53,0.6,0.68);
  leg3->AddEntry(hs3, "Uncorrected multijet ratio");
  leg3->AddEntry(hl3, "Corrected multijet ratio");
  leg3->SetLineWidth(0);
  leg3->Draw("same");

  TLine * l3 = new TLine(20,1,60, 1);
  l3->SetLineStyle(9);
  l3->Draw("same");

  d.drawAll({"p + p #sqrt{s} = 200 GeV"},{
      //"|vz| < 60 cm", 
      Form("|#eta_{jet}| < %1.1f",1.1-ana::JetRs[ir]),
      Form("Jet R=%1.1f",ana::JetRs[ir])
      },.20,.83,21, 700);

  c3->SaveAs(Form("/home/samson72/sphnx/gammajet/pdfs/note/corrected_insitu_multijet_%s_%s.pdf", form, rname));



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
  float xmin = 5;
  float xmax = 80;
  float ymin = 0.9;
  float ymax = 1.15;
  TH1F* frame2 = c2->DrawFrame(xmin, ymin, xmax, ymax);
  frame2->GetXaxis()->SetTitle("jet p_{T} [GeV]");
  frame2->GetYaxis()->SetTitle("#it{in situ} value");
  //c2->GetFrame()->SetLineWidth(0);
  gStyle->SetOptStat(0);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  c2->Modified();
  c2->Update();
  
  
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
  

  //d.drawText("Best func: ", .15,.25,1,20);
  //d.drawText(Form("%.5f + %.5f*pT + %.5f*pT^{2}",minpa,minpb,minpc), .15,.20,1,20);
}
