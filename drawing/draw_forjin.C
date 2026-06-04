#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_forjin() {

  drawer d;

  TH1D * hs[ana::nPtBins]; // scaled
  TH1D * hu[ana::nPtBins]; // unscaled

  for (int i = 0; i < ana::nPtBins; i++) {

    hu[i] = d.get(Form("h2_%i",i),1,-1);
    hu[i]->SetLineColor(kMagenta+1);
    hu[i]->SetMarkerColor(kMagenta+1);
    hu[i]->SetMarkerStyle(20);

    hs[i] = d.get(Form("h1_%i",i),1,-1);
    hs[i]->SetLineColor(kGreen+2);
    hs[i]->SetMarkerColor(kGreen+2);
    hs[i]->SetMarkerStyle(21);
  }

  // ---------------------------------------------------------
  // Histograms for means with unequal pT binning
  // ---------------------------------------------------------

  TH1D * hMeanS = new TH1D(
      "hMeanS",
      ";p_{T};Mean",
      ana::nPtBins,
      ana::ptBins);

  TH1D * hMeanU = new TH1D(
      "hMeanU",
      ";p_{T};Mean",
      ana::nPtBins,
      ana::ptBins);

  TH1D * hRatio = new TH1D(
      "hRatio",
      ";p_{T};Scaled / Unscaled",
      ana::nPtBins,
      ana::ptBins);

  // Fill mean histograms
  for (int i = 0; i < ana::nPtBins; i++) {

    double meanS = hs[i]->GetMean();
    double errS  = hs[i]->GetMeanError();

    double meanU = hu[i]->GetMean();
    double errU  = hu[i]->GetMeanError();

    cout << setprecision(15) << meanS << " " << meanU << " " << meanS/meanU << endl;

    hMeanS->SetBinContent(i+1, meanS);
    hMeanS->SetBinError(i+1, errS);

    hMeanU->SetBinContent(i+1, meanU);
    hMeanU->SetBinError(i+1, errU);

    // Ratio + propagated uncertainty
    double ratio = 0.0;
    double ratioErr = 0.0;

    if (meanU != 0.0) {

      ratio = meanS / meanU;

      double relErrS = (meanS != 0.0) ? errS / meanS : 0.0;
      double relErrU = errU / meanU;

      ratioErr = ratio * std::sqrt(relErrS*relErrS +
                                   relErrU*relErrU);
    }

    hRatio->SetBinContent(i+1, ratio);
    hRatio->SetBinError(i+1, ratioErr);
  }

  // ---------------------------------------------------------
  // Style
  // ---------------------------------------------------------

  hMeanU->SetLineColor(kMagenta+1);
  hMeanU->SetMarkerColor(kMagenta+1);
  hMeanU->SetMarkerStyle(20);

  hMeanS->SetLineColor(kGreen+2);
  hMeanS->SetMarkerColor(kGreen+2);
  hMeanS->SetMarkerStyle(21);

  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20);

  // ---------------------------------------------------------
  // Canvas with two equal pads
  // ---------------------------------------------------------

  TCanvas * c = new TCanvas("c","c",1000,1000);

  TPad * pTop = new TPad("pTop","pTop",0.0,0.5,1.0,1.0);
  TPad * pBot = new TPad("pBot","pBot",0.0,0.0,1.0,0.5);

  pTop->SetBottomMargin(0.02);
  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.15);

  pTop->Draw();
  pBot->Draw();

  // ---------------------------------------------------------
  // Top pad
  // ---------------------------------------------------------

  pTop->cd();

  hMeanS->SetStats(0);
  hMeanS->SetTitle("");

  hMeanS->GetYaxis()->SetTitle("Mean");
  hMeanS->GetYaxis()->SetTitleSize(0.05);
  hMeanS->GetYaxis()->SetLabelSize(0.04);
  hMeanS->GetYaxis()->SetRangeUser(0.5,1.5);

  hMeanS->GetXaxis()->SetLabelSize(0);

  hMeanS->Draw("E1");
  hMeanU->Draw("E1 SAME");

  TLegend * leg = new TLegend(0.55,0.62,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  leg->AddEntry(hMeanS,"Jets scaled by 1/0.975","pl");
  leg->AddEntry(hMeanU,"Jets unscaled","pl");

  leg->Draw();

  d.drawAll({"Pythia 8 #gamma-Jet"},{"Jet R=0.4"},.25,.8,20,500);

  // ---------------------------------------------------------
  // Bottom pad
  // ---------------------------------------------------------

  pBot->cd();

  hRatio->SetStats(0);
  hRatio->SetTitle("");

  hRatio->GetYaxis()->SetTitle("Ratio scaled / unscaled");
  hRatio->GetYaxis()->SetTitleSize(0.06);
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetYaxis()->SetLabelSize(0.05);
  hRatio->GetYaxis()->SetRangeUser(0.95,1.2);

  hRatio->GetXaxis()->SetTitle("p_{T}^{#gamma}");
  hRatio->GetXaxis()->SetTitleSize(0.06);
  hRatio->GetXaxis()->SetLabelSize(0.05);

  hRatio->Draw("E1");

  // Optional unity line
  TLine * line = new TLine(
      ana::ptBins[0],
      1.0,
      ana::ptBins[ana::nPtBins],
      1.0);

  line->SetLineStyle(2);
  line->Draw("SAME");

  TLine * line2 = new TLine(
      ana::ptBins[0],
      1.025,
      ana::ptBins[ana::nPtBins],
      1.025);

  line2->SetLineColor(kRed);
  line2->SetLineStyle(9);
  line2->Draw();
  c->cd();
  c->Update();



  TCanvas * cmany = new TCanvas("cmany","",700*4,700*3);
  const char * cname = "/home/samson72/sphnx/gammajet/pdfs/many_gridtest.pdf";
  gStyle->SetOptStat(0);

  float drawx = 0.77;
  float drawy = 0.92;
  float fontsize = 50;
  TLegend * l2 = new TLegend(drawx,drawy-.35,0.99,drawy-.15);

  int nx = 3;
  int ny = 3;

  TFile * fcorr = TFile::Open("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_nominal_R04.root","READ");
  fcorr->ls();

  int offset = 0;
  TPad * p = new TPad("p","",0,0,1,1);
  p->SetLeftMargin(0.25);
  p->SetBottomMargin(0.25);
  p->Divide(4,3,0,0);
  p->Draw();
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    TH1D * hd = (TH1D*)fcorr->Get(Form("hxj_MC_%i",ipt));
    TH1D * hs = (TH1D*)fcorr->Get(Form("hxj_data_%i",ipt));
    TH1D * hc = (TH1D*)fcorr->Get(Form("hxj_corrected_%i",ipt));
    hd->SetLineColor(kMagenta+1);
    hd->SetMarkerColor(kMagenta+1);
    hd->SetMarkerStyle(20);
    hd->SetMarkerSize(1);
    hd->Scale(1.0/hd->Integral());
    hs->SetLineColor(kGreen+2);
    hs->SetMarkerColor(kGreen+2);
    hs->SetMarkerStyle(21);
    hs->SetMarkerSize(1);
    hs->Scale(1.0/hs->Integral());
    hc->SetLineColor(kBlue);
    hc->SetMarkerColor(kBlue);
    hc->SetMarkerStyle(22);
    hc->SetMarkerSize(1);
    hc->Scale(1.0/hc->Integral());

    hd->GetYaxis()->SetRangeUser(0,0.16);

    int index = ipt + 1 + offset;
    if (index % (nx+1) == 0) {index++; offset++;}
    p->cd(index);
    gPad->SetRightMargin(0.01);
    gPad->SetTicks(1,1);
    if ((ipt) < nx*(ny-1)) gPad->SetBottomMargin(0.04);
    gPad->SetTopMargin(0.01);
    if ((ipt) % nx != 0) gPad->SetLeftMargin(0.05);
    if (index == (ny*(nx+1) - 1)) hd->GetXaxis()->SetTitle("x_{J#gamma}"); // bottom right plot
    else hd->GetXaxis()->SetTitle("");
    hd->GetXaxis()->SetTitleSize(0.10);
    hd->GetXaxis()->SetLabelSize(0.08);
    if (index == (nx+1)*(ny-1)+1) hd->GetXaxis()->SetLabelSize(0.07); // bottom left plot
    if (ipt < nx*(ny-1) ) hd->GetXaxis()->SetLabelSize(0.00);
    hd->GetYaxis()->ChangeLabel(-1,-1,0);
    if (index != (ny*(nx+1) - 1)) hd->GetXaxis()->ChangeLabel(-1,-1,0);
    if (index == 1) hd->GetYaxis()->SetTitle("Normalized Counts"); // top left plot
    else hd->GetYaxis()->SetTitle("");
    hd->GetYaxis()->SetTitleSize(0.10);
    hd->GetYaxis()->SetLabelSize(0.08);
    if (ipt % nx != 0 ) hd->GetYaxis()->SetLabelSize(0.00);
    if (ipt == nx*(ny-1)) hd->GetYaxis()->SetLabelSize(0.06); // bottom left plot
    hd->GetYaxis()->SetLabelOffset(0.04);
    hd->GetYaxis()->SetMaxDigits(3);
    hd->GetYaxis()->SetDecimals(2);

    hd->Draw("hist same");
    hs->Draw("hist same");
    hc->Draw("hist same");

    if (ipt == 0) {
      l2->AddEntry(hd,"Unscaled Jets");
      l2->AddEntry(hs,"Scaled Jets");
      l2->AddEntry(hc,"Corrected Jets");
    }
  }
  p->cd();
  d.drawAll({
      //info1, 
      //info2,
      },
      {
      //"Analysis cuts",
      Form("Jet R=%1.1f",ana::JetRs[2]),
      //calibstring,
      //t1[iabcd1].c_str()
      },
      drawx,drawy,fontsize,c->GetWh());
  l2->SetLineWidth(0);
  l2->Draw();
  cmany->SaveAs(cname);




  // Corrected version of first plot
  TH1D * hc[ana::nPtBins]; // corrected
  TH1D * hd[ana::nPtBins]; // unscaled

  for (int i = 0; i < ana::nPtBins; i++) {

    hd[i] = (TH1D*)fcorr->Get(Form("hxj_default_%i",i));
    hd[i]->SetLineColor(kMagenta+1);
    hd[i]->SetMarkerColor(kMagenta+1);
    hd[i]->SetMarkerStyle(20);

    hc[i] = (TH1D*)fcorr->Get(Form("hxj_corrected_%i",i));;
    hc[i]->SetLineColor(kBlue);
    hc[i]->SetMarkerColor(kBlue);
    hc[i]->SetMarkerStyle(22);
  }

  // ---------------------------------------------------------
  // Histograms for means with unequal pT binning
  // ---------------------------------------------------------

  TH1D * hMeanC = new TH1D(
      "hMeanC",
      ";p_{T};Mean",
      ana::nPtBins,
      ana::ptBins);

  TH1D * hMeanD = new TH1D(
      "hMeanD",
      ";p_{T};Mean",
      ana::nPtBins,
      ana::ptBins);

  TH1D * hRatio2 = new TH1D(
      "hRatio2",
      ";p_{T};Corrected / Unscaled",
      ana::nPtBins,
      ana::ptBins);

  // Fill mean histograms
  for (int i = 0; i < ana::nPtBins; i++) {

    double meanC = hc[i]->GetMean();
    double errC  = hc[i]->GetMeanError();

    double meanD = hd[i]->GetMean();
    double errD  = hd[i]->GetMeanError();

    cout << setprecision(15) << meanC << " " << meanD << " " << meanC/meanD << endl;

    hMeanC->SetBinContent(i+1, meanC);
    hMeanC->SetBinError(i+1, errC);

    hMeanD->SetBinContent(i+1, meanD);
    hMeanD->SetBinError(i+1, errD);

    // Ratio + propagated uncertainty
    double ratio = 0.0;
    double ratioErr = 0.0;

    if (meanD != 0.0) {

      ratio = meanC / meanD;

      double relErrC = (meanC != 0.0) ? errC / meanC : 0.0;
      double relErrD = errD / meanD;

      ratioErr = ratio * std::sqrt(relErrC*relErrC +
                                   relErrD*relErrD);
    }

    hRatio2->SetBinContent(i+1, ratio);
    hRatio2->SetBinError(i+1, ratioErr);
  }

  // ---------------------------------------------------------
  // Style
  // ---------------------------------------------------------

  hMeanD->SetLineColor(kMagenta+1);
  hMeanD->SetMarkerColor(kMagenta+1);
  hMeanD->SetMarkerStyle(20);

  hMeanC->SetLineColor(kBlue);
  hMeanC->SetMarkerColor(kBlue);
  hMeanC->SetMarkerStyle(22);

  hRatio2->SetLineColor(kBlack);
  hRatio2->SetMarkerColor(kBlack);
  hRatio2->SetMarkerStyle(20);

  // ---------------------------------------------------------
  // Canvas with two equal pads
  // ---------------------------------------------------------

  TCanvas * c2 = new TCanvas("c2","c2",1000,1000);

  TPad * pTop2 = new TPad("pTop2","pTop2",0.0,0.5,1.0,1.0);
  TPad * pBot2 = new TPad("pBot2","pBot2",0.0,0.0,1.0,0.5);

  pTop2->SetBottomMargin(0.02);
  pBot2->SetTopMargin(0.02);
  pBot2->SetBottomMargin(0.15);

  pTop2->Draw();
  pBot2->Draw();

  // ---------------------------------------------------------
  // Top pad
  // ---------------------------------------------------------

  pTop2->cd();

  hMeanC->SetStats(0);
  hMeanC->SetTitle("");

  hMeanC->GetYaxis()->SetTitle("Mean");
  hMeanC->GetYaxis()->SetTitleSize(0.05);
  hMeanC->GetYaxis()->SetLabelSize(0.04);
  hMeanC->GetYaxis()->SetRangeUser(0.5,1.5);

  hMeanC->GetXaxis()->SetLabelSize(0);

  hMeanC->Draw("E1");
  hMeanD->Draw("E1 SAME");

  TLegend * leg2 = new TLegend(0.35,0.42,0.88,0.68);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);

  leg2->AddEntry(hMeanC,"Jets corrected with f = 1.028 + -0.00016*pT","pl");
  leg2->AddEntry(hMeanD,"Jets unscaled","pl");

  leg2->Draw();

  d.drawAll({"Pythia 8 #gamma-Jet"},{"Jet R=0.4"},.25,.8,20,500);

  // ---------------------------------------------------------
  // Bottom pad
  // ---------------------------------------------------------

  pBot2->cd();

  hRatio2->SetStats(0);
  hRatio2->SetTitle("");

  hRatio2->GetYaxis()->SetTitle("Ratio scaled / unscaled");
  hRatio2->GetYaxis()->SetTitleSize(0.06);
  hRatio2->GetYaxis()->SetTitleOffset(0.8);
  hRatio2->GetYaxis()->SetLabelSize(0.05);
  hRatio2->GetYaxis()->SetRangeUser(0.95,1.2);

  hRatio2->GetXaxis()->SetTitle("p_{T}^{#gamma}");
  hRatio2->GetXaxis()->SetTitleSize(0.06);
  hRatio2->GetXaxis()->SetLabelSize(0.05);

  hRatio2->Draw("E1");

  // Optional unity line
  TLine * line3 = new TLine(
      ana::ptBins[0],
      1.0,
      ana::ptBins[ana::nPtBins],
      1.0);

  line3->SetLineStyle(2);
  line3->Draw("SAME");

  TLine * line4 = new TLine(
      ana::ptBins[0],
      1.025,
      ana::ptBins[ana::nPtBins],
      1.025);

  line4->SetLineColor(kRed);
  line4->SetLineStyle(9);
  //line4->Draw();
  c->cd();
  c->Update();


}
