#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_forjin() {

  drawer d;

  TH1D * hs[ana::nPtBins]; // scaled
  TH1D * hu[ana::nPtBins]; // unscaled

  for (int i = 0; i < ana::nPtBins; i++) {

    hu[i] = d.get(Form("h1_%i",i),1);
    hu[i]->SetLineColor(kGreen+2);
    hu[i]->SetMarkerColor(kGreen+2);
    hu[i]->SetMarkerStyle(20);

    hs[i] = d.get(Form("h2_%i",i),1);
    hs[i]->SetLineColor(kMagenta+1);
    hs[i]->SetMarkerColor(kMagenta+1);
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

    cout << meanS << " " << meanU << " " << meanS/meanU << endl;

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

  leg->AddEntry(hMeanS,"Jets scaled by 1/0.95","pl");
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
  hRatio->GetYaxis()->SetRangeUser(0.95,1.1);

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
      1.05,
      ana::ptBins[ana::nPtBins],
      1.05);

  line2->SetLineColor(kRed);
  line2->SetLineStyle(9);
  line2->Draw();
  c->cd();
  c->Update();
}
