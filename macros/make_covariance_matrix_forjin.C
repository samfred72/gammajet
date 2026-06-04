void make_covariance_matrix_forjin() {
  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/insitu_fit_linear_forjin.root","READ");
  TH2D * h2 = (TH2D*)f->Get("hchisq2d");
  
  // -----------------------------
  // Compute covariance from TH2D
  // -----------------------------
  double sumw = 0.0;
  double mean_a = 0.0;
  double mean_b = 0.0;
  double minChi = 23.7094;

  // First pass: weighted mean
  for (int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h2->GetNbinsY(); iy++) {

      double chi = h2->GetBinContent(ix, iy);
      if (!std::isfinite(chi)) continue;

      double pa = h2->GetXaxis()->GetBinCenter(ix);
      double pb = h2->GetYaxis()->GetBinCenter(iy);

      // Convert χ² → likelihood weight
      double w = exp(-0.5 * (chi - minChi));  // minChi must be known

      sumw   += w;
      mean_a += w * pa;
      mean_b += w * pb;
    }
  }

  mean_a /= sumw;
  mean_b /= sumw;

  // Second pass: covariance
  double cov_aa = 0.0;
  double cov_bb = 0.0;
  double cov_ab = 0.0;

  for (int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h2->GetNbinsY(); iy++) {

      double chi = h2->GetBinContent(ix, iy);
      if (!std::isfinite(chi)) continue;

      double pa = h2->GetXaxis()->GetBinCenter(ix);
      double pb = h2->GetYaxis()->GetBinCenter(iy);

      double w = exp(-0.5 * (chi - minChi));

      double da = pa - mean_a;
      double db = pb - mean_b;

      cov_aa += w * da * da;
      cov_bb += w * db * db;
      cov_ab += w * da * db;
    }
  }

  cov_aa /= sumw;
  cov_bb /= sumw;
  cov_ab /= sumw;

  // -----------------------------
  // Print result
  // -----------------------------
  cout << "\nCovariance Matrix:\n";
  cout << "[ " << cov_aa << " , " << cov_ab << " ]\n";
  cout << "[ " << cov_ab << " , " << cov_bb << " ]\n";

  cout << "\nUncertainties:\n";
  cout << "sigma_pa = " << sqrt(cov_aa) << endl;
  cout << "sigma_pb = " << sqrt(cov_bb) << endl;
  cout << "corr = " << cov_ab / sqrt(cov_aa * cov_bb) << endl;



  // Assume: TH2D* h2 already exists

  // -----------------------------
  // Style
  // -----------------------------
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  // -----------------------------
  // Find minimum χ² and location
  // -----------------------------
  minChi = 1e30;
  int minBinX = -1;
  int minBinY = -1;

  for (int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h2->GetNbinsY(); iy++) {
      double chi = h2->GetBinContent(ix, iy);
      if (!std::isfinite(chi)) continue;

      if (chi < minChi) {
        minChi = chi;
        minBinX = ix;
        minBinY = iy;
      }
    }
  }

  double best_pa = h2->GetXaxis()->GetBinCenter(minBinX);
  double best_pb = h2->GetYaxis()->GetBinCenter(minBinY);

  // -----------------------------
  // Canvas
  // -----------------------------
  TCanvas *c = new TCanvas("c_chi2","chi2 surface",800,700);
  gPad->SetLeftMargin(.17);
  gPad->SetRightMargin(.15);
  gPad->SetBottomMargin(.15);
  gPad->SetTicks(1,1);

  // Draw heatmap
  h2->SetTitle(";p_{a};p_{b};#chi^{2}");
  h2->Draw("COLZ");

  // -----------------------------
  // Draw best-fit point
  // -----------------------------
  TMarker *m = new TMarker(best_pa, best_pb, 29);
  m->SetMarkerSize(2.0);
  m->SetMarkerColor(kBlack);
  m->Draw("SAME");

  // -----------------------------
  // Draw 1σ contour (Δχ² = 2.30)
  // -----------------------------
  double level = minChi + 2.30;

  // Clone histogram for contour drawing
  TH2D *hcont = (TH2D*)h2->Clone("hcont");
  hcont->SetContour(1, &level);
  hcont->SetLineColor(kBlack);
  hcont->SetLineWidth(3);
  hcont->SetFillStyle(0);

  // IMPORTANT: draw contours on top
  hcont->Draw("CONT3 SAME");

  // -----------------------------
  // Legend
  // -----------------------------
  TLegend *leg = new TLegend(0.20,0.75,0.45,0.88);
  leg->AddEntry(m, "Best fit", "p");
  leg->AddEntry(hcont, "1#sigma contour", "l");
  leg->Draw();

  // -----------------------------
  // Print info
  // -----------------------------
  cout << "Best fit:\n";
  cout << "pa = " << best_pa << endl;
  cout << "pb = " << best_pb << endl;
  cout << "chi2 = " << minChi << endl;
}
