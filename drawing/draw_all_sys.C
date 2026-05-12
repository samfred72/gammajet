void draw_all_sys() {
  TFile * f2 = TFile::Open("/home/samson72/sphnx/gammajet/hists/funcs_insitu_R02.root","READ");
  TFile * f3 = TFile::Open("/home/samson72/sphnx/gammajet/hists/funcs_insitu_R03.root","READ");
  TFile * f4 = TFile::Open("/home/samson72/sphnx/gammajet/hists/funcs_insitu_R04.root","READ");
  TFile * f5 = TFile::Open("/home/samson72/sphnx/gammajet/hists/funcs_insitu_R05.root","READ");
  TFile * f6 = TFile::Open("/home/samson72/sphnx/gammajet/hists/funcs_insitu_R06.root","READ");
  TFile * f7 = TFile::Open("/home/samson72/sphnx/gammajet/hists/funcs_insitu_R07.root","READ");
  TFile * f8 = TFile::Open("/home/samson72/sphnx/gammajet/hists/funcs_insitu_R08.root","READ");
  TF1 * func2 = (TF1*)f2->Get("fhigh");
  TF1 * func3 = (TF1*)f3->Get("fhigh");
  TF1 * func4 = (TF1*)f4->Get("fhigh");
  TF1 * func5 = (TF1*)f5->Get("fhigh");
  TF1 * func6 = (TF1*)f6->Get("fhigh");
  TF1 * func7 = (TF1*)f7->Get("fhigh");
  TF1 * func8 = (TF1*)f8->Get("fhigh");
  func2->GetYaxis()->SetRangeUser(0.9,1.15);
  func2->Draw();
  func3->Draw("same");
  func4->Draw("same");
  func5->Draw("same");
  func6->Draw("same");
  func7->Draw("same");
  func8->Draw("same");
}
