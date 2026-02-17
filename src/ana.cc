#include "/home/samson72/sphnx/gammajet/src/ana.h"
using namespace std;
ana::ana() {
}

Bool_t ana::PassEtaCut(float eta, float vz = 0)
{
  //float loweta = GetShiftedEta(vz,etamin);
  //float higheta = GetShiftedEta(vz,etamax);
  if (eta < etamin || eta > etamax) return false;
  else return true;
}

Double_t ana::GetShiftedEta(float _vz, float _eta)
{
  double theta = 2*atan(exp(-_eta));
  double z = radius / tan(theta);
  double zshifted = z - _vz;
  double thetashifted = atan2(radius,zshifted);
  double etashifted = -log(tan(thetashifted/2.0));
  return etashifted;
}

Int_t ana::findPtBin(double value)
{
  for (int i = 0; i < nPtBins; ++i) {
    if (value >= ptBins[i] && value < ptBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}
Int_t ana::findxjBin(double value)
{
  for (int i = 0; i < nxjBins; ++i) {
    if (value >= xjBins[i] && value < xjBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}

Int_t ana::findabcdBin(double value)
{
  for (int i = 0; i < nabcdbins; ++i) {
    if (value >= abcdbins[i] && value < abcdbins[i + 1]) {
      return  i;
    }
  }
  if (value > abcdbins[nabcdbins-1]) return nabcdbins-1;
  else return -1;
}

float ana::getPurity(float low, float high) {
  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/purity.root");
  TF1 * func = (TF1*)f->Get("func");
  float val = func->Integral(low,high)/(high-low);
  f->Close();
  return val;
}
