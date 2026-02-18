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
Int_t ana::findBdtBin(double value)
{
  for (int i = 0; i < nBdtBins; ++i) {
    if (value >= bdtBins[i] && value < bdtBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}

Int_t ana::findabcdBin(double iso, double bdt, int bin)
{
  int isiso;
  int isbdt;
  if (iso < isoBins[bin]) {
    isiso = 1;
  }
  else if (iso > isoBinsHigh[bin]) {
    isiso = 0;
  }
  else {
    isiso = -1;
  }
  if (bdt > bdtGoodLow[bin] && bdt < bdtGoodHigh[bin]) {
    isbdt = 1;
  }
  else if (bdt > bdtBadLow[bin] && bdt < bdtBadHigh[bin]) {
    isbdt = 0;
  }
  else {
    isbdt = -1;
  }

  if (isiso == -1 || isbdt == -1) {
    return -1;
  }
  else {
    bool b_isiso = isiso;
    bool b_isbdt = isbdt;
    int iabcd = (((b_isiso << 0b1) | b_isbdt) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    return iabcd;
  }
}

Int_t ana::findabcdBin(double iso, int showershape, int bin) {
  int isiso;
  int istight;
  if (iso < isoBins[bin]) {
    isiso = 1;
  }
  else if (iso > isoBinsHigh[bin]) {
    isiso = 0;
  }
  else {
    isiso = -1;
  }
  if (showershape == 2) {
    istight = 1;
  }
  else if (showershape == 1) {
    istight = 0;
  }
  else {
    istight = -1;
  }
  if (isiso == -1 || istight == -1) {
    return -1;
  }
  else {
    bool b_isiso = isiso;
    bool b_istight = istight;
    int iabcd = (((b_isiso << 0b1) | b_istight) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    return iabcd;
  }
}

float ana::getPurity(float low, float high) {
  TFile * f = TFile::Open("/home/samson72/sphnx/gammajet/hists/purity.root");
  TF1 * func = (TF1*)f->Get("func");
  float val = func->Integral(low,high)/(high-low);
  f->Close();
  return val;
}
