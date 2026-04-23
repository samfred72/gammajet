
#include "/home/samson72/sphnx/gammajet/src/drawer.h"
#include "/home/samson72/sphnx/gammajet/src/ana.h"

void draw_TH1D(const char * sim, bool unfolding, const char * histname, int sample = 1, int subsample = -1, const char * opt = "hist") {
  drawer d = drawer(sim);
  if (unfolding){
    d = drawer(unfolding,sim);
  }
  TH1D * h = d.get(histname,sample,subsample);
  h->Draw(opt);
}
void draw_TH2D(const char * sim, bool unfolding, const char * histname, int sample = 1, int subsample = -1, const char * opt = "colz") {
  drawer d = drawer(sim);
  if (unfolding){
    d = drawer(unfolding,sim);
  }
  TH2D * h = d.get2d(histname,sample,subsample);
  h->Draw();
}
