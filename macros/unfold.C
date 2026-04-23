#include "/home/samson72/sphnx/gammajet/src/unfolder.h"
#include "/home/samson72/sphnx/gammajet/src/object.h"
R__LOAD_LIBRARY(libgammajet.so);
R__LOAD_LIBRARY(libRooUnfold.so);

void unfold(string trigger="Photon5") {
  
  bool dodraw = false;
  
  unfolder uf(trigger);
  uf.set_dodraw(dodraw);
  uf.fill_matrix();
  if (!dodraw) {
    uf.unfold();
    uf.end();
  }
}
