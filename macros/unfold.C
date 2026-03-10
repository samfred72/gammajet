#include "/home/samson72/sphnx/gammajet/src/unfolder.h"
#include "/home/samson72/sphnx/gammajet/src/object.h"
R__LOAD_LIBRARY(libgammajet.so);
R__LOAD_LIBRARY(libRooUnfold.so);

void unfold(string trigger) {
  gInterpreter->GenerateDictionary("std::vector<std::vector<float>>", "vector");
  
  unfolder uf(trigger);
  uf.fill_matrix();
  uf.end();
}
