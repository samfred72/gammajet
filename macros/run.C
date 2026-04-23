#include "/home/samson72/sphnx/gammajet/src/histmaker.h"
#include "/home/samson72/sphnx/gammajet/src/object.h"
R__LOAD_LIBRARY(libgammajet.so);

void run(string trigger, string sim = "pythia") {
//  gInterpreter->GenerateDictionary("std::vector<std::vector<float>>", "vector");
  
  histmaker hm(trigger, sim);
  hm.make_hists();
  hm.end();
}
