#for i in {0..6}; do
  i=2
  root -b -l -q "grid_insitu_linear.C(\"nominal\",$i)"
  root -b -l -q "grid_insitu_linear.C(\"bdt\",$i)"
  root -b -l -q "grid_insitu_linear.C(\"iso\",$i)"
  root -b -l -q "grid_insitu_linear.C(\"3jet\",$i)"
  root -b -l -q "grid_insitu_linear.C(\"JERhigh\",$i)"
  root -b -l -q "grid_insitu_linear.C(\"JERlow\",$i)"
  root -b -l -q "grid_insitu_linear.C(\"HERWIG\",$i)"
#done
