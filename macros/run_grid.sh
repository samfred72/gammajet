#for i in {0..6}; do
  i=2
  root -b -l -q "grid_insitu.C(\"nominal\",$i,0)"
  root -b -l -q "grid_insitu.C(\"bdt\",$i,0)"
  root -b -l -q "grid_insitu.C(\"iso\",$i,0)"
  root -b -l -q "grid_insitu.C(\"3jet\",$i,0)"
  root -b -l -q "grid_insitu.C(\"JERhigh\",$i,0)"
  root -b -l -q "grid_insitu.C(\"JERlow\",$i,0)"
  root -b -l -q "grid_insitu.C(\"scalehigh\",$i,0)"
  root -b -l -q "grid_insitu.C(\"scalelow\",$i,0)"
  root -b -l -q "grid_insitu.C(\"HERWIG\",$i,0)"
#done
