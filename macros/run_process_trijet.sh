for i in {0..6}; do
  #i=5
  root -b -l -q "process_trijet.C(\"\",$i)"
  root -b -l -q "process_trijet.C(\"JERHigh\",$i)"
  root -b -l -q "process_trijet.C(\"JERLow\",$i)"
  #root -b -l -q "process_trijet.C(\"JERReco\",$i)"
  root -b -l -q "process_trijet.C(\"HERWIG\",$i)"
done
