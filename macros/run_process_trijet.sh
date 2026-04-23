for i in {0..3}; do
  root -b -l -q "process_trijet.C(\"\",$i)"
  root -b -l -q "process_trijet.C(\"JERHigh\",$i)"
  root -b -l -q "process_trijet.C(\"JERLow\",$i)"
  root -b -l -q "process_trijet.C(\"HERWIG\",$i)"
done
