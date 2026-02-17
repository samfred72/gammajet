#/bin/bash
root -l -q "truthhistmaker.C(\"Photon5\")" &
root -l -q "truthhistmaker.C(\"Photon10\")" &
root -l -q "truthhistmaker.C(\"Photon20\")" &
wait
#root -l -q "histmaker.C(\"Jet5\")" &
#root -l -q "histmaker.C(\"Jet10\")" &
#root -l -q "histmaker.C(\"Jet20\")" &
#root -l -q "histmaker.C(\"Jet30\")" &
#root -l -q "histmaker.C(\"Jet50\")" &
#root -l -q "histmaker.C(\"Jet70\")" &
#wait
echo "All done!"

#root -b -l -q draw_abcd.C
#root -b -l -q draw_axj.C
