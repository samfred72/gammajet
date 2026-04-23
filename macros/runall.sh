#/bin/bash
bash run.sh Data &
bash run.sh Photon5 &
bash run.sh Photon10 &
bash run.sh Photon20 &
bash run.sh Photon5 herwig &
bash run.sh Photon10 herwig &
bash run.sh Photon20 herwig &
wait
bash run.sh Jet5 &
bash run.sh Jet10 &
bash run.sh Jet20 &
bash run.sh Jet30 &
bash run.sh Jet50 &
bash run.sh Jet70 &
wait
echo "All Done!"
