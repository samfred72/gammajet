#/bin/bash
bash run.sh Data &
bash run.sh Photon5 &
bash run.sh Photon10 &
bash run.sh Photon20 &
wait
echo "All Done!"
