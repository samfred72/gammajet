#/bin/bash
bash run_unfold.sh Photon5 &
bash run_unfold.sh Photon10 &
bash run_unfold.sh Photon20 &
wait
echo "All Done!"
