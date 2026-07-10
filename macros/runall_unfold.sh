#/bin/bash
DODATA=0
DOPHOTON=0
DOJET=0

if [[ -z $1 ]]; then
  DOPHOTON=1
  DOJET=1
  DODATA=1
elif [[ $1 == "jet" ]]; then
  DOJET=1
elif [[ $1 == "photon" ]]; then
  DOPHOTON=1
elif [[ $1 == "data" ]]; then
  DODATA=1
fi

if [[ $DODATA == 1 ]]; then
  bash run_unfold.sh Data &
fi
if [[ $DOPHOTON == 1 ]]; then
  bash run_unfold.sh Photon5 &
  bash run_unfold.sh Photon10 &
  bash run_unfold.sh Photon20 &
  #bash run_unfold.sh Photon5 herwig &
  #bash run_unfold.sh Photon10 herwig &
  #bash run_unfold.sh Photon20 herwig &
fi
wait
if [[ $DOJET == 1 ]]; then
  bash run_unfold.sh Jet5 &
  bash run_unfold.sh Jet8 &
  bash run_unfold.sh Jet12 &
  bash run_unfold.sh Jet20 &
  bash run_unfold.sh Jet30 &
  bash run_unfold.sh Jet50 &
  bash run_unfold.sh Jet70 &
fi
wait
echo "All Done!"
