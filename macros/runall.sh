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
  bash run.sh Data &
fi
if [[ $DOPHOTON == 1 ]]; then
  bash run.sh Photon5 &
  bash run.sh Photon10 &
  bash run.sh Photon20 &
  bash run.sh Photon5 herwig &
  bash run.sh Photon10 herwig &
  bash run.sh Photon20 herwig &
fi
wait
if [[ $DOJET == 1 ]]; then
  bash run.sh Jet5 &
  bash run.sh Jet8 &
  bash run.sh Jet12 &
  bash run.sh Jet20 &
  bash run.sh Jet30 &
  bash run.sh Jet50 &
  bash run.sh Jet70 &
fi
wait
echo "All Done!"
