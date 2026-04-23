SIM="pythia"
[ -n "$2" ] && SIM=$2
root -b -l -q "run.C(\"$1\",\"$SIM\")"
